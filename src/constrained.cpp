#include "constrained.h"
#include "mpi_telemetry.h"
#include "igraph.h"

void build_displacements_output_file(long * displacements, long* size_array, int nprocs) {
    displacements[0] = 0;
    for (int i = 1; i < nprocs; i ++) {
        displacements[i] = displacements[i-1] + size_array[i-1];
    }
}

int update_cluster_id_array(long * cluster_ids, int cluster_size, int previous_cluster_id, long* cluster_displacements, int nprocs)  {
  int current_cluster = 0;
  int updated_cluster = previous_cluster_id;
  int count = 1;
  for (int i = 0; i < cluster_size; i++) {
    if (cluster_ids[i] > current_cluster) {
        current_cluster = cluster_ids[i];
        updated_cluster++;
    } else if (i != 0  && nprocs > 1 && i == cluster_displacements[count]) {
        count ++;
        updated_cluster++;
        current_cluster = 0;
        // current_cluster = updated_cluster;
    }
    cluster_ids[i] =  updated_cluster;
  }
  if (cluster_size > 0) {
    updated_cluster++;
  }
  return updated_cluster;
}

long ConstrainedClustering::WriteClusterQueueMPI(std::queue<std::vector<long>>* cluster_queue, uint64_t* opCount, int my_rank, int nprocs) {
    this->WriteToLogFile("========================================", Log::debug, my_rank);
    this->WriteToLogFile("Starting WriteClusterQueueMPI", Log::info, my_rank);
    this->WriteToLogFile("Initial queue size: " + std::to_string(cluster_queue->size()), Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::debug, my_rank);
    
    int current_cluster_id = 0;
    
    // Use heap allocation instead of VLAs
    std::vector<long> node_id_arr;
    std::vector<long> cluster_id_arr;
    long index_count = 0;

    // Write to individual cluster files
    this->WriteToLogFile("Processing local clusters for rank " + std::to_string(my_rank), Log::debug, my_rank);

    while(!cluster_queue->empty()) {
        std::vector<long> current_cluster = cluster_queue->front();
        cluster_queue->pop();
        // this->WriteToLogFile("new cluster size: " + std::to_string(current_cluster.size()), Log::debug, my_rank);
        for(size_t i = 0; i < current_cluster.size(); i ++) {
            // this->WriteToLogFile(VAS(graph, "name", current_cluster[i]), Log::debug, my_rank);
            node_id_arr.push_back(current_cluster[i]);
            cluster_id_arr.push_back(current_cluster_id);
            index_count++;
        }
        current_cluster_id++;
    }
    
    this->WriteToLogFile("Processed " + std::to_string(current_cluster_id) + " local clusters", Log::info, my_rank);
    this->WriteToLogFile("Total nodes in local clusters: " + std::to_string(index_count), Log::info, my_rank);
    
    // Write to aggregate cluster file
    std::vector<long> index_count_arr(nprocs);
    std::vector<long> cluster_displacements(nprocs);
    long node_cluster_id_agg_size;

    this->WriteToLogFile("========================================", Log::debug, my_rank);
    this->WriteToLogFile("Starting MPI aggregation across " + std::to_string(nprocs) + " ranks", Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::debug, my_rank);
    
    // Send index counts to ROOT
    this->WriteToLogFile("Phase 1: Gathering node counts from all ranks", Log::debug, my_rank);
    if (my_rank == 0) {
        this->WriteToLogFile("ROOT: Gathering index counts from all " + std::to_string(nprocs) + " ranks", Log::debug, my_rank);
        MPI_Gather(&index_count, 1, MPI_LONG, index_count_arr.data(), 1, MPI_LONG, 0, MPI_COMM_WORLD, my_rank, -1, 3, opCount);
        
        this->WriteToLogFile("ROOT: Received counts from all ranks", Log::debug, my_rank);
        for (int i = 0; i < nprocs; i++) {
            this->WriteToLogFile("  Rank " + std::to_string(i) + ": " + std::to_string(index_count_arr[i]) + " nodes", Log::debug, my_rank);
        }
        
        build_displacements_output_file(cluster_displacements.data(), index_count_arr.data(), nprocs);
        node_cluster_id_agg_size = cluster_displacements[nprocs-1] + index_count_arr[nprocs-1];
        
        this->WriteToLogFile("ROOT: Total aggregated size: " + std::to_string(node_cluster_id_agg_size) + " nodes", Log::info, my_rank);
    } else {
        this->WriteToLogFile("Sending count (" + std::to_string(index_count) + ") to ROOT", Log::debug, my_rank);
        MPI_Gather(&index_count, 1, MPI_LONG, NULL, 0, MPI_LONG, 0, MPI_COMM_WORLD, my_rank, -1, 3, opCount);
    }
    
    this->WriteToLogFile("Waiting at barrier after count gathering", Log::debug, my_rank);
    MPI_Barrier(my_rank, -1, 5, opCount);
    this->WriteToLogFile("Passed barrier, proceeding to data gathering", Log::debug, my_rank);
    

    // Send node_ids and cluster_ids to ROOT
    this->WriteToLogFile("Phase 2: Gathering node IDs and cluster IDs", Log::debug, my_rank);
    if (my_rank == 0) {
        this->WriteToLogFile("ROOT: Preparing to write aggregated clusters to file", Log::info, my_rank);
        this->WriteToLogFile("ROOT: Output file: " + this->output_file, Log::info, my_rank);
        
        std::ofstream mpi_clustering_output;
        mpi_clustering_output.open(this->output_file);
        
        // Use heap allocation for aggregated arrays
        std::vector<long> node_id_arr_agg(node_cluster_id_agg_size);
        std::vector<long> cluster_id_arr_agg(node_cluster_id_agg_size);
        
        this->WriteToLogFile("ROOT: Gathering node IDs from all ranks", Log::debug, my_rank);
        MPI_Gatherv(node_id_arr.data(), index_count, MPI_LONG, node_id_arr_agg.data(), index_count_arr.data(), cluster_displacements.data(), MPI_LONG, 0, MPI_COMM_WORLD, my_rank, -1, 4, opCount);
        
        this->WriteToLogFile("ROOT: Gathering cluster IDs from all ranks", Log::debug, my_rank);
        MPI_Gatherv(cluster_id_arr.data(), index_count, MPI_LONG, cluster_id_arr_agg.data(), index_count_arr.data(), cluster_displacements.data(), MPI_LONG, 0, MPI_COMM_WORLD, my_rank, -1, 4, opCount);
        
        this->WriteToLogFile("ROOT: Updating cluster IDs for consistency across ranks", Log::debug, my_rank);
        update_cluster_id_array(cluster_id_arr_agg.data(), node_cluster_id_agg_size, 0, cluster_displacements.data(), nprocs);
        
        this->WriteToLogFile("ROOT: Writing " + std::to_string(node_cluster_id_agg_size) + " node-cluster pairs to file", Log::info, my_rank);
        for (long i = 0; i < node_cluster_id_agg_size; i++) {
            mpi_clustering_output << node_id_arr_agg[i] << " " << cluster_id_arr_agg[i] << "\n";
        }
        mpi_clustering_output.close();
        this->WriteToLogFile("ROOT: Successfully wrote output file", Log::info, my_rank);
    } else {
        this->WriteToLogFile("Sending " + std::to_string(index_count) + " node IDs to ROOT", Log::debug, my_rank);
        MPI_Gatherv(node_id_arr.data(), index_count, MPI_LONG, NULL, NULL, NULL, MPI_LONG, 0, MPI_COMM_WORLD, my_rank, -1, 4, opCount);
        
        this->WriteToLogFile("Sending " + std::to_string(index_count) + " cluster IDs to ROOT", Log::debug, my_rank);
        MPI_Gatherv(cluster_id_arr.data(), index_count, MPI_LONG, NULL, NULL, NULL, MPI_LONG, 0, MPI_COMM_WORLD, my_rank, -1, 4, opCount);
    }

    this->WriteToLogFile("========================================", Log::debug, my_rank);
    this->WriteToLogFile("Completed WriteClusterQueueMPI", Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::debug, my_rank);
    
    // Broadcast previous_cluster_id
    // MPI_Bcast(&previous_cluster_id, 1, MPI_LONG, 0, MPI_COMM_WORLD, my_rank, -1, 0, opCount);
    return -1;
}

int ConstrainedClustering::WriteToLogFile(std::string message, Log message_type) {
    if(this->log_level >= message_type) {
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        std::string log_message_prefix;
        if(message_type == Log::info) {
            log_message_prefix = "[INFO]";
        } else if(message_type == Log::debug) {
            log_message_prefix = "[DEBUG]";
        } else if(message_type == Log::error) {
            log_message_prefix = "[ERROR]";
        }
        auto days_elapsed = std::chrono::duration_cast<std::chrono::days>(now - this->start_time);
        auto hours_elapsed = std::chrono::duration_cast<std::chrono::hours>(now - this->start_time - days_elapsed);
        auto minutes_elapsed = std::chrono::duration_cast<std::chrono::minutes>(now - this->start_time - days_elapsed - hours_elapsed);
        auto seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time - days_elapsed - hours_elapsed - minutes_elapsed);
        auto total_seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time);
        log_message_prefix += "[";
        log_message_prefix += std::to_string(days_elapsed.count());
        log_message_prefix += "-";
        log_message_prefix += std::to_string(hours_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(minutes_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(seconds_elapsed.count());
        log_message_prefix += "]";

        log_message_prefix += "(t=";
        log_message_prefix += std::to_string(total_seconds_elapsed.count());
        log_message_prefix += "s)";
        this->log_file_handle << log_message_prefix << " " << message << '\n';

        if(this->num_calls_to_log_write % 10 == 0) {
            std::flush(this->log_file_handle);
        }
        this->num_calls_to_log_write ++;
    }
    return 0;
}

int ConstrainedClustering::WriteToLogFile(std::string message, Log message_type, int my_rank) {
    if(this->log_level >= message_type) {
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        std::string log_message_prefix;
        if(message_type == Log::info) {
            log_message_prefix = "[INFO]";
        } else if(message_type == Log::debug) {
            log_message_prefix = "[DEBUG]";
        } else if(message_type == Log::error) {
            log_message_prefix = "[ERROR]";
        }
        auto days_elapsed = std::chrono::duration_cast<std::chrono::days>(now - this->start_time);
        auto hours_elapsed = std::chrono::duration_cast<std::chrono::hours>(now - this->start_time - days_elapsed);
        auto minutes_elapsed = std::chrono::duration_cast<std::chrono::minutes>(now - this->start_time - days_elapsed - hours_elapsed);
        auto seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time - days_elapsed - hours_elapsed - minutes_elapsed);
        auto total_milliseconds_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - this->start_time);
        log_message_prefix += "[";
        log_message_prefix += std::to_string(days_elapsed.count());
        log_message_prefix += "-";
        log_message_prefix += std::to_string(hours_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(minutes_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(seconds_elapsed.count());
        log_message_prefix += "]";

        log_message_prefix += "(t=";
        log_message_prefix += std::to_string(total_milliseconds_elapsed.count());
        log_message_prefix += "ms)";
        this->log_file_handle <<"my_rank: " << my_rank << " " << log_message_prefix << " " << message << '\n';

        if(this->num_calls_to_log_write % 10 == 0) {
            std::flush(this->log_file_handle);
        }
        this->num_calls_to_log_write ++;
    }
    return 0;
}