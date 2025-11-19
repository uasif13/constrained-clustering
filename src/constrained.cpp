#include "constrained.h"
#include "mpi_telemetry.h"
#include "igraph.h"

void build_displacements_output_file(int * displacements, int* size_array, int nprocs) {
    displacements[0] = 0;
    for (int i = 1; i < nprocs; i ++) {
        displacements[i] = displacements[i-1] + size_array[i-1];
    }
}

int update_cluster_id_array(int * cluster_ids, int cluster_size, int previous_cluster_id, int* cluster_displacements, int nprocs)  {
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

std::map<int, std::string> ConstrainedClustering::InvertMap(const std::map<std::string, int>& original_to_new_id_map) {
    std::map<int, std::string> new_to_original_id_map;
    for(auto const& [original_id, new_id] : original_to_new_id_map) {
        new_to_original_id_map[new_id] = original_id;
    }
    return new_to_original_id_map;
}

std::map<std::string, int> ConstrainedClustering::GetOriginalToNewIdMap(std::string edgelist) {
    std::map<std::string, int> original_to_new_id_map;
    char delimiter = get_delimiter(edgelist);
    std::ifstream edgelist_file(edgelist);
    std::string line;
    int line_no = 0;
    int next_node_id = 0;
    while(std::getline(edgelist_file, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        if(line_no == 0) {
            line_no ++;
            continue;
        }
        std::string source = current_line[0];
        std::string target = current_line[1];
        if (!original_to_new_id_map.contains(source)) {
            original_to_new_id_map[source] = next_node_id;
            next_node_id ++;
        }
        if (!original_to_new_id_map.contains(target)) {
            original_to_new_id_map[target] = next_node_id;
            next_node_id ++;
        }
        line_no ++;
    }
    this->num_edges = line_no - 1;
    return original_to_new_id_map;
}

void ConstrainedClustering::LoadEdgesFromFile(igraph_t* graph, std::string edgelist, const std::map<std::string, int>& original_to_new_id_map) {
    /* igraph_setup(); */
    igraph_add_vertices(graph, original_to_new_id_map.size(), NULL);
    char delimiter = get_delimiter(edgelist);
    std::ifstream edgelist_file(edgelist);
    std::string line;
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, this->num_edges * 2);
    int line_no = 0;
    int edge_index = 0;
    while(std::getline(edgelist_file, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        if(line_no == 0) {
            line_no ++;
            continue;
        }
        std::string source = current_line[0];
        std::string target = current_line[1];
        /* igraph_add_edge(graph, original_to_new_id_map.at(source), original_to_new_id_map.at(target)); */
        VECTOR(edges)[edge_index] = original_to_new_id_map.at(source);
        edge_index ++;
        VECTOR(edges)[edge_index] = original_to_new_id_map.at(target);
        edge_index ++;
        line_no ++;

    }
    igraph_add_edges(graph, &edges, NULL);
    igraph_vector_int_destroy(&edges);
}

void ConstrainedClustering::WriteClusterQueue(std::queue<std::vector<int>>& cluster_queue, igraph_t* graph, const std::map<int, std::string>& new_to_original_id_map) {
    std::ofstream clustering_output(this->output_file);
    int current_cluster_id = 0;
    this->WriteToLogFile("final clusters:", Log::debug);
    clustering_output << "node_id,cluster_id" << '\n';
    while(!cluster_queue.empty()) {
        std::vector<int> current_cluster = cluster_queue.front();
        cluster_queue.pop();
        this->WriteToLogFile("new cluster", Log::debug);
        for(size_t i = 0; i < current_cluster.size(); i ++) {
            this->WriteToLogFile(std::to_string(current_cluster[i]), Log::debug);
            /* clustering_output << VAS(graph, "name", current_cluster[i]) << "," << current_cluster_id << '\n'; */
            clustering_output << new_to_original_id_map.at(current_cluster[i]) << "," << current_cluster_id << '\n';
        }
        current_cluster_id ++;
    }
    clustering_output.close();
}

int ConstrainedClustering::WriteClusterQueueMPI(std::queue<std::vector<int>>* cluster_queue, igraph_t* graph, std::map<int, std::string>* new_to_original_id_map, int cluster_start, int previous_cluster_id, int iteration, uint64_t* opCount) {
    // std::ofstream clustering_output;
    // clustering_output.open(this->output_file, std::ios_base::app);
    int current_cluster_id = cluster_start;

    int v_count = igraph_vcount(graph);
    int node_id_arr[v_count];
    int cluster_id_arr[v_count];
    int index_count = 0;

    // Write to individual cluster files
    this->WriteToLogFile("final clusters:", Log::debug, my_rank);

    while(!cluster_queue->empty()) {
        std::vector<int> current_cluster = cluster_queue->front();
        cluster_queue->pop();
        this->WriteToLogFile("new cluster size: " + std::to_string(current_cluster.size()), Log::debug, my_rank);
        for(size_t i = 0; i < current_cluster.size(); i ++) {
            this->WriteToLogFile(std::to_string(current_cluster[i]), Log::debug, my_rank);
            // clustering_output << VAS(graph, "name", current_cluster[i]) << " " << current_cluster_id << '\n';
            node_id_arr[index_count] = stoi(new_to_original_id_map ->at(current_cluster[i]));
            cluster_id_arr[index_count] = current_cluster_id;
            index_count++;

        }
        current_cluster_id ++;
    }
    // Write to aggregate cluster file
    int index_count_arr[nprocs];
    int cluster_displacements[nprocs];
    int node_cluster_id_agg_size;

    this -> WriteToLogFile("Write to MPI Output", Log::debug, my_rank);
    // Send index counts to ROOT
    if (my_rank == 0) {
      MPI_Gather(&index_count, 1, MPI_INT, index_count_arr, 1, MPI_INT, 0, MPI_COMM_WORLD, my_rank, iteration, 3, opCount);
      //output_arr(index_count_arr, nprocs);
      build_displacements_output_file(cluster_displacements, index_count_arr, nprocs);
      //output_arr(cluster_displacements, nprocs);
      node_cluster_id_agg_size = cluster_displacements[nprocs-1]+index_count_arr[nprocs-1];

    } else {
      MPI_Gather(&index_count, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD, my_rank, iteration, 3, opCount);
    }
    MPI_Barrier(my_rank, iteration, 5, opCount);
    

   // Send node_ids and cluster_ids to ROOT
    if (my_rank == 0) {
        this -> WriteToLogFile("my_rank: " + std::to_string(my_rank) + " Write to cluster mpi file", Log::info, my_rank);
        std::ofstream mpi_clustering_output;
        mpi_clustering_output.open(this->output_file, std::ios_base::app);
        int node_id_arr_agg[node_cluster_id_agg_size];
        int cluster_id_arr_agg[node_cluster_id_agg_size];
        MPI_Gatherv(node_id_arr, index_count, MPI_INT,node_id_arr_agg, index_count_arr, cluster_displacements, MPI_INT, 0, MPI_COMM_WORLD, my_rank, iteration, 4, opCount);
        MPI_Gatherv(cluster_id_arr, index_count, MPI_INT, cluster_id_arr_agg, index_count_arr, cluster_displacements, MPI_INT, 0, MPI_COMM_WORLD, my_rank, iteration, 4, opCount);
        previous_cluster_id = update_cluster_id_array(cluster_id_arr_agg, node_cluster_id_agg_size, previous_cluster_id, cluster_displacements, nprocs);
        for (int i = 0; i < node_cluster_id_agg_size; i++) {
            mpi_clustering_output << node_id_arr_agg[i] << "," << cluster_id_arr_agg[i] << "\n";
        }
        mpi_clustering_output.close();
    } else {
      MPI_Gatherv(node_id_arr, index_count, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD, my_rank, iteration, 4, opCount);
      MPI_Gatherv(cluster_id_arr, index_count, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD, my_rank, iteration, 4, opCount);
    }

    // Broadcast previous_cluster_id
    MPI_Bcast(&previous_cluster_id,1, MPI_INT, 0, MPI_COMM_WORLD, my_rank, iteration, 0, opCount);
    return previous_cluster_id;
}

void ConstrainedClustering::WritePartitionMap(std::map<int, int>& final_partition) {
    std::ofstream clustering_output(this->output_file);
    for(auto const& [node_id, cluster_id]: final_partition) {
        clustering_output << node_id << "," << cluster_id << '\n';
    }
    clustering_output.close();
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
