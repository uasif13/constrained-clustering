#include "constrained.h"
#include "mpi.h"


template <typename T>
void output_vec(std::vector<T> vec) {
    for (auto item: vec) {
        cout << item << " ";
    }
    cout << "\n";
}

template <typename T>
void output_arr(T * arr, int arr_size) {
    for (int i = 0 ; i < arr_size; i++) {
        cout << arr[i] << " ";
    }
    cout << "\n";
}

void build_displacements_output_file(int * displacements, int* size_array, int nprocs) {
    displacements[0] = 0;
    for (int i = 1; i < nprocs; i ++) {
        displacements[i] = displacements[i-1] + size_array[i-1];
    }
}

void ConstrainedClustering::WriteClusterQueue(std::queue<std::vector<int>>& cluster_queue, igraph_t* graph) {
    std::ofstream clustering_output(this->output_file);
    int current_cluster_id = 0;
    this->WriteToLogFile("final clusters:", Log::debug);
    int v_count = igraph_vcount(graph);
    int node_id_arr[v_count];
    int cluster_id_arr[v_count];
    int index_count = 0;
    int index_count_arr[nprocs];
    int cluster_displacements[nprocs];
    int node_cluster_id_agg_size;
    while(!cluster_queue.empty()) {
        std::vector<int> current_cluster = cluster_queue.front();
        cluster_queue.pop();
        output_vec(current_cluster);
        this->WriteToLogFile("new cluster", Log::debug);
        for(size_t i = 0; i < current_cluster.size(); i ++) {
            this->WriteToLogFile(std::to_string(current_cluster[i]), Log::debug);
            clustering_output << VAS(graph, "name", current_cluster[i]) << " " << current_cluster_id << '\n';
            node_id_arr[index_count] = stoi(VAS(graph, "name", current_cluster[i]));
            cluster_id_arr[index_count] = current_cluster_id;
            index_count++;
        }
        current_cluster_id ++;
    }
    clustering_output.close();
    
}

void ConstrainedClustering::WriteClusterQueue(std::queue<std::vector<int>>& cluster_queue, igraph_t* graph, int cluster_start) {
    std::ofstream clustering_output(this->output_file);
    int current_cluster_id = cluster_start;
    this->WriteToLogFile("final clusters:", Log::debug);
    int v_count = igraph_vcount(graph);
    int node_id_arr[v_count];
    int cluster_id_arr[v_count];
    int index_count = 0;
    int index_count_arr[nprocs];
    int cluster_displacements[nprocs];
    int node_cluster_id_agg_size;
    while(!cluster_queue.empty()) {
        std::vector<int> current_cluster = cluster_queue.front();
        cluster_queue.pop();
        output_vec(current_cluster);
        this->WriteToLogFile("new cluster", Log::debug);
        for(size_t i = 0; i < current_cluster.size(); i ++) {
            this->WriteToLogFile(std::to_string(current_cluster[i]), Log::debug);
            clustering_output << VAS(graph, "name", current_cluster[i]) << " " << current_cluster_id << '\n';
            node_id_arr[index_count] = stoi(VAS(graph, "name", current_cluster[i]));
            cluster_id_arr[index_count] = current_cluster_id;
            index_count++;
        }
        current_cluster_id ++;
    }
    clustering_output.close();
    std::ofstream mpi_clustering_output("./output_clusters_mpi.tsv");

    MPI_Allgather(&index_count, 1, MPI_INT, index_count_arr, 1, MPI_INT, MPI_COMM_WORLD);
    output_arr(index_count_arr, nprocs);
    build_displacements_output_file(cluster_displacements, index_count_arr, nprocs);
    output_arr(cluster_displacements, nprocs);
    node_cluster_id_agg_size = cluster_displacements[nprocs-1]+index_count_arr[nprocs-1];
    int node_id_arr_agg[node_cluster_id_agg_size];
    int cluster_id_arr_agg[node_cluster_id_agg_size];
    MPI_Allgatherv(node_id_arr, index_count, MPI_INT, node_id_arr_agg, index_count_arr, cluster_displacements, MPI_INT, MPI_COMM_WORLD );
    MPI_Allgatherv(cluster_id_arr, index_count, MPI_INT, cluster_id_arr_agg, index_count_arr, cluster_displacements, MPI_INT, MPI_COMM_WORLD );
    output_arr(node_id_arr_agg, node_cluster_id_agg_size);
    output_arr(cluster_id_arr_agg, node_cluster_id_agg_size);
    if (my_rank == 0) {
        for (int i = 0; i < node_cluster_id_agg_size; i++) {
            mpi_clustering_output << node_id_arr_agg[i] << " " << cluster_id_arr_agg[i] << "\n";
        }
    }
    mpi_clustering_output.close();
    
}

void ConstrainedClustering::WritePartitionMap(std::map<int, int>& final_partition) {
    std::ofstream clustering_output(this->output_file);
    for(auto const& [node_id, cluster_id]: final_partition) {
        clustering_output << node_id << " " << cluster_id << '\n';
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
