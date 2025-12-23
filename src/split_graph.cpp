#include "split_graph.h"

int SplitGraph::main() {
    std::map<int, int> node_id_to_cluster_id_map;
    
    std::map<int, int> cluster_id_to_new_cluster_id_map;
    std::ifstream existing_clustering_file(this -> existing_clustering);
    
    int node_id = 1;
    int cluster_id = 1;
    int cluster_id_new = 0;
    while (existing_clustering_file >> node_id >> cluster_id) {
        if (!cluster_id_to_new_cluster_id_map.contains(cluster_id)) {
            cluster_id_to_new_cluster_id_map[cluster_id] = cluster_id_new;
            cluster_id_new++;
        }
        node_id_to_cluster_id_map[node_id] = cluster_id_to_new_cluster_id_map[cluster_id];
    }
    this -> WriteToLogFile("Finished reading Communities cluster_count is: " + to_string(cluster_id_new),Log::info);
    int cluster_size = (cluster_id_new)/this -> num_partitions;
    if ((cluster_id_new)%(this -> num_partitions) != 0) {
        cluster_size ++;
    }
    std::ifstream input_graph(this->edgelist);
    int start_node = 1;
    int end_node = 1;
    std::map<int, std::ofstream> output_files;
    for (int i = 0; i < this -> num_partitions; i++) {
        output_files[i] = std::ofstream(this -> output_header + "_" + to_string(i) + ".tsv");
    }
    this -> WriteToLogFile("Finished creating output streams", Log::info);
    int count = 1;
    string line;
    while (input_graph >> start_node >> end_node) {
        if (node_id_to_cluster_id_map.contains(start_node) && node_id_to_cluster_id_map.contains(end_node)
        && node_id_to_cluster_id_map.at(start_node) == node_id_to_cluster_id_map.at(end_node)) {
            int partition = node_id_to_cluster_id_map[start_node] % this -> num_partitions;
            this -> WriteToLogFile("Edge to write: " + to_string(start_node) + "\t" + to_string(end_node) + " Partition: " + to_string(partition), Log::debug);
            output_files[partition] << start_node << "\t" << end_node << "\n";
            // std::cerr << "Intercluster edges written, count=" << count << std::endl;
            count ++;
        }
        else {
            // inter cluster edges and edges with nodes with partial cluster info
            // this -> WriteToLogFile("Edge to write: " + to_string(start_node) + "\t" + to_string(end_node) + " Partition: " + to_string(this -> num_partitions), Log::debug);

            // output_files[this->num_partitions] << start_node << "\t" << end_node << "\n";
        }
    }
    std::cerr << "Outside reading input loop" << std::endl;

    this -> WriteToLogFile("Finished writing to output streams", Log::info);
    
    // for (int i = 0; i < this -> num_partitions; i++) {
    //     output_files[i].close();
    // }
    std::cerr << "after partition deallocations" << std::endl;
    std::cerr.flush();
    return 0;
}

int SplitGraph::WriteToLogFile(std::string message, Log message_type) {
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
        this->log_file_handle << log_message_prefix << " " << message << '\n';

        if(this->num_calls_to_log_write % 10 == 0) {
            std::flush(this->log_file_handle);
        }
        this->num_calls_to_log_write ++;
    }
    return 0;
}

