#include "split_graph.h"

int SplitGraph::main() {
    this->WriteToLogFile("Loading the initial graph" , Log::info);

    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_t graph;
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    MMapGraphLoader::LoadEdgelistMMap(this->edgelist, &graph,&original_to_new_id_unordered_map,false);
    ConstrainedClustering::SetIgraphAllEdgesWeight(&graph, 1.0);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info);
    // /* std::cerr << EAN(&graph, "weight", 0) << std::endl; */
    
    int graph_vcount = igraph_vcount(&graph);
    string edge_count;
    edge_count = "before rice edge_count " + to_string(igraph_ecount(&graph));
    this -> WriteToLogFile(edge_count, Log::info);
    std::unordered_map<long, long> node_id_to_cluster_id_unordered_map;
        // non mpi
    this->WriteToLogFile("Loading the new id to cluster id map" , Log::debug);
    MMapGraphLoader::LoadClusteringMMap(this->existing_clustering, &node_id_to_cluster_id_unordered_map, original_to_new_id_unordered_map);
    this->WriteToLogFile("Finished loading the new id to cluster id map" , Log::debug);    
    // output_map(original_to_new_id_map_nonmpi);
    this->WriteToLogFile("Removing Inter cluster edges vertices: " + std::to_string(igraph_vcount(&graph)) + " edges: " + std::to_string(igraph_ecount(&graph)) , Log::debug);
    ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_unordered_map, this-> num_processors);

    // this->WriteToLogFile("Finished removing Inter cluster edges vertices: " + std::to_string(igraph_vcount(&graph)) + " edges: " + std::to_string(igraph_ecount(&graph)) , Log::debug);

    // /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    std::map<int, std::ofstream> output_files;
    for (int i = 0; i < this -> num_partitions; i++) {
        output_files[i] = std::ofstream(this -> output_header + "_" + to_string(this->num_partitions) + "_" + to_string(i) + ".tsv");
    }
    for (int i = 0; i < connected_components_vector.size(); i++) {
        int partition = i % this -> num_partitions;
        igraph_vector_int_t nodes_to_keep;
        igraph_vector_int_t new_id_to_old_id_vector_map;
        igraph_vector_int_init(&nodes_to_keep, connected_components_vector[i].size());
        for(size_t j = 0; j < connected_components_vector[i].size(); j ++) {
            VECTOR(nodes_to_keep)[j] = connected_components_vector[i][j];
        }
        igraph_t induced_subgraph;
        igraph_vector_int_init(&new_id_to_old_id_vector_map, igraph_vector_int_size(&nodes_to_keep));
        igraph_induced_subgraph_map(&graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_vector_map);
        igraph_eit_t eit;
        igraph_eit_create(&induced_subgraph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
        for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
            long from_node = VECTOR(new_id_to_old_id_vector_map)[IGRAPH_FROM(&induced_subgraph, current_edge)];
            long to_node = VECTOR(new_id_to_old_id_vector_map)[IGRAPH_TO(&induced_subgraph, current_edge)];
            output_files[partition] << from_node << " " << to_node << "\n";
        }
        output_files[partition] << "-\n";
    }
    return 1;
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

