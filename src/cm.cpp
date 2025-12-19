#include "cm.h"
#include <igraph.h>
#include <string>
#include <unordered_map>

int CM::main() {
    this->WriteToLogFile("Loading the initial graph" , Log::info);

    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_t graph;
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    MMapGraphLoader::LoadEdgelistMMap(this->edgelist, &graph,&original_to_new_id_unordered_map,false);
    SetIgraphAllEdgesWeight(&graph, 1.0);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info);
    /* std::cerr << EAN(&graph, "weight", 0) << std::endl; */
    
    int graph_vcount = igraph_vcount(&graph);
    string edge_count;
    edge_count = "before rice edge_count " + to_string(igraph_ecount(&graph));
    this -> WriteToLogFile(edge_count, Log::info);
    std::unordered_map<long, long> node_id_to_cluster_id_unordered_map;
    
    if(this->existing_clustering == "") {
        /* int seed = uni(rng); */
        // not working
        int seed = 0;
        std::map<long, long> node_id_to_cluster_id_map = ConstrainedClustering::GetCommunities("", this->algorithm, seed, this->clustering_parameter, &graph);
        // output_map(node_id_to_cluster_id_map);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_map);
    } else if(this->existing_clustering != "") {
        // non mpi
        this->WriteToLogFile("Loading the new id to cluster id map" , Log::debug);
        MMapGraphLoader::LoadClusteringMMap(this->existing_clustering, &node_id_to_cluster_id_unordered_map, original_to_new_id_unordered_map);
        this->WriteToLogFile("Finished loading the new id to cluster id map" , Log::debug);    
        // output_map(original_to_new_id_map_nonmpi);
        this->WriteToLogFile("Removing Inter cluster edges vertices: " + std::to_string(igraph_vcount(&graph)) + " edges: " + std::to_string(igraph_ecount(&graph)) , Log::debug);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_unordered_map, this-> num_processors);
        this->WriteToLogFile("Finished removing Inter cluster edges vertices: " + std::to_string(igraph_vcount(&graph)) + " edges: " + std::to_string(igraph_ecount(&graph)) , Log::debug);
    }

    /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    // store the results into the queue that each thread pulls from
    
    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    
    int current_components_vector_index = 0;
    int cc_count = connected_components_vector.size();

    /** SECTION MinCutOnceAndCluster Each Connected Component START **/
    this->WriteToLogFile(std::to_string(cc_count) + " [connected components / clusters] to be mincut", Log::debug);
    /* before_mincut_number_of_clusters = CM::to_be_mincut_clusters.size(); */
    /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
    
    /* start the threads */
    while (current_components_vector_index < cc_count) {
        /* start the threads */
        // for(int i = 0; i < this->num_processors; i ++) {
        //     MincutOnly::to_be_mincut_clusters.push({-1});
        // }
        std::vector<std::thread> thread_vector;
        for(int i = 0; i < this->num_processors; i ++) {
            // printf("current_components_vector_index: %d\n", current_components_vector_index);
            // if ( connected_components_vector[current_components_vector_index].size() == 80) {
            //     printf("cluster causes address not mapped\n");
            // }
            // else 
            if (current_components_vector_index < cc_count ) {
                thread_vector.push_back(std::thread(CM::MinCutOrClusterWorkerRecursive,connected_components_vector[current_components_vector_index], &graph, algorithm, 0, clustering_parameter));
            }
            else {
                break;
            }
            current_components_vector_index++;
        }
        /* get the result back from threads */
        /* the results from each thread gets stored in to_be_clustered_clusters */
        for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
            thread_vector[thread_index].join();
        }
    }
    
    this->WriteToLogFile(std::to_string(CM::done_being_clustered_clusters.size() - previous_done_being_clustered_size) + " [connected components / clusters] were found to be well connected", Log::debug);
    previous_done_being_clustered_size = CM::done_being_clustered_clusters.size();
    /** SECTION MinCutOnceAndCluster Each Connected Component END **/

    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteClusterQueue(CM::done_being_clustered_clusters, &graph);
    igraph_destroy(&graph);
    return 0;
}
