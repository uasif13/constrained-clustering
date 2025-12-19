#include "cm.h"

int CM::main() {
    /* std::random_device rd; */
    /* std::mt19937 rng{rd()}; */
    /* std::uniform_int_distribution<int> uni(0, 100000); */
    this->WriteToLogFile("Loading the initial graph" , Log::info);
    FILE* edgelist_file = fopen(this->edgelist.c_str(), "r");
    igraph_t graph;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_ncol(&graph, edgelist_file, NULL, 1, IGRAPH_ADD_WEIGHTS_IF_PRESENT, IGRAPH_UNDIRECTED);
    if(!igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_EDGE, "weight")) {
        SetIgraphAllEdgesWeight(&graph, 1.0);
    }
    /* igraph_read_graph_edgelist(&graph, edgelist_file, 0, false); */
    fclose(edgelist_file);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info);
    /* std::cerr << EAN(&graph, "weight", 0) << std::endl; */


    /* int before_mincut_number_of_clusters = -1; */
    int after_mincut_number_of_clusters = -2;
    int iter_count = 0;

    if(this->existing_clustering == "") {
        /* int seed = uni(rng); */
        int seed = 0;
        std::map<int, int> node_id_to_cluster_id_map = ConstrainedClustering::GetCommunities("", this->algorithm, seed, this->clustering_parameter, &graph);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_map);
    } else if(this->existing_clustering != "") {
        std::map<std::string, int> original_to_new_id_map = ConstrainedClustering::GetOriginalToNewIdMap(&graph);
        std::map<int, int> new_node_id_to_cluster_id_map = ConstrainedClustering::ReadCommunities(original_to_new_id_map, this->existing_clustering);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, new_node_id_to_cluster_id_map);
    }

    /** SECTION Get Connected Components START **/
    std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    // store the results into the queue that each thread pulls from
    
    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    
    int current_components_vector_index = 0;
    int cc_count = connected_components_vector.size();

    this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
    if(iter_count % 10 == 0) {
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
    }

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
