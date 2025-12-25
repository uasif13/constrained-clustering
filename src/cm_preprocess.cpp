#include "cm_preprocess.h"
#include <igraph.h>
#include <string>
#include <unordered_map>

void output(std::vector<std::vector<long>> vec) {
    for (auto const& subvec: vec) {
        for (int i = 0; i < subvec.size(); i+=2) {
            std::cout << subvec[i] << " " << subvec[i+1] << std::endl;
        }
        std::cout << std::endl;
    }
}

void output(std::unordered_map<long, long> map) {
    for (auto const& [key,val]: map) {
        std::cout << key << " " << val << endl;
    }
}

int CMPreprocess::main() {
    this -> WriteToLogFile("Read subgraph edgelist file", Log::info);
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> subgraph_edges_vector = MMapSubgraphLoader::LoadEdgelistMMap(this->edgelist, &original_to_new_id_unordered_map,false);

    this -> WriteToLogFile("Finished reading subgraph edgelist file", Log::info);

    output(subgraph_edges_vector);
    output(original_to_new_id_unordered_map);

    // store the results into the queue that each thread pulls from
    
    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    
    int current_components_vector_index = 0;
    int cc_count = subgraph_edges_vector.size();

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
                thread_vector.push_back(std::thread(CMPreprocess::MinCutOrClusterWorkerRecursive,subgraph_edges_vector[current_components_vector_index], algorithm, 0, clustering_parameter));
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
    
    this->WriteToLogFile(std::to_string(CMPreprocess::done_being_clustered_clusters.size() - previous_done_being_clustered_size) + " [connected components / clusters] were found to be well connected", Log::debug);
    previous_done_being_clustered_size = CMPreprocess::done_being_clustered_clusters.size();
    /** SECTION MinCutOnceAndCluster Each Connected Component END **/

    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteClusterQueue(CMPreprocess::done_being_clustered_clusters);
    return 0;
}