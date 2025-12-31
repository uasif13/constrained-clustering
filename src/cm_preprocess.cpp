#include "cm_preprocess.h"

int CMPreprocess::main() {
    this -> WriteToLogFile("Read subgraph edgelist file", Log::info);
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> subgraph_edges_vector = MMapSubgraphLoader::LoadEdgelistMMap(this->edgelist, &original_to_new_id_unordered_map,false);

    this -> WriteToLogFile("Finished reading subgraph edgelist file", Log::info);

    // store the results into the queue that each thread pulls from
    
    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    
    int current_components_vector_index = 0;
    int cc_count = subgraph_edges_vector.size();

    /** SECTION MinCutOnceAndCluster Each Connected Component START **/
    this->WriteToLogFile(std::to_string(cc_count) + " [connected components / clusters] to be mincut", Log::debug);
    // store the results into the queue that each thread pulls from
    for(size_t i = 0; i < subgraph_edges_vector.size(); i ++) {
        CMPreprocess::to_be_mincut_clusters.push(subgraph_edges_vector[i]);
    }
    
    /** SECTION Get Connected Components END **/
    long iter_count = 0;
    while (true) {
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
        if(iter_count % 10 == 0) {
            this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
        }

        /** SECTION MinCutOnceAndCluster Each Connected Component START **/
        this->WriteToLogFile(std::to_string(CMPreprocess::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::debug);
        /* before_mincut_number_of_clusters = CM::to_be_mincut_clusters.size(); */
        /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */

        /* start the threads */
        std::vector<std::thread> thread_vector;
        for(int i = 0; i < this->num_processors; i ++) {
            /* int seed = uni(rng); */
            int seed = 0;
            thread_vector.push_back(std::thread(CMPreprocess::MinCutOrClusterWorker, this->algorithm, seed, this->clustering_parameter));
        }

        /* get the result back from threads */
        /* the results from each thread gets stored in to_be_clustered_clusters */
        for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
            thread_vector[thread_index].join();
        }
        this->WriteToLogFile(std::to_string(CMPreprocess::to_be_clustered_clusters.size()) + " [connected components / clusters] to be clustered after a round of mincuts", Log::debug);
        this->WriteToLogFile(std::to_string(CMPreprocess::done_being_clustered_clusters.size() - previous_done_being_clustered_size) + " [connected components / clusters] were found to be well connected", Log::debug);
        previous_done_being_clustered_size = CMPreprocess::done_being_clustered_clusters.size();
        /** SECTION MinCutOnceAndCluster Each Connected Component END **/

        /** SECTION Check If All Clusters Are Well-Connected START **/
        long after_mincut_number_of_clusters = CMPreprocess::to_be_mincut_clusters.size();
        if(after_mincut_number_of_clusters == 0) {
            this->WriteToLogFile("all clusters are well-connected", Log::info);
            this->WriteToLogFile("Total number of iterations: " + std::to_string(iter_count + 1), Log::info);
            break;
        }
        /** SECTION Check If All Clusters Are Well-Connected END **/

        while(!CMPreprocess::to_be_clustered_clusters.empty()) {
            CMPreprocess::to_be_mincut_clusters.push(CMPreprocess::to_be_clustered_clusters.front());
            CMPreprocess::to_be_clustered_clusters.pop();
        }

        iter_count ++;
    }


    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteClusterQueue(CMPreprocess::done_being_clustered_clusters);
    return 0;
}
