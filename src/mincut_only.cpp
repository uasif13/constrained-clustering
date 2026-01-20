#include "mincut_only.h"

int MincutOnly::main() {

    std::cerr << "\n========================================" << std::endl;
    std::cerr << "INSTRUMENTED VERSION - DETAILED DIAGNOSTICS ENABLED" << std::endl;
    std::cerr << "========================================\n" << std::endl;

    this->WriteToLogFile("Parsing connectedness criterion" , Log::info);
    
    double connectedness_criterion_c = this->connectedness_criterion_c;
    double connectedness_criterion_x = this->connectedness_criterion_x;
    double pre_computed_log = this->pre_computed_log;
    ConnectednessCriterion current_connectedness_criterion = this->current_connectedness_criterion;
    
    std::cerr << "[MAIN] Connectedness Criterion: C=" << connectedness_criterion_c 
              << " x=" << connectedness_criterion_x 
              << " mode=" << (int)current_connectedness_criterion << std::endl;
    
    this->WriteToLogFile("Loading the initial graph" , Log::info);

    std::map<std::string, int> original_to_new_id_map = this->GetOriginalToNewIdMap(this->edgelist);
    std::map<int, std::string> new_to_originial_id_map = this->InvertMap(original_to_new_id_map);
    igraph_t graph;
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    this->LoadEdgesFromFile(&graph, this->edgelist, original_to_new_id_map);
    
    std::cerr << "[MAIN] Graph loaded: vertices=" << igraph_vcount(&graph) 
              << " edges=" << igraph_ecount(&graph) << std::endl;
    
    this->WriteToLogFile("Finished loading the initial graph vertex count: " +std::to_string(igraph_vcount(&graph)) + " edge count: " + std::to_string(igraph_ecount(&graph))  , Log::info);

    int before_mincut_number_of_clusters = -1;
    int after_mincut_number_of_clusters = -2;
    int iter_count = 0;

    std::map<int, int> new_node_id_to_cluster_id_map = ConstrainedClustering::ReadCommunities(original_to_new_id_map, this->existing_clustering);
    ConstrainedClustering::RemoveInterClusterEdges(&graph, new_node_id_to_cluster_id_map);

    std::cerr << "[MAIN] After inter-cluster edge removal: vertices=" << igraph_vcount(&graph) 
              << " edges=" << igraph_ecount(&graph) << std::endl;

    /** SECTION Get Connected Components START **/
    std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    
    std::cerr << "[MAIN] Initial connected components: " << connected_components_vector.size() << std::endl;
    
    // Calculate initial node count for verification
    uint64_t initial_total_nodes = 0;
    for (const auto& cc : connected_components_vector) {
        initial_total_nodes += cc.size();
    }
    std::cerr << "[MAIN] Total nodes in initial CCs: " << initial_total_nodes << std::endl;
    
    this -> WriteToLogFile("Connected components count: " + std::to_string(connected_components_vector.size()), Log::info);
    /** SECTION Get Connected Components END **/
    
    if(current_connectedness_criterion == ConnectednessCriterion::Simple) {
        std::cerr << "[MAIN] Simple mode: all CCs marked as complete" << std::endl;
        for(size_t i = 0; i < connected_components_vector.size(); i ++) {
            MincutOnly::done_being_mincut_clusters.push(connected_components_vector[i]);
        }
    } else {
        // store the results into the queue that each thread pulls from
        std::cerr << "[MAIN] Queueing " << connected_components_vector.size() << " initial CCs for processing" << std::endl;
        for(size_t i = 0; i < connected_components_vector.size(); i ++) {
            MincutOnly::to_be_mincut_clusters.push(connected_components_vector[i]);
        }
        
        std::cerr << "\n[MAIN] ======== Starting Iterative Processing ========\n" << std::endl;
        
        while (true) {
            std::cerr << "\n========================================" << std::endl;
            std::cerr << "[MAIN] Iteration #" << iter_count << std::endl;
            std::cerr << "========================================" << std::endl;
            
            this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
            if(iter_count % 10000 == 0) {
                this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
                this->WriteToLogFile(std::to_string(MincutOnly::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::info);
            }

            /** SECTION MinCut Each Connected Component START **/
            before_mincut_number_of_clusters = MincutOnly::to_be_mincut_clusters.size();
            
            std::cerr << "[MAIN] Queue state BEFORE iteration:" << std::endl;
            std::cerr << "  to_be_mincut: " << before_mincut_number_of_clusters << std::endl;
            std::cerr << "  done_being_mincut: " << MincutOnly::done_being_mincut_clusters.size() << std::endl;
            std::cerr << "  Global counters:" << std::endl;
            std::cerr << "    dequeued: " << total_clusters_dequeued << std::endl;
            std::cerr << "    enqueued: " << total_clusters_enqueued << std::endl;
            std::cerr << "    completed: " << total_clusters_completed << std::endl;
            std::cerr << "    singletons_discarded: " << total_singleton_partitions_discarded << std::endl;
            std::cerr << "    node_violations: " << total_node_conservation_violations << std::endl;
            
            // INVARIANT CHECK
            uint64_t expected_in_queue = total_clusters_enqueued - total_clusters_dequeued + connected_components_vector.size();
            uint64_t actual_in_queue = before_mincut_number_of_clusters;
            if (expected_in_queue != actual_in_queue) {
                std::cerr << "  **QUEUE ACCOUNTING MISMATCH**" << std::endl;
                std::cerr << "    Expected in queue: " << expected_in_queue << std::endl;
                std::cerr << "    Actual in queue: " << actual_in_queue << std::endl;
            }
            
            this->WriteToLogFile(std::to_string(before_mincut_number_of_clusters) + " [connected components / clusters] to be mincut", Log::debug);
            
            /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
            if(before_mincut_number_of_clusters > 1) {
                std::cerr << "[MAIN] Spawning " << this->num_processors << " worker threads" << std::endl;
                
                /* start the threads */
                for(int i = 0; i < this->num_processors; i ++) {
                    MincutOnly::to_be_mincut_clusters.push({-1});
                }
                std::vector<std::thread> thread_vector;
                for(int i = 0; i < this->num_processors; i ++) {
                    thread_vector.push_back(std::thread(MincutOnly::MinCutWorker, &graph, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, this->mincut_type));
                }
                
                std::cerr << "[MAIN] Waiting for workers..." << std::endl;
                
                /* get the result back from threads */
                for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
                    thread_vector[thread_index].join();
                }
                
                std::cerr << "[MAIN] All workers completed" << std::endl;
            } else if (before_mincut_number_of_clusters == 1) {
                std::cerr << "[MAIN] Single cluster remaining, processing with single worker" << std::endl;
                MincutOnly::to_be_mincut_clusters.push({-1});
                MincutOnly::MinCutWorker(&graph, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, this->mincut_type);
            } else {
                std::cerr << "[MAIN] No clusters in queue, breaking" << std::endl;
            }
            
            after_mincut_number_of_clusters = MincutOnly::to_be_mincut_clusters.size();
            
            std::cerr << "[MAIN] Queue state AFTER iteration:" << std::endl;
            std::cerr << "  to_be_mincut: " << after_mincut_number_of_clusters << std::endl;
            std::cerr << "  done_being_mincut: " << MincutOnly::done_being_mincut_clusters.size() << std::endl;
            
            this->WriteToLogFile(std::to_string(after_mincut_number_of_clusters) + " [connected components / clusters] to be mincut after a round of mincuts", Log::debug);
            /** SECTION MinCut Each Connected Component END **/

            /** SECTION Check If All Clusters Are Well-Connected START **/
            if(after_mincut_number_of_clusters == 0) {
                std::cerr << "\n========================================" << std::endl;
                std::cerr << "[MAIN] ALL CLUSTERS WELL-CONNECTED!" << std::endl;
                std::cerr << "========================================" << std::endl;
                std::cerr << "Total iterations: " << (iter_count + 1) << std::endl;
                std::cerr << "Final cluster count: " << MincutOnly::done_being_mincut_clusters.size() << std::endl;
                
                this->WriteToLogFile("all clusters are (well) connected", Log::info);
                this->WriteToLogFile("Total number of iterations: " + std::to_string(iter_count + 1), Log::info);
                break;
            }
            /** SECTION Check If All Clusters Are Well-Connected END **/

            iter_count ++;
            
            if (iter_count > 100000) {
                std::cerr << "[MAIN] WARNING: Reached 100000 iterations, breaking" << std::endl;
                break;
            }
        }
    }

    std::cerr << "\n[MAIN] ======== FINAL STATISTICS ========" << std::endl;
    std::cerr << "Total clusters dequeued: " << total_clusters_dequeued << std::endl;
    std::cerr << "Total clusters enqueued: " << total_clusters_enqueued << std::endl;
    std::cerr << "Total clusters completed: " << total_clusters_completed << std::endl;
    std::cerr << "Total singleton partitions discarded: " << total_singleton_partitions_discarded << std::endl;
    std::cerr << "Total node conservation violations: " << total_node_conservation_violations << std::endl;
    
    // Calculate final node count
    uint64_t final_total_nodes = 0;
    std::queue<std::vector<int>> temp_queue = MincutOnly::done_being_mincut_clusters;
    while (!temp_queue.empty()) {
        final_total_nodes += temp_queue.front().size();
        temp_queue.pop();
    }
    
    std::cerr << "\nNode accounting:" << std::endl;
    std::cerr << "  Initial nodes: " << initial_total_nodes << std::endl;
    std::cerr << "  Final nodes: " << final_total_nodes << std::endl;
    std::cerr << "  Discarded singletons: " << total_singleton_partitions_discarded << std::endl;
    std::cerr << "  Expected final: " << (initial_total_nodes - total_singleton_partitions_discarded) << std::endl;
    
    if (final_total_nodes != (initial_total_nodes - total_singleton_partitions_discarded)) {
        std::cerr << "  **CRITICAL: NODE COUNT MISMATCH**" << std::endl;
    } else {
        std::cerr << "  Node accounting verified ✓" << std::endl;
    }

    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    std::cerr << "\n[MAIN] Writing " << MincutOnly::done_being_mincut_clusters.size() 
              << " clusters to: " << this->output_file << std::endl;
    
    this->WriteClusterQueue(MincutOnly::done_being_mincut_clusters, &graph, new_to_originial_id_map);
    igraph_destroy(&graph);
    
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "INSTRUMENTED RUN COMPLETED" << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    return 0;
}