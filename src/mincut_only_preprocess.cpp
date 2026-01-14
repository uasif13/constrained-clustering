#include "mincut_only_preprocess.h"
#include "mmap_subgraph_loader.h"
#include <iomanip>


int MincutOnlyPreProcess::main() {

    this->WriteToLogFile("Parsing connectedness criterion" , Log::info);
/* F(n) = C log_x(n), where C and x are parameters specified by the user (our default is C=1 and x=10) */
/* G(n) = C n^x, where C and x are parameterss specified by the user (here, presumably 0<x<2). Note that x=1 makes it linear. */
        /* static inline bool IsWellConnected(std::string connectedness_criterion, int in_partition_size, int out_partition_size, int edge_cut_size) { */
    size_t log_position = this->connectedness_criterion.find("log_");
    size_t n_caret_position = this->connectedness_criterion.find("n^");
    double connectedness_criterion_c = 1;
    double connectedness_criterion_x = 1;
    double pre_computed_log = -1;
    ConnectednessCriterion current_connectedness_criterion = ConnectednessCriterion::Simple;
    if (log_position != std::string::npos) {
        // is Clog_x(n)
        current_connectedness_criterion = ConnectednessCriterion::Logarithimic;
        if (log_position == 0) {
            connectedness_criterion_c = 1;
        } else {
            connectedness_criterion_c = std::stod(this->connectedness_criterion.substr(0, log_position));
        }
        size_t open_parenthesis_position = this->connectedness_criterion.find("(", log_position + 4);
        connectedness_criterion_x = std::stod(this->connectedness_criterion.substr(log_position + 4, open_parenthesis_position));
    } else if (n_caret_position != std::string::npos) {
        // is cN^x
        current_connectedness_criterion = ConnectednessCriterion::Exponential;
        if (n_caret_position == 0) {
            connectedness_criterion_c = 1;
        } else {
            connectedness_criterion_c = std::stod(this->connectedness_criterion.substr(0, n_caret_position));
        }
        connectedness_criterion_x = std::stod(this->connectedness_criterion.substr(n_caret_position + 2));
    } else if (this->connectedness_criterion != "0") {
        // wasn't log or exponent so if it isn't 0 then it's an error
        // exit
        this->WriteToLogFile("Colud not parse connectedness_criterion" , Log::error);
        this->WriteToLogFile("Accepted forms are Clog_x(n), Cn^x, or 0 where C and x are integers" , Log::error);
        return 1;
    }
    if (current_connectedness_criterion == ConnectednessCriterion::Simple) {
        this->WriteToLogFile("Running with CC mode (mincut of each cluster has to be greater than 0)" , Log::info);
    } else if (current_connectedness_criterion == ConnectednessCriterion::Logarithimic) {
        this->WriteToLogFile("Running with WCC mode (mincut of each cluster has to be greater than " + std::to_string(connectedness_criterion_c) + " times log base " + std::to_string(connectedness_criterion_x) + "of n" , Log::info);
        pre_computed_log = connectedness_criterion_c / std::log(connectedness_criterion_x);
    } else if (current_connectedness_criterion == ConnectednessCriterion::Exponential) {
        this->WriteToLogFile("Running with WCC mode (mincut of each cluster has to be greater than " + std::to_string(connectedness_criterion_c) + " times n to the power of " + std::to_string(connectedness_criterion_x), Log::info);
    } else {
        // should not possible to reach
        exit(1);
    }
    this -> WriteToLogFile("Read subgraph edgelist file", Log::info);
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> subgraph_edges_vector = MMapSubgraphLoader::LoadEdgelistMMap(this->edgelist, &original_to_new_id_unordered_map,false);

    this -> WriteToLogFile("Finished reading subgraph edgelist file", Log::info);
    int before_mincut_number_of_clusters = 0;
    int iter_count = 0;
    /** SECTION Get Connected Components END **/

        // store the results into the queue that each thread pulls from
    for(size_t i = 0; i < subgraph_edges_vector.size(); i ++) {
        MincutOnlyPreProcess::to_be_mincut_clusters.push(subgraph_edges_vector[i]);
    }
    
    int initial_clusters = subgraph_edges_vector.size();
    
    while (true) {
        /* std::cerr << "iter num: " << std::to_string(iter_count) << std::endl; */
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
        if(iter_count % 10000 == 0) {
            this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
            this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::info);
        }

        /** SECTION MinCut Each Connected Component START **/
        this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::debug);
        before_mincut_number_of_clusters = MincutOnlyPreProcess::to_be_mincut_clusters.size();
        /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
        /* std::cerr << "num clusters to be processed: " << std::to_string(before_mincut_number_of_clusters) << std::endl; */
        if(before_mincut_number_of_clusters > 1) {
            /* start the threads */
            for(int i = 0; i < this->num_processors; i ++) {
                MincutOnlyPreProcess::to_be_mincut_clusters.push({-1});
            }
            std::vector<std::thread> thread_vector;
            for(int i = 0; i < this->num_processors; i ++) {
                thread_vector.push_back(std::thread(MincutOnlyPreProcess::MinCutWorker, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log));
            }
            /* get the result back from threads */
            /* the results from each thread gets stored in to_be_clustered_clusters */
            for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
                thread_vector[thread_index].join();
            }
        } else {
            MincutOnlyPreProcess::to_be_mincut_clusters.push({-1});
            MincutOnlyPreProcess::MinCutWorker(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log);
        }
        this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut after a round of mincuts", Log::debug);
        /** SECTION MinCut Each Connected Component END **/

        /** SECTION Check If All Clusters Are Well-Connected START **/
        int after_mincut_number_of_clusters = MincutOnlyPreProcess::to_be_mincut_clusters.size();
        if(after_mincut_number_of_clusters == 0) {
            this->WriteToLogFile("all clusters are (well) connected", Log::info);
            this->WriteToLogFile("Total number of iterations: " + std::to_string(iter_count + 1), Log::info);
            break;
        }
        /** SECTION Check If All Clusters Are Well-Connected END **/

        iter_count ++;
    }


    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteClusterQueue(MincutOnlyPreProcess::done_being_mincut_clusters);

    // ========== DIAGNOSTIC SUMMARY ==========
    std::cerr << "\n";
    std::cerr << "=============================================" << std::endl;
    std::cerr << "   DIAGNOSTIC SUMMARY (QUEUE)" << std::endl;
    std::cerr << "=============================================" << std::endl;
    std::cerr << "Initial clusters:         " << initial_clusters << std::endl;
    std::cerr << "Iterations:               " << iter_count << std::endl;
    std::cerr << "Connectedness criterion:  " << this->connectedness_criterion << std::endl;
    if(current_connectedness_criterion == ConnectednessCriterion::Logarithimic) {
        std::cerr << "  (C=" << connectedness_criterion_c << ", x=" << connectedness_criterion_x << ")" << std::endl;
    } else if(current_connectedness_criterion == ConnectednessCriterion::Exponential) {
        std::cerr << "  (C=" << connectedness_criterion_c << ", x=" << connectedness_criterion_x << ")" << std::endl;
    }
    std::cerr << "Clusters dequeued:        " << MincutOnlyPreProcess::g_clusters_dequeued.load() << std::endl;
    std::cerr << "Mincuts computed:         " << MincutOnlyPreProcess::g_clusters_mincut_computed.load() << std::endl;
    std::cerr << "Output (well-connected):  " << MincutOnlyPreProcess::g_clusters_output_done.load() << std::endl;
    std::cerr << "Output (not-well-conn):   " << MincutOnlyPreProcess::g_clusters_output_notdone.load() << std::endl;
    std::cerr << "Both-singleton fixed:     " << MincutOnlyPreProcess::g_both_singleton_fixed.load() << std::endl;
    std::cerr << "Empty vectors found:      " << MincutOnlyPreProcess::g_empty_vectors_found.load() << std::endl;
    std::cerr << "---------------------------------------------" << std::endl;
    std::cerr << "GCCP Analysis:" << std::endl;
    std::cerr << "  GCCP calls:               " << MincutOnlyPreProcess::g_gccp_calls.load() << std::endl;
    std::cerr << "  igraph found components:  " << MincutOnlyPreProcess::g_igraph_components_found.load() << std::endl;
    std::cerr << "  GCCP returned components: " << MincutOnlyPreProcess::g_gccp_components_returned.load() << std::endl;
    std::cerr << "  Zero-edge components:     " << MincutOnlyPreProcess::g_gccp_zero_edge_components.load() << std::endl;
    
    // Calculate averages
    long gccp_calls = MincutOnlyPreProcess::g_gccp_calls.load();
    if(gccp_calls > 0) {
        double avg_igraph = (double)MincutOnlyPreProcess::g_igraph_components_found.load() / gccp_calls;
        double avg_returned = (double)MincutOnlyPreProcess::g_gccp_components_returned.load() / gccp_calls;
        std::cerr << "  Avg components per call:" << std::endl;
        std::cerr << "    igraph finds:  " << std::fixed << std::setprecision(2) << avg_igraph << std::endl;
        std::cerr << "    GCCP returns:  " << std::fixed << std::setprecision(2) << avg_returned << std::endl;
    }
    std::cerr << "=============================================" << std::endl;
    
    int total_output = MincutOnlyPreProcess::g_clusters_output_done.load() + 
                       MincutOnlyPreProcess::g_clusters_output_notdone.load();
    int mincut_loss = MincutOnlyPreProcess::g_clusters_mincut_computed.load() - total_output;
    int dequeue_loss = MincutOnlyPreProcess::g_clusters_dequeued.load() - 
                       MincutOnlyPreProcess::g_clusters_mincut_computed.load();
    
    std::cerr << "\nDiscrepancy Analysis:" << std::endl;
    if(dequeue_loss > 0) {
        std::cerr << "  ⚠ WARNING: Lost " << dequeue_loss 
                  << " clusters between dequeue and mincut!" << std::endl;
        std::cerr << "     (These are likely sentinel values - expected)" << std::endl;
    }
    if(mincut_loss > 0) {
        std::cerr << "  ⚠ WARNING: Lost " << mincut_loss 
                  << " clusters between mincut and output!" << std::endl;
    }
    if(MincutOnlyPreProcess::g_empty_vectors_found.load() > 0) {
        std::cerr << "  ⚠ CRITICAL: Found " << MincutOnlyPreProcess::g_empty_vectors_found.load() 
                  << " empty vectors from GetConnectedComponents!" << std::endl;
        std::cerr << "     This is likely causing significant node loss." << std::endl;
    }
    if(MincutOnlyPreProcess::g_both_singleton_fixed.load() > 0) {
        std::cerr << "  ✓ INFO: Fixed " << MincutOnlyPreProcess::g_both_singleton_fixed.load() 
                  << " clusters where both partitions were singletons" << std::endl;
    }
    
    if(MincutOnlyPreProcess::g_gccp_components_returned.load() != 
       MincutOnlyPreProcess::g_igraph_components_found.load()) {
        long component_diff = MincutOnlyPreProcess::g_gccp_components_returned.load() - 
                              MincutOnlyPreProcess::g_igraph_components_found.load();
        std::cerr << "  ⚠ GCCP MISMATCH: GCCP returned " << component_diff 
                  << " different components than igraph found!" << std::endl;
        std::cerr << "     This could explain differences in splitting behavior." << std::endl;
    }
    
    if(dequeue_loss == 0 && mincut_loss == 0 && 
       MincutOnlyPreProcess::g_empty_vectors_found.load() == 0) {
        std::cerr << "  ✓ All counters balanced - no obvious loss detected" << std::endl;
    }
    std::cerr << "=============================================" << std::endl;
    
    return 0;
}