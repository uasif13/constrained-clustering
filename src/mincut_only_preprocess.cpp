#include "mincut_only_preprocess.h"
#include "mmap_subgraph_loader.h"


int MincutOnlyPreProcess::main() {

    std::cerr << "\n========================================" << std::endl;
    std::cerr << "[MAIN] Starting MincutOnlyPreProcess::main() (QUEUE-BASED)" << std::endl;
    std::cerr << "========================================\n" << std::endl;

    this->WriteToLogFile("Parsing connectedness criterion" , Log::info);
    std::cerr << "[MAIN] Parsing connectedness criterion: " << this->connectedness_criterion << std::endl;
    
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
        std::cerr << "[MAIN] Detected logarithmic connectedness criterion" << std::endl;
        // is Clog_x(n)
        current_connectedness_criterion = ConnectednessCriterion::Logarithimic;
        if (log_position == 0) {
            connectedness_criterion_c = 1;
        } else {
            connectedness_criterion_c = std::stod(this->connectedness_criterion.substr(0, log_position));
        }
        size_t open_parenthesis_position = this->connectedness_criterion.find("(", log_position + 4);
        connectedness_criterion_x = std::stod(this->connectedness_criterion.substr(log_position + 4, open_parenthesis_position));
        std::cerr << "[MAIN] C=" << connectedness_criterion_c << ", x=" << connectedness_criterion_x << std::endl;
    } else if (n_caret_position != std::string::npos) {
        std::cerr << "[MAIN] Detected exponential connectedness criterion" << std::endl;
        // is cN^x
        current_connectedness_criterion = ConnectednessCriterion::Exponential;
        if (n_caret_position == 0) {
            connectedness_criterion_c = 1;
        } else {
            connectedness_criterion_c = std::stod(this->connectedness_criterion.substr(0, n_caret_position));
        }
        connectedness_criterion_x = std::stod(this->connectedness_criterion.substr(n_caret_position + 2));
        std::cerr << "[MAIN] C=" << connectedness_criterion_c << ", x=" << connectedness_criterion_x << std::endl;
    } else if (this->connectedness_criterion != "0") {
        // wasn't log or exponent so if it isn't 0 then it's an error
        // exit
        std::cerr << "[MAIN ERROR] Could not parse connectedness_criterion: " << this->connectedness_criterion << std::endl;
        this->WriteToLogFile("Could not parse connectedness_criterion" , Log::error);
        this->WriteToLogFile("Accepted forms are Clog_x(n), Cn^x, or 0 where C and x are integers" , Log::error);
        return 1;
    }
    
    if (current_connectedness_criterion == ConnectednessCriterion::Simple) {
        this->WriteToLogFile("Running with CC mode (mincut of each cluster has to be greater than 0)" , Log::info);
        std::cerr << "[MAIN] Mode: Simple Connected Components (CC)" << std::endl;
    } else if (current_connectedness_criterion == ConnectednessCriterion::Logarithimic) {
        this->WriteToLogFile("Running with WCC mode (mincut of each cluster has to be greater than " + std::to_string(connectedness_criterion_c) + " times log base " + std::to_string(connectedness_criterion_x) + " of n" , Log::info);
        pre_computed_log = connectedness_criterion_c / std::log(connectedness_criterion_x);
        std::cerr << "[MAIN] Mode: Well-Connected (logarithmic), pre_computed_log=" << pre_computed_log << std::endl;
    } else if (current_connectedness_criterion == ConnectednessCriterion::Exponential) {
        this->WriteToLogFile("Running with WCC mode (mincut of each cluster has to be greater than " + std::to_string(connectedness_criterion_c) + " times n to the power of " + std::to_string(connectedness_criterion_x), Log::info);
        std::cerr << "[MAIN] Mode: Well-Connected (exponential)" << std::endl;
    } else {
        // should not possible to reach
        std::cerr << "[MAIN ERROR] Unexpected connectedness criterion state" << std::endl;
        exit(1);
    }
    
    std::cerr << "\n[MAIN] ======== Loading Graph ========" << std::endl;
    this->WriteToLogFile("Read subgraph edgelist file", Log::info);
    std::cerr << "[MAIN] Reading from: " << this->edgelist << std::endl;
    
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> subgraph_edges_vector = MMapSubgraphLoader::LoadEdgelistMMap(this->edgelist, &original_to_new_id_unordered_map, false);

    std::cerr << "[MAIN] Finished loading graph" << std::endl;
    std::cerr << "[MAIN] Loaded " << subgraph_edges_vector.size() << " connected components" << std::endl;
    std::cerr << "[MAIN] Total unique nodes: " << original_to_new_id_unordered_map.size() << std::endl;

    this->WriteToLogFile("Finished reading subgraph edgelist file", Log::info);
    
    // Validate that we actually got data
    if (subgraph_edges_vector.empty()) {
        std::cerr << "[MAIN ERROR] No components loaded! Vector is empty!" << std::endl;
        this->WriteToLogFile("ERROR: No components loaded from file", Log::error);
        return 1;
    }
    if (subgraph_edges_vector.size() == 1 && subgraph_edges_vector[0].size() == 1 && subgraph_edges_vector[0][0] == -1) {
        std::cerr << "[MAIN ERROR] File loading failed (error marker detected)" << std::endl;
        this->WriteToLogFile("ERROR: File loading returned error marker", Log::error);
        return 1;
    }
    
    int before_mincut_number_of_clusters = 0;
    int iter_count = 0;
    /** SECTION Get Connected Components END **/

    std::cerr << "\n[MAIN] ======== Populating Initial Work Queue ========" << std::endl;
    std::cerr << "[MAIN] Adding " << subgraph_edges_vector.size() << " components to work queue" << std::endl;
    
    // store the results into the queue that each thread pulls from
    for(size_t i = 0; i < subgraph_edges_vector.size(); i ++) {
        if (i < 10 || i % 100 == 0 || i >= subgraph_edges_vector.size() - 10) {
            std::cerr << "[MAIN] Adding component " << i << " (size: " << subgraph_edges_vector[i].size() << " values)" << std::endl;
        }
        MincutOnlyPreProcess::to_be_mincut_clusters.push(subgraph_edges_vector[i]);
    }
    
    std::cerr << "[MAIN] Initial queue populated with " << MincutOnlyPreProcess::to_be_mincut_clusters.size() << " clusters" << std::endl;
    
    std::cerr << "\n[MAIN] ======== Starting Iterative Processing ========" << std::endl;
    
    while (true) {
        std::cerr << "\n[MAIN] ========================================" << std::endl;
        std::cerr << "[MAIN] Iteration #" << iter_count << std::endl;
        std::cerr << "[MAIN] ========================================" << std::endl;
        
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
        if(iter_count % 10 == 0) {
            this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
            this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::info);
        }

        /** SECTION MinCut Each Connected Component START **/
        before_mincut_number_of_clusters = MincutOnlyPreProcess::to_be_mincut_clusters.size();
        
        std::cerr << "[MAIN] Clusters in work queue: " << before_mincut_number_of_clusters << std::endl;
        std::cerr << "[MAIN] Completed clusters: " << MincutOnlyPreProcess::done_being_mincut_clusters.size() << std::endl;
        
        this->WriteToLogFile(std::to_string(before_mincut_number_of_clusters) + " [connected components / clusters] to be mincut", Log::debug);
        
        /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
        if(before_mincut_number_of_clusters > 1) {
            std::cerr << "[MAIN] Spawning " << this->num_processors << " worker threads" << std::endl;
            
            /* start the threads */
            for(int i = 0; i < this->num_processors; i ++) {
                MincutOnlyPreProcess::to_be_mincut_clusters.push({-1});
                std::cerr << "[MAIN] Added termination signal #" << (i+1) << " to queue" << std::endl;
            }
            
            std::cerr << "[MAIN] Queue now has " << MincutOnlyPreProcess::to_be_mincut_clusters.size() 
                     << " items (including " << this->num_processors << " termination signals)" << std::endl;
            
            std::vector<std::thread> thread_vector;
            for(int i = 0; i < this->num_processors; i ++) {
                std::cerr << "[MAIN] Spawning worker thread #" << i << std::endl;
                thread_vector.push_back(std::thread(MincutOnlyPreProcess::MinCutWorker, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log));
                // Notify the thread
                to_be_mincut_condition_variable.notify_one();
            }
            
            std::cerr << "[MAIN] All " << this->num_processors << " worker threads spawned" << std::endl;
            std::cerr << "[MAIN] Waiting for workers to complete..." << std::endl;
            
            /* get the result back from threads */
            /* the results from each thread gets stored in to_be_clustered_clusters */
            for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
                std::cerr << "[MAIN] Joining worker thread #" << thread_index << "..." << std::endl;
                thread_vector[thread_index].join();
                std::cerr << "[MAIN] Worker thread #" << thread_index << " joined" << std::endl;
            }
            
            std::cerr << "[MAIN] All worker threads completed" << std::endl;
        } else if (before_mincut_number_of_clusters == 1) {
            std::cerr << "[MAIN] Only 1 cluster remaining, processing with single worker" << std::endl;
            MincutOnlyPreProcess::to_be_mincut_clusters.push({-1});
            to_be_mincut_condition_variable.notify_one();
            MincutOnlyPreProcess::MinCutWorker(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log);
            std::cerr << "[MAIN] Single worker completed" << std::endl;
        } else {
            std::cerr << "[MAIN] No clusters in queue, skipping worker spawn" << std::endl;
        }
        
        int after_mincut_number_of_clusters = MincutOnlyPreProcess::to_be_mincut_clusters.size();
        
        std::cerr << "[MAIN] After iteration: work queue has " << after_mincut_number_of_clusters << " clusters" << std::endl;
        this->WriteToLogFile(std::to_string(after_mincut_number_of_clusters) + " [connected components / clusters] to be mincut after a round of mincuts", Log::debug);
        /** SECTION MinCut Each Connected Component END **/

        /** SECTION Check If All Clusters Are Well-Connected START **/
        if(after_mincut_number_of_clusters == 0) {
            std::cerr << "\n[MAIN] ========================================" << std::endl;
            std::cerr << "[MAIN] ALL CLUSTERS ARE WELL-CONNECTED!" << std::endl;
            std::cerr << "[MAIN] ========================================" << std::endl;
            std::cerr << "[MAIN] Total iterations: " << (iter_count + 1) << std::endl;
            std::cerr << "[MAIN] Final cluster count: " << MincutOnlyPreProcess::done_being_mincut_clusters.size() << std::endl;
            
            this->WriteToLogFile("all clusters are (well) connected", Log::info);
            this->WriteToLogFile("Total number of iterations: " + std::to_string(iter_count + 1), Log::info);
            break;
        }
        /** SECTION Check If All Clusters Are Well-Connected END **/

        iter_count ++;
        
        // if (iter_count > 10000) {
        //     std::cerr << "[MAIN WARNING] Reached 10000 iterations, possible infinite loop!" << std::endl;
        //     std::cerr << "[MAIN] Breaking out of loop..." << std::endl;
        //     break;
        // }
    }

    std::cerr << "\n[MAIN] ======== Writing Output ========" << std::endl;
    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    std::cerr << "[MAIN] Output file: " << this->output_file << std::endl;
    std::cerr << "[MAIN] Writing " << MincutOnlyPreProcess::done_being_mincut_clusters.size() << " clusters" << std::endl;
    
    this->WriteClusterQueue(MincutOnlyPreProcess::done_being_mincut_clusters);
    
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "[MAIN] MincutOnlyPreProcess completed successfully!" << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    return 0;
}