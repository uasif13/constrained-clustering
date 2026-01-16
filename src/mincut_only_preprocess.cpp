#include "mincut_only_preprocess.h"
#include "mmap_subgraph_loader.h"

int MincutOnlyPreProcess::main() {

    std::cerr << "\n========================================" << std::endl;
    std::cerr << "[MAIN] Starting MincutOnlyPreProcess::main()" << std::endl;
    std::cerr << "========================================\n" << std::endl;

    this->WriteToLogFile("Parsing connectedness criterion" , Log::info);
    std::cerr << "[MAIN] Parsing connectedness criterion: " << this->connectedness_criterion << std::endl;
    
/* F(n) = C log_x(n), where C and x are parameters specified by the user (our default is C=1 and x=10) */
/* G(n) = C n^x, where C and x are parameters specified by the user (here, presumably 0<x<2). Note that x=1 makes it linear. */
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

    this->WriteToLogFile("Finished reading subgraph edgelist file", Log::info);
    long before_mincut_number_of_clusters = 0;
    int iter_count = 0;

    std::cerr << "\n[MAIN] ======== Starting MinCut Processing ========" << std::endl;
    this->WriteToLogFile("Start Mincut", Log::info);

    int previous_done_being_clustered_size = 0;
    int previous_cluster_id = 0;

    this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
    if(iter_count % 10000 == 0) {
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
        this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::info);
    }

    /** SECTION MinCut Each Connected Component START **/
    this->WriteToLogFile(std::to_string(subgraph_edges_vector.size()) + " [connected components / clusters] to be mincut", Log::debug);
    std::cerr << "[MAIN] Processing " << subgraph_edges_vector.size() << " components" << std::endl;
    
    before_mincut_number_of_clusters = subgraph_edges_vector.size();
    
    long current_components_vector_index = 0;
    
    std::cerr << "[MAIN] Starting threaded processing with " << this->num_processors << " threads" << std::endl;
    
    while (current_components_vector_index < before_mincut_number_of_clusters) {
        std::cerr << "\n[MAIN] Thread batch starting at index " << current_components_vector_index 
                 << " / " << before_mincut_number_of_clusters << std::endl;
        
        std::vector<std::thread> thread_vector;
        
        for(int i = 0; i < this->num_processors; i++) {
            if (current_components_vector_index < before_mincut_number_of_clusters) {
                std::cerr << "[MAIN] Spawning thread " << i << " for component " << current_components_vector_index << std::endl;
                std::cerr << "[MAIN]   Component size: " << subgraph_edges_vector[current_components_vector_index].size() 
                         << " edges" << std::endl;
                
                // Safety check
                if (subgraph_edges_vector[current_components_vector_index].empty()) {
                    std::cerr << "[MAIN WARNING] Component " << current_components_vector_index 
                             << " is empty! Skipping..." << std::endl;
                    current_components_vector_index++;
                    continue;
                }
                
                thread_vector.push_back(std::thread(
                    MincutOnlyPreProcess::ComputeMinCutRecursive,
                    subgraph_edges_vector[current_components_vector_index], 
                    current_connectedness_criterion, 
                    connectedness_criterion_c, 
                    connectedness_criterion_x, 
                    pre_computed_log
                ));
                
                current_components_vector_index++;
            } else {
                std::cerr << "[MAIN] No more components to process, stopping thread spawn" << std::endl;
                break;
            }
        }
        
        std::cerr << "[MAIN] Waiting for " << thread_vector.size() << " threads to complete..." << std::endl;
        
        for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index++) {
            std::cerr << "[MAIN] Joining thread " << thread_index << std::endl;
            thread_vector[thread_index].join();
        }
        
        std::cerr << "[MAIN] Thread batch completed" << std::endl;
    }
    
    std::cerr << "\n[MAIN] ======== All Processing Complete ========" << std::endl;
    /** SECTION MinCut Each Connected Component END **/

    /** SECTION Check If All Clusters Are Well-Connected START **/
    this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()) + " [connected components / clusters] mincut after a round of mincuts", Log::debug);
    std::cerr << "[MAIN] Final cluster count: " << MincutOnlyPreProcess::done_being_mincut_clusters.size() << std::endl;
        
    std::cerr << "\n[MAIN] ======== Writing Output ========" << std::endl;
    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    std::cerr << "[MAIN] Output file: " << this->output_file << std::endl;
    
    this->WriteClusterQueue(MincutOnlyPreProcess::done_being_mincut_clusters);
    
    std::cerr << "\n========================================" << std::endl;
    std::cerr << "[MAIN] MincutOnlyPreProcess completed successfully!" << std::endl;
    std::cerr << "========================================\n" << std::endl;
    
    return 0;
}