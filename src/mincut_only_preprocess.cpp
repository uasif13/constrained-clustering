#include "mincut_only_preprocess.h"
#include "mmap_subgraph_loader.h"

int MincutOnlyPreProcess::main() {

    this->WriteToLogFile("========================================", Log::info);
    this->WriteToLogFile("Starting MincutOnlyPreProcess::main()", Log::info);
    this->WriteToLogFile("========================================", Log::info);

    this->WriteToLogFile("Parsing connectedness criterion: " + this->connectedness_criterion, Log::info);
    
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
        this->WriteToLogFile("Detected logarithmic connectedness criterion", Log::info);
        // is Clog_x(n)
        current_connectedness_criterion = ConnectednessCriterion::Logarithimic;
        if (log_position == 0) {
            connectedness_criterion_c = 1;
        } else {
            connectedness_criterion_c = std::stod(this->connectedness_criterion.substr(0, log_position));
        }
        size_t open_parenthesis_position = this->connectedness_criterion.find("(", log_position + 4);
        connectedness_criterion_x = std::stod(this->connectedness_criterion.substr(log_position + 4, open_parenthesis_position));
        this->WriteToLogFile("C=" + std::to_string(connectedness_criterion_c) + ", x=" + std::to_string(connectedness_criterion_x), Log::info);
    } else if (n_caret_position != std::string::npos) {
        this->WriteToLogFile("Detected exponential connectedness criterion", Log::info);
        // is cN^x
        current_connectedness_criterion = ConnectednessCriterion::Exponential;
        if (n_caret_position == 0) {
            connectedness_criterion_c = 1;
        } else {
            connectedness_criterion_c = std::stod(this->connectedness_criterion.substr(0, n_caret_position));
        }
        connectedness_criterion_x = std::stod(this->connectedness_criterion.substr(n_caret_position + 2));
        this->WriteToLogFile("C=" + std::to_string(connectedness_criterion_c) + ", x=" + std::to_string(connectedness_criterion_x), Log::info);
    } else if (this->connectedness_criterion != "0") {
        // wasn't log or exponent so if it isn't 0 then it's an error
        // exit
        this->WriteToLogFile("Could not parse connectedness_criterion: " + this->connectedness_criterion, Log::error);
        this->WriteToLogFile("Could not parse connectedness_criterion" , Log::error);
        this->WriteToLogFile("Accepted forms are Clog_x(n), Cn^x, or 0 where C and x are integers" , Log::error);
        return 1;
    }
    
    if (current_connectedness_criterion == ConnectednessCriterion::Simple) {
        this->WriteToLogFile("Running with CC mode (mincut of each cluster has to be greater than 0)" , Log::info);
    } else if (current_connectedness_criterion == ConnectednessCriterion::Logarithimic) {
        this->WriteToLogFile("Running with WCC mode (mincut of each cluster has to be greater than " + std::to_string(connectedness_criterion_c) + " times log base " + std::to_string(connectedness_criterion_x) + " of n" , Log::info);
        pre_computed_log = connectedness_criterion_c / std::log(connectedness_criterion_x);
        this->WriteToLogFile("Mode: Well-Connected (logarithmic), pre_computed_log=" + std::to_string(pre_computed_log), Log::info);
    } else if (current_connectedness_criterion == ConnectednessCriterion::Exponential) {
        this->WriteToLogFile("Running with WCC mode (mincut of each cluster has to be greater than " + std::to_string(connectedness_criterion_c) + " times n to the power of " + std::to_string(connectedness_criterion_x), Log::info);
    } else {
        // should not possible to reach
        this->WriteToLogFile("Unexpected connectedness criterion state", Log::error);
        exit(1);
    }
    
    this->WriteToLogFile("======== Loading Graph ========", Log::info);
    this->WriteToLogFile("Read subgraph edgelist file", Log::info);
    this->WriteToLogFile("Reading from: " + this->edgelist, Log::info);
    
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> subgraph_edges_vector = MMapSubgraphLoader::LoadEdgelistMMap(this->edgelist, &original_to_new_id_unordered_map, false);

    this->WriteToLogFile("Finished loading graph", Log::info);
    this->WriteToLogFile("Loaded " + std::to_string(subgraph_edges_vector.size()) + " connected components", Log::info);
    this->WriteToLogFile("Total unique nodes: " + std::to_string(original_to_new_id_unordered_map.size()), Log::info);
    
    // Validate that we actually got data
    if (subgraph_edges_vector.empty()) {
        this->WriteToLogFile("ERROR: No components loaded! Vector is empty!", Log::error);
        return 1;
    }
    if (subgraph_edges_vector.size() == 1 && subgraph_edges_vector[0].size() == 1 && subgraph_edges_vector[0][0] == -1) {
        this->WriteToLogFile("ERROR: File loading failed (error marker detected)", Log::error);
        return 1;
    }

    this->WriteToLogFile("Finished reading subgraph edgelist file", Log::info);
    long before_mincut_number_of_clusters = 0;
    int iter_count = 0;

    this->WriteToLogFile("======== Starting MinCut Processing ========", Log::info);
    this->WriteToLogFile("Start Mincut", Log::info);

    int previous_done_being_clustered_size = 0;
    int previous_cluster_id = 0;

    /** SECTION MinCut Each Connected Component START **/
    this->WriteToLogFile(std::to_string(subgraph_edges_vector.size()) + " [connected components / clusters] to be mincut", Log::debug);
    this->WriteToLogFile("Processing " + std::to_string(subgraph_edges_vector.size()) + " components", Log::info);
    
    before_mincut_number_of_clusters = subgraph_edges_vector.size();
    
    long current_components_vector_index = 0;
    
    this->WriteToLogFile("Starting threaded processing with " + std::to_string(this->num_processors) + " threads", Log::info);
    
    while (current_components_vector_index < before_mincut_number_of_clusters) {
        this->WriteToLogFile("Thread batch starting at index " + std::to_string(current_components_vector_index) + 
                 " / " + std::to_string(before_mincut_number_of_clusters), Log::debug);
        
        std::vector<std::thread> thread_vector;
        
        for(int i = 0; i < this->num_processors; i++) {
            if (current_components_vector_index < before_mincut_number_of_clusters) {
                this->WriteToLogFile("Spawning thread " + std::to_string(i) + " for component " + std::to_string(current_components_vector_index), Log::debug);
                this->WriteToLogFile("Component size: " + std::to_string(subgraph_edges_vector[current_components_vector_index].size()) + 
                         " edges", Log::debug);
                
                // Safety check
                if (subgraph_edges_vector[current_components_vector_index].empty()) {
                    this->WriteToLogFile("Component " + std::to_string(current_components_vector_index) + 
                             " is empty! Skipping...", Log::debug);
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
                this->WriteToLogFile("No more components to process, stopping thread spawn", Log::debug);
                break;
            }
        }
        
        this->WriteToLogFile("Waiting for " + std::to_string(thread_vector.size()) + " threads to complete...", Log::debug);
        
        for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index++) {
            thread_vector[thread_index].join();
        }
        
        this->WriteToLogFile("Thread batch completed", Log::debug);
    }
    
    this->WriteToLogFile("======== All Processing Complete ========", Log::info);
    /** SECTION MinCut Each Connected Component END **/

    /** SECTION Check If All Clusters Are Well-Connected START **/
    this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()) + " [connected components / clusters] mincut after a round of mincuts", Log::debug);
    this->WriteToLogFile("Final cluster count: " + std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()), Log::info);
        
    this->WriteToLogFile("======== Writing Output ========", Log::info);
    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteToLogFile("Output file: " + this->output_file, Log::info);
    
    this->WriteClusterQueue(MincutOnlyPreProcess::done_being_mincut_clusters);
    
    this->WriteToLogFile("========================================", Log::info);
    this->WriteToLogFile("MincutOnlyPreProcess completed successfully!", Log::info);
    this->WriteToLogFile("========================================", Log::info);
    
    return 0;
}