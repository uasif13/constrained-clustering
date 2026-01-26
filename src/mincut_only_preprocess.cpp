#include "mincut_only_preprocess.h"
#include "mmap_subgraph_loader.h"


int MincutOnlyPreProcess::main() {

    this->WriteToLogFile("========================================", Log::info);
    this->WriteToLogFile("Starting MincutOnlyPreProcess::main() (QUEUE-BASED)", Log::info);
    this->WriteToLogFile("========================================", Log::info);

    this->WriteToLogFile("Parsing connectedness criterion: " + this->connectedness_criterion, Log::info);
    
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

    this->WriteToLogFile("Finished reading subgraph edgelist file", Log::info);
    
    // Validate that we actually got data
    if (subgraph_edges_vector.empty()) {
        this->WriteToLogFile("ERROR: No components loaded! Vector is empty!", Log::error);
        return 1;
    }
    if (subgraph_edges_vector.size() == 1 && subgraph_edges_vector[0].size() == 1 && subgraph_edges_vector[0][0] == -1) {
        this->WriteToLogFile("ERROR: File loading failed (error marker detected)", Log::error);
        return 1;
    }
    
    int before_mincut_number_of_clusters = 0;
    int iter_count = 0;
    /** SECTION Get Connected Components END **/

    this->WriteToLogFile("======== Populating Initial Work Queue ========", Log::info);
    this->WriteToLogFile("Adding " + std::to_string(subgraph_edges_vector.size()) + " components to work queue", Log::info);
    
    // store the results into the queue that each thread pulls from
    for(size_t i = 0; i < subgraph_edges_vector.size(); i ++) {
        MincutOnlyPreProcess::to_be_mincut_clusters.push(subgraph_edges_vector[i]);
    }
    
    this->WriteToLogFile("Initial queue populated with " + std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + " clusters", Log::info);
    
    this->WriteToLogFile("======== Starting Iterative Processing ========", Log::info);
    
    while (true) {
        this->WriteToLogFile("========================================", Log::debug);
        this->WriteToLogFile("Iteration #" + std::to_string(iter_count), Log::debug);
        this->WriteToLogFile("========================================", Log::debug);
        
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
        if(iter_count % 10 == 0) {
            this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
            this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::info);
        }

        /** SECTION MinCut Each Connected Component START **/
        before_mincut_number_of_clusters = MincutOnlyPreProcess::to_be_mincut_clusters.size();
        
        this->WriteToLogFile("Clusters in work queue: " + std::to_string(before_mincut_number_of_clusters), Log::debug);
        this->WriteToLogFile("Completed clusters: " + std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()), Log::debug);
        
        this->WriteToLogFile(std::to_string(before_mincut_number_of_clusters) + " [connected components / clusters] to be mincut", Log::debug);
        
        /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
        if(before_mincut_number_of_clusters > 1) {
            this->WriteToLogFile("Spawning " + std::to_string(this->num_processors) + " worker threads", Log::debug);
            
            /* start the threads */
            for(int i = 0; i < this->num_processors; i ++) {
                MincutOnlyPreProcess::to_be_mincut_clusters.push({-1});
            }
            
            this->WriteToLogFile("Queue now has " + std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + 
                     " items (including " + std::to_string(this->num_processors) + " termination signals)", Log::debug);
            
            std::vector<std::thread> thread_vector;
            for(int i = 0; i < this->num_processors; i ++) {
                thread_vector.push_back(std::thread(MincutOnlyPreProcess::MinCutWorker, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log));
                // Notify the thread
                to_be_mincut_condition_variable.notify_one();
            }
            
            this->WriteToLogFile("All " + std::to_string(this->num_processors) + " worker threads spawned", Log::debug);
            this->WriteToLogFile("Waiting for workers to complete...", Log::debug);
            
            /* get the result back from threads */
            /* the results from each thread gets stored in to_be_clustered_clusters */
            for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
                thread_vector[thread_index].join();
            }
            
            this->WriteToLogFile("All worker threads completed", Log::debug);
        } else if (before_mincut_number_of_clusters == 1) {
            this->WriteToLogFile("Only 1 cluster remaining, processing with single worker", Log::debug);
            MincutOnlyPreProcess::to_be_mincut_clusters.push({-1});
            to_be_mincut_condition_variable.notify_one();
            MincutOnlyPreProcess::MinCutWorker(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log);
            this->WriteToLogFile("Single worker completed", Log::debug);
        } else {
            this->WriteToLogFile("No clusters in queue, skipping worker spawn", Log::debug);
        }
        
        int after_mincut_number_of_clusters = MincutOnlyPreProcess::to_be_mincut_clusters.size();
        
        this->WriteToLogFile("After iteration: work queue has " + std::to_string(after_mincut_number_of_clusters) + " clusters", Log::debug);
        this->WriteToLogFile(std::to_string(after_mincut_number_of_clusters) + " [connected components / clusters] to be mincut after a round of mincuts", Log::debug);
        /** SECTION MinCut Each Connected Component END **/

        /** SECTION Check If All Clusters Are Well-Connected START **/
        if(after_mincut_number_of_clusters == 0) {
            this->WriteToLogFile("========================================", Log::debug);
            this->WriteToLogFile("ALL CLUSTERS ARE WELL-CONNECTED!", Log::debug);
            this->WriteToLogFile("========================================", Log::debug);
            this->WriteToLogFile("Total iterations: " + std::to_string(iter_count + 1), Log::debug);
            this->WriteToLogFile("Final cluster count: " + std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()), Log::debug);
            
            this->WriteToLogFile("all clusters are (well) connected", Log::info);
            this->WriteToLogFile("Total number of iterations: " + std::to_string(iter_count + 1), Log::info);
            break;
        }
        /** SECTION Check If All Clusters Are Well-Connected END **/

        iter_count ++;
    }

    this->WriteToLogFile("======== Writing Output ========", Log::info);
    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteToLogFile("Output file: " + this->output_file, Log::info);
    this->WriteToLogFile("Writing " + std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()) + " clusters", Log::info);
    
    this->WriteClusterQueue(MincutOnlyPreProcess::done_being_mincut_clusters);
    
    this->WriteToLogFile("========================================", Log::info);
    this->WriteToLogFile("MincutOnlyPreProcess completed successfully!", Log::info);
    this->WriteToLogFile("========================================", Log::info);
    
    return 0;
}