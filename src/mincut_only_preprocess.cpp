#include "mincut_only_preprocess.h"
#include "mmap_subgraph_loader.h"

int MincutOnlyPreProcess::main() {

    this->WriteToLogFile("Parsing connectedness criterion" , Log::info);
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
    long before_mincut_number_of_clusters = 0;
    int iter_count = 0;


    this->WriteToLogFile("Start Mincut" , Log::info);

    int previous_done_being_clustered_size = 0;
    int previous_cluster_id = 0;
    // for (int i = 0; i < nprocs; i++) {
    //     mincut_continue[i] = 1;
    // }
    // store the results into the queue that each thread pulls from
    // for(size_t i = 0; i < connected_components_vector.size(); i ++) {
    //     MincutOnly::to_be_mincut_clusters.push(connected_components_vector[i]);
    // }
    
    /* std::cerr << "iter num: " << std::to_string(iter_count) << std::endl; */
    this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
    if(iter_count % 10000 == 0) {
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
        this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::info);
    }

    /** SECTION MinCut Each Connected Component START **/
    this->WriteToLogFile(std::to_string(subgraph_edges_vector.size()) + " [connected components / clusters] to be mincut", Log::debug);
    before_mincut_number_of_clusters = subgraph_edges_vector.size();
    /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
    /* std::cerr << "num clusters to be processed: " << std::to_string(before_mincut_number_of_clusters) << std::endl; */
    long current_components_vector_index = 0;
    while (current_components_vector_index < before_mincut_number_of_clusters) {
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
            if (current_components_vector_index < before_mincut_number_of_clusters ) {
                thread_vector.push_back(std::thread(MincutOnlyPreProcess::ComputeMinCutRecursive,subgraph_edges_vector[current_components_vector_index], current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log));
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
    // this->WriteToLogFile(std::to_string(MincutOnly::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut after a round of mincuts", Log::debug);
    /** SECTION MinCut Each Connected Component END **/

    /** SECTION Check If All Clusters Are Well-Connected START **/
    // after_mincut_number_of_clusters = MincutOnly::to_be_mincut_clusters.size();

    this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size())+ " [connected components / clusters] mincut after a round of mincuts", Log::debug);
        
    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteClusterQueue(MincutOnlyPreProcess::done_being_mincut_clusters);
    return 0;
}
