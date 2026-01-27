#include "mincut_only_preprocess.h"
#include "mmap_subgraph_loader.h"
#include <mpi_telemetry.h>

bool checkMC(int * arr, int arr_size) {
  for (int i = 0; i < arr_size; i++) {
    if (arr[i] == 1) {
      return true;
    }
  }
  return false;
}

int mincut_continue[100];

int MincutOnlyPreProcess::main(int my_rank, int nprocs, uint64_t * opCount) {

    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("Starting MincutOnlyPreProcess (MPI Recursive)", Log::info, my_rank);
    this->WriteToLogFile("Rank " + std::to_string(my_rank) + " of " + std::to_string(nprocs), Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);
    
    this->WriteToLogFile("Parsing connectedness criterion: " + this->connectedness_criterion, Log::info, my_rank);
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
        this->WriteToLogFile("Detected logarithmic criterion: C=" + std::to_string(connectedness_criterion_c) + ", x=" + std::to_string(connectedness_criterion_x), Log::info, my_rank);
    } else if (n_caret_position != std::string::npos) {
        // is cN^x
        current_connectedness_criterion = ConnectednessCriterion::Exponential;
        if (n_caret_position == 0) {
            connectedness_criterion_c = 1;
        } else {
            connectedness_criterion_c = std::stod(this->connectedness_criterion.substr(0, n_caret_position));
        }
        connectedness_criterion_x = std::stod(this->connectedness_criterion.substr(n_caret_position + 2));
        this->WriteToLogFile("Detected exponential criterion: C=" + std::to_string(connectedness_criterion_c) + ", x=" + std::to_string(connectedness_criterion_x), Log::info, my_rank);
    } else if (this->connectedness_criterion != "0") {
        // wasn't log or exponent so if it isn't 0 then it's an error
        // exit
        this->WriteToLogFile("Could not parse connectedness_criterion: " + this->connectedness_criterion, Log::error, my_rank);
        this->WriteToLogFile("Accepted forms are Clog_x(n), Cn^x, or 0 where C and x are integers" , Log::error, my_rank);
        return 1;
    }
    if (current_connectedness_criterion == ConnectednessCriterion::Simple) {
        this->WriteToLogFile("Running with CC mode (mincut of each cluster has to be greater than 0)" , Log::info, my_rank);
    } else if (current_connectedness_criterion == ConnectednessCriterion::Logarithimic) {
        this->WriteToLogFile("Running with WCC mode (mincut of each cluster has to be greater than " + std::to_string(connectedness_criterion_c) + " times log base " + std::to_string(connectedness_criterion_x) + " of n" , Log::info, my_rank);
        pre_computed_log = connectedness_criterion_c / std::log(connectedness_criterion_x);
        this->WriteToLogFile("Pre-computed log factor: " + std::to_string(pre_computed_log), Log::info, my_rank);
    } else if (current_connectedness_criterion == ConnectednessCriterion::Exponential) {
        this->WriteToLogFile("Running with WCC mode (mincut of each cluster has to be greater than " + std::to_string(connectedness_criterion_c) + " times n to the power of " + std::to_string(connectedness_criterion_x), Log::info, my_rank);
    } else {
        // should not possible to reach
        this->WriteToLogFile("Unexpected connectedness criterion state", Log::error, my_rank);
        exit(1);
    }
    
    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("Loading graph data (recursive approach)", Log::info, my_rank);
    this->WriteToLogFile("Input file: " + this->edgelist, Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);
    
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> subgraph_edges_vector = MMapSubgraphLoader::LoadEdgelistMMap(this->edgelist, &original_to_new_id_unordered_map,false);

    this->WriteToLogFile("Finished loading graph", Log::info, my_rank);
    this->WriteToLogFile("Connected components loaded: " + std::to_string(subgraph_edges_vector.size()), Log::info, my_rank);
    this->WriteToLogFile("Unique nodes: " + std::to_string(original_to_new_id_unordered_map.size()), Log::info, my_rank);
    int before_mincut_number_of_clusters = 0;
    int after_mincut_number_of_clusters = 0;
    int iter_count = 0;     
    /** SECTION Get Connected Components END **/
    this->WriteToLogFile("Finished Getting Connected Components size: " + std::to_string(subgraph_edges_vector.size()) , Log::info, my_rank);
    
    int current_components_vector_index = 0;

    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("Initializing recursive mincut processing", Log::info, my_rank);
    this->WriteToLogFile("Number of threads per rank: " + std::to_string(this->num_processors), Log::info, my_rank);
    this->WriteToLogFile("Total components to process: " + std::to_string(subgraph_edges_vector.size()), Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);

    int previous_done_being_clustered_size = 0;
    int previous_cluster_id = 0;
    
    // Initialize mincut_continue array
    for (int i = 0; i < nprocs; i++) {
        mincut_continue[i] = 1;
    }
    this->WriteToLogFile("Initialized mincut_continue array for " + std::to_string(nprocs) + " ranks", Log::debug, my_rank);
    
    this->WriteToLogFile("Starting recursive processing (single pass, no iteration)", Log::info, my_rank);
    this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug, my_rank);
    
    if(iter_count % 10000 == 0) {
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info, my_rank);
        this->WriteToLogFile(std::to_string(subgraph_edges_vector.size()) + " [connected components / clusters] to be mincut", Log::info, my_rank);
    }

    /** SECTION MinCut Each Connected Component START **/
    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("Starting recursive mincut on all components", Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);

    this->WriteToLogFile(std::to_string(subgraph_edges_vector.size()) + " [connected components / clusters] to be mincut", Log::debug, my_rank);
    before_mincut_number_of_clusters = subgraph_edges_vector.size();
    this->WriteToLogFile("Total clusters to process: " + std::to_string(before_mincut_number_of_clusters), Log::info, my_rank);
    /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
    int batch_count = 0;
    this->WriteToLogFile("Starting batch processing with " + std::to_string(this->num_processors) + " threads per batch", Log::info, my_rank);
    
    while (current_components_vector_index < before_mincut_number_of_clusters) {
        batch_count++;
        this->WriteToLogFile("----------------------------------------", Log::debug, my_rank);
        this->WriteToLogFile("Processing batch " + std::to_string(batch_count), Log::debug, my_rank);
        this->WriteToLogFile("Starting component index: " + std::to_string(current_components_vector_index) + " / " + std::to_string(before_mincut_number_of_clusters), Log::debug, my_rank);
        this->WriteToLogFile("----------------------------------------", Log::debug, my_rank);
        
        /* start the threads */
        std::vector<std::thread> thread_vector;
        int threads_spawned = 0;
        
        for(int i = 0; i < this->num_processors; i ++) {
            if (current_components_vector_index < before_mincut_number_of_clusters) {
                this->WriteToLogFile("Spawning thread " + std::to_string(i) + " for component " + std::to_string(current_components_vector_index) + 
                                    " (size: " + std::to_string(subgraph_edges_vector[current_components_vector_index].size()) + " edges)", Log::debug, my_rank);
                
                thread_vector.push_back(std::thread(MincutOnlyPreProcess::ComputeMinCutRecursive, 
                                                    subgraph_edges_vector[current_components_vector_index], 
                                                    current_connectedness_criterion, 
                                                    connectedness_criterion_c, 
                                                    connectedness_criterion_x, 
                                                    pre_computed_log));
                threads_spawned++;
                current_components_vector_index++;
            } else {
                this->WriteToLogFile("No more components to process, stopping thread spawn at " + std::to_string(i) + " threads", Log::debug, my_rank);
                break;
            }
        }
        
        this->WriteToLogFile("Spawned " + std::to_string(threads_spawned) + " threads for batch " + std::to_string(batch_count), Log::debug, my_rank);
        this->WriteToLogFile("Waiting for threads to complete...", Log::debug, my_rank);
        
        /* get the result back from threads */
        /* the results from each thread gets stored in to_be_clustered_clusters */
        for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
            thread_vector[thread_index].join();
        }
        
        this->WriteToLogFile("Batch " + std::to_string(batch_count) + " completed", Log::debug, my_rank);
        this->WriteToLogFile("Progress: " + std::to_string(current_components_vector_index) + " / " + std::to_string(before_mincut_number_of_clusters) + " components processed", Log::debug, my_rank);
        
        if(batch_count % 10 == 0) {
            this->WriteToLogFile("========================================", Log::info, my_rank);
            this->WriteToLogFile("Progress Update - Batch " + std::to_string(batch_count), Log::info, my_rank);
            this->WriteToLogFile("Components processed: " + std::to_string(current_components_vector_index) + " / " + std::to_string(before_mincut_number_of_clusters), Log::info, my_rank);
            this->WriteToLogFile("Clusters completed so far: " + std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()), Log::info, my_rank);
            this->WriteToLogFile("========================================", Log::info, my_rank);
        }
    }
    
    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("All local components processed", Log::info, my_rank);
    this->WriteToLogFile("Total batches processed: " + std::to_string(batch_count), Log::info, my_rank);
    this->WriteToLogFile("Local clusters completed: " + std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()), Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);
    
    this->WriteToLogFile("Waiting at MPI barrier for all ranks to complete", Log::info, my_rank);
    MPI_Barrier(my_rank, iter_count, 5, opCount);
    this->WriteToLogFile("All ranks have completed processing", Log::info, my_rank);

    this->WriteToLogFile(std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()) + " [connected components / clusters] completed after recursive mincut", Log::debug, my_rank);

    /** SECTION MinCut Each Connected Component END **/

    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("Starting MPI output aggregation", Log::info, my_rank);
    this->WriteToLogFile("Output file: " + this->output_file, Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);
    
    previous_cluster_id = this->WriteClusterQueueMPI(&MincutOnlyPreProcess::done_being_mincut_clusters, opCount, my_rank, nprocs);
    
    iter_count++;

    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("MincutOnlyPreProcess (Recursive) completed successfully", Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);

    return 0;
}