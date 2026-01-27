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
    this->WriteToLogFile("Starting MincutOnlyPreProcess (MPI Distributed)", Log::info, my_rank);
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
    this->WriteToLogFile("Loading graph data", Log::info, my_rank);
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

    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("Initializing mincut processing", Log::info, my_rank);
    this->WriteToLogFile("Number of threads per rank: " + std::to_string(this->num_processors), Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);

    int previous_done_being_clustered_size = 0;
    int previous_cluster_id = 0;
    
    // Initialize mincut_continue array
    for (int i = 0; i < nprocs; i++) {
        mincut_continue[i] = 1;
    }
    this->WriteToLogFile("Initialized mincut_continue array for " + std::to_string(nprocs) + " ranks", Log::debug, my_rank);
    
    // store the results into the queue that each thread pulls from
    this->WriteToLogFile("Populating initial work queue with " + std::to_string(subgraph_edges_vector.size()) + " components", Log::info, my_rank);
    for(size_t i = 0; i < subgraph_edges_vector.size(); i ++) {
        MincutOnlyPreProcess::to_be_mincut_clusters.push(subgraph_edges_vector[i]);
    }
    this->WriteToLogFile("Work queue populated successfully", Log::info, my_rank);
    
    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("Starting iterative mincut processing", Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);
    while (checkMC(mincut_continue, nprocs)) {
        this->WriteToLogFile("----------------------------------------", Log::debug, my_rank);
        this->WriteToLogFile("Iteration " + std::to_string(iter_count), Log::debug, my_rank);
        this->WriteToLogFile("----------------------------------------", Log::debug, my_rank);
        
        if(iter_count % 10000 == 0) {
            this->WriteToLogFile("========================================", Log::info, my_rank);
            this->WriteToLogFile("Progress Update - Iteration " + std::to_string(iter_count), Log::info, my_rank);
            this->WriteToLogFile("Clusters remaining in queue: " + std::to_string(MincutOnlyPreProcess::to_be_mincut_clusters.size()), Log::info, my_rank);
            this->WriteToLogFile("Clusters completed: " + std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()), Log::info, my_rank);
            this->WriteToLogFile("========================================", Log::info, my_rank);
        }

        /** SECTION MinCut Each Connected Component START **/
        before_mincut_number_of_clusters = MincutOnlyPreProcess::to_be_mincut_clusters.size();
        this->WriteToLogFile("Clusters to process this iteration: " + std::to_string(before_mincut_number_of_clusters), Log::debug, my_rank);
        
        /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
        if(before_mincut_number_of_clusters > 1) {
            this->WriteToLogFile("Spawning " + std::to_string(this->num_processors) + " worker threads", Log::debug, my_rank);
            
            /* start the threads */
            for(int i = 0; i < this->num_processors; i ++) {
                MincutOnlyPreProcess::to_be_mincut_clusters.push({-1});
            }
            
            std::vector<std::thread> thread_vector;
            for(int i = 0; i < this->num_processors; i ++) {
                thread_vector.push_back(std::thread(MincutOnlyPreProcess::MinCutWorker, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log));
            }
            
            this->WriteToLogFile("Waiting for " + std::to_string(this->num_processors) + " threads to complete", Log::debug, my_rank);
            
            /* get the result back from threads */
            /* the results from each thread gets stored in to_be_clustered_clusters */
            for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
                thread_vector[thread_index].join();
            }
            
            this->WriteToLogFile("All threads completed", Log::debug, my_rank);
        } else {
            this->WriteToLogFile("Only 1 cluster, processing with single worker", Log::debug, my_rank);
            MincutOnlyPreProcess::to_be_mincut_clusters.push({-1});
            MincutOnlyPreProcess::MinCutWorker(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log);
        }
        after_mincut_number_of_clusters = MincutOnlyPreProcess::to_be_mincut_clusters.size();
        this->WriteToLogFile("Clusters remaining after processing: " + std::to_string(after_mincut_number_of_clusters), Log::debug, my_rank);
        /** SECTION MinCut Each Connected Component END **/

        /** SECTION Check If All Clusters Are Well-Connected START **/
        if(after_mincut_number_of_clusters == 0) {
            this->WriteToLogFile("========================================", Log::info, my_rank);
            this->WriteToLogFile("All local clusters are well-connected", Log::info, my_rank);
            this->WriteToLogFile("Iterations completed: " + std::to_string(iter_count + 1), Log::info, my_rank);
            this->WriteToLogFile("========================================", Log::info, my_rank);
            mincut_continue[my_rank] = 0;
        } else {
            mincut_continue[my_rank] = 1;
        }
        /** SECTION Check If All Clusters Are Well-Connected END **/
        
        this->WriteToLogFile("Preparing for MPI synchronization", Log::debug, my_rank);
        int mincut_continue_mr = !MincutOnlyPreProcess::to_be_mincut_clusters.empty();
        
        this->WriteToLogFile("Broadcasting convergence status across all ranks", Log::debug, my_rank);
        MPI_Allgather(&mincut_continue_mr, 1, MPI_INT, mincut_continue, 1, MPI_INT, MPI_COMM_WORLD);
        
        // Log convergence status from all ranks
        int still_running = 0;
        for (int i = 0; i < nprocs; i++) {
            if (mincut_continue[i]) still_running++;
        }
        this->WriteToLogFile("Ranks still processing: " + std::to_string(still_running) + " / " + std::to_string(nprocs), Log::debug, my_rank);
        
        iter_count++;
    }
    
    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("All ranks completed mincut processing", Log::info, my_rank);
    this->WriteToLogFile("Total iterations: " + std::to_string(iter_count), Log::info, my_rank);
    this->WriteToLogFile("Final clusters found: " + std::to_string(MincutOnlyPreProcess::done_being_mincut_clusters.size()), Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);
    
    this->WriteToLogFile("Writing final output via MPI aggregation", Log::info, my_rank);
    this->WriteClusterQueueMPI(&MincutOnlyPreProcess::done_being_mincut_clusters, opCount, my_rank, nprocs);

    this->WriteToLogFile("========================================", Log::info, my_rank);
    this->WriteToLogFile("MincutOnlyPreProcess completed successfully", Log::info, my_rank);
    this->WriteToLogFile("========================================", Log::info, my_rank);

    return 0;
}