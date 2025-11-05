#include "mincut_only.h"
#include "mpi_telemetry.h"

int mincut_continue_wcc[100];

bool checkMc(int * arr, int arr_size) {
  for (int i = 0; i < arr_size; i++) {
    if (arr[i] == 1) {
      return true;
    }
  }
  return false;
}

int MincutOnly::main(int my_rank, int nprocs, uint64_t* opCount) {
    this->WriteToLogFile("Loading the initial graph" , Log::info, my_rank);
    FILE* edgelist_file = fopen(this->edgelist.c_str(), "r");
    igraph_t graph;
    /* igraph_read_graph_edgelist(&graph, edgelist_file, 0, false); */
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_ncol(&graph, edgelist_file, NULL, 1, IGRAPH_ADD_WEIGHTS_IF_PRESENT, IGRAPH_UNDIRECTED);
    /* if(!igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_EDGE, "weight")) { */
    /*     SetIgraphAllEdgesWeight(&graph, 1.0); */
    /* } */
    fclose(edgelist_file);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info);

    int before_mincut_number_of_clusters = -1;
    int after_mincut_number_of_clusters = -2;
    int iter_count = 0;

    ConstrainedClustering::initializeSlice(&graph);
    std::map<std::string, int> original_to_new_id_map = ConstrainedClustering::GetOriginalToNewIdMapDistributed(&graph);
    std::map<int, int> new_node_id_to_cluster_id_map = ConstrainedClustering::ReadCommunities(original_to_new_id_map, this->existing_clustering);
    ConstrainedClustering::RemoveInterClusterEdgesDistributed(&graph, new_node_id_to_cluster_id_map);

    /** SECTION Get Connected Components START **/
    std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    /** SECTION Get Connected Components END **/
       int cc_count = connected_components_vector.size();
    // for (int i = 0 ; i < cc_count; i++) {
    //     output_vec(connected_components_vector[i]);
    // }
    int cc_my_work = cc_count/nprocs;
    if (cc_count%nprocs)
        cc_my_work++;
    int cc_start = my_rank*cc_my_work;
    int cc_end = (my_rank+1)*cc_my_work;
    if (my_rank == nprocs-1)
        cc_end = cc_count;
    for (int i = 0; i < nprocs; i++) {
      mincut_continue_wcc[i] = 1;
    }
    int previous_done_being_clustered_size = 0;
    int previous_cluster_id = 0;
    if(this->connectedness_criterion == ConnectednessCriterion::Simple) {
        for(size_t i = cc_start; i < cc_end; i ++) {
            MincutOnly::done_being_mincut_clusters.push(connected_components_vector[i]);
        }
    } else if(this->connectedness_criterion == ConnectednessCriterion::Well) {
        // store the results into the queue that each thread pulls from
        for(size_t i = cc_start; i < cc_end; i ++) {
            MincutOnly::to_be_mincut_clusters.push(connected_components_vector[i]);
        }
        while (checkMc(mincut_continue_wcc, nprocs)) {
            this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
            if(iter_count % 10000 == 0) {
                this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
                this->WriteToLogFile(std::to_string(MincutOnly::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::info);
            }

            /** SECTION MinCut Each Connected Component START **/
            this->WriteToLogFile(std::to_string(MincutOnly::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::debug);
            before_mincut_number_of_clusters = MincutOnly::to_be_mincut_clusters.size();
            /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
            for(int i = 0; i < this->num_processors; i ++) {
                MincutOnly::to_be_mincut_clusters.push({-1});
            }
            if(before_mincut_number_of_clusters > 1) {
                /* start the threads */
                std::vector<std::thread> thread_vector;
                for(int i = 0; i < this->num_processors; i ++) {
                    thread_vector.push_back(std::thread(MincutOnly::MinCutWorker, &graph, this->connectedness_criterion));
                }
                /* get the result back from threads */
                /* the results from each thread gets stored in to_be_clustered_clusters */
                for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
                    thread_vector[thread_index].join();
                }
            } else {
                MincutOnly::MinCutWorker(&graph, this->connectedness_criterion);
            }
            previous_done_being_clustered_size = MincutOnly::done_being_mincut_clusters.size();
            this->WriteToLogFile(std::to_string(MincutOnly::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut after a round of mincuts", Log::debug);
            /** SECTION MinCut Each Connected Component END **/

            /** SECTION Check If All Clusters Are Well-Connected START **/
            after_mincut_number_of_clusters = MincutOnly::to_be_mincut_clusters.size();
            if(after_mincut_number_of_clusters == 0) {
                this->WriteToLogFile("all clusters are (well) connected", Log::info);
                this->WriteToLogFile("Total number of iterations: " + std::to_string(iter_count + 1), Log::info);
                mincut_continue_wcc[my_rank] = 0;
            } else { 
                mincut_continue_wcc[my_rank] = 1;
            }
            /** SECTION Check If All Clusters Are Well-Connected END **/
            this->WriteToLogFile("my_rank: " + to_string(my_rank) + " Writing output to: " + this->output_file, Log::info, my_rank);
            previous_cluster_id = this->WriteClusterQueueMPI(&MincutOnly::done_being_mincut_clusters, &graph, cc_start, previous_cluster_id, iter_count, opCount);
            int mincut_continue_mr = !MincutOnly::to_be_mincut_clusters.empty();
            MPI_Allgather(&mincut_continue_mr, 1, MPI_INT, mincut_continue_wcc, 1, MPI_INT, MPI_COMM_WORLD);
            iter_count ++;
        }
    }


    //this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    //this->WriteClusterQueue(MincutOnly::done_being_mincut_clusters, &graph, cc_start);
    igraph_destroy(&graph);
    return 0;
}
