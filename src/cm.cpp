#include "cm.h"
#include <mpi_telemetry.h>
#include <igraph.h>
#include <string>

using namespace std;
#define ROOT 0

int mincut_continue[100];

template <typename T>
void output_map(std::map<T, int> map) {
    for (auto const& [key, val]: map) {
        cout << key << " " << val << "\n";
    }

}
template <typename T>
void output_arr(T * arr, int arr_size) {
    for (int i = 0; i< arr_size; i++) {
        cout << arr[i] << " ";
    }
    cout << "\n";
}
template <typename T>
void output_vec(std::vector<T> vec) {
    for (auto item: vec) {
        cout << item << " ";
    }
    cout << "\n";
}

void build_displacements(int * displacements, int* size_array, int nprocs) {
    displacements[0] = 0;
    for (int i = 1; i < nprocs; i ++) {
        displacements[i] = displacements[i-1] + size_array[i-1];
    }
}

void convert_vec_to_arr(igraph_vector_int_t vec, int* arr, int vec_size) {
    for (int i = 0 ; i < vec_size; i++) {
        arr[i] = VECTOR(vec)[i];
    }
}

void convert_arr_to_vec(igraph_vector_int_t vec, int* arr, int arr_size) {
    for (int i= 0; i < arr_size; i++) {
        VECTOR(vec)[i] = arr[i];
    }
}

bool checkMC(int * arr, int arr_size) {
  for (int i = 0; i < arr_size; i++) {
    if (arr[i] == 1) {
      return true;
    }
  }
  return false;
}

int CM::main(int my_rank, int nprocs, uint64_t* opCount) {
    /* std::random_device rd; */
    /* std::mt19937 rng{rd()}; */
    /* std::uniform_int_distribution<int> uni(0, 100000); */
    this->WriteToLogFile("Loading the distributed graph" , Log::info, my_rank);

    printf("my_rank: %d load graph\n", my_rank);
    FILE* edgelist_file = fopen((output_header + "_" + to_string(my_rank) + ".tsv" ).c_str(), "r");
    igraph_t graph;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_ncol(&graph, edgelist_file, NULL, 1, IGRAPH_ADD_WEIGHTS_IF_PRESENT, IGRAPH_UNDIRECTED);
    if(!igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_EDGE, "weight")) {
        SetIgraphAllEdgesWeight(&graph, 1.0);
    }
    /* igraph_read_graph_edgelist(&graph, edgelist_file, 0, false); */
    fclose(edgelist_file);
    printf("my_rank: %d finished loading graph\n", my_rank);
    this->WriteToLogFile("Finished loading the initial graph with num_vertices: " + to_string(igraph_vcount(&graph)) + " num_edges: " + to_string(igraph_vcount(&graph)) , Log::info, my_rank);
    /* std::cerr << EAN(&graph, "weight", 0) << std::endl; */
    
    int graph_vcount = igraph_vcount(&graph);
    int after_mincut_number_of_clusters = -2;
    int iter_count = 0;

    /** SECTION Get Connected Components START **/
    std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    // store the results into the queue that each thread pulls from
    
    int cc_count = connected_components_vector.size();
    int cc_count_arr[nprocs];
    int cc_displacements_arr[nprocs];
    MPI_Allgather(&cc_count, 1, MPI_INT, cc_count_arr, 1, MPI_INT, MPI_COMM_WORLD);
    build_displacements(cc_displacements_arr, cc_count_arr, nprocs);

    // for (int i = 0 ; i < cc_count; i++) {
    //     output_vec(connected_components_vector[i]);
    // }
    // int cc_my_work = cc_count/nprocs;
    // if (cc_count%nprocs)
    //     cc_my_work++;
    // int cc_start = my_rank*cc_my_work;
    // int cc_end = (my_rank+1)*cc_my_work;
    // if (my_rank == nprocs-1)
    //     cc_end = cc_count;
    int cc_start = cc_displacements_arr[my_rank];
    int cc_end = cc_start + cc_count;
    //this->WriteToLogFile("CC_count: " + std::to_string(cc_count) + " cc_start: " + std::to_string(cc_start) + " cc_end: " + std::to_string(cc_end), Log::info, my_rank);

    for(size_t i = 0; i < cc_count; i ++) {
        CM::to_be_mincut_clusters.push(connected_components_vector[i]);
    }
    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    int previous_cluster_id = 0;
    for (int i = 0; i < nprocs; i++) {
      mincut_continue[i] = 1;
    }
    //    while (!CM::to_be_mincut_clusters.empty()) {
    while (checkMC(mincut_continue, nprocs)) {
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug, my_rank);
        if(iter_count % 10 == 0) {
            this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info, my_rank);
        }

        /** SECTION MinCutOnceAndCluster Each Connected Component START **/
        this->WriteToLogFile(std::to_string(CM::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::debug, my_rank);
        /* before_mincut_number_of_clusters = CM::to_be_mincut_clusters.size(); */
        /* if a thread gets a cluster {-1}, then they know processing is done and they can stop working */
        for(int i = 0; i < this->num_processors; i ++) {
            CM::to_be_mincut_clusters.push({-1});
        }
        /* start the threads */
        std::vector<std::thread> thread_vector;
        for(int i = 0; i < this->num_processors; i ++) {
            /* int seed = uni(rng); */
            int seed = 0;
            thread_vector.push_back(std::thread(CM::MinCutOrClusterWorker, &graph, this->algorithm, seed, this->clustering_parameter));
        }
        /* get the result back from threads */
        /* the results from each thread gets stored in to_be_clustered_clusters */
        for(size_t thread_index = 0; thread_index < thread_vector.size(); thread_index ++) {
            thread_vector[thread_index].join();
        }
        this->WriteToLogFile(std::to_string(CM::to_be_clustered_clusters.size()) + " [connected components / clusters] to be clustered after a round of mincuts", Log::debug, my_rank);
        this->WriteToLogFile(std::to_string(CM::done_being_clustered_clusters.size() - previous_done_being_clustered_size) + " [connected components / clusters] were found to be well connected", Log::debug, my_rank);
        previous_done_being_clustered_size = CM::done_being_clustered_clusters.size();
        /** SECTION MinCutOnceAndCluster Each Connected Component END **/

        /** SECTION Check If All Clusters Are Well-Connected START **/
        after_mincut_number_of_clusters = CM::to_be_clustered_clusters.size();
        if(after_mincut_number_of_clusters == 0) {
            this->WriteToLogFile("all clusters are well-connected", Log::info, my_rank);
            this->WriteToLogFile("Total number of iterations: " + std::to_string(iter_count + 1), Log::info, my_rank);
            mincut_continue[my_rank] = 0;
        }
        else {
            mincut_continue[my_rank] = 1;
        }
        /** SECTION Check If All Clusters Are Well-Connected END **/

        while(!CM::to_be_clustered_clusters.empty()) {
            CM::to_be_mincut_clusters.push(CM::to_be_clustered_clusters.front());
            CM::to_be_clustered_clusters.pop();
        }
        this->WriteToLogFile("my_rank: " + to_string(my_rank) + " Writing output to: " + this->output_file, Log::info, my_rank);
        previous_cluster_id = this->WriteClusterQueueMPI(&CM::done_being_clustered_clusters, &graph, cc_start, previous_cluster_id, iter_count, opCount);
        int mincut_continue_mr = !CM::to_be_mincut_clusters.empty();
        MPI_Allgather(&mincut_continue_mr, 1, MPI_INT, mincut_continue, 1, MPI_INT, MPI_COMM_WORLD);
        iter_count ++;
    }
    igraph_destroy(&graph);
    return 0;
}
