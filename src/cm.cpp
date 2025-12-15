#include "cm.h"
#include <mpi_telemetry.h>
#include <igraph.h>
#include <string>
#include <unordered_map>

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
    this->WriteToLogFile("Loading the initial graph" , Log::info, my_rank);

    printf("my_rank: %d load graph\n", my_rank);
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_t graph;
    std::unordered_map<long long, long long> original_to_new_id_unordered_map;
    MMapGraphLoader::LoadEdgelistMMap(this->edgelist, &graph,&original_to_new_id_unordered_map,false);
    SetIgraphAllEdgesWeight(&graph, 1.0);
    printf("my_rank: %d finished loading graph\n", my_rank);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info, my_rank);
    /* std::cerr << EAN(&graph, "weight", 0) << std::endl; */
    
    int graph_vcount = igraph_vcount(&graph);
    string edge_count;
    edge_count = "my_rank: %d before rice edge_count " + to_string(igraph_ecount(&graph));
    this -> WriteToLogFile(edge_count, Log::info, my_rank);
    printf("my_rank: %d before rice edge_count: %d\n", my_rank, igraph_ecount(&graph));
    std::unordered_map<long long, long long> node_id_to_cluster_id_unordered_map;
    
    if(this->existing_clustering == "") {
        /* int seed = uni(rng); */
        // not working
        printf("my_rank: %d no clustering data\n", my_rank);
        int seed = 0;
        std::map<long long, long long> node_id_to_cluster_id_map = ConstrainedClustering::GetCommunities("", this->algorithm, seed, this->clustering_parameter, &graph);
        // output_map(node_id_to_cluster_id_map);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_map);
    } else if(this->existing_clustering != "") {
        // non mpi
        this->WriteToLogFile("Loading the new id to cluster id map" , Log::debug, my_rank);
        MMapGraphLoader::LoadClusteringMMap(this->existing_clustering, &node_id_to_cluster_id_unordered_map, original_to_new_id_unordered_map);
        this->WriteToLogFile("Finished loading the new id to cluster id map" , Log::debug, my_rank);    
        // output_map(original_to_new_id_map_nonmpi);
        this->WriteToLogFile("Removing Inter cluster edges vertices: " + std::to_string(igraph_vcount(&graph)) + " edges: " + std::to_string(igraph_ecount(&graph)) , Log::debug, my_rank);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_unordered_map, this-> num_processors);
        this->WriteToLogFile("Finished removing Inter cluster edges vertices: " + std::to_string(igraph_vcount(&graph)) + " edges: " + std::to_string(igraph_ecount(&graph)) , Log::debug, my_rank);
    }

    edge_count = "my_rank: %d after rice edge_count " + to_string(igraph_ecount(&graph));
    this -> WriteToLogFile(edge_count, Log::info, my_rank);

    printf("my_rank: %d after rice edge_count: %d\n", my_rank, igraph_ecount(&graph));

    /** SECTION Get Connected Components START **/
    // std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    //printf("my_rank: %d connected components start\n", my_rank);
    this->WriteToLogFile("Getting all the connected components" , Log::debug, my_rank);
    std::vector<std::vector<long long>> connected_components_vector = ConstrainedClustering::GetConnectedComponentsDistributed(&graph, node_id_to_cluster_id_unordered_map, my_rank, nprocs);

    this->WriteToLogFile("Finished Getting all the connected components" , Log::debug, my_rank);
    // store the results into the queue that each thread pulls from
    
    int cc_count = connected_components_vector.size();

    int cc_start = 0;
    int cc_end = cc_count;
    this->WriteToLogFile("CC_count: " + std::to_string(cc_count) + " cc_start: " + std::to_string(cc_start) + " cc_end: " + std::to_string(cc_end), Log::info, my_rank);

    // non mpi
    // for(size_t i = 0; i < connected_components_vector.size(); i ++) {
    //     CM::to_be_mincut_clusters.push(connected_components_vector[i]);
    // }
    // igraph_t subgraph;
    // igraph_vector_int_t nodes_to_keep;
    // igraph_vector_int_t new_id_to_old_id_vector_map;
    // igraph_vector_int_init(&nodes_to_keep, my_work);
    // for(size_t i = cc_start; i < cc_end; i ++) {
    //     for(int j = 0; j < connected_components_vector[i].size();j++)
    //     VECTOR(nodes_to_keep)[i] = connected_components_vector[i][j];
    // }
    // igraph_induced_subgraph_map(graph, &subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_vector_map);
    for(size_t i = cc_start; i < cc_end; i ++) {
        CM::to_be_mincut_clusters.push(connected_components_vector[i]);
    }
    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    long long previous_cluster_id = 0;
    for (int i = 0; i < nprocs; i++) {
      mincut_continue[i] = 1;
    }
    //    while (!CM::to_be_mincut_clusters.empty()) {
    int iter_count = 0;
    int after_mincut_number_of_clusters = 0;
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
        int mincut_continue_mr = !CM::to_be_mincut_clusters.empty();
        MPI_Allgather(&mincut_continue_mr, 1, MPI_INT, mincut_continue, 1, MPI_INT, MPI_COMM_WORLD);
        iter_count ++;
    }

    this->WriteToLogFile("my_rank: " + to_string(my_rank) + " Writing output to: " + this->output_file, Log::info, my_rank);
    previous_cluster_id = this->WriteClusterQueueMPI(&CM::done_being_clustered_clusters, &graph, cc_start, previous_cluster_id, iter_count, opCount);
        


    //this->WriteToLogFile("Writing output to: " + this->output_file, Log::info, my_rank);
    //this->WriteClusterQueue(CM::done_being_clustered_clusters, &graph, cc_start);
    igraph_destroy(&graph);
    return 0;
}
