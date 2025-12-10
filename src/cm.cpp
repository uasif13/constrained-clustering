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

    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_t graph;
    std::unordered_map<int, int> original_to_new_id_map_nonmpi;
    MMapGraphLoader::LoadEdgelistMMap(this->edgelist, &graph,&original_to_new_id_map_nonmpi,false);
    SetIgraphAllEdgesWeight(&graph, 1.0);

    this->WriteToLogFile("Finished loading the initial graph" , Log::info, my_rank);
    /* std::cerr << EAN(&graph, "weight", 0) << std::endl; */
    
    int graph_vcount = igraph_vcount(&graph);
    string edge_count;
    edge_count = "my_rank: %d before rice edge_count " + to_string(igraph_ecount(&graph));
    // this -> WriteToLogFile(edge_count, Log::info, my_rank);
    printf("my_rank: %d before rice edge_count: %d\n", my_rank, igraph_ecount(&graph));
    ClusteringData data;
    std::unordered_map<int, int> node_id_to_cluster_id_unordered_map;
    
    if(this->existing_clustering == "") {
        /* int seed = uni(rng); */
        // not working
        printf("my_rank: %d no clustering data\n", my_rank);
        int seed = 0;
        this->WriteToLogFile("Loading the communities" , Log::debug, my_rank);

        std::map<int, int> node_id_to_cluster_id_map = ConstrainedClustering::GetCommunities("", this->algorithm, seed, this->clustering_parameter, &graph);
        // output_map(node_id_to_cluster_id_map);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_map);
    } else if(this->existing_clustering != "") {
        // non mpis

        this->WriteToLogFile("Loading the new id to cluster id map" , Log::debug, my_rank);

        // output_map(original_to_new_id_map_nonmpi);
        MMapGraphLoader::LoadClusteringMMap(this->existing_clustering, &node_id_to_cluster_id_unordered_map, original_to_new_id_map_nonmpi);
        this->WriteToLogFile("Finished loading the new id to cluster id map" , Log::debug, my_rank);

        this->WriteToLogFile("Removing Inter cluster edges vertices: " + std::to_string(igraph_vcount(&graph)) + " edges: " + std::to_string(igraph_ecount(&graph)) , Log::debug, my_rank);

        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_unordered_map, this->num_processors);

        this->WriteToLogFile("Finished removing Inter cluster edges vertices: " + std::to_string(igraph_vcount(&graph)) + " edges: " + std::to_string(igraph_ecount(&graph)) , Log::debug, my_rank);

    }

    edge_count = "my_rank: %d after rice edge_count " + to_string(igraph_ecount(&graph));
    // this -> WriteToLogFile(edge_count, Log::info, my_rank);

    printf("my_rank: %d after rice edge_count: %d\n", my_rank, igraph_ecount(&graph));

    this->WriteToLogFile("Loading the node id to cluster id map" , Log::debug, my_rank);

    this->WriteToLogFile("Finished Loading the node id to cluster id map" , Log::debug, my_rank);


    /** SECTION Get Connected Components START **/
    this->WriteToLogFile("Getting all the connected components" , Log::debug, my_rank);

    std::vector<std::vector<int>> connected_components_vector;
    connected_components_vector = ConstrainedClustering::GetConnectedComponentsDistributed(&graph, node_id_to_cluster_id_unordered_map, my_rank, nprocs);

    this->WriteToLogFile("Finished Getting all the connected components" , Log::debug, my_rank);

    int cc_count = connected_components_vector.size();
    int cc_start = 0;
    int cc_end = cc_count;
    // this->WriteToLogFile("CC_count: " + std::to_string(cc_count) + " cc_start: " + std::to_string(cc_start) + " cc_end: " + std::to_string(cc_end), Log::info, my_rank);

    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    int previous_cluster_id = 0;

    int current_components_vector_index = 0;

    /** SECTION MinCutOnceAndCluster Each Connected Component START **/
    this->WriteToLogFile(std::to_string(cc_count) + " [connected components / clusters] to be mincut", Log::debug, my_rank);

    /* start the threads */
    while (current_components_vector_index < cc_count) {
        /* start the threads */
        std::vector<std::thread> thread_vector;
        for(int i = 0; i < this->num_processors; i ++) {
            if (current_components_vector_index < cc_count ) {
                thread_vector.push_back(std::thread(CM::MinCutOrClusterWorkerRecursive,connected_components_vector[current_components_vector_index], &graph, algorithm, 0, clustering_parameter));
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
    this->WriteToLogFile(std::to_string(CM::done_being_clustered_clusters.size())+ " [connected components / clusters] mincut after a round of mincuts", Log::debug);

    
    this->WriteToLogFile("my_rank: " + to_string(my_rank) + " Writing output to: " + this->output_file, Log::info, my_rank);
    // this->WriteClusterQueueMPI(&CM::done_being_clustered_clusters, &graph, cc_start, previous_cluster_id, 1, opCount);
    previous_cluster_id = MMapWriter::WriteDistributed(&CM::done_being_clustered_clusters, &graph, this->output_file, cc_start, previous_cluster_id, my_rank, nprocs);    
    igraph_destroy(&graph);
    return 0;
}
