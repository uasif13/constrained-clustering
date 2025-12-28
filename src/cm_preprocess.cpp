#include "cm_preprocess.h"
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

int CMPreprocess::main(int my_rank, int nprocs, uint64_t* opCount) {
    this -> WriteToLogFile("Read subgraph edgelist file", Log::info);
    std::unordered_map<long, long> original_to_new_id_unordered_map;
    /** SECTION Get Connected Components START **/
    std::vector<std::vector<long>> subgraph_edges_vector = MMapSubgraphLoader::LoadEdgelistMMap(this->edgelist, &original_to_new_id_unordered_map,false);

    this -> WriteToLogFile("Finished reading subgraph edgelist file", Log::info);

    // store the results into the queue that each thread pulls from
    
    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    
    int current_components_vector_index = 0;
    int cc_count = subgraph_edges_vector.size();

    /** SECTION MinCutOnceAndCluster Each Connected Component START **/
    this->WriteToLogFile(std::to_string(cc_count) + " [connected components / clusters] to be mincut", Log::debug);

    /* start the threads */
    while (current_components_vector_index < cc_count) {
        /* start the threads */
        std::vector<std::thread> thread_vector;
        for(int i = 0; i < this->num_processors; i ++) {
            if (current_components_vector_index < cc_count ) {
                thread_vector.push_back(std::thread(CMPreprocess::MinCutOrClusterWorkerRecursive,subgraph_edges_vector[current_components_vector_index], algorithm, 0, clustering_parameter));
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
    this->WriteToLogFile(std::to_string(CMPreprocess::done_being_clustered_clusters.size())+ " [connected components / clusters] mincut after a round of mincuts", Log::debug);

    
    this->WriteToLogFile("my_rank: " + to_string(my_rank) + " Writing output to: " + this->output_file, Log::info, my_rank);
    this->WriteClusterQueueMPI(&CMPreprocess::done_being_clustered_clusters, 0, 0, 1, opCount);
    // previous_cluster_id = MMapWriter::WriteDistributed(&CM::done_being_clustered_clusters, &graph, this->output_file, cc_start, previous_cluster_id, my_rank, nprocs);    
    return 0;
}
