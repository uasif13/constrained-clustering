#include "cm.h"
#include <mpi.h>
#include <igraph.h>

using namespace std;
#define ROOT 0

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

void build_displacements(int * displacements, int* size_array, int nprocs) {
    displacements[0] = 0;
    for (int i = 1; i < nprocs; i ++) {
        displacements[i] = displacements[i-1] + size_array[i];
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

int CM::main(int my_rank, int nprocs) {
    /* std::random_device rd; */
    /* std::mt19937 rng{rd()}; */
    /* std::uniform_int_distribution<int> uni(0, 100000); */
    this->WriteToLogFile("Loading the initial graph" , Log::info);

    printf("my_rank: %d load graph\n", my_rank);
    FILE* edgelist_file = fopen(this->edgelist.c_str(), "r");
    igraph_t graph;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_ncol(&graph, edgelist_file, NULL, 1, IGRAPH_ADD_WEIGHTS_IF_PRESENT, IGRAPH_UNDIRECTED);
    if(!igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_EDGE, "weight")) {
        SetIgraphAllEdgesWeight(&graph, 1.0);
    }
    /* igraph_read_graph_edgelist(&graph, edgelist_file, 0, false); */
    fclose(edgelist_file);
    printf("my_rank: %d finished loading graph\n", my_rank);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info);
    /* std::cerr << EAN(&graph, "weight", 0) << std::endl; */

    int graph_vcount = igraph_vcount(&graph);
    int my_work = graph_vcount/nprocs;
    if (graph_vcount % nprocs != 0)
        my_work ++;
    int start_vertex = my_rank*my_work;
    int end_vertex = (my_rank+1)*my_work;
    /* int before_mincut_number_of_clusters = -1; */
    int after_mincut_number_of_clusters = -2;
    int iter_count = 0;
    igraph_vector_int_t rice_vec;
    igraph_es_t rice_agg_es;
    igraph_vector_int_t rice_agg_vec;
    igraph_integer_t rice_agg_size;

    int rice_size;
    int rice_size_arr[nprocs];
    int rice_displacements[nprocs];
    
    if(this->existing_clustering == "") {
        /* int seed = uni(rng); */
        // not working
        printf("my_rank: %d no clustering data\n", my_rank);
        int seed = 0;
        std::map<int, int> node_id_to_cluster_id_map = ConstrainedClustering::GetCommunities("", this->algorithm, seed, this->clustering_parameter, &graph);
        // output_map(node_id_to_cluster_id_map);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_map);
    } else if(this->existing_clustering != "") {
        // non mpi
        std::map<std::string, int> original_to_new_id_map = ConstrainedClustering::GetOriginalToNewIdMap(&graph);
        std::map<int, int> new_node_id_to_cluster_id_map = ConstrainedClustering::ReadCommunities(original_to_new_id_map, this->existing_clustering);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, new_node_id_to_cluster_id_map);


        // MPI
        // printf("my_rank: %d clustering data\n", my_rank);
        // std::map<std::string, int> original_to_new_id_map = ConstrainedClustering::GetOriginalToNewIdMapDistributed(&graph, start_vertex, end_vertex);
        // printf("my_rank: %d read orig to new id map\n", my_rank);
        // // printf("my_rank: %d clustering file %s\n",this -> existing_clustering);
        // ConstrainedClustering::initializeSlice(&graph, original_to_new_id_map);
        // // output_map(original_to_new_id_map);
        
        // std::map<int, int> new_node_id_to_cluster_id_map = ConstrainedClustering::ReadCommunitiesDistributed(original_to_new_id_map, this->existing_clustering);
        // printf("my_rank: %d read communities\n", my_rank);
        // // output_map(new_node_id_to_cluster_id_map);
        // // get edges to delete
        // // printf("my_rank: %d enter rice\n" ,  my_rank);
        // // ConstrainedClustering::RemoveInterClusterEdges(&graph, new_node_id_to_cluster_id_map);
        // printf("my_rank: %d enter rice_distributed\n", my_rank);
        // rice_vec = ConstrainedClustering::RemoveInterClusterEdgesDistributed(&graph, new_node_id_to_cluster_id_map);
        // igraph_vector_int_print(&rice_vec);
        // // aggregates edges to delete`
        // rice_size = igraph_vector_int_size(&rice_vec);
        // printf("my_rank: %d rice_size: %d\n", my_rank, rice_size);
        // MPI_Allgather(&rice_size, 1, MPI_INT, rice_size_arr, 1, MPI_INT, MPI_COMM_WORLD);
        // // output_arr(rice_size_arr, nprocs);
        // build_displacements(rice_displacements, rice_size_arr, nprocs);
        // // output_arr(rice_displacements, nprocs);
        // int rice_arr[rice_size];
        // rice_agg_size = rice_displacements[nprocs-1]+rice_size;
        // int rice_agg[rice_agg_size];
        // convert_vec_to_arr(rice_vec, rice_arr, rice_size);
        // MPI_Allgatherv(rice_arr, rice_size, MPI_INT, rice_agg, rice_size_arr, rice_displacements, MPI_INT, MPI_COMM_WORLD );
        // output_arr(rice_agg, rice_agg_size);
        // // delete edges from graph
        // igraph_vector_int_init(&rice_agg_vec, rice_agg_size);
        // convert_arr_to_vec(rice_agg_vec, rice_agg, rice_agg_size);
        // igraph_es_vector_copy(&rice_agg_es, &rice_agg_vec);
        // igraph_delete_edges(&graph, rice_agg_es);
    }
    

    /** SECTION Get Connected Components START **/
    std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    // store the results into the queue that each thread pulls from
    int cc_count = connected_components_vector.size();
    int cc_my_work = cc_count/nprocs;
    if (cc_count%nprocs)
        cc_count++;
    int cc_start = my_rank*cc_my_work;
    int cc_end = (my_rank+1)*cc_my_work;
    // non mpi
    for(size_t i = 0; i < connected_components_vector.size(); i ++) {
        CM::to_be_mincut_clusters.push(connected_components_vector[i]);
    }
    // igraph_t subgraph;
    // igraph_vector_int_t nodes_to_keep;
    // igraph_vector_int_t new_id_to_old_id_vector_map;
    // igraph_vector_int_init(&nodes_to_keep, my_work);
    // for(size_t i = cc_start; i < cc_end; i ++) {
    //     for(int j = 0; j < connected_components_vector[i].size();j++)
    //     VECTOR(nodes_to_keep)[i] = connected_components_vector[i][j];
    // }
    // igraph_induced_subgraph_map(graph, &subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_vector_map);
    // for(size_t i = cc_start; i < cc_end; i ++) {
    //     CM::to_be_mincut_clusters.push(connected_components_vector[i]);
    // }
    /** SECTION Get Connected Components END **/
    int previous_done_being_clustered_size = 0;
    while (true) {
        this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::debug);
        if(iter_count % 10 == 0) {
            this->WriteToLogFile("Iteration number: " + std::to_string(iter_count), Log::info);
        }

        /** SECTION MinCutOnceAndCluster Each Connected Component START **/
        this->WriteToLogFile(std::to_string(CM::to_be_mincut_clusters.size()) + " [connected components / clusters] to be mincut", Log::debug);
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
        this->WriteToLogFile(std::to_string(CM::to_be_clustered_clusters.size()) + " [connected components / clusters] to be clustered after a round of mincuts", Log::debug);
        this->WriteToLogFile(std::to_string(CM::done_being_clustered_clusters.size() - previous_done_being_clustered_size) + " [connected components / clusters] were found to be well connected", Log::debug);
        previous_done_being_clustered_size = CM::done_being_clustered_clusters.size();
        /** SECTION MinCutOnceAndCluster Each Connected Component END **/

        /** SECTION Check If All Clusters Are Well-Connected START **/
        after_mincut_number_of_clusters = CM::to_be_clustered_clusters.size();
        if(after_mincut_number_of_clusters == 0) {
            this->WriteToLogFile("all clusters are well-connected", Log::info);
            this->WriteToLogFile("Total number of iterations: " + std::to_string(iter_count + 1), Log::info);
            break;
        }
        /** SECTION Check If All Clusters Are Well-Connected END **/

        while(!CM::to_be_clustered_clusters.empty()) {
            CM::to_be_mincut_clusters.push(CM::to_be_clustered_clusters.front());
            CM::to_be_clustered_clusters.pop();
        }

        iter_count ++;
    }


    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteClusterQueue(CM::done_being_clustered_clusters, &graph);
    igraph_destroy(&graph);
    return 0;
}
