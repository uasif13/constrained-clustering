#include "cm.h"
#include "mpi.h"

using namespace std;
#define ROOT 0

int * build_displacements(int* size_array, int nprocs) {
    int * displacements;
    displacements[nprocs];
    displacements[0] = 0;
    for (int i = 1; i < nprocs; i ++) {
        displacements[i] = displacements[i-1] + size_array[i];
    }
    return displacements;
}

int CM::main(int my_rank, int nprocs) {
    /* std::random_device rd; */
    /* std::mt19937 rng{rd()}; */
    /* std::uniform_int_distribution<int> uni(0, 100000); */
    this->WriteToLogFile("Loading the initial graph" , Log::info);
    FILE* edgelist_file = fopen(this->edgelist.c_str(), "r");
    igraph_t graph;
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_ncol(&graph, edgelist_file, NULL, 1, IGRAPH_ADD_WEIGHTS_IF_PRESENT, IGRAPH_UNDIRECTED);
    if(!igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_EDGE, "weight")) {
        SetIgraphAllEdgesWeight(&graph, 1.0);
    }
    /* igraph_read_graph_edgelist(&graph, edgelist_file, 0, false); */
    fclose(edgelist_file);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info);
    /* std::cerr << EAN(&graph, "weight", 0) << std::endl; */


    /* int before_mincut_number_of_clusters = -1; */
    int after_mincut_number_of_clusters = -2;
    int iter_count = 0;
    ConstrainedClustering::initializeSlice(&graph, my_rank, nprocs);
    igraph_vector_int_t rice_vec;
    int * rice;
    int * rice_agg;
    int rice_agg_size;

    int rice_size;
    int * rice_size_arr;

    int * rice_displacements;

    rice_size_arr[nprocs];

    if(this->existing_clustering == "") {
        /* int seed = uni(rng); */
        // not working
        int seed = 0;
        std::map<int, int> node_id_to_cluster_id_map = ConstrainedClustering::GetCommunities("", this->algorithm, seed, this->clustering_parameter, &graph);
        rice = ConstrainedClustering::RemoveInterClusterEdgesArray(&graph, node_id_to_cluster_id_map);
    } else if(this->existing_clustering != "") {
        std::map<std::string, int> original_to_new_id_map = ConstrainedClustering::GetOriginalToNewIdMap(&graph);
        std::map<int, int> new_node_id_to_cluster_id_map = ConstrainedClustering::ReadCommunities(original_to_new_id_map, this->existing_clustering);
        rice = ConstrainedClustering::RemoveInterClusterEdgesArray(&graph, new_node_id_to_cluster_id_map);
    }
    if (my_rank == ROOT) {
        rice_size = igraph_vector_int_size(&rice);
        MPI_Gather(&rice_size, 1, MPI_INT, rice_size_arr, ROOT, MPI_COMM_WORLD);
        rice_displacements = build_displacements(rice_size_arr, nprocs);
        MPI_Gatherv(&rice, rice_size, MPI_INT, rice_agg, &rice_size, rice_displacements, MPI_INT, ROOT, MPI_COMM_WORLD );
        rice_agg_size = nprocs-1+rice_size;
    }
    MPI_Bcast(&rice_agg_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(rice_agg, rice_agg_size, MPI_INT, ROOT, MPI_COMM_WORLD);
    printf("my_rank : %d rice_agg_size: %d \n", my_rank, rice_agg_size);

    /** SECTION Get Connected Components START **/
    std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    // store the results into the queue that each thread pulls from
    for(size_t i = 0; i < connected_components_vector.size(); i ++) {
        CM::to_be_mincut_clusters.push(connected_components_vector[i]);
    }
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
