#ifndef MINCUT_ONLY_H
#define MINCUT_ONLY_H
#include "constrained.h"
#include "mincut_recursive.h"



class MincutOnly : public ConstrainedClustering
{
public:
    MincutOnly(std::string edgelist, std::string existing_clustering, int num_processors, std::string output_file, std::string log_file, std::string connectedness_criterion, int log_level, int my_rank, int nprocs, int thread_coarsening) : ConstrainedClustering(edgelist, "", -1, existing_clustering, num_processors, output_file, log_file, log_level, my_rank, nprocs), connectedness_criterion(connectedness_criterion), thread_coarsening(thread_coarsening) {
                                                                                                                                                                                                                                               };
    int main(int my_rank, int nprocs, uint64_t *opCount) override;

    static inline std::vector<std::vector<int>> GetConnectedComponentsOnPartition(const igraph_t *graph, std::vector<int> &partition)
    {
        std::vector<std::vector<int>> cluster_vectors;
        igraph_vector_int_t nodes_to_keep;
        igraph_vector_int_t new_id_to_old_id_map;
        igraph_vector_int_init(&new_id_to_old_id_map, partition.size());
        igraph_vector_int_init(&nodes_to_keep, partition.size());
        for (size_t i = 0; i < partition.size(); i++)
        {
            VECTOR(nodes_to_keep)
            [i] = partition[i];
        }
        igraph_t induced_subgraph;
        igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
        std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&induced_subgraph);
        for (size_t i = 0; i < connected_components_vector.size(); i++)
        {
            std::vector<int> translated_cluster_vector;
            for (size_t j = 0; j < connected_components_vector.at(i).size(); j++)
            {
                int new_id = connected_components_vector.at(i).at(j);
                translated_cluster_vector.push_back(VECTOR(new_id_to_old_id_map)[new_id]);
            }
            cluster_vectors.push_back(translated_cluster_vector);
        }
        igraph_vector_int_destroy(&nodes_to_keep);
        igraph_vector_int_destroy(&new_id_to_old_id_map);
        igraph_destroy(&induced_subgraph);
        return cluster_vectors;
    }

    static inline void ComputeMinCutRecursive(const vector<int> current_cluster, const igraph_t *original_graph, ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log, int depth)
    {
        // if (current_cluster.size() == 1 || current_cluster[0] == -1)
        // {
        //     // done with work!
        //     /* std::cerr << "thread done" << std::endl; */
        //     return;
        // }
        igraph_vector_int_t nodes_to_keep;
        igraph_vector_int_t new_id_to_old_id_map;
        igraph_vector_int_init(&new_id_to_old_id_map, current_cluster.size());
        igraph_vector_int_init(&nodes_to_keep, current_cluster.size());
        for (size_t i = 0; i < current_cluster.size(); i++)
        {
            /* std::cerr << "node to keep[i]=" << std::to_string(current_cluster.at(i)) << std::endl; */
            VECTOR(nodes_to_keep)
            [i] = current_cluster[i];
        }
        igraph_t induced_subgraph;
        // technically could just pass in the nodes and edges info directly by iterating through the edges and checking if it's inter vs intracluster
        // likely not too much of a memory or time overhead
        // if (2*current_cluster.size() < igraph_vcount(graph)) {
        igraph_induced_subgraph_map(original_graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
        // }
        // else {
        //     igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_COPY_AND_DELETE, NULL, &new_id_to_old_id_map);

        // }
        /* std::cerr << "induced subgraph" << std::endl; */

        igraph_vector_int_destroy(&nodes_to_keep);

        // printf("depth: %d induced subgraph no_of_vertices: %d\n", depth, igraph_vcount(&induced_subgraph));
        MinCutCustom mcc(&induced_subgraph);
        int edge_cut_size = mcc.ComputeMinCut();
        std::vector<int> in_partition = mcc.GetInPartition();
        std::vector<int> out_partition = mcc.GetOutPartition();
        // printf("mincut results edge_cut_size: %d in_partition_size: %d out_partition_size: %d\n", edge_cut_size, in_partition.size(), out_partition.size());
        bool current_criterion = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition.size(), out_partition.size(), edge_cut_size);

        // printf("current_criterion: %d\n", current_criterion);

        // base case
        if (current_criterion)
        {
            igraph_destroy(&induced_subgraph);
            igraph_vector_int_destroy(&new_id_to_old_id_map);
            // printf("cluster is well connected size: %d\n", current_cluster.size());
            // MincutOnly::done_being_mincut_clusters.push(current_cluster);

            // std::unique_lock<std::mutex> lk(cv_m);
            // std::cerr << "Waiting... \n";
            // cv.wait(lk, []{ return i == 1; });
            // std::cerr << "...finished waiting. i == 1\n";
           
            std::unique_lock<std::mutex> done_being_mincut_lock{MincutOnly::done_being_mincut_mutex};
            MincutOnly::done_being_mincut_clusters.push(current_cluster);

            done_being_mincut_lock.unlock();

        }
        else
        {
            // printf("cluster is not well connected depth %d\n",depth);
            if (in_partition.size() > 1)
            {
                // printf("in_partition depth %d\n",depth);

                std::vector<std::vector<int>> in_clusters = GetConnectedComponentsOnPartition(&induced_subgraph, in_partition);
                // igraph_destroy(&induced_subgraph);
                for (size_t i = 0; i < in_clusters.size(); i++)
                {
                    std::vector<int> translated_in_clusters;
                    for (size_t j = 0; j < in_clusters[i].size(); j++)
                    {
                        translated_in_clusters.push_back(VECTOR(new_id_to_old_id_map)[in_clusters[i][j]]);
                    }
                    MincutOnly::ComputeMinCutRecursive(translated_in_clusters, original_graph, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, depth+1);
                }
                // igraph_vector_int_destroy(&new_id_to_old_id_map);
            }
            if (out_partition.size() > 1)
            {
                // printf("out_partition depth %d\n",depth);

                std::vector<std::vector<int>> out_clusters = GetConnectedComponentsOnPartition(&induced_subgraph, out_partition);
                // igraph_destroy(&induced_subgraph);
                for (size_t i = 0; i < out_clusters.size(); i++)
                {
                    std::vector<int> translated_out_clusters;
                    for (size_t j = 0; j < out_clusters[i].size(); j++)
                    {
                        translated_out_clusters.push_back(VECTOR(new_id_to_old_id_map)[out_clusters[i][j]]);
                    }
                    MincutOnly::ComputeMinCutRecursive(translated_out_clusters, original_graph, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, depth+1);

                }
                // igraph_vector_int_destroy(&new_id_to_old_id_map);
            }
            // CRITICAL FIX: Destroy resources to prevent memory leak
            igraph_destroy(&induced_subgraph);
            igraph_vector_int_destroy(&new_id_to_old_id_map);
        }
    }

    static inline void MinCutWorker(igraph_t *graph, ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log, int thread_coarsening)
    {
        while (true)
        {
            std::unique_lock<std::mutex> to_be_mincut_lock{to_be_mincut_mutex};
            to_be_mincut_condition_variable.wait(to_be_mincut_lock, []()
                                                 { return !MincutOnly::to_be_mincut_clusters.empty(); });

            std::vector<std::vector<int>> current_cluster_queue;
            for (int i = 0; i < thread_coarsening; i++)
            {
                if (MincutOnly::to_be_mincut_clusters.size() > 0)
                {
                    current_cluster_queue.push_back(MincutOnly::to_be_mincut_clusters.front());
                    MincutOnly::to_be_mincut_clusters.pop();
                }
            }
            // current_cluster_queue.push_back(MincutOnly::to_be_mincut_clusters.front());
            // MincutOnly::to_be_mincut_clusters.pop();
            to_be_mincut_lock.unlock();
            printf("processing cluster_queue of size: %d\n", current_cluster_queue.size());

            /* std::cerr << "processing cluster of size:" << std::to_string(current_cluster.size()) << std::endl; */
            /* std::cerr << "current cluster size: " << current_cluster.size() << std::endl; */
            for (int i = 0; i < current_cluster_queue.size(); i++)
            {
                vector<int> current_cluster = current_cluster_queue[i];

                if (current_cluster.size() == 1 || current_cluster[0] == -1)
                {
                    // done with work!
                    /* std::cerr << "thread done" << std::endl; */
                    return;
                }
                ComputeMinCutRecursive(current_cluster, graph, current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log,1);
                printf("done with mincut recursive \n");

                // std::thread::id this_id = std::this_thread::get_id();
                // printf("thread id %d current_cluster_size: %d\n", this_id, current_cluster.size());
                
                // Mincut Recursive
                // MinCutRecursive mcr(graph, current_connectedness_criterion, connectedness_criterion_c,connectedness_criterion_x, pre_computed_log);
                

                // Mincut Nonrecursive

                // igraph_vector_int_t nodes_to_keep;
                // igraph_vector_int_t new_id_to_old_id_map;
                // igraph_vector_int_init(&new_id_to_old_id_map, current_cluster.size());
                // igraph_vector_int_init(&nodes_to_keep, current_cluster.size());
                // for(size_t i = 0; i < current_cluster.size(); i ++) {
                //     /* std::cerr << "node to keep[i]=" << std::to_string(current_cluster.at(i)) << std::endl; */
                //     VECTOR(nodes_to_keep)[i] = current_cluster[i];
                // }
                // igraph_t induced_subgraph;
                // // technically could just pass in the nodes and edges info directly by iterating through the edges and checking if it's inter vs intracluster
                // // likely not too much of a memory or time overhead
                // igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
                // /* std::cerr << "induced subgraph" << std::endl; */

                // /* igraph_vector_int_t node_degrees; */
                // /* igraph_vector_int_init(&node_degrees, 0); */
                // /* igraph_degree(induced_subgraph, node_degrees, igraph_vss_vector(&nodes_to_keep), IGRAPH_ALL, IGRAPH_NO_LOOPS); */
                // /* for(size_t i = 0; i < ; i ++) { */
                // /*     VECTOR(nodes_to_keep)[i] = current_cluster[i]; */
                // /* } */
                // /* while(!ConstrainedClustering::IsWellConnected */

                // /* SetIgraphAllEdgesWeight(&induced_subgraph, 1.0); */
                // MinCutCustom mcc(&induced_subgraph);
                // int edge_cut_size = mcc.ComputeMinCut();
                // // mcc.ComputeMinCutRecursive(&induced_subgraph, current_cluster, igraph_vss_all(), igraph_ess_all(IGRAPH_EDGEORDER_ID), &new_id_to_old_id_map);
                // std::vector<int> in_partition = mcc.GetInPartition();
                // std::vector<int> out_partition = mcc.GetOutPartition();
                // /* std::cerr << "got the cuts into " << std::to_string(in_partition.size()) << " and " << std::to_string(out_partition.size()) << std::endl; */
                // bool current_criterion = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition.size(), out_partition.size(), edge_cut_size);
                // /* if(connectedness_criterion == ConnectednessCriterion::Simple) { */
                // /*     current_criterion = ConstrainedClustering::IsConnected(edge_cut_size); */
                // /* } else if(connectedness_criterion == ConnectednessCriterion::Well) { */
                // /*     current_criterion = ConstrainedClustering::IsWellConnected(in_partition, out_partition, edge_cut_size, &induced_subgraph); */
                // /* } */

                // if(!current_criterion) {
                //     /* for(size_t i = 0; i < in_partition.size(); i ++) { */
                //     /*     in_partition[i] = VECTOR(new_id_to_old_id_map)[in_partition[i]]; */
                //     /* } */
                //     /* for(size_t i = 0; i < out_partition.size(); i ++) { */
                //     /*     out_partition[i] = VECTOR(new_id_to_old_id_map)[out_partition[i]]; */
                //     /* } */
                //     printf("cluster is not well connected\n");
                //     if(in_partition.size() > 1) {
                //         std::vector<std::vector<int>> in_clusters = GetConnectedComponentsOnPartition(&induced_subgraph, in_partition);
                //         for(size_t i = 0; i < in_clusters.size(); i ++) {
                //             std::vector<int> translated_in_clusters;
                //             for(size_t j = 0; j < in_clusters[i].size(); j ++) {
                //                 translated_in_clusters.push_back(VECTOR(new_id_to_old_id_map)[in_clusters[i][j]]);
                //             }
                //             {
                //                 std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnly::to_be_mincut_mutex);
                //                 MincutOnly::to_be_mincut_clusters.push(translated_in_clusters);
                //             }
                //         }
                //     }
                //     if(out_partition.size() > 1) {
                //         std::vector<std::vector<int>> out_clusters = GetConnectedComponentsOnPartition(&induced_subgraph, out_partition);
                //         for(size_t i = 0; i < out_clusters.size(); i ++) {
                //             std::vector<int> translated_out_clusters;
                //             for(size_t j = 0; j < out_clusters[i].size(); j ++) {
                //                 translated_out_clusters.push_back(VECTOR(new_id_to_old_id_map)[out_clusters[i][j]]);
                //             }
                //             {
                //                 std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnly::to_be_mincut_mutex);
                //                 MincutOnly::to_be_mincut_clusters.push(translated_out_clusters);
                //             }
                //         }
                //     }
                // } else {
                //     std::unique_lock<std::mutex> done_being_mincut_lock{MincutOnly::done_being_mincut_mutex};
                //     MincutOnly::done_being_mincut_clusters.push(current_cluster);
                //     done_being_mincut_lock.unlock();
                // }
                // igraph_vector_int_destroy(&nodes_to_keep);
                // igraph_vector_int_destroy(&new_id_to_old_id_map);
                // igraph_destroy(&induced_subgraph);
            }
        }
    }

private:
    static inline std::mutex to_be_mincut_mutex;
    static inline std::condition_variable to_be_mincut_condition_variable;
    static inline std::queue<std::vector<int>> to_be_mincut_clusters;
    static inline std::mutex done_being_mincut_mutex;
    static inline std::queue<std::vector<int>> done_being_mincut_clusters;
    std::string connectedness_criterion;
    int thread_coarsening;
    static inline std::condition_variable cv;
    static inline std::mutex cv_m; 
    static inline int waiting = 0;
};

#endif
