#ifndef MINCUT_RECURSIVE_H
#define MINCUT_RECURSIVE_H

#include <algorithms/global_mincut/cactus/cactus_mincut.h>
#include <igraph/igraph.h>
#include "constrained.h"
#include "mincut_only.h"

class MinCutRecursive
{
public:
    MinCutRecursive(
        const igraph_t *original_graph,
        ConnectednessCriterion current_connectedness_criterion,
        double connectedness_criterion_c, double connectedness_criterion_x,
        double pre_computed_log) : original_graph(original_graph), current_connectedness_criterion(current_connectedness_criterion),
                                   connectedness_criterion_c(connectedness_criterion_c), connectedness_criterion_x(connectedness_criterion_x),
                                   pre_computed_log(pre_computed_log) {

                                   };
    static inline std::vector<std::vector<int>> GetConnectedComponentsOnPartition(const igraph_t* graph, std::vector<int>& partition) {
            std::vector<std::vector<int>> cluster_vectors;
            igraph_vector_int_t nodes_to_keep;
            igraph_vector_int_t new_id_to_old_id_map;
            igraph_vector_int_init(&new_id_to_old_id_map, partition.size());
            igraph_vector_int_init(&nodes_to_keep, partition.size());
            for(size_t i = 0; i < partition.size(); i ++) {
                VECTOR(nodes_to_keep)[i] = partition[i];
            }
            igraph_t induced_subgraph;
            igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
            std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&induced_subgraph);
            for(size_t i = 0; i < connected_components_vector.size(); i ++) {
                std::vector<int> translated_cluster_vector;
                for(size_t j = 0; j < connected_components_vector.at(i).size(); j ++) {
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

    /* int ComputeAllMinCuts(); */
    const std::vector<int> &GetInPartition() const;
    const std::vector<int> &GetOutPartition() const;
    const static inline std::vector<std::vector<int>> &GetDoneBeingClusteredClusters();
    const static inline int &SetDoneBeingClusteredClusters();

    /* const std::vector<std::vector<int>>& GetAllPartitions() const; */
private:
    std::vector<int> in_partition;
    std::vector<int> out_partition;
    /* std::vector<std::vector<int>> all_partitions; */
    const igraph_t *original_graph;
    ConnectednessCriterion current_connectedness_criterion;
    double connectedness_criterion_c;
    double connectedness_criterion_x;
    double pre_computed_log;
    static inline std::vector<std::vector<int>> done_being_clustered_clusters;
};

#endif
