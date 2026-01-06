#ifndef MINCUT_H
#define MINCUT_H

#include <algorithms/global_mincut/cactus/cactus_mincut.h>
#include <igraph/igraph.h>


class MinCutCustom {
    public:
        MinCutCustom(const igraph_t* graph) : graph(graph) {
        };
        long ComputeMinCut();
        /* int ComputeAllMinCuts(); */
        const std::vector<long>& GetInPartition() const;
        const std::vector<long>& GetOutPartition() const;
        /* const std::vector<std::vector<int>>& GetAllPartitions() const; */
    private:
        std::vector<long> in_partition;
        std::vector<long> out_partition;
        /* std::vector<std::vector<int>> all_partitions; */
        std::map<int, int> new_to_old_node_id_map;
        const igraph_t* graph;
};

#endif
