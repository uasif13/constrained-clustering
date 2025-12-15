#ifndef MINCUT_H
#define MINCUT_H

#include <algorithms/global_mincut/cactus/cactus_mincut.h>
#include <igraph/igraph.h>


class MinCutCustom {
    public:
        MinCutCustom(const igraph_t* graph) : graph(graph) {
        };
        int ComputeMinCut();
        /* int ComputeAllMinCuts(); */
        const std::vector<long long>& GetInPartition() const;
        const std::vector<long long>& GetOutPartition() const;
        /* const std::vector<std::vector<int>>& GetAllPartitions() const; */
    private:
        std::vector<long long> in_partition;
        std::vector<long long> out_partition;
        /* std::vector<std::vector<int>> all_partitions; */
        std::map<long long, long long> new_to_old_node_id_map;
        const igraph_t* graph;
};

#endif
