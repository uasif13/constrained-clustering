#include "mincut_recursive.h"
#include "mincut_only.h"





const std::vector<int> &MinCutRecursive::GetInPartition() const
{
    return this->in_partition;
}
const std::vector<int> &MinCutRecursive::GetOutPartition() const
{
    return this->out_partition;
}
const std::vector<std::vector<int>> &MinCutRecursive::GetDoneBeingClusteredClusters() 
{
    return MinCutRecursive::done_being_clustered_clusters;
}
const int &MinCutRecursive::SetDoneBeingClusteredClusters() 
{
    MinCutRecursive::done_being_clustered_clusters = {};
    return 1;
}