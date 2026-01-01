#ifndef CMPREPROCESS_H
#define CMPREPROCESS_H
#include "constrained.h"
#include "mmap_subgraph_loader.h"

class CMPreprocess : public ConstrainedClustering {
    public:
        CMPreprocess(std::string subgraph_file, std::string algorithm, double clustering_parameter, int num_processors, std::string output_file, std::string log_file, int log_level) : ConstrainedClustering(subgraph_file, algorithm, clustering_parameter, "", num_processors, output_file, log_file, log_level) {
        };
        int main();

        static inline void output_vec_long(std::vector<std::vector<long>> cc, igraph_t* graph) {
            for (auto const& cluster: cc) {
                for (auto const& item: cluster) {
                    cout << VAS(graph, "name", item) << " ";
                }
                cout << endl;
            }
        }

        static inline std::vector<std::vector<long>> RunClusterOnPartition(const igraph_t* graph, std::unordered_map<long,long> prev_new_id_to_old_id_map, std::string algorithm, int seed, double clustering_parameter, std::vector<long>& partition) {
            std::vector<std::vector<long>> cluster_vectors;
            std::map<long, std::vector<long>> cluster_map;
            std::unordered_map<long, long> new_id_to_old_id_unordered_map;
            igraph_vector_int_t nodes_to_keep;
            igraph_vector_int_t new_id_to_old_id_map;
            igraph_vector_int_init(&new_id_to_old_id_map, partition.size());
            igraph_vector_int_init(&nodes_to_keep, partition.size());
            for(size_t i = 0; i < partition.size(); i ++) {
                VECTOR(nodes_to_keep)[i] = partition[i];
            }
            // printf("runclusteronpartition nodes_to_keep: %d\n", partition.size());
            igraph_t induced_subgraph;
            // igraph_empty(&induced_subgraph, partition.size(), IGRAPH_UNDIRECTED);
            // igraph_vector_int_t edges;
            // vector<long> edges_to_add;
            // igraph_eit_t eit;
            // igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            // long next_node_id = 0;
            // for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            //     igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
            //     long from = IGRAPH_FROM(graph, current_edge);
            //     long to = IGRAPH_TO(graph, current_edge);
            //     if (prev_new_id_to_old_id_map[current_edge] != )
            //     if (new_id_to_old_id_unordered_map.find(from) == new_id_to_old_id_unordered_map.end())
            //     {
            //         new_id_to_old_id_unordered_map[from] = next_node_id;
            //         next_node_id++;
            //     }
            //     if (new_id_to_old_id_unordered_map.find(to) == new_id_to_old_id_unordered_map.end())
            //     {
            //         new_id_to_old_id_unordered_map[to] = next_node_id;
            //         next_node_id++;
            //     }

            // }
            // igraph_add_edges(graph, &edge_vector, NULL);
            igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
            // printf("runclusteronpartition created subgraph\n");

            // igraph_eit_t eit;
            // igraph_eit_create(&induced_subgraph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            // for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
            //     igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
            //     long from = VECTOR(new_id_to_old_id_map)[IGRAPH_FROM(graph, current_edge)];
            //     long to = VECTOR(new_id_to_old_id_map)[IGRAPH_TO(graph, current_edge)];
            //     cout << prev_new_id_to_old_id_map[from] << " " << prev_new_id_to_old_id_map[to];
            // }
            // cout << endl;
            
            std::map<long, long> partition_map = ConstrainedClustering::GetCommunities("", algorithm, seed, clustering_parameter, &induced_subgraph);
            // printf("runclusteronpartition partition map\n");

            ConstrainedClustering::RemoveInterClusterEdges(&induced_subgraph, partition_map);
            // printf("runclusteronpartition remove inter cluster edges\n");

            std::vector<std::vector<long>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&induced_subgraph);
            // printf("runclusteronpartition get connected components %d\n", connected_components_vector.size());
            // // for(size_t i = 0; i < connected_components_vector.size(); i ++) {
            // //     std::vector<long> translated_cluster_vector;
            // //     for(size_t j = 0; j < connected_components_vector.at(i).size(); j ++) {
            // //         long new_id = connected_components_vector.at(i).at(j);
            // //         translated_cluster_vector.push_back(VECTOR(new_id_to_old_id_map)[new_id]);
            // //     }
            // //     cluster_vectors.push_back(translated_cluster_vector);
            // // }
            for (int i = 0; i < connected_components_vector.size(); i++) {
                std::vector<long> translated_cluster_vector;

                igraph_vector_int_t sub_nodes_to_keep;
                igraph_vector_int_t sub_new_id_to_old_id_vector_map;
                // std::unordered_map<long, long> sub_old_id_to_new_id_map;
                // std::unordered_map<long, long> sub_new_id_to_old_id_map;
                igraph_vector_int_init(&sub_nodes_to_keep, connected_components_vector[i].size());
                for(size_t j = 0; j < connected_components_vector[i].size(); j ++) {
                    VECTOR(sub_nodes_to_keep)[j] = connected_components_vector[i][j];
                }
                igraph_t sub_subgraph;
                igraph_vector_int_init(&sub_new_id_to_old_id_vector_map, igraph_vector_int_size(&sub_nodes_to_keep));
                igraph_induced_subgraph_map(&induced_subgraph, &sub_subgraph, igraph_vss_vector(&sub_nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &sub_new_id_to_old_id_vector_map);
                igraph_eit_t eit;
                igraph_eit_create(&sub_subgraph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
                for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                    igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                    long from_node = prev_new_id_to_old_id_map[VECTOR(new_id_to_old_id_map)[VECTOR(sub_new_id_to_old_id_vector_map)[IGRAPH_FROM(&sub_subgraph, current_edge)]]];
                    long to_node = prev_new_id_to_old_id_map[VECTOR(new_id_to_old_id_map)[VECTOR(sub_new_id_to_old_id_vector_map)[IGRAPH_TO(&sub_subgraph, current_edge)]]];
                    // cout << from_node << " " << to_node << " ";
                    // if (current_cluster_set.find(from_node) != current_cluster_set.end() && current_cluster_set.find(to_node) != current_cluster_set.end()) {
                    translated_cluster_vector.push_back(from_node);   
                    translated_cluster_vector.push_back(to_node);  
                    // }
                        //    
      
                }
                cluster_vectors.push_back(translated_cluster_vector);
            }
            // igraph_vector_int_destroy(&nodes_to_keep);
            // igraph_vector_int_destroy(&new_id_to_old_id_map);
            // igraph_destroy(&induced_subgraph);
            return cluster_vectors;
        }

        static inline void MinCutOrClusterWorker(std::string algorithm, int seed, double clustering_parameter) {
	  std::vector<long> current_cluster_edges;
	    std::unique_lock<std::mutex> to_be_mincut_lock{to_be_mincut_mutex};
            if (CMPreprocess::to_be_mincut_clusters.empty()) return;
            else current_cluster_edges = CMPreprocess::to_be_mincut_clusters.front();
            CMPreprocess::to_be_mincut_clusters.pop();
            to_be_mincut_lock.unlock();
            printf("inside Mincutorcluster worker to_be_mincut_cluster_size: %d\n", CMPreprocess::to_be_mincut_clusters.size() );

            if(current_cluster_edges.size() == 1 && current_cluster_edges[0] == -1) {
                // done with work!
                return;
            }            
            std::unordered_map<long, long> new_id_to_old_id_map;
            std::unordered_map<long, long> newnew_id_to_old_id_map;
            std::unordered_map<long, long> old_id_to_new_id_map;

            igraph_vector_int_t edges;
            igraph_vector_int_init(&edges, current_cluster_edges.size());
            std::set<long> current_cluster_set;

            long next_node_id = 0;
            for (int i = 0; i < current_cluster_edges.size(); i++) {
                if (old_id_to_new_id_map.find(current_cluster_edges[i]) == old_id_to_new_id_map.end())
                {
                    current_cluster_set.insert(next_node_id);
                    old_id_to_new_id_map[current_cluster_edges[i]] = next_node_id;
                    new_id_to_old_id_map[next_node_id] = current_cluster_edges[i];
                    // cout << current_cluster_edges[i] << " " << next_node_id << " ";

                    next_node_id++;
                    
                }
                VECTOR(edges)[i] = old_id_to_new_id_map[current_cluster_edges[i]];
                // std::cout << VECTOR(edges)[i] << " ";
            }
            if(current_cluster_set.size() == 1) {
                // done with work!
                return;
            }
            // std::cout << endl;
            igraph_t subgraph;
            igraph_empty(&subgraph, next_node_id, IGRAPH_UNDIRECTED);
            igraph_add_edges(&subgraph, &edges, NULL);          
            std::vector<long> in_partition;
            std::vector<long> out_partition;
            bool is_well_connected = false;
            bool is_non_trivial_cut = false;
            int edge_cut_size = -1;
            // igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_vector_map);
            // printf("before trivial subgraph vcount %d ecount %d\n", igraph_vcount(&subgraph), igraph_ecount(&subgraph));
            // for (;IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
            //     cout << IGRAPH_VIT_GET(vit) << " ";
            // }

            // remove trivial cuts 
            while(!(is_non_trivial_cut || is_well_connected)) {
                MinCutCustom mcc(&subgraph);
                edge_cut_size = mcc.ComputeMinCut();
                std::vector<long> in_partition_candidate = mcc.GetInPartition();
                std::vector<long> out_partition_candidate = mcc.GetOutPartition();
                is_well_connected = ConstrainedClustering::IsWellConnected(in_partition_candidate, out_partition_candidate, edge_cut_size, &subgraph);
                is_non_trivial_cut = in_partition_candidate.size() > 1 && out_partition_candidate.size() > 1;
                if(is_well_connected || is_non_trivial_cut) {
                    in_partition = in_partition_candidate;
                    out_partition = out_partition_candidate;
                } else {
                    long node_to_remove = -1;
                    if(in_partition_candidate.size() == 1) {
                        igraph_vector_int_t newnew_id_to_new_id_map;
                        igraph_vector_int_init(&newnew_id_to_new_id_map, out_partition_candidate.size());
                        node_to_remove = in_partition_candidate.at(0);
                        current_cluster_set.erase(new_id_to_old_id_map[node_to_remove]);
                        igraph_delete_vertices_idx(&subgraph, igraph_vss_1(node_to_remove), NULL, &newnew_id_to_new_id_map);
                        for(long i = 0; i < igraph_vector_int_size(&newnew_id_to_new_id_map); i ++) {
                            long newnew_id = i;
                            long new_id = VECTOR(newnew_id_to_new_id_map)[newnew_id];
                            long old_id = new_id_to_old_id_map[new_id];
                            old_id_to_new_id_map[old_id] = newnew_id;
                            newnew_id_to_old_id_map[newnew_id] = old_id;
                        }
                        new_id_to_old_id_map = newnew_id_to_old_id_map;
                        igraph_vector_int_destroy(&newnew_id_to_new_id_map);
                    }
                    if(out_partition_candidate.size() == 1) {
                        igraph_vector_int_t newnew_id_to_new_id_map;
                        igraph_vector_int_init(&newnew_id_to_new_id_map, out_partition_candidate.size());
                        node_to_remove = out_partition_candidate.at(0);
                        current_cluster_set.erase(new_id_to_old_id_map[node_to_remove]);
                        igraph_delete_vertices_idx(&subgraph, igraph_vss_1(node_to_remove), NULL, &newnew_id_to_new_id_map);
                        for(long i = 0; i < igraph_vector_int_size(&newnew_id_to_new_id_map); i ++) {
                            long newnew_id = i;
                            long new_id = VECTOR(newnew_id_to_new_id_map)[newnew_id];
                            long old_id = new_id_to_old_id_map[new_id];
                            old_id_to_new_id_map[old_id] = newnew_id;
                            newnew_id_to_old_id_map[newnew_id] = old_id;
                        }
                        new_id_to_old_id_map = newnew_id_to_old_id_map;
                        igraph_vector_int_destroy(&newnew_id_to_new_id_map);
                    }
                    if(igraph_vcount(&subgraph) == 0) {
                        // this cluster completely disintegrated somehow
                        // i don't think it'll ever reach 1 before it reaches 2 and does 1 node in both partitions
                        break;
                    }
                }
            }
            // printf("after trivial subgraph vcount %d ecount %d in_partition size %d out_partition_size %d well_connected %d\n", igraph_vcount(&subgraph), igraph_ecount(&subgraph), in_partition.size(), out_partition.size(), is_well_connected);

                if(is_well_connected) {
                    /* std::cerr << "cluster size " << std::to_string(current_cluster.size()) <<  " well connected after performing " << std::to_string(num_nodes_removed) << " trivial mincuts " << std::endl; */
                    // std::vector<long> translated_current_cluster = std::vector<long>(current_cluster_set.begin(), current_cluster_set.end());
                    // for (int i = 0; i < translated_current_cluster.size(); i++) {
                    //     translated_current_cluster[i] = new_id_to_old_id_map[translated_current_cluster[i]];
                    // }
                    // {
                    //     std::lock_guard<std::mutex> done_being_clustered_guard(CMPreprocess::done_being_clustered_mutex);
                    //     CMPreprocess::done_being_clustered_clusters.push(translated_current_cluster);
                    // }
                    std::vector<long> current_cluster;
                    for (int i = 0; i < igraph_vcount(&subgraph); i++) {
                        current_cluster.push_back(new_id_to_old_id_map[i]);
                    }
                    // printf("after trivial current_cluster_size %d \n", current_cluster.size());

                    {
                        std::lock_guard<std::mutex> done_being_clustered_guard(CMPreprocess::done_being_clustered_mutex);
                        CMPreprocess::done_being_clustered_clusters.push(current_cluster);
                    }
                } 
                else if(is_non_trivial_cut) {
                    // printf("after trivial non trivial cut subgraph vcount %d ecount %d\n", igraph_vcount(&subgraph), igraph_ecount(&subgraph));

                    /* std::cerr << "cluster mincut into " << std::to_string(in_partition.size()) << ":" << std::to_string(out_partition.size()) << " after performing " << std::to_string(num_nodes_removed) << " trivial mincuts " << std::endl; */
                    if(in_partition.size() > 1) {
                        std::vector<std::vector<long>> in_clusters = CMPreprocess::RunClusterOnPartition(&subgraph, new_id_to_old_id_map, algorithm, seed, clustering_parameter, in_partition);
                        for(size_t i = 0; i < in_clusters.size(); i ++) {
                            {
                                std::lock_guard<std::mutex> to_be_clustered_guard(CMPreprocess::to_be_clustered_mutex);
                                CMPreprocess::to_be_mincut_clusters.push(in_clusters[i]);
                            }
                        }
                    }
                    if(out_partition.size() > 1) {
                        std::vector<std::vector<long>> out_clusters = CMPreprocess::RunClusterOnPartition(&subgraph, new_id_to_old_id_map, algorithm, seed, clustering_parameter, out_partition);
                        for(size_t i = 0; i < out_clusters.size(); i ++) {
                            {
                                std::lock_guard<std::mutex> to_be_clustered_guard(CMPreprocess::to_be_clustered_mutex);
                                CMPreprocess::to_be_mincut_clusters.push(out_clusters[i]);
                            }
                        }
                    }
                }
            igraph_destroy(&subgraph);
            }
    private:
        static inline std::mutex to_be_mincut_mutex;
        static inline std::condition_variable to_be_mincut_condition_variable;
        static inline std::queue<std::vector<long>> to_be_mincut_clusters;
        static inline std::mutex to_be_clustered_mutex;
        static inline std::queue<std::vector<long>> to_be_clustered_clusters;
        static inline std::mutex done_being_clustered_mutex;
        static inline std::queue<std::vector<long>> done_being_clustered_clusters;
};

#endif
