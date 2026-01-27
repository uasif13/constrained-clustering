#ifndef MINCUT_ONLY_PREPROCESS_H
#define MINCUT_ONLY_PREPROCESS_H
#include "constrained.h"



class MincutOnlyPreProcess : public ConstrainedClustering {
    public:
        MincutOnlyPreProcess(std::string edgelist, int num_processors, std::string output_file, std::string log_file, std::string connectedness_criterion, int log_level) : ConstrainedClustering(edgelist, "", -1, "", num_processors, output_file, log_file, log_level), connectedness_criterion(connectedness_criterion) {
        };
        int main() override;

        virtual ~MincutOnlyPreProcess() {
        }

        static inline std::vector<std::vector<long>> GetConnectedComponentsOnPartition(const igraph_t* graph, std::vector<long>& partition, std::unordered_map<long, long> prev_new_id_to_old_id_map, int worker_id) {
            
            std::vector<std::vector<long>> cluster_vectors;
            igraph_vector_int_t nodes_to_keep;
            igraph_vector_int_t new_id_to_old_id_map;
            
            igraph_vector_int_init(&new_id_to_old_id_map, partition.size());
            igraph_vector_int_init(&nodes_to_keep, partition.size());
            for(size_t i = 0; i < partition.size(); i ++) {
                VECTOR(nodes_to_keep)[i] = partition[i];
            }
            
            igraph_t induced_subgraph;
            igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
            
            std::vector<std::vector<long>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&induced_subgraph);
            
            for (size_t i = 0; i < connected_components_vector.size(); i++) {
                
                std::vector<long> translated_cluster_vector;

                igraph_vector_int_t sub_nodes_to_keep;
                igraph_vector_int_t sub_new_id_to_old_id_vector_map;
                
                igraph_vector_int_init(&sub_nodes_to_keep, connected_components_vector[i].size());
                for(size_t j = 0; j < connected_components_vector[i].size(); j ++) {
                    VECTOR(sub_nodes_to_keep)[j] = connected_components_vector[i][j];
                }
                
                igraph_t sub_subgraph;
                igraph_vector_int_init(&sub_new_id_to_old_id_vector_map, igraph_vector_int_size(&sub_nodes_to_keep));
                igraph_induced_subgraph_map(&induced_subgraph, &sub_subgraph, igraph_vss_vector(&sub_nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &sub_new_id_to_old_id_vector_map);
                
                igraph_eit_t eit;
                igraph_eit_create(&sub_subgraph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
                
                int edge_count = 0;
                for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                    igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                    long from_node = prev_new_id_to_old_id_map[VECTOR(new_id_to_old_id_map)[VECTOR(sub_new_id_to_old_id_vector_map)[IGRAPH_FROM(&sub_subgraph, current_edge)]]];
                    long to_node = prev_new_id_to_old_id_map[VECTOR(new_id_to_old_id_map)[VECTOR(sub_new_id_to_old_id_vector_map)[IGRAPH_TO(&sub_subgraph, current_edge)]]];
                    
                    translated_cluster_vector.push_back(from_node);   
                    translated_cluster_vector.push_back(to_node);
                    edge_count++;
                }
                igraph_eit_destroy(&eit);
                
                igraph_vector_int_destroy(&sub_nodes_to_keep);
                igraph_vector_int_destroy(&sub_new_id_to_old_id_vector_map);
                igraph_destroy(&sub_subgraph);
                
                cluster_vectors.push_back(translated_cluster_vector);
            }
            
            igraph_vector_int_destroy(&nodes_to_keep);
            igraph_vector_int_destroy(&new_id_to_old_id_map);
            igraph_destroy(&induced_subgraph);
            
            return cluster_vectors;
        }

        static inline void MinCutWorker(ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log) {
            // Assign worker ID
            static std::atomic<int> worker_id_counter{0};
            int worker_id = worker_id_counter++;
            
            int clusters_processed = 0;
            
            while (true) {
                std::unique_lock<std::mutex> to_be_mincut_lock{to_be_mincut_mutex};
                
                to_be_mincut_condition_variable.wait(to_be_mincut_lock, []() {
                    return !MincutOnlyPreProcess::to_be_mincut_clusters.empty();
                });
                
                std::vector<long> current_cluster = MincutOnlyPreProcess::to_be_mincut_clusters.front();
                MincutOnlyPreProcess::to_be_mincut_clusters.pop();
                int remaining_in_queue = MincutOnlyPreProcess::to_be_mincut_clusters.size();
                to_be_mincut_lock.unlock();
                
                if(current_cluster.size() == 1 && current_cluster[0] == -1) {
                    return;
                }
                
                clusters_processed++;
                
                // Build ID mappings
                std::unordered_map<long, long> new_id_to_old_id_map;
                std::unordered_map<long, long> old_id_to_new_id_map;

                igraph_vector_int_t edges;
                igraph_vector_int_init(&edges, current_cluster.size());
                std::set<long> current_cluster_set;

                long next_node_id = 0;
                for (size_t i = 0; i < current_cluster.size(); i++) {
                    if (old_id_to_new_id_map.find(current_cluster[i]) == old_id_to_new_id_map.end()) {
                        current_cluster_set.insert(next_node_id);
                        old_id_to_new_id_map[current_cluster[i]] = next_node_id;
                        new_id_to_old_id_map[next_node_id] = current_cluster[i];
                        next_node_id++;
                    }
                    VECTOR(edges)[i] = old_id_to_new_id_map[current_cluster[i]];
                }
                
                igraph_t subgraph;
                igraph_empty(&subgraph, next_node_id, IGRAPH_UNDIRECTED);
                
                igraph_add_edges(&subgraph, &edges, NULL);

                MinCutCustom mcc(&subgraph);
                long edge_cut_size = mcc.ComputeMinCut();
                
                std::vector<long> in_partition = mcc.GetInPartition();
                std::vector<long> out_partition = mcc.GetOutPartition();
                
                bool current_criterion = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition.size(), out_partition.size(), edge_cut_size);

                if(!current_criterion) {
                    
                    int new_clusters_added = 0;
                    
                    if(in_partition.size() > 1) {
                        std::vector<std::vector<long>> in_clusters = GetConnectedComponentsOnPartition(&subgraph, in_partition, new_id_to_old_id_map, worker_id);
                        
                        for(size_t i = 0; i < in_clusters.size(); i ++) {
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnlyPreProcess::to_be_mincut_mutex);
                                MincutOnlyPreProcess::to_be_mincut_clusters.push(in_clusters[i]);
                                new_clusters_added++;
                            }
                        }
                        // Notify other workers
                        to_be_mincut_condition_variable.notify_all();
                    }
                    
                    if(out_partition.size() > 1) {
                        std::vector<std::vector<long>> out_clusters = GetConnectedComponentsOnPartition(&subgraph, out_partition, new_id_to_old_id_map, worker_id);
                        
                        for(size_t i = 0; i < out_clusters.size(); i ++) {
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnlyPreProcess::to_be_mincut_mutex);
                                MincutOnlyPreProcess::to_be_mincut_clusters.push(out_clusters[i]);
                                new_clusters_added++;
                            }
                        }
                        // Notify other workers
                        to_be_mincut_condition_variable.notify_all();
                    }
                    
                } else {
                    
                    std::vector<long> cc;
                    for (int i = 0; i < igraph_vcount(&subgraph); i++) {
                        cc.push_back(new_id_to_old_id_map[i]);
                    }
                    
                    std::unique_lock<std::mutex> done_being_mincut_lock{MincutOnlyPreProcess::done_being_mincut_mutex};
                    MincutOnlyPreProcess::done_being_mincut_clusters.push(cc);
                    int total_done = MincutOnlyPreProcess::done_being_mincut_clusters.size();
                    done_being_mincut_lock.unlock();
                    
                }
                
                igraph_vector_int_destroy(&edges);
                igraph_destroy(&subgraph);
            }
        }
        
    private:
        static inline std::mutex to_be_mincut_mutex;
        static inline std::condition_variable to_be_mincut_condition_variable;
        static inline std::queue<std::vector<long>> to_be_mincut_clusters;
        static inline std::mutex done_being_mincut_mutex;
        static inline std::queue<std::vector<long>> done_being_mincut_clusters;
        std::string connectedness_criterion;
};

#endif