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

        static inline std::vector<std::vector<long>> GetConnectedComponentsOnPartition(const igraph_t* graph, std::vector<long>& partition, std::unordered_map<long, long> prev_new_id_to_old_id_map) {
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

            igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
            
            std::vector<std::vector<long>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&induced_subgraph);
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
                
                // ========== DIAGNOSTIC: Check for empty vectors ==========
                if(translated_cluster_vector.empty()) {
                    static std::mutex diag_mutex;
                    std::lock_guard<std::mutex> guard(diag_mutex);
                    std::cerr << "[DIAG_EMPTY] GetConnectedComponents returned EMPTY vector! "
                              << "Partition had " << connected_components_vector[i].size() << " nodes, "
                              << "but 0 edges in subgraph" << std::endl;
                    g_empty_vectors_found.fetch_add(1);
                }
                // ========================================================
                
                cluster_vectors.push_back(translated_cluster_vector);
            }
            igraph_vector_int_destroy(&nodes_to_keep);
            igraph_vector_int_destroy(&new_id_to_old_id_map);
            igraph_destroy(&induced_subgraph);
            return cluster_vectors;
        }

        static inline void MinCutWorker(ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log) {
            while (true) {
                std::unique_lock<std::mutex> to_be_mincut_lock{to_be_mincut_mutex};
                to_be_mincut_condition_variable.wait(to_be_mincut_lock, []() {
                    return !MincutOnlyPreProcess::to_be_mincut_clusters.empty();
                });
                std::vector<long> current_cluster = MincutOnlyPreProcess::to_be_mincut_clusters.front();
                MincutOnlyPreProcess::to_be_mincut_clusters.pop();
                to_be_mincut_lock.unlock();
                
                // ========== DIAGNOSTIC: Track dequeue ==========
                g_clusters_dequeued.fetch_add(1);
                
                // DIAGNOSTIC: Check for empty input
                if(current_cluster.empty()) {
                    std::cerr << "[DIAG_ERROR] Dequeued EMPTY cluster!" << std::endl;
                    continue;
                }
                // ================================================
                
                if(current_cluster.size() == 1 || current_cluster[0] == -1) {
                    // done with work!
                    /* std::cerr << "thread done" << std::endl; */
                    return;
                }
                /* std::cerr << "processing cluster of size:" << std::to_string(current_cluster.size()) << std::endl; */
                /* std::cerr << "current cluster size: " << current_cluster.size() << std::endl; */
                std::unordered_map<long, long> new_id_to_old_id_map;
                std::unordered_map<long, long> newnew_id_to_old_id_map;
                std::unordered_map<long, long> old_id_to_new_id_map;

                igraph_vector_int_t edges;
                igraph_vector_int_init(&edges, current_cluster.size());
                std::set<long> current_cluster_set;

                long next_node_id = 0;
                for (int i = 0; i < current_cluster.size(); i++) {
                    if (old_id_to_new_id_map.find(current_cluster[i]) == old_id_to_new_id_map.end())
                    {
                        current_cluster_set.insert(next_node_id);
                        old_id_to_new_id_map[current_cluster[i]] = next_node_id;
                        new_id_to_old_id_map[next_node_id] = current_cluster[i];
                        // cout << current_cluster_edges[i] << " " << next_node_id << " ";

                        next_node_id++;
                        
                    }
                    VECTOR(edges)[i] = old_id_to_new_id_map[current_cluster[i]];
                    // std::cout << VECTOR(edges)[i] << " ";
                }
                
                // ========== DIAGNOSTIC: Log small clusters being processed ==========
                if(next_node_id <= 5) {  // Only log clusters with 5 or fewer nodes
                    static std::mutex diag_mutex;
                    std::lock_guard<std::mutex> guard(diag_mutex);
                    std::cerr << "[DIAG_PROC] Processing cluster: edges=" << current_cluster.size()/2 
                              << " unique_nodes=" << next_node_id << std::endl;
                }
                // ===================================================================
                
                igraph_t subgraph;
                igraph_empty(&subgraph, next_node_id, IGRAPH_UNDIRECTED);
                igraph_add_edges(&subgraph, &edges, NULL); 
                // technically could just pass in the nodes and edges info directly by iterating through the edges and checking if it's inter vs intracluster
                // likely not too much of a memory or time overhead
                /* std::cerr << "induced subgraph" << std::endl; */

                /* igraph_vector_int_t node_degrees; */
                /* igraph_vector_int_init(&node_degrees, 0); */
                /* igraph_degree(induced_subgraph, node_degrees, igraph_vss_vector(&nodes_to_keep), IGRAPH_ALL, IGRAPH_NO_LOOPS); */
                /* for(size_t i = 0; i < ; i ++) { */
                /*     VECTOR(nodes_to_keep)[i] = current_cluster[i]; */
                /* } */
                /* while(!ConstrainedClustering::IsWellConnected */

                /* SetIgraphAllEdgesWeight(&induced_subgraph, 1.0); */
                MinCutCustom mcc(&subgraph);
                long edge_cut_size = mcc.ComputeMinCut();
                std::vector<long> in_partition = mcc.GetInPartition();
                std::vector<long> out_partition = mcc.GetOutPartition();
                
                // ========== DIAGNOSTIC: Track mincut completion ==========
                g_clusters_mincut_computed.fetch_add(1);
                
                // DIAGNOSTIC: Log small partitions
                if(in_partition.size() <= 2 || out_partition.size() <= 2) {
                    static std::mutex diag_mutex;
                    std::lock_guard<std::mutex> guard(diag_mutex);
                    std::cerr << "[DIAG_PART] Small partition: in=" << in_partition.size() 
                              << " out=" << out_partition.size() 
                              << " cut=" << edge_cut_size
                              << " total_nodes=" << igraph_vcount(&subgraph) << std::endl;
                }
                // =========================================================
                
                /* std::cerr << "got the cuts into " << std::to_string(in_partition.size()) << " and " << std::to_string(out_partition.size()) << std::endl; */
                bool current_criterion = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition.size(), out_partition.size(), edge_cut_size);
                
                // ========== DIAGNOSTIC: Log criterion result ==========
                if(!current_criterion && igraph_vcount(&subgraph) <= 5) {
                    static std::mutex diag_mutex;
                    std::lock_guard<std::mutex> guard(diag_mutex);
                    double criterion_needed = 0;
                    if(current_connectedness_criterion == ConnectednessCriterion::Logarithimic) {
                        criterion_needed = connectedness_criterion_c * (std::log(in_partition.size() + out_partition.size()) / pre_computed_log);
                    } else if(current_connectedness_criterion == ConnectednessCriterion::Exponential) {
                        criterion_needed = connectedness_criterion_c * std::pow(in_partition.size() + out_partition.size(), connectedness_criterion_x);
                    }
                    std::cerr << "[DIAG_CRIT] NOT well-connected: in=" << in_partition.size() 
                              << " out=" << out_partition.size()
                              << " cut=" << edge_cut_size 
                              << " criterion_needed=" << criterion_needed << std::endl;
                }
                // ======================================================
                
                /* if(connectedness_criterion == ConnectednessCriterion::Simple) { */
                /*     current_criterion = ConstrainedClustering::IsConnected(edge_cut_size); */
                /* } else if(connectedness_criterion == ConnectednessCriterion::Well) { */
                /*     current_criterion = ConstrainedClustering::IsWellConnected(in_partition, out_partition, edge_cut_size, &induced_subgraph); */
                /* } */

                if(!current_criterion) {
                    /* for(size_t i = 0; i < in_partition.size(); i ++) { */
                    /*     in_partition[i] = VECTOR(new_id_to_old_id_map)[in_partition[i]]; */
                    /* } */
                    /* for(size_t i = 0; i < out_partition.size(); i ++) { */
                    /*     out_partition[i] = VECTOR(new_id_to_old_id_map)[out_partition[i]]; */
                    /* } */
                    if(in_partition.size() == 1 || out_partition.size() == 1) {
                        static std::mutex diag_mutex;
                        std::lock_guard<std::mutex> guard(diag_mutex);
                        std::cerr << "[DROP] in=" << in_partition.size() 
                                << " out=" << out_partition.size() 
                                << " cut=" << edge_cut_size << std::endl;
                    }
                    if(in_partition.size() == 1 && out_partition.size() == 1) {
                        // ========== DIAGNOSTIC: Track both-singleton case ==========
                        g_both_singleton_fixed.fetch_add(1);
                        // ===========================================================
                        
                        // Output original cluster unsplit
                        std::vector<long> cc;
                        for (int i = 0; i < igraph_vcount(&subgraph); i++) {
                            cc.push_back(new_id_to_old_id_map[i]);
                        }
                        {
                            std::unique_lock<std::mutex> lock{MincutOnlyPreProcess::done_being_mincut_mutex};
                            MincutOnlyPreProcess::done_being_mincut_clusters.push(cc);
                        }
                        
                        // ========== DIAGNOSTIC: Track output ==========
                        g_clusters_output_notdone.fetch_add(1);
                        // ==============================================
                        
                        igraph_destroy(&subgraph);
                        return;  // Exit early
                    }
                    
                    // ========== DIAGNOSTIC: Track not-well-connected path ==========
                    g_clusters_output_notdone.fetch_add(1);
                    // ===============================================================
                    
                    if(in_partition.size() > 1) {
                        std::vector<std::vector<long>> in_clusters = GetConnectedComponentsOnPartition(&subgraph, in_partition, new_id_to_old_id_map);
                        
                        // ========== DIAGNOSTIC: Check for empty returns ==========
                        int empty_count = 0;
                        for(size_t i = 0; i < in_clusters.size(); i++) {
                            if(in_clusters[i].empty()) {
                                empty_count++;
                            }
                        }
                        if(empty_count > 0) {
                            static std::mutex diag_mutex;
                            std::lock_guard<std::mutex> guard(diag_mutex);
                            std::cerr << "[DIAG_EMPTY_RETURN] GetConnectedComponents returned " << empty_count 
                                      << " EMPTY vectors for in_partition (size=" << in_partition.size() << ")!" << std::endl;
                        }
                        // =========================================================
                        
                        for(size_t i = 0; i < in_clusters.size(); i ++) {
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnlyPreProcess::to_be_mincut_mutex);
                                MincutOnlyPreProcess::to_be_mincut_clusters.push(in_clusters[i]);
                            }
                        }
                    }
                    if(out_partition.size() > 1) {
                        std::vector<std::vector<long>> out_clusters = GetConnectedComponentsOnPartition(&subgraph, out_partition, new_id_to_old_id_map);
                        
                        // ========== DIAGNOSTIC: Check for empty returns ==========
                        int empty_count = 0;
                        for(size_t i = 0; i < out_clusters.size(); i++) {
                            if(out_clusters[i].empty()) {
                                empty_count++;
                            }
                        }
                        if(empty_count > 0) {
                            static std::mutex diag_mutex;
                            std::lock_guard<std::mutex> guard(diag_mutex);
                            std::cerr << "[DIAG_EMPTY_RETURN] GetConnectedComponents returned " << empty_count 
                                      << " EMPTY vectors for out_partition (size=" << out_partition.size() << ")!" << std::endl;
                        }
                        // =========================================================
                        
                        for(size_t i = 0; i < out_clusters.size(); i ++) {
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnlyPreProcess::to_be_mincut_mutex);
                                MincutOnlyPreProcess::to_be_mincut_clusters.push(out_clusters[i]);
                            }
                        }
                    }
                } else {
                    std::vector<long> cc;
                    for (int i = 0; i < igraph_vcount(&subgraph); i++) {
                        cc.push_back(new_id_to_old_id_map[i]);
                    }
                    // printf("after trivial current_cluster_size %d \n", current_cluster.size());
                    std::unique_lock<std::mutex> done_being_mincut_lock{MincutOnlyPreProcess::done_being_mincut_mutex};
                    MincutOnlyPreProcess::done_being_mincut_clusters.push(cc);
                    done_being_mincut_lock.unlock();
                    
                    // ========== DIAGNOSTIC: Track well-connected output ==========
                    g_clusters_output_done.fetch_add(1);
                    // =============================================================
                }
                igraph_destroy(&subgraph);
            }
        }
        
        // ========== DIAGNOSTIC: Global counters ==========
        static inline std::atomic<int> g_clusters_dequeued{0};
        static inline std::atomic<int> g_clusters_mincut_computed{0};
        static inline std::atomic<int> g_clusters_output_done{0};
        static inline std::atomic<int> g_clusters_output_notdone{0};
        static inline std::atomic<int> g_empty_vectors_found{0};
        static inline std::atomic<int> g_both_singleton_fixed{0};
        // =================================================
        
    private:
        static inline std::mutex to_be_mincut_mutex;
        static inline std::condition_variable to_be_mincut_condition_variable;
        static inline std::queue<std::vector<long>> to_be_mincut_clusters;
        static inline std::mutex done_being_mincut_mutex;
        static inline std::queue<std::vector<long>> done_being_mincut_clusters;
        std::string connectedness_criterion;
};

#endif