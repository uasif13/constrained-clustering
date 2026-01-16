#ifndef MINCUT_ONLY_PREPROCESS_H
#define MINCUT_ONLY_PREPROCESS_H
#include "constrained.h"



class MincutOnlyPreProcess : public ConstrainedClustering {
    public:
        MincutOnlyPreProcess(std::string edgelist, int num_processors, std::string output_file, std::string log_file, std::string connectedness_criterion, int log_level) : ConstrainedClustering(edgelist, "", -1, "", num_processors, output_file, log_file, log_level), connectedness_criterion(connectedness_criterion) {
            std::cerr << "[CONSTRUCTOR] MincutOnlyPreProcess created (QUEUE-BASED)" << std::endl;
            std::cerr << "[CONSTRUCTOR]   edgelist: " << edgelist << std::endl;
            std::cerr << "[CONSTRUCTOR]   num_processors: " << num_processors << std::endl;
            std::cerr << "[CONSTRUCTOR]   output_file: " << output_file << std::endl;
            std::cerr << "[CONSTRUCTOR]   connectedness_criterion: " << connectedness_criterion << std::endl;
        };
        int main() override;

        virtual ~MincutOnlyPreProcess() {
            std::cerr << "[DESTRUCTOR] MincutOnlyPreProcess destroyed" << std::endl;
        }

        static inline std::vector<std::vector<long>> GetConnectedComponentsOnPartition(const igraph_t* graph, std::vector<long>& partition, std::unordered_map<long, long> prev_new_id_to_old_id_map, int worker_id) {
            std::cerr << "[WORKER-" << worker_id << "][GETCC] Starting GetConnectedComponentsOnPartition" << std::endl;
            std::cerr << "[WORKER-" << worker_id << "][GETCC] Partition size: " << partition.size() << " nodes" << std::endl;
            
            std::vector<std::vector<long>> cluster_vectors;
            igraph_vector_int_t nodes_to_keep;
            igraph_vector_int_t new_id_to_old_id_map;
            
            igraph_vector_int_init(&new_id_to_old_id_map, partition.size());
            igraph_vector_int_init(&nodes_to_keep, partition.size());
            for(size_t i = 0; i < partition.size(); i ++) {
                VECTOR(nodes_to_keep)[i] = partition[i];
            }
            
            std::cerr << "[WORKER-" << worker_id << "][GETCC] Creating induced subgraph..." << std::endl;
            igraph_t induced_subgraph;
            igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
            
            std::cerr << "[WORKER-" << worker_id << "][GETCC] Induced subgraph: " << igraph_vcount(&induced_subgraph) 
                     << " vertices, " << igraph_ecount(&induced_subgraph) << " edges" << std::endl;
            
            std::cerr << "[WORKER-" << worker_id << "][GETCC] Finding connected components..." << std::endl;
            std::vector<std::vector<long>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&induced_subgraph);
            std::cerr << "[WORKER-" << worker_id << "][GETCC] Found " << connected_components_vector.size() << " connected components" << std::endl;
            
            for (size_t i = 0; i < connected_components_vector.size(); i++) {
                std::cerr << "[WORKER-" << worker_id << "][GETCC] Processing component " << i 
                         << " with " << connected_components_vector[i].size() << " nodes" << std::endl;
                
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
                
                std::cerr << "[WORKER-" << worker_id << "][GETCC]   Translated " << edge_count << " edges" << std::endl;
                
                igraph_vector_int_destroy(&sub_nodes_to_keep);
                igraph_vector_int_destroy(&sub_new_id_to_old_id_vector_map);
                igraph_destroy(&sub_subgraph);
                
                cluster_vectors.push_back(translated_cluster_vector);
            }
            
            igraph_vector_int_destroy(&nodes_to_keep);
            igraph_vector_int_destroy(&new_id_to_old_id_map);
            igraph_destroy(&induced_subgraph);
            
            std::cerr << "[WORKER-" << worker_id << "][GETCC] Complete, returning " << cluster_vectors.size() << " clusters" << std::endl;
            return cluster_vectors;
        }

        static inline void MinCutWorker(ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log) {
            // Assign worker ID
            static std::atomic<int> worker_id_counter{0};
            int worker_id = worker_id_counter++;
            
            std::cerr << "\n[WORKER-" << worker_id << "] ======== Worker Thread Started ========" << std::endl;
            
            int clusters_processed = 0;
            
            while (true) {
                std::unique_lock<std::mutex> to_be_mincut_lock{to_be_mincut_mutex};
                
                std::cerr << "[WORKER-" << worker_id << "] Waiting for work... (queue size: " 
                         << MincutOnlyPreProcess::to_be_mincut_clusters.size() << ")" << std::endl;
                
                to_be_mincut_condition_variable.wait(to_be_mincut_lock, []() {
                    return !MincutOnlyPreProcess::to_be_mincut_clusters.empty();
                });
                
                std::vector<long> current_cluster = MincutOnlyPreProcess::to_be_mincut_clusters.front();
                MincutOnlyPreProcess::to_be_mincut_clusters.pop();
                int remaining_in_queue = MincutOnlyPreProcess::to_be_mincut_clusters.size();
                to_be_mincut_lock.unlock();
                
                if(current_cluster.size() == 1 && current_cluster[0] == -1) {
                    std::cerr << "[WORKER-" << worker_id << "] Termination signal received, shutting down" << std::endl;
                    std::cerr << "[WORKER-" << worker_id << "] Total clusters processed: " << clusters_processed << std::endl;
                    std::cerr << "[WORKER-" << worker_id << "] ======== Worker Thread Terminated ========\n" << std::endl;
                    return;
                }
                
                clusters_processed++;
                std::cerr << "\n[WORKER-" << worker_id << "] ======== Processing Cluster #" << clusters_processed << " ========" << std::endl;
                std::cerr << "[WORKER-" << worker_id << "] Cluster size: " << current_cluster.size() << " values" << std::endl;
                std::cerr << "[WORKER-" << worker_id << "] Queue remaining: " << remaining_in_queue << " clusters" << std::endl;
                
                // Build ID mappings
                std::cerr << "[WORKER-" << worker_id << "] Building node ID mappings..." << std::endl;
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
                
                std::cerr << "[WORKER-" << worker_id << "] Unique nodes: " << next_node_id << std::endl;
                std::cerr << "[WORKER-" << worker_id << "] Total edges: " << (current_cluster.size() / 2) << std::endl;
                
                std::cerr << "[WORKER-" << worker_id << "] Creating igraph subgraph..." << std::endl;
                igraph_t subgraph;
                igraph_empty(&subgraph, next_node_id, IGRAPH_UNDIRECTED);
                
                std::cerr << "[WORKER-" << worker_id << "] Adding edges to subgraph..." << std::endl;
                igraph_add_edges(&subgraph, &edges, NULL);
                
                std::cerr << "[WORKER-" << worker_id << "] Subgraph created: " << igraph_vcount(&subgraph) 
                         << " vertices, " << igraph_ecount(&subgraph) << " edges" << std::endl;

                std::cerr << "[WORKER-" << worker_id << "] Computing mincut..." << std::endl;
                MinCutCustom mcc(&subgraph);
                long edge_cut_size = mcc.ComputeMinCut();
                
                std::cerr << "[WORKER-" << worker_id << "] Mincut computed, edge_cut_size=" << edge_cut_size << std::endl;
                
                std::vector<long> in_partition = mcc.GetInPartition();
                std::vector<long> out_partition = mcc.GetOutPartition();
                
                std::cerr << "[WORKER-" << worker_id << "] Partition sizes: in=" << in_partition.size() 
                         << ", out=" << out_partition.size() << std::endl;
                
                std::cerr << "[WORKER-" << worker_id << "] Checking connectedness criterion..." << std::endl;
                bool current_criterion = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition.size(), out_partition.size(), edge_cut_size);
                
                std::cerr << "[WORKER-" << worker_id << "] Is well-connected: " << (current_criterion ? "YES" : "NO") << std::endl;

                if(!current_criterion) {
                    std::cerr << "[WORKER-" << worker_id << "] Not well-connected, splitting partitions" << std::endl;
                    
                    int new_clusters_added = 0;
                    
                    if(in_partition.size() > 1) {
                        std::cerr << "[WORKER-" << worker_id << "] Processing IN partition (" << in_partition.size() << " nodes)" << std::endl;
                        std::vector<std::vector<long>> in_clusters = GetConnectedComponentsOnPartition(&subgraph, in_partition, new_id_to_old_id_map, worker_id);
                        std::cerr << "[WORKER-" << worker_id << "] IN partition split into " << in_clusters.size() << " clusters" << std::endl;
                        
                        for(size_t i = 0; i < in_clusters.size(); i ++) {
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnlyPreProcess::to_be_mincut_mutex);
                                MincutOnlyPreProcess::to_be_mincut_clusters.push(in_clusters[i]);
                                new_clusters_added++;
                            }
                        }
                        // Notify other workers
                        to_be_mincut_condition_variable.notify_all();
                    } else if (in_partition.size() == 1) {
                        std::cerr << "[WORKER-" << worker_id << "] IN partition is singleton, skipping" << std::endl;
                    }
                    
                    if(out_partition.size() > 1) {
                        std::cerr << "[WORKER-" << worker_id << "] Processing OUT partition (" << out_partition.size() << " nodes)" << std::endl;
                        std::vector<std::vector<long>> out_clusters = GetConnectedComponentsOnPartition(&subgraph, out_partition, new_id_to_old_id_map, worker_id);
                        std::cerr << "[WORKER-" << worker_id << "] OUT partition split into " << out_clusters.size() << " clusters" << std::endl;
                        
                        for(size_t i = 0; i < out_clusters.size(); i ++) {
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnlyPreProcess::to_be_mincut_mutex);
                                MincutOnlyPreProcess::to_be_mincut_clusters.push(out_clusters[i]);
                                new_clusters_added++;
                            }
                        }
                        // Notify other workers
                        to_be_mincut_condition_variable.notify_all();
                    } else if (out_partition.size() == 1) {
                        std::cerr << "[WORKER-" << worker_id << "] OUT partition is singleton, skipping" << std::endl;
                    }
                    
                    std::cerr << "[WORKER-" << worker_id << "] Added " << new_clusters_added << " new clusters to queue" << std::endl;
                } else {
                    std::cerr << "[WORKER-" << worker_id << "] Well-connected! Storing as final cluster" << std::endl;
                    
                    std::vector<long> cc;
                    for (int i = 0; i < igraph_vcount(&subgraph); i++) {
                        cc.push_back(new_id_to_old_id_map[i]);
                    }
                    
                    std::cerr << "[WORKER-" << worker_id << "] Final cluster has " << cc.size() << " nodes" << std::endl;
                    
                    std::unique_lock<std::mutex> done_being_mincut_lock{MincutOnlyPreProcess::done_being_mincut_mutex};
                    MincutOnlyPreProcess::done_being_mincut_clusters.push(cc);
                    int total_done = MincutOnlyPreProcess::done_being_mincut_clusters.size();
                    done_being_mincut_lock.unlock();
                    
                    std::cerr << "[WORKER-" << worker_id << "] Total completed clusters: " << total_done << std::endl;
                }
                
                std::cerr << "[WORKER-" << worker_id << "] Cleaning up..." << std::endl;
                igraph_vector_int_destroy(&edges);
                igraph_destroy(&subgraph);
                
                std::cerr << "[WORKER-" << worker_id << "] ======== Cluster #" << clusters_processed << " Complete ========\n" << std::endl;
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