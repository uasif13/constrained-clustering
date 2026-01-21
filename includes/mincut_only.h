#ifndef MINCUT_ONLY_H
#define MINCUT_ONLY_H
#include "constrained.h"
#include <atomic>
#include <iomanip>
#include <sstream>

class MincutOnly : public ConstrainedClustering {
    public:
        MincutOnly(std::string edgelist, std::string existing_clustering, int num_processors, std::string output_file, std::string log_file, int log_level, std::string connectedness_criterion, std::string mincut_type) : ConstrainedClustering(edgelist, "", -1, existing_clustering, num_processors, output_file, log_file, "", log_level, connectedness_criterion, mincut_type) {
        };
        int main() override;

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
            
            // DIAGNOSTIC: Log CC extraction from partition
            {
                std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                std::cerr << "  [CC_EXTRACT] Partition size=" << partition.size() 
                         << " -> " << connected_components_vector.size() << " components" << std::endl;
            }
            
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

        static inline void MinCutWorker(igraph_t* graph, ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log, std::string mincut_type = "cactus") {
            // Thread-local statistics
            static thread_local int thread_id = -1;
            if (thread_id == -1) {
                static std::atomic<int> next_thread_id{0};
                thread_id = next_thread_id.fetch_add(1);
            }
            
            int clusters_processed = 0;
            int clusters_completed = 0;
            int clusters_split = 0;
            int singleton_in_count = 0;
            int singleton_out_count = 0;
            int64_t total_nodes_processed = 0;
            int64_t total_nodes_completed = 0;
            int64_t total_nodes_split = 0;
            
            {
                std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                std::cerr << "[WORKER-" << thread_id << "] Starting" << std::endl;
            }
            
            while (true) {
                std::unique_lock<std::mutex> to_be_mincut_lock{to_be_mincut_mutex};
                to_be_mincut_condition_variable.wait(to_be_mincut_lock, []() {
                    return !MincutOnly::to_be_mincut_clusters.empty();
                });
                
                std::vector<int> current_cluster = MincutOnly::to_be_mincut_clusters.front();
                MincutOnly::to_be_mincut_clusters.pop();
                
                // DIAGNOSTIC: Track dequeue operation
                total_clusters_dequeued++;
                size_t queue_size_after_dequeue = MincutOnly::to_be_mincut_clusters.size();
                
                to_be_mincut_lock.unlock();
                
                // Check for termination signal
                if(current_cluster.size() == 1 && current_cluster[0] == -1) {
                    std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                    std::cerr << "[WORKER-" << thread_id << "] TERMINATING" << std::endl;
                    std::cerr << "    Processed: " << clusters_processed
                              << " | Completed: " << clusters_completed 
                              << " | Split: " << clusters_split << std::endl;
                    std::cerr << "    Singleton IN: " << singleton_in_count
                              << " | Singleton OUT: " << singleton_out_count << std::endl;
                    std::cerr << "    Nodes: Processed=" << total_nodes_processed
                              << " Completed=" << total_nodes_completed
                              << " Split=" << total_nodes_split << std::endl;
                    return;
                }
                
                // DIAGNOSTIC: Log cluster dequeue
                clusters_processed++;
                total_nodes_processed += current_cluster.size();
                
                if (clusters_processed <= 5 || clusters_processed % 100 == 0) {
                    std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                    std::cerr << "[WORKER-" << thread_id << "] #" << clusters_processed
                              << " Dequeued cluster: size=" << current_cluster.size()
                              << " queue_remaining=" << queue_size_after_dequeue << std::endl;
                }
                
                // Validate cluster
                if (current_cluster.empty()) {
                    std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                    std::cerr << "[WORKER-" << thread_id << "] ERROR: Empty cluster dequeued!" << std::endl;
                    continue;
                }
                
                // Build induced subgraph
                igraph_vector_int_t nodes_to_keep;
                igraph_vector_int_t new_id_to_old_id_map;
                igraph_vector_int_init(&new_id_to_old_id_map, current_cluster.size());
                igraph_vector_int_init(&nodes_to_keep, current_cluster.size());
                for(size_t i = 0; i < current_cluster.size(); i ++) {
                    VECTOR(nodes_to_keep)[i] = current_cluster[i];
                }
                igraph_t induced_subgraph;
                igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);

                // Compute mincut
                MinCutCustom mcc(&induced_subgraph, mincut_type);
                int edge_cut_size = mcc.ComputeMinCut();
                std::vector<int> in_partition = mcc.GetInPartition();
                std::vector<int> out_partition = mcc.GetOutPartition();
                
                // DIAGNOSTIC: Log mincut results
                bool is_singleton_in = (in_partition.size() == 1);
                bool is_singleton_out = (out_partition.size() == 1);
                if (is_singleton_in) singleton_in_count++;
                if (is_singleton_out) singleton_out_count++;
                
                if (clusters_processed <= 5 || is_singleton_in || is_singleton_out || clusters_processed % 100 == 0) {
                    std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                    std::cerr << "[WORKER-" << thread_id << "] #" << clusters_processed << " MINCUT"
                              << " IN=" << in_partition.size() 
                              << " OUT=" << out_partition.size()
                              << " CUT=" << edge_cut_size;
                    if (is_singleton_in || is_singleton_out) {
                        std::cerr << " **SINGLETON**";
                    }
                    std::cerr << std::endl;
                }
                
                bool current_criterion = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition.size(), out_partition.size(), edge_cut_size);

                if(!current_criterion) {
                    // NOT well-connected: split the cluster
                    clusters_split++;
                    total_nodes_split += current_cluster.size();
                    
                    int components_added = 0;
                    int nodes_in_components = 0;
                    
                    // Process IN partition
                    if(in_partition.size() > 1) {
                        std::vector<std::vector<int>> in_clusters = GetConnectedComponentsOnPartition(&induced_subgraph, in_partition);
                        for(size_t i = 0; i < in_clusters.size(); i ++) {
                            std::vector<int> translated_in_clusters;
                            for(size_t j = 0; j < in_clusters[i].size(); j ++) {
                                translated_in_clusters.push_back(VECTOR(new_id_to_old_id_map)[in_clusters[i][j]]);
                            }
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnly::to_be_mincut_mutex);
                                MincutOnly::to_be_mincut_clusters.push(translated_in_clusters);
                                total_clusters_enqueued++;
                            }
                            components_added++;
                            nodes_in_components += translated_in_clusters.size();
                        }
                    } else if (in_partition.size() == 1) {
                        // DIAGNOSTIC: Track singleton handling
                        total_singleton_partitions_discarded++;
                        std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                        std::cerr << "[WORKER-" << thread_id << "] #" << clusters_processed 
                                  << " SINGLETON_IN discarded (size=1)" << std::endl;
                    }
                    
                    // Process OUT partition
                    if(out_partition.size() > 1) {
                        std::vector<std::vector<int>> out_clusters = GetConnectedComponentsOnPartition(&induced_subgraph, out_partition);
                        for(size_t i = 0; i < out_clusters.size(); i ++) {
                            std::vector<int> translated_out_clusters;
                            for(size_t j = 0; j < out_clusters[i].size(); j ++) {
                                translated_out_clusters.push_back(VECTOR(new_id_to_old_id_map)[out_clusters[i][j]]);
                            }
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnly::to_be_mincut_mutex);
                                MincutOnly::to_be_mincut_clusters.push(translated_out_clusters);
                                total_clusters_enqueued++;
                            }
                            components_added++;
                            nodes_in_components += translated_out_clusters.size();
                        }
                    } else if (out_partition.size() == 1) {
                        // DIAGNOSTIC: Track singleton handling
                        total_singleton_partitions_discarded++;
                        std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                        std::cerr << "[WORKER-" << thread_id << "] #" << clusters_processed 
                                  << " SINGLETON_OUT discarded (size=1)" << std::endl;
                    }
                    
                    // DIAGNOSTIC: Log split operation
                    if (clusters_processed <= 5 || clusters_processed % 100 == 0) {
                        std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                        std::cerr << "[WORKER-" << thread_id << "] #" << clusters_processed << " SPLIT"
                                  << " original_size=" << current_cluster.size()
                                  << " -> " << components_added << " components"
                                  << " total_nodes=" << nodes_in_components << std::endl;
                    }
                    
                    // CRITICAL CHECK: Node conservation
                    int expected_nodes = current_cluster.size();
                    int discarded_singletons = 0;
                    if (in_partition.size() == 1) discarded_singletons++;
                    if (out_partition.size() == 1) discarded_singletons++;
                    int expected_nodes_in_components = expected_nodes - discarded_singletons;
                    
                    if (nodes_in_components != expected_nodes_in_components) {
                        std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                        std::cerr << "[WORKER-" << thread_id << "] #" << clusters_processed 
                                  << " **NODE MISMATCH** Expected=" << expected_nodes_in_components
                                  << " Got=" << nodes_in_components
                                  << " (Original=" << expected_nodes 
                                  << " Discarded=" << discarded_singletons << ")" << std::endl;
                        total_node_conservation_violations++;
                    }
                    
                } else {
                    // Well-connected: mark as complete
                    clusters_completed++;
                    total_nodes_completed += current_cluster.size();
                    
                    std::unique_lock<std::mutex> done_being_mincut_lock{MincutOnly::done_being_mincut_mutex};
                    MincutOnly::done_being_mincut_clusters.push(current_cluster);
                    total_clusters_completed++;
                    done_being_mincut_lock.unlock();
                    
                    if (clusters_processed <= 5 || clusters_processed % 100 == 0) {
                        std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
                        std::cerr << "[WORKER-" << thread_id << "] #" << clusters_processed 
                                  << " COMPLETED size=" << current_cluster.size() << std::endl;
                    }
                }
                
                igraph_vector_int_destroy(&nodes_to_keep);
                igraph_vector_int_destroy(&new_id_to_old_id_map);
                igraph_destroy(&induced_subgraph);
            }
        }
        
    private:
        static inline std::mutex to_be_mincut_mutex;
        static inline std::condition_variable to_be_mincut_condition_variable;
        static inline std::queue<std::vector<int>> to_be_mincut_clusters;
        static inline std::mutex done_being_mincut_mutex;
        static inline std::queue<std::vector<int>> done_being_mincut_clusters;
        
        // DIAGNOSTIC: Global counters
        static inline std::atomic<uint64_t> total_clusters_dequeued{0};
        static inline std::atomic<uint64_t> total_clusters_enqueued{0};
        static inline std::atomic<uint64_t> total_clusters_completed{0};
        static inline std::atomic<uint64_t> total_singleton_partitions_discarded{0};
        static inline std::atomic<uint64_t> total_node_conservation_violations{0};
        
        // DIAGNOSTIC: Mutex for thread-safe stderr logging
        static inline std::mutex diagnostic_log_mutex;
        
        // Helper for thread-safe logging
        static inline void DiagnosticLog(const std::string& message) {
            std::lock_guard<std::mutex> guard(diagnostic_log_mutex);
            std::cerr << message << std::endl;
        }
};

#endif