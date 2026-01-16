#ifndef MINCUT_ONLY_PREPROCESS_H
#define MINCUT_ONLY_PREPROCESS_H
#include "constrained.h"



class MincutOnlyPreProcess : public ConstrainedClustering {
    public:
        // ✅ FIX: Added my_rank and nprocs parameters to match base class
        MincutOnlyPreProcess(std::string edgelist, int num_processors, std::string output_file, std::string log_file, std::string connectedness_criterion, int log_level) 
            : ConstrainedClustering(edgelist, "", -1, "", num_processors, output_file, log_file, log_level), 
              connectedness_criterion(connectedness_criterion) {
            std::cerr << "[CONSTRUCTOR] MincutOnlyPreProcess created" << std::endl;
            std::cerr << "[CONSTRUCTOR]   edgelist: " << edgelist << std::endl;
            std::cerr << "[CONSTRUCTOR]   num_processors: " << num_processors << std::endl;
            std::cerr << "[CONSTRUCTOR]   output_file: " << output_file << std::endl;
            std::cerr << "[CONSTRUCTOR]   connectedness_criterion: " << connectedness_criterion << std::endl;
        };
        
        // ✅ FIX: Changed signature to match base class virtual function
        int main() override;

        virtual ~MincutOnlyPreProcess() {
            std::cerr << "[DESTRUCTOR] MincutOnlyPreProcess destroyed" << std::endl;
        }

        static inline std::vector<std::vector<long>> GetConnectedComponentsOnPartition(const igraph_t* graph, std::vector<long>& partition, std::unordered_map<long, long> prev_new_id_to_old_id_map) {
            std::cerr << "[GETCC] Starting GetConnectedComponentsOnPartition" << std::endl;
            std::cerr << "[GETCC] Partition size: " << partition.size() << " nodes" << std::endl;
            
            std::vector<std::vector<long>> cluster_vectors;
            std::map<long, std::vector<long>> cluster_map;
            std::unordered_map<long, long> new_id_to_old_id_unordered_map;
            igraph_vector_int_t nodes_to_keep;
            igraph_vector_int_t new_id_to_old_id_map;
            
            std::cerr << "[GETCC] Initializing igraph vectors..." << std::endl;
            igraph_vector_int_init(&new_id_to_old_id_map, partition.size());
            igraph_vector_int_init(&nodes_to_keep, partition.size());
            
            for(size_t i = 0; i < partition.size(); i++) {
                VECTOR(nodes_to_keep)[i] = partition[i];
            }
            
            std::cerr << "[GETCC] Creating induced subgraph..." << std::endl;
            igraph_t induced_subgraph;

            igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
            
            std::cerr << "[GETCC] Induced subgraph created with " << igraph_vcount(&induced_subgraph) 
                     << " vertices and " << igraph_ecount(&induced_subgraph) << " edges" << std::endl;
            
            std::cerr << "[GETCC] Finding connected components in induced subgraph..." << std::endl;
            std::vector<std::vector<long>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&induced_subgraph);
            std::cerr << "[GETCC] Found " << connected_components_vector.size() << " connected components" << std::endl;
            
            for (int i = 0; i < connected_components_vector.size(); i++) {
                std::cerr << "[GETCC] Processing component " << i << " with " 
                         << connected_components_vector[i].size() << " nodes" << std::endl;
                
                std::vector<long> translated_cluster_vector;

                igraph_vector_int_t sub_nodes_to_keep;
                igraph_vector_int_t sub_new_id_to_old_id_vector_map;
                
                igraph_vector_int_init(&sub_nodes_to_keep, connected_components_vector[i].size());
                for(size_t j = 0; j < connected_components_vector[i].size(); j++) {
                    VECTOR(sub_nodes_to_keep)[j] = connected_components_vector[i][j];
                }
                
                std::cerr << "[GETCC]   Creating sub-subgraph..." << std::endl;
                igraph_t sub_subgraph;
                igraph_vector_int_init(&sub_new_id_to_old_id_vector_map, igraph_vector_int_size(&sub_nodes_to_keep));
                igraph_induced_subgraph_map(&induced_subgraph, &sub_subgraph, igraph_vss_vector(&sub_nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &sub_new_id_to_old_id_vector_map);
                
                std::cerr << "[GETCC]   Sub-subgraph: " << igraph_vcount(&sub_subgraph) 
                         << " vertices, " << igraph_ecount(&sub_subgraph) << " edges" << std::endl;
                
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
                
                std::cerr << "[GETCC]   Translated " << edge_count << " edges to original node IDs" << std::endl;
                
                igraph_vector_int_destroy(&sub_nodes_to_keep);
                igraph_vector_int_destroy(&sub_new_id_to_old_id_vector_map);
                igraph_destroy(&sub_subgraph);
                
                cluster_vectors.push_back(translated_cluster_vector);
            }
            
            igraph_vector_int_destroy(&nodes_to_keep);
            igraph_vector_int_destroy(&new_id_to_old_id_map);
            igraph_destroy(&induced_subgraph);
            
            std::cerr << "[GETCC] GetConnectedComponentsOnPartition complete, returning " 
                     << cluster_vectors.size() << " clusters" << std::endl;
            return cluster_vectors;
        }

        static inline void ComputeMinCutRecursive(std::vector<long> current_cluster, ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log) {
            static std::atomic<int> call_counter{0};
            int call_id = call_counter++;
            
            std::cerr << "\n[RECURSIVE #" << call_id << "] ======== ComputeMinCutRecursive START ========" << std::endl;
            std::cerr << "[RECURSIVE #" << call_id << "] Cluster size: " << current_cluster.size() << " values" << std::endl;
            
            if(current_cluster.size() == 1 || current_cluster[0] == -1) {
                std::cerr << "[RECURSIVE #" << call_id << "] Termination condition met (size=1 or marker=-1)" << std::endl;
                return;
            }
            
            std::cerr << "[RECURSIVE #" << call_id << "] Building node ID mappings..." << std::endl;
            std::unordered_map<long, long> new_id_to_old_id_map;
            std::unordered_map<long, long> newnew_id_to_old_id_map;
            std::unordered_map<long, long> old_id_to_new_id_map;

            igraph_vector_int_t edges;
            igraph_vector_int_init(&edges, current_cluster.size());
            std::set<long> current_cluster_set;

            long next_node_id = 0;
            for (int i = 0; i < current_cluster.size(); i++) {
                if (old_id_to_new_id_map.find(current_cluster[i]) == old_id_to_new_id_map.end()) {
                    current_cluster_set.insert(next_node_id);
                    old_id_to_new_id_map[current_cluster[i]] = next_node_id;
                    new_id_to_old_id_map[next_node_id] = current_cluster[i];
                    next_node_id++;
                }
                VECTOR(edges)[i] = old_id_to_new_id_map[current_cluster[i]];
            }
            
            std::cerr << "[RECURSIVE #" << call_id << "] Unique nodes: " << next_node_id << std::endl;
            std::cerr << "[RECURSIVE #" << call_id << "] Total edges: " << (current_cluster.size() / 2) << std::endl;
            
            std::cerr << "[RECURSIVE #" << call_id << "] Creating igraph subgraph..." << std::endl;
            igraph_t subgraph;
            igraph_empty(&subgraph, next_node_id, IGRAPH_UNDIRECTED);
            
            std::cerr << "[RECURSIVE #" << call_id << "] Adding edges to subgraph..." << std::endl;
            igraph_add_edges(&subgraph, &edges, NULL);
            
            std::cerr << "[RECURSIVE #" << call_id << "] Subgraph created: " << igraph_vcount(&subgraph) 
                     << " vertices, " << igraph_ecount(&subgraph) << " edges" << std::endl;

            std::cerr << "[RECURSIVE #" << call_id << "] Computing mincut..." << std::endl;
            MinCutCustom mcc(&subgraph);
            long edge_cut_size = mcc.ComputeMinCut();
            
            std::cerr << "[RECURSIVE #" << call_id << "] Mincut computed, edge_cut_size=" << edge_cut_size << std::endl;
            
            std::vector<long> in_partition = mcc.GetInPartition();
            std::vector<long> out_partition = mcc.GetOutPartition();
            
            std::cerr << "[RECURSIVE #" << call_id << "] Partition sizes: in=" << in_partition.size() 
                     << ", out=" << out_partition.size() << std::endl;
            
            std::cerr << "[RECURSIVE #" << call_id << "] Checking connectedness criterion..." << std::endl;
            bool current_criterion = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition.size(), out_partition.size(), edge_cut_size);
            
            std::cerr << "[RECURSIVE #" << call_id << "] Is well-connected: " << (current_criterion ? "YES" : "NO") << std::endl;

            if(!current_criterion) {
                std::cerr << "[RECURSIVE #" << call_id << "] Not well-connected, need to recurse on partitions" << std::endl;
                
                if(in_partition.size() > 1) {
                    std::cerr << "[RECURSIVE #" << call_id << "] Processing IN partition (" << in_partition.size() << " nodes)" << std::endl;
                    std::vector<std::vector<long>> in_clusters = GetConnectedComponentsOnPartition(&subgraph, in_partition, new_id_to_old_id_map);
                    std::cerr << "[RECURSIVE #" << call_id << "] IN partition split into " << in_clusters.size() << " clusters" << std::endl;
                    
                    for(size_t i = 0; i < in_clusters.size(); i++) {
                        std::cerr << "[RECURSIVE #" << call_id << "] Recursing on IN cluster " << i << std::endl;
                        MincutOnlyPreProcess::ComputeMinCutRecursive(in_clusters[i], current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log);
                    }
                } else {
                    std::cerr << "[RECURSIVE #" << call_id << "] IN partition too small (size=" << in_partition.size() << "), skipping" << std::endl;
                }
                
                if(out_partition.size() > 1) {
                    std::cerr << "[RECURSIVE #" << call_id << "] Processing OUT partition (" << out_partition.size() << " nodes)" << std::endl;
                    std::vector<std::vector<long>> out_clusters = GetConnectedComponentsOnPartition(&subgraph, out_partition, new_id_to_old_id_map);
                    std::cerr << "[RECURSIVE #" << call_id << "] OUT partition split into " << out_clusters.size() << " clusters" << std::endl;
                    
                    for(size_t i = 0; i < out_clusters.size(); i++) {
                        std::cerr << "[RECURSIVE #" << call_id << "] Recursing on OUT cluster " << i << std::endl;
                        MincutOnlyPreProcess::ComputeMinCutRecursive(out_clusters[i], current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log);
                    }
                } else {
                    std::cerr << "[RECURSIVE #" << call_id << "] OUT partition too small (size=" << out_partition.size() << "), skipping" << std::endl;
                }
            } else {
                std::cerr << "[RECURSIVE #" << call_id << "] Well-connected! Adding to done_being_mincut_clusters" << std::endl;
                
                std::vector<long> cc;
                for (int i = 0; i < igraph_vcount(&subgraph); i++) {
                    cc.push_back(new_id_to_old_id_map[i]);
                }
                
                std::cerr << "[RECURSIVE #" << call_id << "] Final cluster has " << cc.size() << " nodes" << std::endl;
                
                std::unique_lock<std::mutex> done_being_mincut_lock{MincutOnlyPreProcess::done_being_mincut_mutex};
                MincutOnlyPreProcess::done_being_mincut_clusters.push(cc);
                int total_done = MincutOnlyPreProcess::done_being_mincut_clusters.size();
                done_being_mincut_lock.unlock();
                
                std::cerr << "[RECURSIVE #" << call_id << "] Total completed clusters: " << total_done << std::endl;
            }
            
            std::cerr << "[RECURSIVE #" << call_id << "] Cleaning up..." << std::endl;
            igraph_destroy(&subgraph);
            igraph_vector_int_destroy(&edges);
            
            std::cerr << "[RECURSIVE #" << call_id << "] ======== ComputeMinCutRecursive END ========\n" << std::endl;
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