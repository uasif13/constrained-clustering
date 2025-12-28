#ifndef CONSTRAINED_H
#define CONSTRAINED_H
// mincut_custom.h needs to be included first because of string override
#include "mincut_custom.h"

#include <cmath>
#include <chrono>
#include <condition_variable>
#include <random>
#include <thread>
#include <map>
#include <fstream>
#include <string>

#include <GraphHelper.h>
#include <libleidenalg/Optimiser.h>
#include <libleidenalg/CPMVertexPartition.h>
#include <libleidenalg/ModularityVertexPartition.h>

using namespace std;

enum Log {info, debug, error = -1};

template <typename T>
void output_set(set<T> input) {
    for (T val : input) {
        cout << val << " ";
    }
    cout << "\n";
}
template <typename T>
void output_vec_cc(std::vector<T> vec, igraph_t * graph) {
    for (auto item: vec) {
        cout << VAS(graph, "name",item) << " ";
    }
    cout << "\n";
}

class ConstrainedClustering {
    public:
        ConstrainedClustering(std::string edgelist, std::string algorithm, double clustering_parameter, std::string existing_clustering, int num_processors, std::string output_file, std::string log_file, int log_level, int my_rank, int nprocs) : edgelist(edgelist), algorithm(algorithm), clustering_parameter(clustering_parameter), existing_clustering(existing_clustering), num_processors(num_processors), output_file(output_file), log_file(log_file), log_level(log_level), my_rank(my_rank), nprocs(nprocs){
            if(this->log_level > -1) {
                this->start_time = std::chrono::steady_clock::now();
                this->log_file_handle.open(this->log_file);
            }
            this->num_calls_to_log_write = 0;
        };

        virtual ~ConstrainedClustering() {
            if(this->log_level > -1) {
                this->log_file_handle.close();
            }
        }

        virtual int main(int my_rank, int nprocs, uint64_t* opCount) = 0;
        int WriteToLogFile(std::string message, Log message_type, int my_rank = -1);
        void WritePartitionMap(std::map<long,long>& final_partition);
        void WriteClusterQueue(std::queue<std::vector<long>>& to_be_clustered_clusters, igraph_t* graph);
        long WriteClusterQueueMPI(std::queue<std::vector<long>>* cluster_queue, long cluster_start, long previous_cluster_id, int iteration, uint64_t* opCount);
        long WriteClusterQueueMPI(std::queue<std::vector<long>>* cluster_queue, igraph_t* graph, long cluster_start, long previous_cluster_id, int iteration, uint64_t* opCount);

        int initializeSlice(igraph_t * graph){
            this -> vertex_count = igraph_vcount(graph);
            long my_work = vertex_count/nprocs;
            if (vertex_count % nprocs != 0)
                my_work ++;
            this -> start_vertex = my_rank*my_work;
            this -> end_vertex = (my_rank+1)*my_work;
            printf("my_rank: %d start: %ld end: %ld\n", my_rank, this -> start_vertex, this -> end_vertex);
            igraph_vector_int_init(&(this-> vertex_vec),0);
            igraph_eit_t eit;
            igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                long from_node = IGRAPH_FROM(graph, current_edge);
                long to_node = IGRAPH_TO(graph, current_edge);
                if (from_node >= this -> start_vertex && from_node < this -> end_vertex || to_node >= this -> start_vertex && to_node < this -> end_vertex) {
                    igraph_vector_int_push_back(&(this -> vertex_vec),from_node);
                    igraph_vector_int_push_back(&(this -> vertex_vec), to_node);
                    this -> vertex_set.insert(from_node);
                    this -> vertex_set.insert(to_node);
                }
            }
            return 0;
        }
            
        static inline std::map<std::string, long> GetOriginalToNewIdMap(igraph_t* graph) {
            std::map<std::string, long> original_to_new_id_map;
            igraph_vit_t vit;
            igraph_vit_create(graph, igraph_vss_all(), &vit);
            for(; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
                 igraph_integer_t current_node = IGRAPH_VIT_GET(vit);
                 original_to_new_id_map[VAS(graph, "name", current_node)] = current_node;
            }
            igraph_vit_destroy(&vit);
            return original_to_new_id_map;
        }

        std::map<std::string, long> GetOriginalToNewIdMapDistributed(igraph_t* graph) {
            std::map<std::string, long> original_to_new_id_map;
            igraph_vit_t vertex_iter;
            igraph_vs_t vertex_select = igraph_vss_vector(&(this -> vertex_vec));
            igraph_vit_create(graph, vertex_select, &vertex_iter);
            for(;!IGRAPH_VIT_END(vertex_iter); IGRAPH_VIT_NEXT(vertex_iter)) {
                long node_id = IGRAPH_VIT_GET(vertex_iter);
                original_to_new_id_map[VAS(graph, "name", node_id)] = node_id;
            }
            return original_to_new_id_map;
        }

        static inline std::map<long, long> ReadCommunities(const std::map<std::string, long>& original_to_new_id_map, std::string existing_clustering) {
            std::map<long, long> partition_map;
            std::ifstream existing_clustering_file(existing_clustering);
            std::string node_id;
            long cluster_id = -1;
            while (existing_clustering_file >> node_id >> cluster_id) {
                if (original_to_new_id_map.contains(node_id)) {
                    long new_node_id = original_to_new_id_map.at(node_id);
                    partition_map[new_node_id] = cluster_id;
                }
            }
            return partition_map;
        }

        static inline std::map<long, long> ReadCommunities(std::string existing_clustering) {
            std::map<long, long> partition_map;
            std::ifstream existing_clustering_file(existing_clustering);
            long node_id = -1;
            long cluster_id = -1;
            while (existing_clustering_file >> node_id >> cluster_id) {
                partition_map[node_id] = cluster_id;
            }
            return partition_map;
        }

        std::map<long, long> ReadCommunitiesDistributed(std::string existing_clustering) {
            std::map<long, long> partition_map;
            std::ifstream existing_clustering_file(existing_clustering);
            long node_id = -1;
            long cluster_id = -1;
            while (existing_clustering_file >> node_id >> cluster_id) {
                if (this -> vertex_set.contains(node_id)) {
                    partition_map[node_id] = cluster_id;
                }
            }
            return partition_map;
        }

        static inline void ClusterQueueToMap(std::queue<std::vector<long>>& cluster_queue, std::map<long, long>& node_id_to_cluster_id_map) {
            long cluster_id = 0;
            while(!cluster_queue.empty()) {
                std::vector<long> current_cluster = cluster_queue.front();
                cluster_queue.pop();
                for(size_t i = 0; i < current_cluster.size(); i ++) {
                    node_id_to_cluster_id_map[current_cluster[i]] = cluster_id;
                }
                cluster_id ++;
            }
        }
        
        static inline void RemoveInterClusterEdges(igraph_t* graph, const std::map<long,long>& node_id_to_cluster_id_map) {
            igraph_vector_int_t edges_to_keep;
            igraph_vector_int_init(&edges_to_keep, 0);
            igraph_eit_t eit;
            igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                long from_node = IGRAPH_FROM(graph, current_edge);
                long to_node = IGRAPH_TO(graph, current_edge);
                if(node_id_to_cluster_id_map.contains(from_node) && node_id_to_cluster_id_map.contains(to_node)
                    && (node_id_to_cluster_id_map.at(from_node) == node_id_to_cluster_id_map.at(to_node))) {
                    // keep the edge
                    igraph_vector_int_push_back(&edges_to_keep, current_edge);
                } 
            }
            igraph_t new_graph;
            igraph_es_t es;
            igraph_es_vector_copy(&es, &edges_to_keep);
            igraph_subgraph_from_edges(graph, &new_graph, es, false);
            igraph_destroy(graph);
            *graph = new_graph;
            igraph_eit_destroy(&eit);
            igraph_es_destroy(&es);
            igraph_vector_int_destroy(&edges_to_keep);
        }

        // currently keeps only those edges that go from within these clusters defined in the map
        static inline void RemoveInterClusterEdges(igraph_t* graph, const std::unordered_map<long,long>& node_id_to_cluster_id_map, int num_threads) {
            const long num_nodes = igraph_vcount(graph);
            const long total_edges = igraph_ecount(graph);
            
            std::vector<long> node_to_cluster(num_nodes, -1);
            for (const auto& [node, cluster]: node_id_to_cluster_id_map) {
                if (node < num_nodes) {
                    node_to_cluster[node] = cluster;
                }
            }

            long thread_edges_size = total_edges / num_threads / 2;
            
            std::vector<std::vector<long>> thread_local_edges(num_threads);
            for (int i = 0; i < num_threads; i++) {
                thread_local_edges[i].reserve(thread_edges_size);
            }

            auto filter_edges = [&](long thread_id, long start_edge, long end_edge) {
                igraph_eit_t eit;
                igraph_es_t es;
                igraph_es_range(&es, start_edge, end_edge);
                igraph_eit_create(graph, es, &eit);

                auto& local_edges = thread_local_edges[thread_id];

                for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                    igraph_integer_t edge_id = IGRAPH_EIT_GET(eit);
                    long from_node = IGRAPH_FROM(graph, edge_id);
                    long to_node = IGRAPH_TO(graph, edge_id);

                    long from_cluster = node_to_cluster[from_node];
                    long to_cluster = node_to_cluster[to_node];

                    if (from_cluster != -1 && to_cluster != -1 && from_cluster == to_cluster) {
                        local_edges.push_back(edge_id);
                    }
                }
                igraph_eit_destroy(&eit);
                igraph_es_destroy(&es);
            };

            std::vector<std::thread> threads;
            const long chunk_size = (total_edges + num_threads - 1) / num_threads;

            for (int i = 0; i < num_threads; i++) {
                long start = static_cast<long>(i) * chunk_size;
                long end = std::min(start + chunk_size, total_edges);
                if (start < end) {
                    threads.emplace_back(filter_edges, i, start, end);
                }
            }

            for (auto& t: threads) {
                t.join();
            }

            long total_kept = 0;
            for (const auto& local: thread_local_edges) {
                total_kept += local.size();
            }

            igraph_vector_int_t edges_to_keep;
            igraph_vector_int_init(&edges_to_keep, total_kept);

            long idx = 0;
            for (const auto& local: thread_local_edges) {
                for (long edge: local) {
                    VECTOR(edges_to_keep)[idx++] = edge;
                }
            }
            
            igraph_t new_graph;
            igraph_es_t es;
            igraph_es_vector_copy(&es, &edges_to_keep);
            igraph_subgraph_from_edges(graph, &new_graph, es, false);
            igraph_destroy(graph);
            *graph = new_graph;
            igraph_es_destroy(&es);
            igraph_vector_int_destroy(&edges_to_keep);
        }

        // currently keeps only those edges that go from within these clusters defined in the map
        static inline igraph_vector_int_t RemoveInterClusterEdgesDistributed(igraph_t* graph, const std::map<long,long>& node_id_to_cluster_id_map) {
            igraph_vector_int_t edges_to_remove;
            igraph_vector_int_init(&edges_to_remove, 0);
            igraph_eit_t eit;
            igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                long from_node = IGRAPH_FROM(graph, current_edge);
                long to_node = IGRAPH_TO(graph, current_edge);
                if(node_id_to_cluster_id_map.contains(from_node) && node_id_to_cluster_id_map.contains(to_node)) {
                    if (node_id_to_cluster_id_map.at(from_node) != node_id_to_cluster_id_map.at(to_node)) {
                        igraph_vector_int_push_back(&edges_to_remove, IGRAPH_EIT_GET(eit));
                    }
                }
            }
            return edges_to_remove;
        }

        static inline void RestoreInterClusterEdges(const igraph_t* original_graph, igraph_t* graph, const std::map<long,long>& node_id_to_cluster_id_map) {
            struct ClusterEdgeStruct {
                long cluster_id;
                long other_cluster_id;
                double intercluster_num_edges;
                double cluster_sum_edges;
                bool operator<(const struct ClusterEdgeStruct& other) const {
                    return (this->intercluster_num_edges / this->cluster_sum_edges) < (other.intercluster_num_edges / other.cluster_sum_edges);
                }
            };

            std::map<long, std::vector<long>> cluster_id_to_node_vec_map;
            std::map<long, std::set<long>> cluster_id_to_node_set_map;
            std::map<long, long> cluster_id_to_num_edges_map;
            std::map<long, std::map<long, long>> cluster_id_to_cluster_id_to_num_edges_map;
            std::priority_queue<ClusterEdgeStruct> pq;
            std::map<long, std::vector<long>> clusters_to_be_merged;
            std::set<long> clusters_already_chosen_to_be_merged;
            for(auto const& [node_id, cluster_id]: node_id_to_cluster_id_map) {
                cluster_id_to_node_vec_map[cluster_id].push_back(node_id);
                cluster_id_to_node_set_map[cluster_id].insert(node_id);
            }

            igraph_eit_t eit;
            igraph_eit_create(original_graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                long from_node = IGRAPH_FROM(original_graph, current_edge);
                long to_node = IGRAPH_TO(original_graph, current_edge);
                for(auto const& [cluster_id, node_id_set]: cluster_id_to_node_set_map) {
                    if(node_id_set.contains(from_node) && node_id_set.contains(to_node)) {
                        cluster_id_to_num_edges_map[cluster_id] += 1;
                    }
                    for(auto const& [other_cluster_id, other_node_id_set]: cluster_id_to_node_set_map) {
                        if((node_id_set.contains(from_node) && other_node_id_set.contains(to_node)) || (node_id_set.contains(to_node) && other_node_id_set.contains(from_node))) {
                            cluster_id_to_cluster_id_to_num_edges_map[cluster_id][other_cluster_id] += 1;
                        }
                    }
                }
            }
            igraph_eit_destroy(&eit);

            for(auto const& [cluster_id, other_cluster_id_to_num_edges_map]: cluster_id_to_cluster_id_to_num_edges_map) {
                for(auto const& [other_cluster_id, num_edges_from_cluster_id_to_other_cluster_id]: other_cluster_id_to_num_edges_map) {
                    long current_cluster_id_num_edges = cluster_id_to_num_edges_map[cluster_id];
                    long other_cluster_id_num_edges = cluster_id_to_num_edges_map[other_cluster_id];
                    double intercluster_num_edges = num_edges_from_cluster_id_to_other_cluster_id;
                    double cluster_sum_edges = current_cluster_id_num_edges + other_cluster_id_num_edges;
                    double threshold_value = 4;
                    if((current_cluster_id_num_edges / threshold_value > intercluster_num_edges) && (other_cluster_id_num_edges / threshold_value > intercluster_num_edges)) {
                        pq.push(ClusterEdgeStruct{.cluster_id=cluster_id, .other_cluster_id=other_cluster_id, .intercluster_num_edges=intercluster_num_edges, .cluster_sum_edges=cluster_sum_edges});
                    }
                }
            }
            while(!pq.empty()) {
                ClusterEdgeStruct current_cluster_edge_struct = pq.top();
                pq.pop();
                long cluster_id = current_cluster_edge_struct.cluster_id;
                long other_cluster_id = current_cluster_edge_struct.other_cluster_id;
                if(clusters_already_chosen_to_be_merged.contains(cluster_id) || clusters_already_chosen_to_be_merged.contains(other_cluster_id)) {
                } else {
                    clusters_already_chosen_to_be_merged.insert(cluster_id);
                    clusters_already_chosen_to_be_merged.insert(other_cluster_id);
                }
                ConstrainedClustering::RestoreInterClusterEdges(original_graph, graph, node_id_to_cluster_id_map, cluster_id, other_cluster_id);
            }
        }
        
        static inline void RestoreInterClusterEdges(const igraph_t* original_graph, igraph_t* graph, const std::map<long,long>& node_id_to_cluster_id_map, long cluster_id, long other_cluster_id) {
            igraph_eit_t eit;
            igraph_eit_create(original_graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                long from_node = IGRAPH_FROM(original_graph, current_edge);
                long to_node = IGRAPH_TO(original_graph, current_edge);
                if(node_id_to_cluster_id_map.at(from_node) == cluster_id &&  node_id_to_cluster_id_map.at(to_node) == other_cluster_id) {
                    igraph_add_edge(graph, from_node, to_node);
                }
                if(node_id_to_cluster_id_map.at(to_node) == cluster_id &&  node_id_to_cluster_id_map.at(from_node) == other_cluster_id) {
                    igraph_add_edge(graph, to_node, from_node);
                }
            }
            igraph_eit_destroy(&eit);
        }

        static inline void RunLouvainAndUpdatePartition(std::map<long, long>& partition_map, int seed, igraph_t* graph) {
            igraph_vector_int_t membership;
            igraph_vector_int_init(&membership, 0);
            igraph_rng_seed(igraph_rng_default(), seed);
            igraph_community_multilevel(graph, 0, 1, &membership, 0, 0);

            igraph_eit_t eit;
            igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            std::set<long> visited;
            for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                long from_node = IGRAPH_FROM(graph, current_edge);
                long to_node = IGRAPH_TO(graph, current_edge);
                if(!visited.contains(from_node)) {
                    visited.insert(from_node);
                    partition_map[from_node] = VECTOR(membership)[from_node];
                }
                if(!visited.contains(to_node)) {
                    visited.insert(to_node);
                    partition_map[to_node] = VECTOR(membership)[to_node];
                }
            }
            igraph_eit_destroy(&eit);
            igraph_vector_int_destroy(&membership);
        }

        static inline void RunLeidenAndUpdatePartition(std::map<long, long>& partition_map, MutableVertexPartition* partition, int seed, igraph_t* graph, int num_iter= 2) {
            Optimiser o;
            o.set_rng_seed(seed);
            for(int i = 0; i < num_iter; i ++) {
                o.optimise_partition(partition);
            }
            igraph_eit_t eit;
            igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            std::set<long> visited;
            for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                igraph_integer_t current_edge = IGRAPH_EIT_GET(eit);
                long from_node = IGRAPH_FROM(graph, current_edge);
                long to_node = IGRAPH_TO(graph, current_edge);
                if(!visited.contains(from_node)) {
                    visited.insert(from_node);
                    partition_map[from_node] = partition->membership(from_node);
                }
                if(!visited.contains(to_node)) {
                    visited.insert(to_node);
                    partition_map[to_node] = partition->membership(to_node);
                }
            }
            igraph_eit_destroy(&eit);
        }

        static inline void SetIgraphAllEdgesWeight(igraph_t* graph, double weight) {
            igraph_eit_t eit;
            igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
            for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                SETEAN(graph, "weight", IGRAPH_EIT_GET(eit), 1);
            }
            igraph_eit_destroy(&eit);
        }

        static inline std::map<long, long> GetCommunities(std::string edgelist, std::string algorithm, int seed, double clustering_parameter, igraph_t* graph_ptr) {
            std::map<long, long> partition_map;
            igraph_t graph;
            if(graph_ptr == NULL) {
                FILE* edgelist_file = fopen(edgelist.c_str(), "r");
                igraph_set_attribute_table(&igraph_cattribute_table);
                igraph_read_graph_ncol(&graph, edgelist_file, NULL, 1, IGRAPH_ADD_WEIGHTS_IF_PRESENT, IGRAPH_UNDIRECTED);
                fclose(edgelist_file);
            } else {
                graph = *graph_ptr;
            }

            if(algorithm == "louvain") {
                ConstrainedClustering::RunLouvainAndUpdatePartition(partition_map, seed, &graph);
            } else if(algorithm == "leiden-cpm") {
                std::vector<double> edge_weights;
                igraph_eit_t eit;
                igraph_eit_create(&graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
                for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                    double current_edge_weight = 1.0;
                    edge_weights.push_back(current_edge_weight);
                }
                igraph_eit_destroy(&eit);
                Graph* leiden_graph = Graph::GraphFromEdgeWeights(&graph, edge_weights);
                CPMVertexPartition partition(leiden_graph, clustering_parameter);
                ConstrainedClustering::RunLeidenAndUpdatePartition(partition_map, &partition, seed, &graph);
            } else if(algorithm == "leiden-mod") {
                std::vector<double> edge_weights;
                igraph_eit_t eit;
                igraph_eit_create(&graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
                for(; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
                    double current_edge_weight = EAN(&graph, "weight", IGRAPH_EIT_GET(eit));
                    edge_weights.push_back(current_edge_weight);
                }
                igraph_eit_destroy(&eit);
                Graph* leiden_graph = Graph::GraphFromEdgeWeights(&graph, edge_weights);
                ModularityVertexPartition partition(leiden_graph);
                ConstrainedClustering::RunLeidenAndUpdatePartition(partition_map, &partition, seed, &graph);
            } else {
                throw std::invalid_argument("GetCommunities(): Unsupported algorithm");
            }

            if(graph_ptr == NULL) {
                igraph_destroy(&graph);
            }
            return partition_map;
        }

        static inline std::vector<std::vector<long>> GetConnectedComponents(igraph_t* graph_ptr) {
            std::vector<std::vector<long>> connected_components_vector;
            std::map<long, std::vector<long>> component_id_to_member_vector_map;
            igraph_vector_int_t component_id_vector;
            igraph_vector_int_init(&component_id_vector, 0);
            igraph_vector_int_t membership_size_vector;
            igraph_vector_int_init(&membership_size_vector, 0);
            igraph_integer_t number_of_components;
            igraph_connected_components(graph_ptr, &component_id_vector, &membership_size_vector, &number_of_components, IGRAPH_WEAK);
            for(long node_id = 0; node_id < igraph_vcount(graph_ptr); node_id ++) {
                long current_component_id = VECTOR(component_id_vector)[node_id];
                if(VECTOR(membership_size_vector)[current_component_id] > 1) {
                    component_id_to_member_vector_map[current_component_id].push_back(node_id);
                }
            }
            igraph_vector_int_destroy(&component_id_vector);
            igraph_vector_int_destroy(&membership_size_vector);
            for(auto const& [component_id, member_vector] : component_id_to_member_vector_map) {
                connected_components_vector.push_back(member_vector);
            }
            return connected_components_vector;
        }

        std::vector<std::vector<long>> GetConnectedComponentsDistributed(igraph_t* graph_ptr, const std::unordered_map<long,long>& node_id_to_cluster_id_map, int my_rank, int nprocs) {
            std::vector<std::vector<long>> connected_components_vector;
            std::unordered_map<long, std::vector<long>> component_id_to_member_vector_map;
            igraph_vector_int_t component_id_vector;
            igraph_vector_int_init(&component_id_vector, 0);
            igraph_vector_int_t membership_size_vector;
            igraph_vector_int_init(&membership_size_vector, 0);

            printf("my_rank: %d connected components algorithm start\n", my_rank);
            igraph_integer_t number_of_components;
            igraph_bitset_t already_added;
            igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
            igraph_vector_int_t neis = IGRAPH_VECTOR_NULL;
            long no_of_vertices = igraph_vcount(graph_ptr);
            igraph_vector_int_resize(&component_id_vector, no_of_vertices);
            for (long vertex = 0; vertex < no_of_vertices; ++vertex) {
                VECTOR(component_id_vector)[vertex] = -1;
            }
            igraph_bitset_init(&already_added, no_of_vertices);
            igraph_dqueue_int_init(&q, no_of_vertices > 100000 ? 10000: no_of_vertices);
            igraph_vector_int_init(&neis, 0);
            long no_of_components = 0;

            std::vector<long> node_to_original_id(igraph_vcount(graph_ptr));
            // for(long i = 0; i < igraph_vcount(graph_ptr); i++) {
            //     node_to_original_id[i] = stol(VAS(graph_ptr, "name", i));
            // }
            //printf("my_rank: %d connected components algorithm initialized\n");
            for (long vertex = 0; vertex < no_of_vertices; ++vertex) {
                if (!node_id_to_cluster_id_map.contains(vertex)) {
                    continue;
                }
                if (node_id_to_cluster_id_map.at(vertex) % nprocs != my_rank) {
                    continue;
                }
                //printf("vertex to check connection: %ld\n", vertex);
                long act_component_size;
                if (IGRAPH_BIT_TEST(already_added, vertex)) {
                    continue;
                }
                IGRAPH_BIT_SET(already_added, vertex);
                act_component_size = 1;
                VECTOR(component_id_vector)[vertex] = no_of_components;
                igraph_dqueue_int_push(&q, vertex);

                // printf("my_rank: %d connected components queue start\n", my_rank);
                while (!igraph_dqueue_int_empty(&q)) {
                    long act_node = igraph_dqueue_int_pop(&q);
                    // printf("act_node inside queue: %ld\n", act_node);

                    igraph_neighbors(graph_ptr, &neis, act_node, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE);
                    long nei_count = igraph_vector_int_size(&neis);
                    for (long i = 0; i < nei_count; i++) {
                        long neighbor = VECTOR(neis)[i];
                        // printf("neighbor inside queue: %ld\n", neighbor);
                        if (!node_id_to_cluster_id_map.contains(neighbor)) {
                            continue;
                        }
                        if (node_id_to_cluster_id_map.at(neighbor) % nprocs != my_rank) {
                            continue;
                        }
                        if (IGRAPH_BIT_TEST(already_added, neighbor)) {
                            continue;
                        }
                        igraph_dqueue_int_push(&q, neighbor);
                        IGRAPH_BIT_SET(already_added, neighbor);
                        act_component_size++;
                        VECTOR(component_id_vector)[neighbor] = no_of_components;
                    }
                }
                no_of_components++;
                igraph_vector_int_push_back(&membership_size_vector, act_component_size);
            }
            // igraph_vector_int_print(&component_id_vector);
            /* std::cerr << "num con comp: " << number_of_components << std::endl; */
            for (long node_id = 0; node_id < no_of_vertices; node_id++) {
                if (!node_id_to_cluster_id_map.contains(node_id)) {
                    continue;
                }
                if (node_id_to_cluster_id_map.at(node_id) % nprocs != my_rank) {
                    continue;
                }
                long current_component_id = VECTOR(component_id_vector)[node_id];
                /* std::cerr << "component id: " << current_component_id << std::endl; */
                /* std::cerr << "component size: " << VECTOR(membership_size_vector)[current_component_id] << std::endl; */
                /* std::cerr << "graph node id: " << node_id << std::endl; */
                /* std::cerr << "original node id: " << VAS(graph_ptr, "name", node_id) << std::endl; */
                if (VECTOR(membership_size_vector)[current_component_id] > 1) {
                    component_id_to_member_vector_map[current_component_id].push_back(node_id);
                }
            }
            igraph_vector_int_destroy(&component_id_vector);
            igraph_vector_int_destroy(&membership_size_vector);
            for (auto const& [component_id, member_vector] : component_id_to_member_vector_map) {
                if (component_id != -1) {
                    connected_components_vector.push_back(member_vector);
                }
            }
            return connected_components_vector;
        }
        static inline bool IsConnected(int edge_cut_size) {
            return edge_cut_size >= 1;
        }

        static inline bool IsWellConnected(long in_partition_size, long out_partition_size, int edge_cut_size) {
            bool node_connectivity = log10(in_partition_size + out_partition_size) < edge_cut_size;
            return node_connectivity;
        }

        static inline bool IsWellConnected(const std::vector<long>& in_partition, const std::vector<long>& out_partition, int edge_cut_size, const igraph_t* induced_subgraph) {
            if(edge_cut_size == 0) {
                return false;
            }
            bool node_connectivity = log10(in_partition.size() + out_partition.size()) < edge_cut_size;
            return node_connectivity;
        }

    protected:
        std::string edgelist;
        std::string algorithm;
        double clustering_parameter;
        std::string existing_clustering;
        int num_processors;
        int my_rank;
        int nprocs;
        std::string output_file;
        std::string log_file;
        std::chrono::steady_clock::time_point start_time;
        std::ofstream log_file_handle;
        int log_level;
        int num_calls_to_log_write;
        long vertex_count;
        long start_vertex;
        long end_vertex;
        igraph_vector_int_t edge_slice;
        set<long> vertices_in_slice;
        igraph_vector_int_t vertex_vec;
        set<long> vertex_set;
};

#endif