#ifndef MMAP_GRAPH_LOADER_H
#define MMAP_GRAPH_LOADER_H
#include <igraph.h>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>

class MMapGraphLoader
{
public:
    struct Edge
    {
        long from;
        long to;
        double weight;
    };
    
    static int LoadEdgelistMMap(const std::string &filename,
                                igraph_t *graph,
                                std::unordered_map<long,long>* node_id_to_new_node_id_map,
                                bool has_weights)
    {
        int fd = open(filename.c_str(), O_RDONLY);
        if (fd == -1)
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return -1;
        }
        struct stat sb;
        if (fstat(fd, &sb) == -1)
        {
            std::cerr << "Failed to get file stats" << std::endl;
            close(fd);
            return -1;
        }

        size_t file_size = sb.st_size;
        char *mapped = (char *)mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapped == MAP_FAILED)
        {
            std::cerr << "mmap failed" << std::endl;
            close(fd);
            return -1;
        }
        madvise(mapped, file_size, MADV_SEQUENTIAL | MADV_WILLNEED);

        const char *ptr = mapped;
        const char *end = mapped + file_size;

        if (ptr < end && (*ptr < '0' || *ptr > '9'))
        {
            SkipToNextLine(ptr, end);
        }

        size_t estimated_edges = EstimateEdgeCount(mapped, file_size);
        printf("estimated_edge count of edgefile: %ld\n", estimated_edges);
        std::vector<Edge> edges;
        edges.reserve(estimated_edges);

        long next_node_id = 0;

        while (ptr < end)
        {
            while (ptr < end && (*ptr == '\n' || *ptr == '\r'))
                ptr++;
            if (ptr >= end)
                break;

            long from = ParseLong(ptr, end);

            while (ptr < end && (*ptr == '\t' || *ptr == ' '))
                ptr++;

            long to = ParseLong(ptr, end);

            double weight = 1.0;

            if (node_id_to_new_node_id_map->find(from) == node_id_to_new_node_id_map->end())
            {
                node_id_to_new_node_id_map->insert({from, next_node_id++});
            }
            if (node_id_to_new_node_id_map->find(to) == node_id_to_new_node_id_map->end())
            {
                node_id_to_new_node_id_map->insert({to, next_node_id++});
            }

            edges.push_back({node_id_to_new_node_id_map->at(from), node_id_to_new_node_id_map->at(to), weight});

            SkipToNextLine(ptr, end);
        }

        igraph_vector_int_t edge_vector;
        igraph_vector_int_init(&edge_vector, edges.size() * 2);
        for (size_t i = 0; i < edges.size(); i++)
        {
            VECTOR(edge_vector)[2 * i] = edges[i].from;
            VECTOR(edge_vector)[2 * i + 1] = edges[i].to;
        }

        igraph_empty(graph, next_node_id, IGRAPH_UNDIRECTED);
        igraph_add_edges(graph, &edge_vector, NULL);
        igraph_vector_int_destroy(&edge_vector);

        for (const auto& pair: *node_id_to_new_node_id_map) {
            SETVAS(graph, "name", pair.second, std::to_string(pair.first).c_str());
        }

        munmap(mapped, file_size);
        close(fd);
        return 0;
    }

    static int LoadClusteringMMap(const std::string &filename,
                                  std::unordered_map<long, long>* node_to_cluster,
                                  std::unordered_map<long, long> &node_id_to_new_node_id_cluster)
    {
        int fd = open(filename.c_str(), O_RDONLY);
        if (fd == -1) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return -1;
        }
        struct stat sb;
        if (fstat(fd, &sb) == -1) {
            std::cerr << "Failed to get file stats" << std::endl;
            close(fd);
            return -1;
        }
        size_t file_size = sb.st_size;
        char *mapped = (char *)mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapped == MAP_FAILED) {
            std::cerr << "mmap failed" << std::endl;
            close(fd);
            return -1;
        }
        madvise(mapped, file_size, MADV_SEQUENTIAL | MADV_WILLNEED);

        const char *ptr = mapped;
        const char *end = mapped + file_size;

        if (ptr < end && (*ptr < '0' || *ptr > '9'))
        {
            SkipToNextLine(ptr, end);
        }
        size_t estimated_clusters = EstimateEdgeCount(mapped, file_size);
        std::unordered_map<long, long> cluster_id_to_new_cluster_id_map;

        long cluster_id_new = 0;
        while (ptr < end) 
        {
            while (ptr < end && (*ptr == '\n' || *ptr == '\r'))
                ptr++;
            if (ptr >= end)
                break;

            long node = ParseLong(ptr, end);

            while (ptr < end && (*ptr == '\t' || *ptr == ' '))
                ptr++;

            long cluster = ParseLong(ptr, end);

            if (cluster_id_to_new_cluster_id_map.find(cluster) == cluster_id_to_new_cluster_id_map.end()) {
                cluster_id_to_new_cluster_id_map[cluster] = cluster_id_new++;
            }

            node_to_cluster->insert({node_id_to_new_node_id_cluster[node], cluster_id_to_new_cluster_id_map[cluster]});

            SkipToNextLine(ptr, end);
        }

        munmap(mapped, file_size);
        close(fd);
        return 0;
    }

private:
    static inline long ParseLong(const char *&ptr, const char *end)
    {
        long result = 0;
        while (ptr < end && (*ptr == ' ' || *ptr == '\t'))
            ptr++;

        while (ptr < end && *ptr >= '0' && *ptr <= '9')
        {
            result = result * 10 + (*ptr - '0');
            ptr++;
        }
        return result;
    }
    
    static inline void SkipToNextLine(const char *&ptr, const char *end)
    {
        while (ptr < end && *ptr != '\n')
            ptr++;
        if (ptr < end)
            ptr++;
    }
    
    static size_t EstimateEdgeCount(const char *data, size_t size);
};
#endif