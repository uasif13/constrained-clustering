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
        int from;
        int to;
        double weight;
    };
    static int LoadEdgelistMMap(const std::string &filename,
                                igraph_t *graph,
                                //   std::unordered_map<int,int> node_id_to_new_node_id_map,
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
        std::vector<Edge> edges;
        edges.reserve(estimated_edges);

        int next_node_id = 0;
        std::unordered_map<int, int> node_id_to_new_node_id_map;

        while (ptr < end)
        {
            while (ptr < end && (*ptr == '\n' || *ptr == '\r'))
                ptr++;
            if (ptr >= end)
                break;

            int from = ParseInt(ptr, end);

            while (ptr < end && (*ptr == '\t' || *ptr == ' '))
                ptr++;

            int to = ParseInt(ptr, end);

            double weight = 1.0;

            if (node_id_to_new_node_id_map.find(from) == node_id_to_new_node_id_map.end())
            {
                node_id_to_new_node_id_map[from] = next_node_id++;
            }
            if (node_id_to_new_node_id_map.find(to) == node_id_to_new_node_id_map.end())
            {
                node_id_to_new_node_id_map[to] = next_node_id++;
            }

            edges.push_back({node_id_to_new_node_id_map[from], node_id_to_new_node_id_map[to], weight});

            SkipToNextLine(ptr, end);
        }

        igraph_vector_int_t edge_vector;
        igraph_vector_int_init(&edge_vector, edges.size() * 2);
        for (size_t i = 0; i < edges.size(); i++)
        {
            VECTOR(edge_vector)
            [2 * i] = edges[i].from;
            VECTOR(edge_vector)
            [2 * i + 1] = edges[i].to;
        }

        igraph_empty(graph, next_node_id, IGRAPH_UNDIRECTED);
        igraph_add_edges(graph, &edge_vector, NULL);
        igraph_vector_int_destroy(&edge_vector);

        munmap(mapped, file_size);
        close(fd);
        return 0;
    }

    static int LoadClusteringMMap(const std::string &filename,
                                  std::unordered_map<int, int> &node_to_cluster);

private:
    static inline int ParseInt(const char *&ptr, const char *end)
    {
        int result = 0;
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