#ifndef MMAP_WRITER_H
#define MMAP_WRITER_H

#include <igraph.h>
#include "mpi_telemetry.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <queue>
#include <vector>
#include <cstring>
#include <iostream>

class MMapWriter
{
public:
    static int WriteDistributed(std::queue<std::vector<int>> *cluster_queue,
                                igraph_t *graph,
                                const std::string &output_file,
                                int cluster_id_start,
                                int previous_cluster_id,
                                int my_rank,
                                int nprocs)
    {
        std::vector<std::pair<int, int>> my_entries;
        PrepareEntries(cluster_queue, graph, cluster_id_start, my_entries);

        int my_cluster_count = cluster_id_start > 0 ? (my_entries.empty() ? 0 : my_entries.back().second - cluster_id_start + 1) : 0;

        size_t my_bytes = CalculateBytesNeeded(my_entries);

        std::vector<size_t> all_bytes(nprocs);
        MPI_Allgather(&my_bytes, sizeof(size_t), MPI_BYTE, all_bytes.data(), sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

        std::vector<int> all_cluster_counts(nprocs);
        MPI_Allgather(&my_cluster_count, 1, MPI_INT, all_cluster_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        size_t my_offset = 0;
        int my_cluster_offset = previous_cluster_id;
        size_t total_bytes = 0;

        for (int i = 0; i < nprocs; i++)
        {
            if (i < my_rank)
            {
                my_offset += all_bytes[i];
                my_cluster_offset += all_cluster_counts[i];
            }
            total_bytes += all_bytes[i];
        }
        for (auto &entry : my_entries)
        {
            entry.second += my_cluster_offset;
        }

        if (my_rank == 0)
        {
            int fd = open(output_file.c_str(), O_RDWR | O_CREAT, 0644);
            if (fd == -1)
            {
                std::cerr << "Rank 0: Failed to open/create file: " << output_file << std::endl;
                PMPI_Abort(MPI_COMM_WORLD, 1);
            }
            struct stat st;
            fstat(fd, &st);
            size_t current_size = st.st_size;

            size_t new_size = current_size + total_bytes;
            if (ftruncate(fd, new_size) == -1)
            {
                std::cerr << "Rank 0: Failed to extend file" << std::endl;
                close(fd);
                PMPI_Abort(MPI_COMM_WORLD, 1);
            }
            close(fd);
        }

        PMPI_Barrier(MPI_COMM_WORLD);

        int final_cluster_id = my_cluster_offset;

        if (my_bytes > 0)
        {
            int fd = open(output_file.c_str(), O_RDWR);
            if (fd == -1)
            {
                std::cerr << "Rank " << my_rank << ": Failed to open file for writing" << std::endl;
                PMPI_Abort(MPI_COMM_WORLD, 1);
            }

            struct stat st;
            fstat(fd, &st);
            size_t file_base = st.st_size - total_bytes;
            size_t actual_offset = file_base + my_offset;

            char *mapped = (char *)mmap(NULL, my_bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, actual_offset);

            if (mapped == MAP_FAILED)
            {
                std::cerr << "Rank " << my_rank << ": mmap failed" << std::endl;
                close(fd);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            close(fd);

            char *ptr = mapped;
            for (const auto &entry : my_entries)
            {
                ptr = WriteEntry(ptr, entry.first, entry.second);
            }

            msync(mapped, my_bytes, MS_SYNC);
            munmap(mapped, my_bytes);

            if (!my_entries.empty())
                final_cluster_id = my_entries.back().second + 1;
        }

        PMPI_Barrier(MPI_COMM_WORLD);

        int global_final_cluster_id;
        if (my_rank == nprocs - 1)
        {
            global_final_cluster_id = final_cluster_id;
        }
        PMPI_Bcast(&global_final_cluster_id, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);

        return global_final_cluster_id;
    };

private:
    static size_t CalculateBytesNeeded(const std::vector<std::pair<int, int>> &entries);

    static int IntToString(int value, char *buffer);

    static char *WriteEntry(char *ptr, int node_id, int cluster_id);

    static void PrepareEntries(std::queue<std::vector<int>> *cluster_queue,
                                           igraph_t *graph,
                                           int cluster_id_start,
                                           std::vector<std::pair<int, int>> &entries)
    {
        int current_cluster_id = cluster_id_start;
        while (!cluster_queue->empty())
        {
            std::vector<int> current_cluster = cluster_queue->front();
            cluster_queue->pop();

            for (size_t i = 0; i < current_cluster.size(); i++)
            {
                int node_id = std::stoi(VAS(graph, "name", current_cluster[i]));
                entries.push_back({node_id, current_cluster_id});
            }
            current_cluster_id++;
        }
    }
};
#endif