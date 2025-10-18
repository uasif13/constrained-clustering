#ifndef MPI_TELEMETRY_H
#define MPI_TELEMETRY_H
#include <mpi.h>
#include <cstdint>



int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int emitter_rank, MPI_Comm communicator, int my_rank, int iteration, int message_id, uint64_t* opCount);
int MPI_Allgather(const void *buffer_send, int count_send, MPI_Datatype datatype_send, void *buffer_recv, int count_recv, MPI_Datatype datatype_recv, MPI_Comm communicator,int my_rank, int iteration, int message_id, uint64_t* opCount);
int MPI_Allgatherv(const void *buffer_send, int count_send, MPI_Datatype datatype_send, void *buffer_recv, const int* count_recv,const int* displacements, MPI_Datatype datatype_recv, MPI_Comm communicator,int my_rank, int iteration, int message_id, uint64_t* opCount);
int MPI_Gather(const void *buffer_send, int count_send, MPI_Datatype datatype_send, void *buffer_recv, int count_recv, MPI_Datatype datatype_recv,int root, MPI_Comm communicator,int my_rank, int iteration, int message_id, uint64_t* opCount);
int MPI_Gatherv(const void *buffer_send, int count_send, MPI_Datatype datatype_send, void *buffer_recv, const int* count_recv, const int* displacements, MPI_Datatype datatype_recv, int root, MPI_Comm communicator,int my_rank, int iteration, int message_id, uint64_t* opCount);
int MPI_Barrier(int my_rank, int iteration, int message_id, uint64_t* opCount);
int MPI_Finalize(int my_rank, int nprocs, uint64_t * opCount, std::string mpi_log_file_path);



#endif