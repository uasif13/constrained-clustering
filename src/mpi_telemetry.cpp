#include <mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdint>
#include <climits>

#define MAX_BUF_SIZE 1 << 25
uint64_t nVVBuf[1000];
uint64_t checkBuf[1000];

uint64_t my_rankOps[MAX_BUF_SIZE];
uint64_t iterationOps[MAX_BUF_SIZE];
uint64_t messageOps[MAX_BUF_SIZE];
double totalTimeOps[MAX_BUF_SIZE];
uint64_t totalBytesOps[MAX_BUF_SIZE];

int MPI_Bcast(void *buffer, long long count, MPI_Datatype datatype, int emitter_rank, MPI_Comm communicator, int my_rank, int iteration, int message_id, uint64_t* opCount)
{
  if (count > INT_MAX) {
    std::cerr << "Error: MPI_Bcast count " << count << " exceeds INT_MAX" << std::endl;
    return MPI_ERR_COUNT;
  }
  
  double tstart = MPI_Wtime();
  int size;
  int result = PMPI_Bcast(buffer, (int)count, datatype, emitter_rank, communicator);
  double totalTime = MPI_Wtime() - tstart;
  MPI_Type_size(datatype, &size);
  my_rankOps[opCount[my_rank]] = my_rank;
  iterationOps[opCount[my_rank]] = iteration;
  messageOps[opCount[my_rank]] = message_id;
  totalTimeOps[opCount[my_rank]] = totalTime;
  totalBytesOps[opCount[my_rank]] = (uint64_t)count * size;
  opCount[my_rank]++;
  return result;
}

int MPI_Allgather(const void *buffer_send, long long count_send, MPI_Datatype datatype_send, void *buffer_recv, long long count_recv, MPI_Datatype datatype_recv, MPI_Comm communicator, int my_rank, int iteration, int message_id, uint64_t* opCount)
{
  if (count_send > INT_MAX || count_recv > INT_MAX) {
    std::cerr << "Error: MPI_Allgather count exceeds INT_MAX" << std::endl;
    return MPI_ERR_COUNT;
  }
  
  double tstart = MPI_Wtime();
  int size;
  int result = PMPI_Allgather(buffer_send, (int)count_send, datatype_send, buffer_recv, (int)count_recv, datatype_recv, communicator);
  double totalTime = MPI_Wtime() - tstart;
  MPI_Type_size(datatype_send, &size);
  my_rankOps[opCount[my_rank]] = my_rank;
  iterationOps[opCount[my_rank]] = iteration;
  messageOps[opCount[my_rank]] = message_id;
  totalTimeOps[opCount[my_rank]] = totalTime;
  totalBytesOps[opCount[my_rank]] = (uint64_t)count_send * size;
  opCount[my_rank]++;
  return result;
}

int MPI_Allgatherv(const void *buffer_send, long long count_send, MPI_Datatype datatype_send, void *buffer_recv, const long long* count_recv, const long long* displacements, MPI_Datatype datatype_recv, MPI_Comm communicator, int my_rank, int iteration, int message_id, uint64_t* opCount)
{
  if (count_send > INT_MAX) {
    std::cerr << "Error: MPI_Allgatherv count_send exceeds INT_MAX" << std::endl;
    return MPI_ERR_COUNT;
  }
  
  double tstart = MPI_Wtime();
  int size;
  int result;
  
  // Handle NULL case
  if (count_recv == NULL || displacements == NULL) {
    result = PMPI_Allgatherv(buffer_send, (int)count_send, datatype_send, buffer_recv, NULL, NULL, datatype_recv, communicator);
  } else {
    // Need to convert long long arrays to int arrays for PMPI call
    int nprocs;
    PMPI_Comm_size(communicator, &nprocs);
    
    int* count_recv_int = new int[nprocs];
    int* displacements_int = new int[nprocs];
    
    for (int i = 0; i < nprocs; i++) {
      if (count_recv[i] > INT_MAX || displacements[i] > INT_MAX) {
        std::cerr << "Error: MPI_Allgatherv array values exceed INT_MAX" << std::endl;
        delete[] count_recv_int;
        delete[] displacements_int;
        return MPI_ERR_COUNT;
      }
      count_recv_int[i] = (int)count_recv[i];
      displacements_int[i] = (int)displacements[i];
    }
    
    result = PMPI_Allgatherv(buffer_send, (int)count_send, datatype_send, buffer_recv, count_recv_int, displacements_int, datatype_recv, communicator);
    
    delete[] count_recv_int;
    delete[] displacements_int;
  }
  
  double totalTime = MPI_Wtime() - tstart;
  MPI_Type_size(datatype_send, &size);
  my_rankOps[opCount[my_rank]] = my_rank;
  iterationOps[opCount[my_rank]] = iteration;
  messageOps[opCount[my_rank]] = message_id;
  totalTimeOps[opCount[my_rank]] = totalTime;
  totalBytesOps[opCount[my_rank]] = (uint64_t)count_send * size;
  opCount[my_rank]++;
  
  return result;
}

int MPI_Gather(const void *buffer_send, long long count_send, MPI_Datatype datatype_send, void *buffer_recv, long long count_recv, MPI_Datatype datatype_recv, int root, MPI_Comm communicator, int my_rank, int iteration, int message_id, uint64_t* opCount)
{
  if (count_send > INT_MAX || count_recv > INT_MAX) {
    std::cerr << "Error: MPI_Gather count exceeds INT_MAX" << std::endl;
    return MPI_ERR_COUNT;
  }
  
  double tstart = MPI_Wtime();
  int size;
  int result = PMPI_Gather(buffer_send, (int)count_send, datatype_send, buffer_recv, (int)count_recv, datatype_recv, root, communicator);
  double totalTime = MPI_Wtime() - tstart;
  MPI_Type_size(datatype_send, &size);
  my_rankOps[opCount[my_rank]] = my_rank;
  iterationOps[opCount[my_rank]] = iteration;
  messageOps[opCount[my_rank]] = message_id;
  totalTimeOps[opCount[my_rank]] = totalTime;
  totalBytesOps[opCount[my_rank]] = (uint64_t)count_send * size;
  opCount[my_rank]++;
  return result;
}

int MPI_Gatherv(const void *buffer_send, long long count_send, MPI_Datatype datatype_send, void *buffer_recv, const long long* count_recv, const long long* displacements, MPI_Datatype datatype_recv, int root, MPI_Comm communicator, int my_rank, int iteration, int message_id, uint64_t* opCount)
{
  if (count_send > INT_MAX) {
    std::cerr << "Error: MPI_Gatherv count_send exceeds INT_MAX" << std::endl;
    return MPI_ERR_COUNT;
  }
  
  double tstart = MPI_Wtime();
  int size;
  int result;
  
  // Handle NULL case (non-root ranks)
  if (count_recv == NULL || displacements == NULL) {
    result = PMPI_Gatherv(buffer_send, (int)count_send, datatype_send, buffer_recv, NULL, NULL, datatype_recv, root, communicator);
  } else {
    // Need to convert long long arrays to int arrays for PMPI call
    int nprocs;
    PMPI_Comm_size(communicator, &nprocs);
    
    int* count_recv_int = new int[nprocs];
    int* displacements_int = new int[nprocs];
    
    for (int i = 0; i < nprocs; i++) {
      if (count_recv[i] > INT_MAX || displacements[i] > INT_MAX) {
        std::cerr << "Error: MPI_Gatherv array values exceed INT_MAX" << std::endl;
        delete[] count_recv_int;
        delete[] displacements_int;
        return MPI_ERR_COUNT;
      }
      count_recv_int[i] = (int)count_recv[i];
      displacements_int[i] = (int)displacements[i];
    }
    
    result = PMPI_Gatherv(buffer_send, (int)count_send, datatype_send, buffer_recv, count_recv_int, displacements_int, datatype_recv, root, communicator);
    
    delete[] count_recv_int;
    delete[] displacements_int;
  }
  
  double totalTime = MPI_Wtime() - tstart;
  MPI_Type_size(datatype_send, &size);
  my_rankOps[opCount[my_rank]] = my_rank;
  iterationOps[opCount[my_rank]] = iteration;
  messageOps[opCount[my_rank]] = message_id;
  totalTimeOps[opCount[my_rank]] = totalTime;
  totalBytesOps[opCount[my_rank]] = (uint64_t)count_send * size;
  opCount[my_rank]++;
  
  return result;
}

int MPI_Barrier(int my_rank, int iteration, int message_id, uint64_t* opCount)
{
  double tstart = MPI_Wtime();
  int result = PMPI_Barrier(MPI_COMM_WORLD);
  double totalTime = MPI_Wtime() - tstart;
  my_rankOps[opCount[my_rank]] = my_rank;
  iterationOps[opCount[my_rank]] = iteration;
  messageOps[opCount[my_rank]] = message_id;
  totalTimeOps[opCount[my_rank]] = totalTime;
  totalBytesOps[opCount[my_rank]] = 0;
  opCount[my_rank]++;
  return result;
}

int MPI_Finalize(int my_rank, int nprocs, uint64_t* opCount, std::string mpi_log_file_path)
{
  std::string opInfoArr[6] = {"MPI_Bcast", "MPI_Allgather", "MPI_Allgatherv", "MPI_Gather", "MPI_Gatherv", "MPI_Barrier"};
  // aggregate wall time, bytes transferred, op type
  uint64_t opCountArr[nprocs];
  MPI_Allgather(&opCount[my_rank], 1, MPI_UINT64_T, opCountArr, 1, MPI_UINT64_T, MPI_COMM_WORLD, my_rank, -1, 1, opCount);
  uint64_t max_ops = 0;
  for (int i = 0; i < nprocs; i++) {
    if (opCountArr[i] > max_ops)
      max_ops = opCountArr[i];
  }

  uint64_t* my_rankOpsAgg = new uint64_t[max_ops * nprocs];
  int* iterationOpsAgg = new int[max_ops * nprocs];
  uint64_t* messageOpsAgg = new uint64_t[max_ops * nprocs];
  double* totalTimeOpsAgg = new double[max_ops * nprocs];
  uint64_t* totalBytesOpsAgg = new uint64_t[max_ops * nprocs];

  PMPI_Gather(my_rankOps, max_ops, MPI_UINT64_T, my_rankOpsAgg, max_ops, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  PMPI_Gather(iterationOps, max_ops, MPI_INT, iterationOpsAgg, max_ops, MPI_INT, 0, MPI_COMM_WORLD);
  PMPI_Gather(messageOps, max_ops, MPI_UINT64_T, messageOpsAgg, max_ops, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  PMPI_Gather(totalTimeOps, max_ops, MPI_DOUBLE, totalTimeOpsAgg, max_ops, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  PMPI_Gather(totalBytesOps, max_ops, MPI_UINT64_T, totalBytesOpsAgg, max_ops, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    double totalMpiTime[6] = {0, 0, 0, 0, 0, 0};
    uint64_t totalBytesTransferredTime[6] = {0, 0, 0, 0, 0, 0};
    std::ofstream mpi_log_file(mpi_log_file_path);
    uint64_t log_index;
    std::string log_line;
    for (int i = 0; i < nprocs; i++) {
      for (uint64_t j = 0; j < opCountArr[i]; j++) {
        log_index = i * max_ops + j;
        log_line = "my_rank: " + std::to_string(my_rankOpsAgg[log_index]) + " iteration: " + std::to_string(iterationOpsAgg[log_index]) + " MPI_Routine: " + opInfoArr[messageOpsAgg[log_index]] + " WallTime: " + std::to_string(totalTimeOpsAgg[log_index]) + " Bytes Transferred: " + std::to_string(totalBytesOpsAgg[log_index]) + "\n";
        mpi_log_file << log_line;
        if (messageOpsAgg[log_index] == 0) {
          totalMpiTime[0] += totalTimeOpsAgg[log_index];
          totalBytesTransferredTime[0] += totalBytesOpsAgg[log_index];
        }
        else if (messageOpsAgg[log_index] == 1) {
          totalMpiTime[1] += totalTimeOpsAgg[log_index];
          totalBytesTransferredTime[1] += totalBytesOpsAgg[log_index];
        }
        else if (messageOpsAgg[log_index] == 2) {
          totalMpiTime[2] += totalTimeOpsAgg[log_index];
          totalBytesTransferredTime[2] += totalBytesOpsAgg[log_index];
        }
        else if (messageOpsAgg[log_index] == 3) {
          totalMpiTime[3] += totalTimeOpsAgg[log_index];
          totalBytesTransferredTime[3] += totalBytesOpsAgg[log_index];
        }
        else if (messageOpsAgg[log_index] == 4) {
          totalMpiTime[4] += totalTimeOpsAgg[log_index];
          totalBytesTransferredTime[4] += totalBytesOpsAgg[log_index];
        }
        else {
          totalMpiTime[5] += totalTimeOpsAgg[log_index];
          totalBytesTransferredTime[5] += totalBytesOpsAgg[log_index];
        }
      }
    }
    mpi_log_file << "---------------------Cumulative Time and Bytes Transferred by Routine Type--------------------\n";
    mpi_log_file << "Wall Time\n";
    log_line = "MPI_Bcast: " + std::to_string(totalMpiTime[0]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Allgather: " + std::to_string(totalMpiTime[1]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Allgatherv: " + std::to_string(totalMpiTime[2]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Gather: " + std::to_string(totalMpiTime[3]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Gatherv: " + std::to_string(totalMpiTime[4]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Barrier: " + std::to_string(totalMpiTime[5]) + "\n";
    mpi_log_file << log_line;
    mpi_log_file << "Bytes Transferred\n";
    log_line = "MPI_Bcast: " + std::to_string(totalBytesTransferredTime[0]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Allgather: " + std::to_string(totalBytesTransferredTime[1]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Allgatherv: " + std::to_string(totalBytesTransferredTime[2]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Gather: " + std::to_string(totalBytesTransferredTime[3]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Gatherv: " + std::to_string(totalBytesTransferredTime[4]) + "\n";
    mpi_log_file << log_line;
    log_line = "MPI_Barrier: " + std::to_string(totalBytesTransferredTime[5]) + "\n";
    mpi_log_file << log_line;
    mpi_log_file.close();
  }
  
  delete[] my_rankOpsAgg;
  delete[] iterationOpsAgg;
  delete[] messageOpsAgg;
  delete[] totalTimeOpsAgg;
  delete[] totalBytesOpsAgg;
  
  int result = PMPI_Finalize();
  return result;
}