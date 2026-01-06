#include "constrained.h"
#include "igraph.h"

void build_displacements_output_file(int * displacements, int* size_array, int nprocs) {
    displacements[0] = 0;
    for (int i = 1; i < nprocs; i ++) {
        displacements[i] = displacements[i-1] + size_array[i-1];
    }
}

int update_cluster_id_array(int * cluster_ids, int cluster_size, int previous_cluster_id, int* cluster_displacements, int nprocs)  {
  int current_cluster = 0;
  int updated_cluster = previous_cluster_id;
  int count = 1;
  for (int i = 0; i < cluster_size; i++) {
    if (cluster_ids[i] > current_cluster) {
        current_cluster = cluster_ids[i];
        updated_cluster++;
    } else if (i != 0  && nprocs > 1 && i == cluster_displacements[count]) {
        count ++;
        updated_cluster++;
        current_cluster = 0;
        // current_cluster = updated_cluster;
    }
    cluster_ids[i] =  updated_cluster;
  }
  if (cluster_size > 0) {
    updated_cluster++;
  }
  return updated_cluster;
}

void ConstrainedClustering::WriteClusterQueue(std::queue<std::vector<long>>& cluster_queue) {
    std::ofstream clustering_output(this->output_file);
    int current_cluster_id = 0;
    this->WriteToLogFile("final clusters:", Log::debug);

    while(!cluster_queue.empty()) {
        std::vector<long> current_cluster = cluster_queue.front();
        cluster_queue.pop();
        this->WriteToLogFile("new cluster size: " + to_string(current_cluster.size()), Log::debug);
        for(size_t i = 0; i < current_cluster.size(); i ++) {
            this->WriteToLogFile(to_string(current_cluster[i]), Log::debug);
            clustering_output << current_cluster[i] << " " << current_cluster_id << '\n';
        }
        current_cluster_id ++;
    }
    clustering_output.close();
    
}

int ConstrainedClustering::WriteToLogFile(std::string message, Log message_type) {
    if(this->log_level >= message_type) {
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        std::string log_message_prefix;
        if(message_type == Log::info) {
            log_message_prefix = "[INFO]";
        } else if(message_type == Log::debug) {
            log_message_prefix = "[DEBUG]";
        } else if(message_type == Log::error) {
            log_message_prefix = "[ERROR]";
        }
        auto days_elapsed = std::chrono::duration_cast<std::chrono::days>(now - this->start_time);
        auto hours_elapsed = std::chrono::duration_cast<std::chrono::hours>(now - this->start_time - days_elapsed);
        auto minutes_elapsed = std::chrono::duration_cast<std::chrono::minutes>(now - this->start_time - days_elapsed - hours_elapsed);
        auto seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time - days_elapsed - hours_elapsed - minutes_elapsed);
        auto total_seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time);
        log_message_prefix += "[";
        log_message_prefix += std::to_string(days_elapsed.count());
        log_message_prefix += "-";
        log_message_prefix += std::to_string(hours_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(minutes_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(seconds_elapsed.count());
        log_message_prefix += "]";

        log_message_prefix += "(t=";
        log_message_prefix += std::to_string(total_seconds_elapsed.count());
        log_message_prefix += "s)";
        this->log_file_handle << log_message_prefix << " " << message << '\n';

        if(this->num_calls_to_log_write % 10 == 0) {
            std::flush(this->log_file_handle);
        }
        this->num_calls_to_log_write ++;
    }
    return 0;
}
