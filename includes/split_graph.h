#ifndef SPLIT_GRAPH_H
#define SPLIT_GRAPH_H

#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <map>
#include <constrained.h>
#include <mmap_graph_loader.h>

class SplitGraph{
    public:
        SplitGraph(std::string edgelist, std::string existing_clustering, int num_partitions, int num_processors, std::string output_header, std::string log_file, int log_level) : edgelist(edgelist), existing_clustering(existing_clustering), num_partitions(num_partitions), num_processors(num_processors), output_header(output_header), log_file(log_file), log_level(log_level){
            if(this->log_level > -1) {
                this->start_time = std::chrono::steady_clock::now();
                this->log_file_handle.open(this->log_file);
            }
            this->num_calls_to_log_write = 0;
        };
        int main();
        virtual ~SplitGraph() {
            if (this -> log_level > -1) {
                this->log_file_handle.close();
            }
        }
        int WriteToLogFile(std::string message, Log message_id);
    protected:
        std::string edgelist;
        std::string existing_clustering;
        int num_partitions;
        std::string output_header;
        std::string log_file;
        int log_level;
        std::chrono::steady_clock::time_point start_time;
        std::ofstream log_file_handle;
        int num_calls_to_log_write;
        int num_processors;
};
    
#endif