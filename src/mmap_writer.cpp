#include "mmap_writer.h"
#include <iostream>
#include <algorithm>
#include <cstdio>

int MMapWriter::IntToString(int value, char* buffer) {
    if (value == 0) {
        buffer[0] = '0';
        return 1;
    }
    int length = 0;
    int temp = value;

    while (temp > 0) {
        length++;
        temp /= 10;
    }

    int pos = length - 1;
    while (value > 0) {
        buffer[pos--] = '0' + (value%10);
        value /= 10;
    }
    return length;
}

char* MMapWriter::WriteEntry(char* ptr, int node_id, int cluster_id) {
    ptr += IntToString(node_id, ptr);

    *ptr++ = ' ';

    ptr += IntToString(cluster_id, ptr);

    *ptr++ = '\n';

    return ptr;
}

size_t MMapWriter::CalculateBytesNeeded(const std::vector<std::pair <int, int>>& entries) {
    size_t total = 0;
    char temp_buffer[32];

    for (const auto& entry: entries) {
        int len1 = IntToString(entry.first, temp_buffer);
        int len2 = IntToString(entry.second, temp_buffer);
        total += len1 + 1 + len2 + 1;
    }
    return total;
}


                                