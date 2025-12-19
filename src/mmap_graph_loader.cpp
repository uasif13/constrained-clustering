#include "mmap_graph_loader.h"
#include <iostream>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <vector>

size_t MMapGraphLoader::EstimateEdgeCount(const char *data, size_t size)
{
    size_t line_count = 0;
    const char *ptr = data;
    const char *end = data + size;

    size_t sample_size = std::min(size, (size_t)(1048576));
    const char *sample_end = data + sample_size;

    while (ptr < sample_end)
    {
        if (*ptr == '\n')
            line_count++;
        ptr++;
    }

    return (line_count * size) / sample_size;
}