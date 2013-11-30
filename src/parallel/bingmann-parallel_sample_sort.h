#ifndef BINGMANN_PARALLEL_SAMPLE_SORT_H_
#define BINGMANN_PARALLEL_SAMPLE_SORT_H_

#include "../tools/stringtools.h"

namespace bingmann_parallel_sample_sort {

void parallel_sample_sort_numa(
    stringtools::string * strings, size_t n, size_t depth,
    int numaNode, int numberOfThreads
    );

} // namespace bingmann_parallel_sample_sort

#endif // BINGMANN_PARALLEL_SAMPLE_SORT_H_
