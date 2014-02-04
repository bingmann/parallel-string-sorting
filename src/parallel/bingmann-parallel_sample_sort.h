#ifndef BINGMANN_PARALLEL_SAMPLE_SORT_H_
#define BINGMANN_PARALLEL_SAMPLE_SORT_H_

#include "../tools/stringtools.h"

namespace bingmann_parallel_sample_sort_lcp {

using namespace stringtools;

void parallel_sample_sort_numa(string *strings, size_t n,
                               int numaNode, int numberOfThreads,
                               const LcpCacheStringPtr& output);

} // namespace bingmann_parallel_sample_sort_lcp

#endif // BINGMANN_PARALLEL_SAMPLE_SORT_H_
