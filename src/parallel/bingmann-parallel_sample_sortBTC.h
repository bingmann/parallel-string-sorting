#ifndef BINGMANN_PARALLEL_SAMPLE_SORTBTC_H_
#define BINGMANN_PARALLEL_SAMPLE_SORTBTC_H_

#include "../tools/stringtools.h"

namespace bingmann_parallel_sample_sortBTC {

void parallel_sample_sortBTC_numa(
    stringtools::string * strings, size_t n, size_t depth,
    int numaNode, int numberOfThreads
    );

} // namespace bingmann_parallel_sample_sortBTC

#endif // BINGMANN_PARALLEL_SAMPLE_SORTBTC_H_
