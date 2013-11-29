#ifndef BINGMANN_PARALLEL_SAMPLE_SORTBTCE_H_
#define BINGMANN_PARALLEL_SAMPLE_SORTBTCE_H_

#include "../tools/stringtools.h"

namespace bingmann_parallel_sample_sortBTCE {

void parallel_sample_sortBTCEU_numa(stringtools::string * strings, size_t n, int numaNode,
		int numberOfThreads);

void parallel_sample_sortBTCEU1(stringtools::string* strings, size_t n);

} // namespace bingmann_parallel_sample_sortBTCE

#endif // BINGMANN_PARALLEL_SAMPLE_SORTBTCE_H_
