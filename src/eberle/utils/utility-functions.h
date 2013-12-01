#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <stdlib.h>

namespace eberle_utils {

static inline
void calculateRanges(std::pair<size_t, size_t>* ranges, unsigned numberOfSplits,
		size_t lengthToSplit) {
	const size_t split = lengthToSplit / numberOfSplits;
	for (unsigned i = 0; i < numberOfSplits - 1; ++i) {
		ranges[i] = std::make_pair(i * split, split);
	}
	ranges[numberOfSplits - 1] = std::make_pair((numberOfSplits - 1) * split,
			lengthToSplit - (numberOfSplits - 1) * split);
}

}
// namespace eberle_utils

#endif // UTILITY_FUNCTIONS_H_
