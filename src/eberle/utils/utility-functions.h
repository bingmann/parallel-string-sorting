#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

namespace eberle_utils {

template<unsigned K>
static inline
void calculateRanges(boost::array<std::pair<size_t, size_t>, K> &ranges, size_t length) {
	const size_t split = length / K;
	for (unsigned i = 0; i < K - 1; ++i) {
		ranges[i] = std::make_pair(i * split, split);
	}
	ranges[K - 1] = std::make_pair((K - 1) * split, length - (K - 1) * split);
}

}
// namespace eberle_utils

#endif // UTILITY_FUNCTIONS_H_
