#include <iostream>
#include "../utils/string-losertree.h"

using namespace eberle_utils;

namespace eberle_mergesort {

typedef unsigned char* string;

template<unsigned K>
static inline
void eberle_mergesort_losertree_kway(string *strings, string *tmp,
		size_t length) {

	if (length <= K) {
		return eberle_mergesort(strings, tmp, length);
	}

	//create ranges of the parts
	const size_t split = size_t(double(length) / double(K));
	boost::array<std::pair<size_t, size_t>, K> ranges;
	for (unsigned i = 0; i < K - 1; ++i) {
		ranges[i] = std::make_pair(i * split, split);
	}
	ranges[K - 1] = std::make_pair((K - 1) * split, length - (K - 1) * split);

	// execute mergesorts for parts
	for (unsigned i = 0; i < K; i++) {
		const size_t offset = ranges[i].first;
		eberle_mergesort_losertree_kway<K>(strings + offset, tmp + offset,
				ranges[i].second);
	}

	//merge
	StringLoserTree<K> *loserTree = new StringLoserTree<K>(strings, ranges);
	loserTree->writeElementsToStream(tmp, length);
	delete loserTree;

	memcpy(strings, tmp, sizeof(string) * length);
}

template<unsigned K> // K must be a power of two
static inline
void eberle_mergesort_losertree_kway(string *strings, size_t n) {
	string *tmp = static_cast<string *>(malloc(n * sizeof(string)));

	eberle_mergesort_losertree_kway<K>(strings, tmp, n);

	free(tmp);
}

void eberle_mergesort_losertree_4way(string *strings, size_t n) {
	eberle_mergesort_losertree_kway<4>(strings, n);
}

void eberle_mergesort_losertree_16way(string *strings, size_t n) {
	eberle_mergesort_losertree_kway<16>(strings, n);
}

void eberle_mergesort_losertree_32way(string *strings, size_t n) {
	eberle_mergesort_losertree_kway<32>(strings, n);
}

void eberle_mergesort_losertree_64way(string *strings, size_t n) {
	eberle_mergesort_losertree_kway<64>(strings, n);
}

void eberle_mergesort_losertree_128way(string *strings, size_t n) {
	eberle_mergesort_losertree_kway<128>(strings, n);
}

void eberle_mergesort_losertree_512way(string *strings, size_t n) {
	eberle_mergesort_losertree_kway<512>(strings, n);
}

void eberle_mergesort_losertree_1024way(string *strings, size_t n) {
	eberle_mergesort_losertree_kway<1024>(strings, n);
}

CONTESTANT_REGISTER(eberle_mergesort_losertree_4way, "eberle/mergesort_losertree_4way", "Mergesort with Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_16way, "eberle/mergesort_losertree_16way", "Mergesort with Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_32way, "eberle/mergesort_losertree_32way", "Mergesort with Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_64way, "eberle/mergesort_losertree_64way", "Mergesort with Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_128way, "eberle/mergesort_losertree_128way", "Mergesort with Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_512way, "eberle/mergesort_losertree_512way", "Mergesort with Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_1024way, "eberle/mergesort_losertree_1024way", "Mergesort with Losertree by Andreas Eberle")

}
 // namespace eberle_mergesort
