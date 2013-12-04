#ifndef EBERLE_LCP_LOSERTREE_MERGESORT_H_
#define EBERLE_LCP_LOSERTREE_MERGESORT_H_

#include <iostream>
#include "../utils/types.h"
#include "../utils/utility-functions.h"
#include "../utils/lcp-string-losertree.h"
#include "../utils/eberle-inssort-lcp.h"
#include "../utils/verification-functions.h"
#include "eberle-lcp-mergesort.h"

//#define EBERLE_LCP_LOSERTREE_MERGESORT_CHECK_LCPS

namespace eberle_mergesort {

using namespace eberle_lcp_utils;
using namespace types;

typedef unsigned char* string;
typedef unsigned int UINT;

template<unsigned K>
static inline
void eberle_mergesort_losertree_lcp_kway(string* strings, AS* tmp, AS* output,
		size_t length) {

	if (length <= 2 * K) {
		//	return eberle_mergesort_lcp::eberle_lcp_mergesort(strings, tmp, output, length);
		return eberle_inssort_lcp::inssort_lcp(strings, output, length);
	}

	//create ranges of the parts
	pair < size_t, size_t > ranges[K];
	eberle_utils::calculateRanges(ranges, K, length);

	// execute mergesorts for parts
	for (unsigned i = 0; i < K; i++) {
		const size_t offset = ranges[i].first;
		eberle_mergesort_losertree_lcp_kway<K>(strings + offset,
				output + offset, tmp + offset, ranges[i].second);
	}

	//merge
	LcpStringLoserTree<K> *loserTree = new LcpStringLoserTree<K>(tmp, ranges);
	loserTree->writeElementsToStream(output, length);
	delete loserTree;
}

template<unsigned K> // K must be a power of two
static inline
void eberle_mergesort_losertree_lcp_kway(string *strings, size_t n) {
	AS *tmp = static_cast<AS *>(malloc(n * sizeof(AS)));
	AS *output = static_cast<AS *>(malloc(n * sizeof(AS)));

	eberle_mergesort_losertree_lcp_kway<K>(strings, tmp, output, n);

#ifdef EBERLE_LCP_LOSERTREE_MERGESORT_CHECK_LCPS
	//check lcps
	eberle_utils::checkLcps(output, n, 0);
#endif //EBERLE_LCP_LOSERTREE_MERGESORT_CHECK_LCPS

	for (size_t i = 0; i < n; i++) {
		strings[i] = output[i].text;
	}

	free(tmp);
	free(output);
}

void eberle_mergesort_losertree_lcp_4way(string *strings, size_t n) {
	eberle_mergesort_losertree_lcp_kway<4>(strings, n);
}

void eberle_mergesort_losertree_lcp_16way(string *strings, size_t n) {
	eberle_mergesort_losertree_lcp_kway<16>(strings, n);
}

void eberle_mergesort_losertree_lcp_32way(string *strings, size_t n) {
	eberle_mergesort_losertree_lcp_kway<32>(strings, n);
}

void eberle_mergesort_losertree_lcp_64way(string *strings, size_t n) {
	eberle_mergesort_losertree_lcp_kway<64>(strings, n);
}

void eberle_mergesort_losertree_lcp_128way(string *strings, size_t n) {
	eberle_mergesort_losertree_lcp_kway<128>(strings, n);
}

void eberle_mergesort_losertree_lcp_512way(string *strings, size_t n) {
	eberle_mergesort_losertree_lcp_kway<512>(strings, n);
}

void eberle_mergesort_losertree_lcp_1024way(string *strings, size_t n) {
	eberle_mergesort_losertree_lcp_kway<1024>(strings, n);
}

CONTESTANT_REGISTER(eberle_mergesort_losertree_lcp_4way, "eberle/mergesort_losertree_lcp_4way", "Mergesort with lcp aware Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_lcp_16way, "eberle/mergesort_losertree_lcp_16way", "Mergesort with lcp aware Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_lcp_32way, "eberle/mergesort_losertree_lcp_32way", "Mergesort with lcp aware Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_lcp_64way, "eberle/mergesort_losertree_lcp_64way", "Mergesort with lcp aware Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_lcp_128way, "eberle/mergesort_losertree_lcp_128way", "Mergesort with lcp aware Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_lcp_512way, "eberle/mergesort_losertree_lcp_512way", "Mergesort with lcp aware Losertree by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_losertree_lcp_1024way, "eberle/mergesort_losertree_lcp_1024way", "Mergesort with lcp aware Losertree by Andreas Eberle")

}
 // namespace eberle_mergesort

#endif // EBERLE_LCP_LOSERTREE_MERGESORT_H_
