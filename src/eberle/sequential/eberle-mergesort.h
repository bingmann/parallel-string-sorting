#ifndef EBERLE_MERGESORT_H_
#define EBERLE_MERGESORT_H_

namespace eberle_mergesort {

typedef unsigned char* string;

static inline
int scmp(string s1, string s2) {
	while (*s1 != '\0' && *s1 == *s2)
		s1++, s2++;
	return (*s1 - *s2);
}

static inline
void eberle_merge(string *strings, string *tmp, size_t length1,
		size_t length2) {

	size_t idx1 = 0;
	const size_t end1 = length1;
	size_t idx2 = end1;
	const size_t end2 = length1 + length2;

	size_t mergedIdx = 0;

	for (; idx1 < end1 && idx2 < end2; mergedIdx++) {
		if (scmp(strings[idx1], strings[idx2]) < 0) {
			tmp[mergedIdx] = strings[idx1];
			idx1++;
		} else {
			tmp[mergedIdx] = strings[idx2];
			idx2++;
		}
	}

	if (idx1 < end1) { // add the rest of part 1 to the end.
		memcpy(strings + mergedIdx, strings + idx1,
				(end1 - idx1) * sizeof(string));
	}

	// copy the merged part from tmp to strings
	memcpy(strings, tmp, mergedIdx * sizeof(string));
}

static inline
void eberle_mergesort(string *strings, string *tmp, size_t length) {

	if (length > 1) {
		const size_t length1 = length / 2;
		const size_t length2 = length - length1;

		eberle_mergesort(strings, tmp, length1);
		eberle_mergesort(strings + length1, tmp + length1, length2);

		eberle_merge(strings, tmp, length1, length2);
	}
}

void eberle_mergesort(string *strings, size_t n) {
	string *tmp = static_cast<string *>(malloc(n * sizeof(string)));

	eberle_mergesort(strings, tmp, n);

	free(tmp);
}

CONTESTANT_REGISTER(eberle_mergesort, "eberle/mergesort", "Mergesort by Andreas Eberle")

template<unsigned K>
static inline
void eberle_mergesort_kway(string *strings, string *tmp, size_t length) {

	if (length <= K) {
		return eberle_mergesort(strings, tmp, length);
	}

	//create ranges of the parts
	pair < size_t, size_t > ranges[K];
	calculateRanges(ranges, K, length);

	// execute mergesorts for parts
	for (unsigned i = 0; i < K; i++) {
		const size_t offset = ranges[i].first;
		eberle_mergesort_kway<K>(strings + offset, tmp + offset, ranges[i].second);
	}

	//merge treewise
	unsigned stepSize = 1;
	while (stepSize < K) {
		for (unsigned i = 0; i < K; i += 2 * stepSize) {
			//			std::cout << "mergeloop: stepSize: " << stepSize << " i: " << i << std::endl;

			std::pair < size_t, size_t > part1 = ranges[i];
			std::pair < size_t, size_t > part2 = ranges[i + stepSize];

			//		std::cout << "merge called with: "<< part1.first <<" " << part1.second <<" " << part2.first <<" " << part2.second << std::endl;

			eberle_merge(strings + part1.first, tmp + part1.first, part1.second, part2.second);

			ranges[i].second = part1.second + part2.second;// increase the length for the next run
		}
		stepSize *= 2;
	}
}

template<unsigned K> // K must be a power of two
static inline
void eberle_mergesort_kway(string *strings, size_t n) {
	string *tmp = static_cast<string *>(malloc(n * sizeof(string)));

	eberle_mergesort_kway<K>(strings, tmp, n);

	free(tmp);
}

void eberle_mergesort_4way(string *strings, size_t n) {
	eberle_mergesort_kway<4>(strings, n);
}

void eberle_mergesort_16way(string *strings, size_t n) {
	eberle_mergesort_kway<16>(strings, n);
}

void eberle_mergesort_32way(string *strings, size_t n) {
	eberle_mergesort_kway<32>(strings, n);
}

void eberle_mergesort_1024way(string *strings, size_t n) {
	eberle_mergesort_kway<1024>(strings, n);
}

CONTESTANT_REGISTER(eberle_mergesort_4way, "eberle/mergesort_4way", "Mergesort by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_16way, "eberle/mergesort_16way", "Mergesort by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_32way, "eberle/mergesort_32way", "Mergesort by Andreas Eberle")
CONTESTANT_REGISTER(eberle_mergesort_1024way, "eberle/mergesort_1024way", "Mergesort by Andreas Eberle")

}
 // namespace eberle_mergesort

#endif // EBERLE_MERGESORT_H_
