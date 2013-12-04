#ifndef LCP_MERGESORT_H_
#define LCP_MERGESORT_H_

#include "../utils/types.h"

namespace eberle_mergesort_lcp {

using namespace types;

//typedefs
typedef unsigned char* string;
typedef unsigned int UINT;

//implementation follows

static inline
void eberle_lcp_merge(AS* input1, size_t length1, AS* input2, size_t length2,
		AS* output) {
	const AS* end1 = input1 + length1;
	const AS* end2 = input2 + length2;
	AS* outStream = output;

	//do the merge
	while (input1 < end1 && input2 < end2) {
		AS* a = input1;
		AS* b = input2;

		if (a->lcp == b->lcp) { // CASE 1 lcps are equal
			string s1 = a->text + a->lcp;
			string s2 = b->text + a->lcp;

			// check the strings starting after lcp and calculate new lcp
			while (*s1 != '\0' && *s1 == *s2)
				s1++, s2++;

			const UINT lcp = s1 - a->text;

			if (*s1 <= *s2) { 	// CASE 1.1: a <= b
				*outStream = *a;
				input1++;
				b->lcp = lcp;

			} else { 			// CASE 1.2: a > b
				*outStream = *b;
				input2++;
				a->lcp = lcp;
			}

		} else if (a->lcp < b->lcp) { // CASE 2: a > b
			*outStream = *b;
			input2++;

		} else { // CASE 3: a < b
			*outStream = *a;
			input1++;
		}

		outStream++;
	}

	if (input1 < end1) { // if there are remaining elements in stream1, copy them to the end
		memcpy(outStream, input1, (end1 - input1) * sizeof(AS));
	} else {
		memcpy(outStream, input2, (end2 - input2) * sizeof(AS));
	}
}

static inline
void eberle_lcp_merge(AS* input1, size_t length1, AS* input2, size_t length2,
		string* output) {
	const AS* end1 = input1 + length1;
	const AS* end2 = input2 + length2;

	//do the merge
	while (input1 < end1 && input2 < end2) {
		AS* a = input1;
		AS* b = input2;

		if (a->lcp == b->lcp) { // CASE 1 lcps are equal
			string s1 = a->text + a->lcp;
			string s2 = b->text + a->lcp;

			// check the strings starting after lcp and calculate new lcp
			while (*s1 != '\0' && *s1 == *s2)
				s1++, s2++;

			const UINT lcp = s1 - a->text;

			if (*s1 <= *s2) { 	// CASE 1.1: a <= b
				*output = a->text;
				input1++;
				b->lcp = lcp;

			} else { 			// CASE 1.2: a > b
				*output = b->text;
				input2++;
				a->lcp = lcp;
			}

		} else if (a->lcp < b->lcp) { // CASE 2: a > b
			*output = b->text;
			input2++;

		} else { // CASE 3: a < b
			*output = a->text;
			input1++;
		}

		output++;
	}

	if (input1 < end1) { // if there are remaining elements in stream1, copy them to the end
		//memcpy(output, input1, (end1 - input1) * sizeof(AS));
		for (; input1 < end1; input1++, output++) {
			*output = input1->text;
		}
	} else {
		//memcpy(output, input2, (end2 - input2) * sizeof(AS));
		for (; input2 < end2; input2++, output++) {
			*output = input2->text;
		}
	}
}

static inline
void eberle_lcp_mergesort(string *strings, AS *tmp, AS *output, size_t length) {

	if (length <= 1) {
		output->lcp = 0;
		output->text = *strings;
		return;
	}

	const size_t length1 = length / 2;
	const size_t length2 = length - length1;
	AS * tmp2 = tmp + length1;

	eberle_lcp_mergesort(strings, output, tmp, length1);
	eberle_lcp_mergesort(strings + length1, output + length1, tmp2, length2);

	eberle_lcp_merge(tmp, length1, tmp2, length2, output);
}

void eberle_lcp_mergesort(string *strings, size_t n) {
	//allocate memory for annotated strings
	AS *tmp = static_cast<AS *>(malloc(n * sizeof(AS)));
	AS *output = static_cast<AS *>(malloc(n * sizeof(AS)));

	eberle_lcp_mergesort(strings, tmp, output, n);

	for (size_t i = 0; i < n; i++) {
		strings[i] = output[i].text;
	}

	free(tmp);
	free(output);
}

CONTESTANT_REGISTER(eberle_lcp_mergesort, "eberle/mergesort_lcp_binary", "Binary Mergesort with LCP-usage by Andreas Eberle")

}
		// namespace eberle_lcp_mergesort

#endif // LCP_MERGESORT_H_

