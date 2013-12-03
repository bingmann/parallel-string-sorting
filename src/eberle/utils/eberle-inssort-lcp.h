#ifndef EBERLE_INSSORT_H_
#define EBERLE_INSSORT_H_

#include <stdlib.h>
#include "types.h"
#include "verification-functions.h"

namespace eberle_inssort_lcp {

using namespace types;

//typedefs
typedef unsigned char* string;
typedef unsigned int UINT;

//implementations

static inline
void inssort_lcp(string* strings, AS* output, size_t length) {

	for (size_t n = 0; n < length; n++) {
		const string candidateText = strings[n];
		unsigned candidateLcp = 0;

		size_t insIdx = 0;
		for (; insIdx < n; insIdx++) {
			AS* curr = output + insIdx;

			if (candidateLcp == curr->lcp) { // CASE 1 lcps are equal
				string s1 = candidateText + candidateLcp;
				string s2 = curr->text + candidateLcp;

				// check the strings starting after lcp and calculate new lcp
				while (*s1 != '\0' && *s1 == *s2)
					s1++, s2++;

				const UINT lcp = s1 - candidateText;

				if (*s1 <= *s2) { 	// CASE 1.1: candidate <= curr => insert it
					curr->lcp = lcp;
					break;
				} else { 			// CASE 1.2: candidate > curr
					candidateLcp = lcp;
				}

			} else if (candidateLcp > curr->lcp) { // CASE 2: candidate < curr => insert it
				break;

			} //  CASE 3: candidate > curr => nothing to do
		}

		// move the tail one back
		for (size_t i = n; i > insIdx; i--) {
			output[i] = output[i - 1];
		}

		// insert the new element
		output[insIdx].lcp = candidateLcp;
		output[insIdx].text = candidateText;
	}
}

static inline
void eberle_lcp_inssort(string *strings, size_t n) {
	AS *output = static_cast<AS *>(malloc(n * sizeof(AS)));

	inssort_lcp(strings, output, n);

	//check lcps
	eberle_utils::checkLcps(output, n, 0);

	for (size_t i = 0; i < n; i++) {
		strings[i] = output[i].text;
	}

	free(output);
}

CONTESTANT_REGISTER(eberle_lcp_inssort, "eberle/lcp_inssort", "LCP aware inssertion sort by Andreas Eberle")

}
 // namespace eberle_inssort_lcp

#endif // EBERLE_INSSORT_H_
