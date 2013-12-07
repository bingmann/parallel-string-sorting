#ifndef VERIFICATION_FUNCTIONS_H_
#define VERIFICATION_FUNCTIONS_H_

#include "types.h"
#include "utility-functions.h"

namespace eberle_utils {

using namespace types;
using namespace stringtools;

static inline
bool checkLcps(AS* output, size_t n, unsigned expectedFirstLcp) {
	bool allValid = true;

	if (output[0].lcp != expectedFirstLcp) {
		std::cout << "output[0].lcp " << output[0].lcp << " but should be "
				<< expectedFirstLcp << std::endl;
		allValid = false;
	}

	for (size_t i = 1; i < n; ++i) {
		string s1 = output[i - 1].text, s2 = output[i].text;
		size_t lcp = 0;
		while (*s1 != 0 && *s1 == *s2)
			++lcp, ++s1, ++s2;

		if (lcp != output[i].lcp) {
			std::cout << "output[" << i << "].lcp mismatch " << lcp << " != "
					<< output[i].lcp << std::endl;
			allValid = false;
		}
	}

	if (allValid) {
		std::cout << "All LCPs valid!" << std::endl;
	} else {
		std::cout << "Found invalid LCPS!" << std::endl;
	}

	return allValid;
}

static inline void checkSorting(AS* stream, size_t length) {
	for (size_t i = 1; i < length; i++) {
		if (scmp(stream[i - 1].text, stream[i].text) > 0) {
			std::cout << "SORT ERROR! ( " << stream[i - 1].text << " | "
					<< stream[i].text << " )" << std::endl;
		}
	}
}

} // namespace eberle_utils

#endif // VERIFICATION_FUNCTIONS_H_
