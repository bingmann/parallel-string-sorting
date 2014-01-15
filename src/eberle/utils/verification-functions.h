#ifndef VERIFICATION_FUNCTIONS_H_
#define VERIFICATION_FUNCTIONS_H_

#include "types.h"

namespace eberle_utils
{

using namespace types;

static inline
bool
checkLcps(LcpStringPtr output, size_t n, lcp_t expectedFirstLcp)
{
    bool allValid = true;

    if (output.lcp() != expectedFirstLcp)
    {
        std::cout << "output.lcp() " << output.lcp() << " but should be " << expectedFirstLcp << std::endl;
        allValid = false;
    }

    string lastText = output.str();
    ++output;

    for (size_t i = 1; i < n; ++i, ++output)
    {
        string s1 = lastText, s2 = output.str();
        lcp_t lcp = 0;
        while (*s1 != 0 && *s1 == *s2)
            ++lcp, ++s1, ++s2;

        if (lcp != output.lcp())
        {
            std::cout << "output[" << i << "].lcp mismatch " << lcp << " != " << output.lcp() << std::endl;
            allValid = false;
        }

        lastText = output.str();
    }

    if (allValid)
    {
        std::cout << "All LCPs valid!" << std::endl;
    }
    else
    {
        std::cout << "Found invalid LCPS!" << std::endl;
    }

    return allValid;
}

} // namespace eberle_utils

#endif // VERIFICATION_FUNCTIONS_H_
