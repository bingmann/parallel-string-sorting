#ifndef VERIFICATION_FUNCTIONS_H_
#define VERIFICATION_FUNCTIONS_H_

#include "types.h"
#include "../../tools/stringtools.h"

namespace eberle_utils
{

using namespace types;

static inline
bool
checkLcps(AS* output, size_t n, lcp_t expectedFirstLcp)
{
    bool allValid = true;

    if (output[0].lcp != expectedFirstLcp)
    {
        std::cout << "output[0].lcp " << output[0].lcp << " but should be " << expectedFirstLcp << std::endl;
        allValid = false;
    }

    for (size_t i = 1; i < n; ++i)
    {
        string s1 = output[i - 1].text, s2 = output[i].text;
        lcp_t lcp = 0;
        while (*s1 != 0 && *s1 == *s2)
            ++lcp, ++s1, ++s2;

        if (lcp != output[i].lcp)
        {
            std::cout << "output[" << i << "].lcp mismatch " << lcp << " != " << output[i].lcp << std::endl;
            allValid = false;
        }
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

static inline void
checkSorting(AS* stream, size_t length)
{
    for (size_t i = 1; i < length; i++)
    {
        if (stringtools::scmp(stream[i - 1].text, stream[i].text) > 0)
        {
            std::cout << "SORT ERROR! ( " << stream[i - 1].text << " | " << stream[i].text << " )" << std::endl;
        }
    }
}

} // namespace eberle_utils

#endif // VERIFICATION_FUNCTIONS_H_
