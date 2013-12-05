#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <stdlib.h>
#include "types.h"

namespace eberle_utils
{

using namespace types;

static inline
void
calculateRanges(std::pair<size_t, size_t>* ranges, unsigned numberOfSplits, size_t lengthToSplit)
{
    const size_t split = lengthToSplit / numberOfSplits;
    for (unsigned i = 0; i < numberOfSplits - 1; ++i)
    {
        ranges[i] = std::make_pair(i * split, split);
    }
    ranges[numberOfSplits - 1] = std::make_pair((numberOfSplits - 1) * split, lengthToSplit - (numberOfSplits - 1) * split);
}

static inline
int
scmp(string s1, string s2)
{
    while (*s1 != '\0' && *s1 == *s2)
        s1++, s2++;
    return (*s1 - *s2);
}

static inline
unsigned
calculateLcp(string s1, string s2)
{
    unsigned lcp = 0;
    while (*s1 != '\0' && *s1 == *s2)
        s1++, s2++, lcp++;

    return lcp;
}

static inline
unsigned
getNextHigherPowerOfTwo(unsigned v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

}
// namespace eberle_utils

#endif // UTILITY_FUNCTIONS_H_
