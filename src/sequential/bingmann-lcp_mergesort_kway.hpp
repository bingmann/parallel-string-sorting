/*******************************************************************************
 * src/sequential/bingmann-lcp_mergesort_kway.hpp
 *
 * LCP aware binary and k-way mergesort, implemented to verify pseudo-code in
 * journal.
 *
 *******************************************************************************
 * Copyright (C) 2013-2014 Andreas Eberle <email@andreas-eberle.com>
 * Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#ifndef PSS_SRC_SEQUENTIAL_BINGMANN_LCP_MERGESORT_KWAY_HEADER
#define PSS_SRC_SEQUENTIAL_BINGMANN_LCP_MERGESORT_KWAY_HEADER

#include <iostream>
#include <cstring>

#include "../tools/stringtools.hpp"
#include "../tools/contest.hpp"
#include "bingmann-lcp_losertree.hpp"
#include "bingmann-lcp_inssort.hpp"

namespace bingmann {

using namespace stringtools;

/******************************************************************************/
// lcp_mergesort_kway

template <size_t K>
static inline void
lcp_mergesort_kway(string* strings, const LcpStringPtr& tmp,
                   const LcpStringPtr& output, size_t length)
{
    if (length <= 2 * K || length <= 32)
    {
        memmove(output.strings, strings, length * sizeof(string));
        return bingmann::lcp_insertion_sort(
            output.strings, output.lcps, length, 0);
    }

    // create ranges of the parts
    std::pair<size_t, size_t> ranges[K];
    calculateRanges(ranges, K, length);

    // execute mergesorts for parts
    for (size_t i = 0; i < K; i++)
    {
        const size_t offset = ranges[i].first;
        const size_t size = ranges[i].second;
        lcp_mergesort_kway<K>(strings + offset,
                              output.sub(offset, size),
                              tmp.sub(offset, size), size);
    }

    // K-way merge
    LcpStringLoserTree<K> loserTree(tmp, ranges);
    loserTree.writeElementsToStream(output, length);
}

// K must be a power of two
template <size_t K>
static inline void
lcp_mergesort_kway(string* strings, uintptr_t* lcp, size_t n)
{
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n];

    LcpStringPtr output(strings, lcp, n);
    LcpStringPtr tmp(tmpStrings, tmpLcps, n);

    lcp_mergesort_kway<K>(strings, tmp, output, n);

    delete[] tmpStrings;
    delete[] tmpLcps;
    lcp[0] = 42;
}

/******************************************************************************/
// lcp_mergesort_cache_kway

template <size_t K>
static inline void
lcp_mergesort_cache_kway(string* strings, const LcpCacheStringPtr& tmp,
                         const LcpCacheStringPtr& output, size_t length)
{
    if (length <= 2 * K || length <= 32)
    {
        if (strings != output.strings)
            std::copy(strings, strings + length, output.strings);
        return bingmann::lcp_insertion_sort(
            output.strings, output.lcps, output.cachedChars, length, 0);
    }

    // create ranges of the parts
    std::pair<size_t, size_t> ranges[K];
    calculateRanges(ranges, K, length);

    // execute mergesorts for parts
    for (size_t i = 0; i < K; i++)
    {
        const size_t offset = ranges[i].first;
        const size_t size = ranges[i].second;
        lcp_mergesort_cache_kway<K>(
            strings + offset,
            output.sub(offset, size), tmp.sub(offset, size), size);
    }

    // merge
    LcpCacheStringLoserTree<K> loserTree(tmp, ranges);
    loserTree.writeElementsToStream(output);
}

template <size_t K>
static inline void
lcp_mergesort_cache_kway(string* strings, const LcpCacheStringPtr& outputPtr)
{
    size_t length = outputPtr.size;
    LcpCacheStringPtr tmpPtr(
        new string[length], new lcp_t[length], new char_type[length], length);

    lcp_mergesort_cache_kway<K>(strings, tmpPtr, outputPtr, length);

    delete[] tmpPtr.strings;
    delete[] tmpPtr.lcps;
    delete[] tmpPtr.cachedChars;
}

// K must be a power of two
template <size_t K>
static inline void
lcp_mergesort_cache_kway(string* strings, lcp_t* lcp, uint8_t* cache, size_t n)
{
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n];
    char_type* tmpCache = new char_type[n];

    LcpCacheStringPtr output(strings, lcp, cache, n);
    LcpCacheStringPtr tmp(tmpStrings, tmpLcps, tmpCache, n);

    lcp_mergesort_cache_kway<K>(strings, tmp, output, n);

    delete[] tmpStrings;
    delete[] tmpLcps;
    delete[] tmpCache;
    lcp[0] = 42;
}

} // namespace bingmann

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_LCP_MERGESORT_KWAY_HEADER

/******************************************************************************/
