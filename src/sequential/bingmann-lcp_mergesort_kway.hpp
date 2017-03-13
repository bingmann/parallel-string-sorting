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
#include "../tools/eberle-utilities.hpp"
#include "../tools/contest.hpp"
#include "bingmann-lcp_losertree.hpp"
#include "bingmann-lcp_inssort.hpp"

namespace bingmann {

using namespace stringtools;

/******************************************************************************/
// lcp_mergesort_kway

template <unsigned K>
static inline void
lcp_mergesort_kway(string* strings, const LcpStringPtr& tmp,
                   const LcpStringPtr& output, size_t length)
{
    if (length <= 2 * K)
    {
        memmove(output.strings, strings, length * sizeof(string));
        return bingmann::lcp_insertion_sort(
            output.strings, output.lcps, length, 0);
    }

    // create ranges of the parts
    std::pair<size_t, size_t> ranges[K];

    const size_t split = length / K;
    for (unsigned i = 0; i < K - 1; ++i)
    {
        ranges[i] = std::make_pair(i * split, split);
    }
    ranges[K - 1] = std::make_pair((K - 1) * split, length - (K - 1) * split);

    // execute mergesorts for parts
    for (unsigned i = 0; i < K; i++)
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
template <unsigned K>
static inline void
lcp_mergesort_kway(string* strings, uintptr_t* lcp, size_t n)
{
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n + 1];

    LcpStringPtr output(strings, lcp, n);
    LcpStringPtr tmp(tmpStrings, tmpLcps, n);

    lcp_mergesort_kway<K>(strings, tmp, output, n);

    delete[] tmpStrings;
    delete[] tmpLcps;
    lcp[0] = 42;
}

static inline void
lcp_mergesort_4way(string* strings, uintptr_t* lcp, size_t n)
{
    lcp_mergesort_kway<4>(strings, lcp, n);
}

static inline void
lcp_mergesort_8way(string* strings, uintptr_t* lcp, size_t n)
{
    lcp_mergesort_kway<8>(strings, lcp, n);
}

static inline void
lcp_mergesort_16way(string* strings, uintptr_t* lcp, size_t n)
{
    lcp_mergesort_kway<16>(strings, lcp, n);
}

static inline void
lcp_mergesort_64way(string* strings, uintptr_t* lcp, size_t n)
{
    lcp_mergesort_kway<64>(strings, lcp, n);
}

PSS_CONTESTANT(lcp_mergesort_4way, "bingmann/lcp_mergesort_4way",
               "4-way LCP-Mergesort by Timo Bingmann and Andreas Eberle")
PSS_CONTESTANT(lcp_mergesort_8way, "bingmann/lcp_mergesort_8way",
               "8-way LCP-Mergesort by Timo Bingmann and Andreas Eberle")
PSS_CONTESTANT(lcp_mergesort_16way, "bingmann/lcp_mergesort_16way",
               "16-way LCP-Mergesort by Timo Bingmann and Andreas Eberle")
PSS_CONTESTANT(lcp_mergesort_64way, "bingmann/lcp_mergesort_64way",
               "64-way LCP-Mergesort by Timo Bingmann and Andreas Eberle")

/******************************************************************************/
// lcp_mergesort_cache_kway

template <unsigned K>
static inline void
lcp_mergesort_cache_kway(string* strings, const LcpCacheStringPtr& tmp,
                         const LcpCacheStringPtr& output, size_t length)
{
    if (length <= 2 * K)
    {
        // return eberle_inssort_lcp::inssort_lcp(strings, output, length);
        if (strings != output.strings)
            std::copy(strings, strings + length, output.strings);
        return bingmann::lcp_insertion_sort(
            output.strings, output.lcps, output.cachedChars, length, 0);
    }

    //create ranges of the parts
    std::pair<size_t, size_t> ranges[K];
    eberle_utils::calculateRanges(ranges, K, length);

    // execute mergesorts for parts
    for (unsigned i = 0; i < K; i++)
    {
        const size_t offset = ranges[i].first;
        const size_t size = ranges[i].second;
        lcp_mergesort_cache_kway<K>(
            strings + offset,
            output.sub(offset, size), tmp.sub(offset, size), size);
    }

    //merge
    LcpCacheStringLoserTree<K> loserTree(tmp, ranges);
    loserTree.writeElementsToStream(output);
}

template <unsigned K>
static inline void
lcp_mergesort_cache_kway(string* strings, const LcpCacheStringPtr& outputPtr)
{
    size_t length = outputPtr.size;
    LcpCacheStringPtr tmpPtr(new string[length], new lcp_t[length], new char_type[length], length);

    lcp_mergesort_cache_kway<K>(strings, tmpPtr, outputPtr, length);

    delete[] tmpPtr.strings;
    delete[] tmpPtr.lcps;
    delete[] tmpPtr.cachedChars;
}

// K must be a power of two
template <unsigned K>
static inline void
lcp_mergesort_cache_kway(string* strings, size_t n)
{
    lcp_t* outputLcps = new lcp_t[n];
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n];

    char_type* outputCache = new char_type[n];
    char_type* tmpCache = new char_type[n];

    LcpCacheStringPtr output(strings, outputLcps, outputCache, n);
    LcpCacheStringPtr tmp(tmpStrings, tmpLcps, tmpCache, n);

    lcp_mergesort_cache_kway<K>(strings, tmp, output, n);

#ifdef EBERLE_LCP_LOSERTREE_MERGESORT_CHECK_LCPS
    //check lcps
    stringtools::verify_lcp_cache(output.strings, output.lcps, output.cachedChars, output.size, 0);
#endif  //EBERLE_LCP_LOSERTREE_MERGESORT_CHECK_LCPS

    delete[] outputLcps;
    delete[] tmpStrings;
    delete[] tmpLcps;

    delete[] outputCache;
    delete[] tmpCache;
}

void
lcp_mergesort_cache_4way(string* strings, size_t n)
{
    lcp_mergesort_cache_kway<4>(strings, n);
}

void
lcp_mergesort_cache_8way(string* strings, size_t n)
{
    lcp_mergesort_cache_kway<8>(strings, n);
}

void
lcp_mergesort_cache_16way(string* strings, size_t n)
{
    lcp_mergesort_cache_kway<16>(strings, n);
}

void
lcp_mergesort_cache_64way(string* strings, size_t n)
{
    lcp_mergesort_cache_kway<64>(strings, n);
}

PSS_CONTESTANT(lcp_mergesort_cache_4way, "bingmann/lcp_mergesort_cache_4way",
               "4-way LCP-Mergesort by Timo Bingmann and Andreas Eberle")
PSS_CONTESTANT(lcp_mergesort_cache_8way, "bingmann/lcp_mergesort_cache_8way",
               "8-way LCP-Mergesort by Timo Bingmann and Andreas Eberle")
PSS_CONTESTANT(lcp_mergesort_cache_16way, "bingmann/lcp_mergesort_cache_16way",
               "16-way LCP-Mergesort by Timo Bingmann and Andreas Eberle")
PSS_CONTESTANT(lcp_mergesort_cache_64way, "bingmann/lcp_mergesort_cache_64way",
               "64-way LCP-Mergesort by Timo Bingmann and Andreas Eberle")

} // namespace bingmann

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_LCP_MERGESORT_KWAY_HEADER

/******************************************************************************/
