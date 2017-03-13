/*******************************************************************************
 * src/sequential/bingmann-lcp-mergesort.cpp
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

#ifndef BINGMANN_LCP_MERGESORT_H_
#define BINGMANN_LCP_MERGESORT_H_

#include <iostream>
#include <cstring>

#include "../tools/stringtools.hpp"
#include "../tools/contest.hpp"
#include "bingmann-lcp_inssort.hpp"

namespace bingmann_lcp_mergesort {

using namespace stringtools;

//! archetypal lcp-aware comparison method
static inline void
lcp_compare(unsigned int a, string inputA, lcp_t lcpA,
            unsigned int b, string inputB, lcp_t lcpB,
            unsigned int& outSmaller, lcp_t& outLcpSmaller,
            unsigned int& outLarger, lcp_t& outLcpAB)
{
    if (lcpA == lcpB)
    {
        // CASE 1 lcps are equal => do string comparision starting at lcp+1st
        // position
        string sA = inputA + lcpA;
        string sB = inputB + lcpA;

        // check the strings starting after lcp and calculate new lcp
        while (*sA != 0 && *sA == *sB)
            sA++, sB++;

        const lcp_t h = sA - inputA;

        if (*sA <= *sB) // CASE 1.1: s1 <= s2
        {
            outSmaller = a;
            outLcpSmaller = lcpA;
            outLarger = b;
            outLcpAB = h;
        }
        else // CASE 1.2: s1 > s2
        {
            outSmaller = b;
            outLcpSmaller = lcpB;
            outLarger = a;
            outLcpAB = h;
        }
    }
    else if (lcpA > lcpB) // CASE 2: lcp1 < lcp2 -> s_1 > s_2
    {
        outSmaller = a;
        outLcpSmaller = lcpA;
        outLarger = b;
        outLcpAB = lcpB;
    }
    else // CASE 3: lcp1 > lcp2 -> s_1 < s_2
    {
        outSmaller = b;
        outLcpSmaller = lcpB;
        outLarger = a;
        outLcpAB = lcpA;
    }

    assert(calc_lcp(inputA, inputB) == outLcpAB);
}

/******************************************************************************/
// bingmann/lcp_mergesort_binary

static inline void
lcp_merge_binary(const string* input1, const lcp_t* lcps1, size_t length1,
                 const string* input2, const lcp_t* lcps2, size_t length2,
                 string* output, lcp_t* output_lcp)
{
    const string* end1 = input1 + length1;
    const string* end2 = input2 + length2;

    string prev = (string)""; // sentinel
    (void)prev;

    lcp_t lcp1 = *lcps1;
    lcp_t lcp2 = *lcps2;

    // do the merge
    while (input1 < end1 && input2 < end2)
    {
        unsigned int cmpSmaller, cmpLarger;
        lcp_t cmpLcp, cmpLcpSmaller;

        //std::cout << "prev = " << prev << " - input1 = " << *input1 << "\n";
        assert(calc_lcp(prev, *input1) == lcp1);
        assert(calc_lcp(prev, *input2) == lcp2);

        lcp_compare(0, *input1, lcp1, 1, *input2, lcp2,
                    cmpSmaller, cmpLcpSmaller, cmpLarger, cmpLcp);

        if (cmpSmaller == 0)
        {
            prev = *input1;
            *output++ = *input1;
            *output_lcp++ = lcp1;

            ++input1;
            ++lcps1;

            lcp1 = input1 < end1 ? *lcps1 : 0;
            lcp2 = cmpLcp;
        }
        else
        {
            prev = *input2;
            *output++ = *input2;
            *output_lcp++ = lcp2;

            ++input2;
            ++lcps2;

            lcp1 = cmpLcp;
            lcp2 = input2 < end2 ? *lcps2 : 0;
        }
    }

    if (input1 < end1)
    {
        // if there are remaining elements in stream1, copy them to the end
        memcpy(output, input1, (end1 - input1) * sizeof(string));
        memcpy(output_lcp, lcps1, (end1 - input1) * sizeof(lcp_t));
        *output_lcp = lcp1;
    }
    else
    {
        memcpy(output, input2, (end2 - input2) * sizeof(string));
        memcpy(output_lcp, lcps2, (end2 - input2) * sizeof(lcp_t));
        *output_lcp = lcp2;
    }
}

static inline void
lcp_mergesort_binary(string* strings, const LcpStringPtr& tmp,
                     const LcpStringPtr& out, size_t length)
{
    if (length == 0) {
        return;
    }
    else if (length == 1) {
        out.setFirst(*strings, 0);
        return;
    }

    size_t length1 = length / 2;
    size_t length2 = length - length1;

    LcpStringPtr out1 = out.sub(0, length1);
    LcpStringPtr out2 = out.sub(length1, length2);

    LcpStringPtr tmp1 = tmp.sub(0, length1);
    LcpStringPtr tmp2 = tmp.sub(length1, length2);

    lcp_mergesort_binary(strings, out1, tmp1, length1);
    lcp_mergesort_binary(strings + length1, out2, tmp2, length2);

    lcp_merge_binary(tmp1.strings, tmp1.lcps, tmp1.size,
                     tmp2.strings, tmp2.lcps, tmp2.size,
                     out.strings, out.lcps);
}

static inline void
lcp_mergesort_binary(string* strings, uintptr_t* lcp, size_t n)
{
    // Allocate memory for LCPs and temporary string array
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n];

    LcpStringPtr output(strings, lcp, n);
    LcpStringPtr tmp(tmpStrings, tmpLcps, n);

    // execute lcp mergesort
    lcp_mergesort_binary(strings, tmp, output, n);

    delete[] tmpStrings;
    delete[] tmpLcps;
    lcp[0] = 42;
}

PSS_CONTESTANT(lcp_mergesort_binary, "bingmann/lcp_mergesort_binary",
               "Binary Mergesort with LCP-merge by Andreas Eberle and Timo Bingmann")

/******************************************************************************/

static inline void
lcp_merge_binary_opt(const string* input1, const lcp_t* lcps1, size_t length1,
                     const string* input2, const lcp_t* lcps2, size_t length2,
                     string* output, lcp_t* output_lcp)
{
    const string* end1 = input1 + length1;
    const string* end2 = input2 + length2;

    lcp_t lcp1 = *lcps1;
    lcp_t lcp2 = *lcps2;

    // do the merge
    while (input1 < end1 && input2 < end2)
    {
        if (lcp1 == lcp2)
        {
            // CASE 1 lcps are equal => do string comparision starting at
            // lcp+1st position
            string s1 = *input1 + lcp1;
            string s2 = *input2 + lcp1;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != '\0' && *s1 == *s2)
                s1++, s2++;

            const lcp_t lcp = s1 - *input1;

            if (*s1 <= *s2)
            {
                // CASE 1.1: input1 <= input2
                *output = *input1;
                *output_lcp = lcp1;
                ++input1;
                ++lcps1;
                lcp1 = input1 < end1 ? *lcps1 : 0;
                lcp2 = lcp;
            }
            else
            {
                // CASE 1.2: input1 > input2
                *output = *input2;
                *output_lcp = lcp2;
                ++input2;
                ++lcps2;
                lcp1 = lcp;
                lcp2 = input2 < end2 ? *lcps2 : 0;
            }
        }
        else if (lcp1 < lcp2)
        {
            // CASE 2: input1 > input2
            *output = *input2;
            *output_lcp = lcp2;
            ++input2;
            ++lcps2;
            lcp2 = input2 < end2 ? *lcps2 : 0;
        }
        else
        {
            // CASE 3: input1 < input2
            *output = *input1;
            *output_lcp = lcp1;
            ++input1;
            ++lcps1;
            lcp1 = input1 < end1 ? *lcps1 : 0;
        }

        ++output;
        ++output_lcp;
    }

    if (input1 < end1)
    {
        // if there are remaining elements in stream1, copy them to the end
        memcpy(output, input1, (end1 - input1) * sizeof(string));
        memcpy(output_lcp, lcps1, (end1 - input1) * sizeof(lcp_t));
        *output_lcp = lcp1;
    }
    else
    {
        memcpy(output, input2, (end2 - input2) * sizeof(string));
        memcpy(output_lcp, lcps2, (end2 - input2) * sizeof(lcp_t));
        *output_lcp = lcp2;
    }
}

static inline void
lcp_merge_binary_opt(const LcpStringPtr& input1, const LcpStringPtr& input2,
                     const LcpStringPtr& output)
{
    lcp_merge_binary_opt(input1.strings, input1.lcps, input1.size,
                         input2.strings, input2.lcps, input2.size,
                         output.strings, output.lcps);
}

static inline void
lcp_mergesort_binary_opt(
    string* strings, size_t length, const LcpStringPtr& tmp,
    const LcpStringPtr& out)
{
    if (length == 0) {
        return;
    }
    else if (length == 1) {
        out.setFirst(*strings, 0);
        return;
    }

    const size_t length1 = length / 2;
    const size_t length2 = length - length1;

    const LcpStringPtr out1 = out.sub(0, length1);
    const LcpStringPtr out2 = out.sub(length1, length2);

    const LcpStringPtr tmp1 = tmp.sub(0, length1);
    const LcpStringPtr tmp2 = tmp.sub(length1, length2);

    lcp_mergesort_binary_opt(strings, length1, out1, tmp1);
    lcp_mergesort_binary_opt(strings + length1, length2, out2, tmp2);

    lcp_merge_binary_opt(tmp1, tmp2, out);
}

void
lcp_mergesort_binary_opt(string* strings, uintptr_t* lcp, size_t n)
{
    // Allocate memory for LCPs and temporary string array
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n];

    LcpStringPtr output(strings, lcp, n);
    LcpStringPtr tmp(tmpStrings, tmpLcps, n);

    // execute lcp mergesort
    lcp_mergesort_binary_opt(strings, n, tmp, output);

    delete[] tmpStrings;
    delete[] tmpLcps;
    lcp[0] = 42;
}

PSS_CONTESTANT(lcp_mergesort_binary_opt, "eberle/mergesort_lcp_binary_opt",
               "Binary Mergesort with LCP-usage by Andreas Eberle and Timo Bingmann")

/******************************************************************************/

template <bool OutputLcp, bool OutputCache>
static inline void
lcp_merge_binary_cache_opt(
    const string* input1, const lcp_t* lcps1, char_type* cache1, size_t length1,
    const string* input2, const lcp_t* lcps2, char_type* cache2, size_t length2,
    string* output, lcp_t* output_lcp, char_type* output_cache)
{
    const string* end1 = input1 + length1;
    const string* end2 = input2 + length2;

    lcp_t lcp1 = *lcps1;
    lcp_t lcp2 = *lcps2;

    //do the merge
    while (input1 < end1 && input2 < end2)
    {
        if (lcp1 == lcp2)
        {
            // CASE 1 lcps are equal => do string comparision starting at lcp+1st position
            lcp_t lcp = lcp1;

            char_type ch1 = *cache1;
            char_type ch2 = *cache2;

            assert((*input1)[lcp] == ch1);
            assert((*input2)[lcp] == ch2);

            // check the strings starting after lcp and calculate new lcp
            while (ch1 != 0 && ch1 == ch2)
            {
                ++lcp;
                ch1 = (*input1)[lcp];
                ch2 = (*input2)[lcp];
            }

            if (ch1 < ch2)
            {
                // CASE 1.1: lcp1 < lcp2
                *output = *input1;
                if (OutputLcp)
                    *output_lcp = lcp1;
                if (OutputCache)
                    *output_cache = ch1;
                ++input1, ++lcps1, ++cache1;
                lcp1 = *lcps1;
                lcp2 = lcp;

                *cache2 = ch2;
            }
            else
            {
                // CASE 1.2: lcp1 >= lcp2
                *output = *input2;
                if (OutputLcp)
                    *output_lcp = lcp2;
                if (OutputCache)
                    *output_cache = ch2;
                ++input2, ++lcps2, ++cache2;
                lcp1 = lcp;
                lcp2 = *lcps2;

                *cache1 = ch1;
            }
        }
        else if (lcp1 < lcp2)
        {
            // CASE 2: lcp1 > lcp2
            *output = *input2;
            if (OutputLcp)
                *output_lcp = lcp2;
            if (OutputCache)
                *output_cache = *cache2;
            ++input2, ++lcps2, ++cache2;
            lcp2 = *lcps2;
        }
        else
        {
            // CASE 3: lcp1 < lcp2
            *output = *input1;
            if (OutputLcp)
                *output_lcp = lcp1;
            if (OutputCache)
                *output_cache = *cache1;
            ++input1, ++lcps1, ++cache1;
            lcp1 = *lcps1;
        }

        ++output;
    }

    if (input1 < end1)
    {
        // if there are remaining elements in stream1, copy them to the end
        memcpy(output, input1, (end1 - input1) * sizeof(string));

        if (OutputLcp) {
            memcpy(output_lcp, lcps1, (end1 - input1) * sizeof(lcp_t));
            *output_lcp = lcp1;
        }
        if (OutputCache) {
            memcpy(output_cache, cache1, (end1 - input1) * sizeof(uint8_t));
            //*output_cache = ch1;
        }
    }
    else
    {
        memcpy(output, input2, (end2 - input2) * sizeof(string));

        if (OutputLcp) {
            memcpy(output_lcp, lcps2, (end2 - input2) * sizeof(lcp_t));
            *output_lcp = lcp2;
        }
        if (OutputCache) {
            memcpy(output_cache, cache2, (end2 - input2) * sizeof(uint8_t));
            //*output_cache = ch2;
        }
    }
}

static inline void
lcp_merge_binary_cache_opt(
    const LcpCacheStringPtr& input1, const LcpCacheStringPtr& input2,
    string* output)
{
    if (0)
    {
        stringtools::verify_lcp_cache(
            input1.strings, input1.lcps, input1.cachedChars,
            input1.size, -1);

        stringtools::verify_lcp_cache(
            input2.strings, input2.lcps, input2.cachedChars,
            input2.size, -1);
    }

    lcp_merge_binary_cache_opt</* OutputLcp */ false, /* OutputCache */ false>(
        input1.strings, input1.lcps, input1.cachedChars, input1.size,
        input2.strings, input2.lcps, input2.cachedChars, input2.size,
        output, nullptr, nullptr);
}

static inline void
lcp_merge_binary_cache_opt(
    const LcpCacheStringPtr& input1, const LcpCacheStringPtr& input2,
    const LcpCacheStringPtr& output)
{
    if (0)
    {
        stringtools::verify_lcp_cache(
            input1.strings, input1.lcps, input1.cachedChars,
            input1.size, -1);

        stringtools::verify_lcp_cache(
            input2.strings, input2.lcps, input2.cachedChars,
            input2.size, -1);
    }

    lcp_merge_binary_cache_opt</* OutputLcp */ true, /* OutputCache */ true>(
        input1.strings, input1.lcps, input1.cachedChars, input1.size,
        input2.strings, input2.lcps, input2.cachedChars, input2.size,
        output.strings, output.lcps, output.cachedChars);
}

} // namespace bingmann_lcp_mergesort

#endif // BINGMANN_LCP_MERGESORT_H_

/******************************************************************************/
