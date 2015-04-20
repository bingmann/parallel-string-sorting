/*******************************************************************************
 * src/sequential/eberle-mergesort-lcp.h
 *
 * LCP aware binary mergesort.
 *
 *******************************************************************************
 * Copyright (C) 2013-2014 Andreas Eberle <email@andreas-eberle.com>
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

#ifndef PSS_SRC_SEQUENTIAL_EBERLE_MERGESORT_LCP_HEADER
#define PSS_SRC_SEQUENTIAL_EBERLE_MERGESORT_LCP_HEADER

#include "../tools/stringtools.h"
#include "../tools/stats_writer.h"

namespace eberle_mergesort_lcp {

using namespace stringtools;

static inline void
eberle_lcp_merge(string* input1, lcp_t* lcps1, size_t length1, string* input2, lcp_t* lcps2, size_t length2, string* output, lcp_t* outputLcps)
{
    const string* end1 = input1 + length1;
    const string* end2 = input2 + length2;

    lcp_t lcp1 = *lcps1;
    lcp_t lcp2 = *lcps2;

    //do the merge
    while (input1 < end1 && input2 < end2)
    {
        if (lcp1 == lcp2)
        {       // CASE 1 lcps are equal => do string comparision starting at lcp+1st position
            string s1 = *input1 + lcp1;
            string s2 = *input2 + lcp1;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != '\0' && *s1 == *s2)
                s1++, s2++;

            const lcp_t lcp = s1 - *input1;

            if (*s1 <= *s2)
            {       // CASE 1.1: input1 <= input2
                *output = *input1;
                *outputLcps = lcp1;
                ++input1;
                ++lcps1;
                lcp1 = *lcps1;
                lcp2 = lcp;
            }
            else
            {       // CASE 1.2: input1 > input2
                *output = *input2;
                *outputLcps = lcp2;
                ++input2;
                ++lcps2;
                lcp1 = lcp;
                lcp2 = *lcps2;
            }
        }
        else if (lcp1 < lcp2)
        {       // CASE 2: input1 > input2
            *output = *input2;
            *outputLcps = lcp2;
            ++input2;
            ++lcps2;
            lcp2 = *lcps2;
        }
        else
        {       // CASE 3: input1 < input2
            *output = *input1;
            *outputLcps = lcp1;
            ++input1;
            ++lcps1;
            lcp1 = *lcps1;
        }

        ++output;
        ++outputLcps;
    }

    if (input1 < end1)
    {       // if there are remaining elements in stream1, copy them to the end
        memcpy(output, input1, (end1 - input1) * sizeof(string));
        memcpy(outputLcps, lcps1, (end1 - input1) * sizeof(lcp_t));
        *outputLcps = lcp1;
    }
    else
    {
        memcpy(output, input2, (end2 - input2) * sizeof(string));
        memcpy(outputLcps, lcps2, (end2 - input2) * sizeof(lcp_t));
        *outputLcps = lcp2;
    }
}

static inline void
eberle_lcp_merge(const LcpStringPtr& input1, const LcpStringPtr& input2, const LcpStringPtr& output)
{
    eberle_lcp_merge(input1.strings, input1.lcps, input1.size,
                     input2.strings, input2.lcps, input2.size,
                     output.strings, output.lcps);
}

static inline void
eberle_lcp_mergesort(string* strings, size_t length, const LcpStringPtr& tmp, const LcpStringPtr& out)
{
    if (length <= 1)
    {
        out.setFirst(*strings, 0);
        return;
    }

    const size_t length1 = length / 2;
    const size_t length2 = length - length1;

    const LcpStringPtr out1 = out.sub(0, length1);
    const LcpStringPtr out2 = out.sub(length1, length2);

    const LcpStringPtr tmp1 = tmp.sub(0, length1);
    const LcpStringPtr tmp2 = tmp.sub(length1, length2);

    eberle_lcp_mergesort(strings, length1, out1, tmp1);
    eberle_lcp_mergesort(strings + length1, length2, out2, tmp2);

    eberle_lcp_merge(tmp1, tmp2, out);
}

void
eberle_lcp_mergesort(string* strings, size_t n)
{
    // Allocate memory for LCPs and temporary string array
    lcp_t* outputLcps = new lcp_t[n];
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n];

    LcpStringPtr output(strings, outputLcps, n);
    LcpStringPtr tmp(tmpStrings, tmpLcps, n);

    // execute lcp mergesort
    eberle_lcp_mergesort(strings, n, tmp, output);

    delete[] outputLcps;
    delete[] tmpStrings;
    delete[] tmpLcps;
}

PSS_CONTESTANT(eberle_lcp_mergesort, "eberle/mergesort_lcp_binary", "Binary Mergesort with LCP-usage by Andreas Eberle")

static inline void
eberle_lcp_merge(string* input1, lcp_t* lcps1, size_t length1,
                 string* input2, lcp_t* lcps2, size_t length2,
                 string* output)
{
    const string* end1 = input1 + length1;
    const string* end2 = input2 + length2;

    lcp_t lcp1 = *lcps1;
    lcp_t lcp2 = *lcps2;

    //do the merge
    while (input1 < end1 && input2 < end2)
    {
        if (lcp1 == lcp2)
        {       // CASE 1 lcps are equal => do string comparision starting at lcp+1st position
            string s1 = *input1 + lcp1;
            string s2 = *input2 + lcp1;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != '\0' && *s1 == *s2)
                s1++, s2++;

            const lcp_t lcp = s1 - *input1;

            if (*s1 <= *s2)
            {       // CASE 1.1: lcp1 <= lcp2
                *output = *input1;
                ++input1;
                ++lcps1;
                lcp1 = *lcps1;
                lcp2 = lcp;
            }
            else
            {       // CASE 1.2: lcp1 > lcp2
                *output = *input2;
                ++input2;
                ++lcps2;
                lcp1 = lcp;
                lcp2 = *lcps2;
            }
        }
        else if (lcp1 < lcp2)
        {       // CASE 2: lcp1 > lcp2
            *output = *input2;
            ++input2;
            ++lcps2;
            lcp2 = *lcps2;
        }
        else
        {       // CASE 3: lcp1 < lcp2
            *output = *input1;
            ++input1;
            ++lcps1;
            lcp1 = *lcps1;
        }

        ++output;
    }

    if (input1 < end1)
    {       // if there are remaining elements in stream1, copy them to the end
        memcpy(output, input1, (end1 - input1) * sizeof(string));
    }
    else
    {
        memcpy(output, input2, (end2 - input2) * sizeof(string));
    }
}

static inline void
eberle_lcp_merge(const LcpStringPtr& input1, const LcpStringPtr& input2, string* output)
{
    eberle_lcp_merge(input1.strings, input1.lcps, input1.size,
                     input2.strings, input2.lcps, input2.size, output);
}

static inline void
eberle_lcp_merge(string* input1, lcp_t* lcps1, char_type* cache1, size_t length1,
                 string* input2, lcp_t* lcps2, char_type* cache2, size_t length2,
                 string* output)
{
    const string* end1 = input1 + length1;
    const string* end2 = input2 + length2;

    lcp_t lcp1 = *lcps1;
    lcp_t lcp2 = *lcps2;

    //do the merge
    while (input1 < end1 && input2 < end2)
    {
        if (lcp1 == lcp2)
        {       // CASE 1 lcps are equal => do string comparision starting at lcp+1st position
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
            {       // CASE 1.1: lcp1 < lcp2
                *output = *input1;
                ++input1, ++lcps1, ++cache1;
                lcp1 = *lcps1;
                lcp2 = lcp;

                *cache2 = ch2;
            }
            else
            {       // CASE 1.2: lcp1 >= lcp2
                *output = *input2;
                ++input2, ++lcps2, ++cache2;
                lcp1 = lcp;
                lcp2 = *lcps2;

                *cache1 = ch1;
            }
        }
        else if (lcp1 < lcp2)
        {       // CASE 2: lcp1 > lcp2
            *output = *input2;
            ++input2, ++lcps2, ++cache2;
            lcp2 = *lcps2;
        }
        else
        {       // CASE 3: lcp1 < lcp2
            *output = *input1;
            ++input1, ++lcps1, ++cache1;
            lcp1 = *lcps1;
        }

        ++output;
    }

    if (input1 < end1)
    {       // if there are remaining elements in stream1, copy them to the end
        memcpy(output, input1, (end1 - input1) * sizeof(string));
    }
    else
    {
        memcpy(output, input2, (end2 - input2) * sizeof(string));
    }
}

static inline void
eberle_lcp_merge(const LcpCacheStringPtr& input1, const LcpCacheStringPtr& input2, string* output)
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

#if 1
    eberle_lcp_merge(input1.strings, input1.lcps, input1.cachedChars, input1.size,
                     input2.strings, input2.lcps, input2.cachedChars, input2.size,
                     output);
#else
    eberle_lcp_merge(input1.strings, input1.lcps, input1.size,
                     input2.strings, input2.lcps, input2.size,
                     output);
#endif
}

} // namespace eberle_lcp_mergesort

#endif // !PSS_SRC_SEQUENTIAL_EBERLE_MERGESORT_LCP_HEADER

/******************************************************************************/
