/******************************************************************************
 * src/sequential/bingmann-lcp-mergesort.h
 *
 * LCP aware binary mergesort, implemented to verify pseudo-code in journal.
 * Not necessarily the fastest implementation, and just binary.
 *
 ******************************************************************************
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
 *****************************************************************************/

#ifndef BINGMANN_LCP_MERGESORT_H_
#define BINGMANN_LCP_MERGESORT_H_

#include <iostream>
#include <cstring>

#include "../tools/stringtools.h"
#include "../tools/contest.h"

namespace bingmann_lcp_mergesort
{

using namespace stringtools;

static inline void
lcp_compare(int a, string inputA, lcp_t lcpA,
            int b, string inputB, lcp_t lcpB,
            int& output, lcp_t& outputLcp)
{
    if (lcpA == lcpB)
    { // CASE 1 lcps are equal => do string comparision starting at lcp+1st position
        string sA = inputA + lcpA;
        string sB = inputB + lcpA;

        // check the strings starting after lcp and calculate new lcp
        while (*sA != 0 && *sA == *sB)
            sA++, sB++;

        const lcp_t h = sA - inputA;

        if (*sA <= *sB) // CASE 1.1: s1 <= s2
        { 	
            output = a;
            outputLcp = h;
        }
        else // CASE 1.2: s1 > s2
        {
            output = b;
            outputLcp = h;
        }
    }
    else if (lcpA > lcpB) // CASE 2: lcp1 < lcp2 -> s_1 > s_2
    {
        output = a;
        outputLcp = lcpB;
    }
    else // CASE 3: lcp1 > lcp2 -> s_1 < s_2
    {
        output = b;
        outputLcp = lcpA;
    }
}

static inline void
lcp_merge_binary(string* input1, lcp_t* lcps1, size_t length1,
                 string* input2, lcp_t* lcps2, size_t length2,
                 string* output, lcp_t* outputLcps)
{
    const string* end1 = input1 + length1;
    const string* end2 = input2 + length2;

    string prev = (string)""; // sentinel

    lcp_t lcp1 = *lcps1;
    lcp_t lcp2 = *lcps2;

    //do the merge
    while (input1 < end1 && input2 < end2)
    {
        int cmp;
        lcp_t cmpLcp;

        //std::cout << "prev = " << prev << " - input1 = " << *input1 << "\n";
        assert( calc_lcp(prev, *input1) == lcp1 );
        assert( calc_lcp(prev, *input2) == lcp2 );

        lcp_compare(1, *input1, lcp1, 2, *input2, lcp2, cmp, cmpLcp);

        assert( calc_lcp(*input1, *input2) == cmpLcp );

        if (cmp == 1)
        {
            prev = *input1;
            *output++ = *input1;
            *outputLcps++ = lcp1;

            ++input1;
            ++lcps1;

            lcp1 = input1 < end1 ? *lcps1 : 0;
            lcp2 = cmpLcp;
        }
        else
        {
            prev = *input2;
            *output++ = *input2;
            *outputLcps++ = lcp2;

            ++input2;
            ++lcps2;

            lcp1 = cmpLcp;
            lcp2 = input2 < end2 ? *lcps2 : 0;
        }
    }

    if (input1 < end1)
    {   // if there are remaining elements in stream1, copy them to the end
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
lcp_mergesort_binary(string *strings, const LcpStringPtr& tmp, const LcpStringPtr& output, size_t length)
{
    if (length <= 1)
    {
        output.set(*strings, 0);
        return;
    }

    const size_t length1 = length / 2;
    const size_t length2 = length - length1;
    const LcpStringPtr tmp2 = tmp + length1;

    lcp_mergesort_binary(strings, output, tmp, length1);
    lcp_mergesort_binary(strings + length1, output + length1, tmp2, length2);

    lcp_merge_binary(tmp.strings, tmp.lcps, length1,
                     tmp2.strings, tmp2.lcps, length2,
                     output.strings, output.lcps);
}

void
lcp_mergesort_binary(string *strings, size_t n)
{
    // Allocate memory for LCPs and temporary string array
    lcp_t* outputLcps = new lcp_t[n];
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n];

    LcpStringPtr output(strings, outputLcps);
    LcpStringPtr tmp(tmpStrings, tmpLcps);

    // execute lcp mergesort
    lcp_mergesort_binary(strings, tmp, output, n);

    // verify result
    stringtools::verify_lcp(strings, output.lcps, n, 0);

    delete[] outputLcps;
    delete[] tmpStrings;
    delete[] tmpLcps;
}

CONTESTANT_REGISTER(lcp_mergesort_binary, "bingmann/lcp_mergesort_binary", "Binary Mergesort with LCP-merge by Andreas Eberle and Timo Bingmann")

}
// namespace bingmann_lcp_mergesort

#endif // BINGMANN_LCP_MERGESORT_H_

