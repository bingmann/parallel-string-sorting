/******************************************************************************
 * src/sequential/bingmann-lcp-mergesort.h
 *
 * LCP aware binary and k-way mergesort, implemented to verify pseudo-code in
 * journal.  Not necessarily the fastest implementations.
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
#include "bingmann-lcp_inssort.h"

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

static inline void
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

CONTESTANT_REGISTER(lcp_mergesort_binary, "bingmann/lcp_mergesort_binary",
    "Binary Mergesort with LCP-merge by Andreas Eberle and Timo Bingmann")

////////////////////////////////////////////////////////////////////////////////

template <unsigned K>
class LcpLoserTree
{
    struct STREAM
    {
        LcpStringPtr elements;
        unsigned length;
        bool isEmpty;
    };

private:
    STREAM streams[K];
    unsigned nodes[K];
    lcp_t lcps[K];

    /*
     * Returns the winner of all games.
     */
    inline unsigned
    updateNode(unsigned &defenderIdx, unsigned contenderIdx)
    {
        const STREAM* defenderStream = streams + defenderIdx;

        if (defenderStream->isEmpty)
            return contenderIdx;

        const STREAM* contenderStream = streams + contenderIdx;

        lcp_t* contenderLcp = lcps + contenderIdx;
        lcp_t* defenderLcp = lcps + defenderIdx;

        if (contenderStream->isEmpty || *defenderLcp > *contenderLcp)
        { // CASE 2: curr->lcp > contender->lcp => curr < contender
            std::swap(defenderIdx, contenderIdx);

        }
        else if (*defenderLcp == *contenderLcp)
        { // CASE 1: curr.lcp == contender.lcp
            lcp_t lcp = *defenderLcp;
            string s1 = defenderStream->elements.str() + lcp;
            string s2 = contenderStream->elements.str() + lcp;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != 0 && *s1 == *s2)
                s1++, s2++, lcp++;

            if (*s1 < *s2)
            { 	// CASE 1.1: curr < contender
                *contenderLcp = lcp;
                std::swap(defenderIdx, contenderIdx);
            }
            else
            {	// CASE 1.2: curr >= contender
                *defenderLcp = lcp;
            }
        } // else // CASE 3: curr->lcp < contender->lcp => contender < curr  => nothing to do

        return contenderIdx;
    }

    inline void
    initTree(lcp_t knownCommonLcp)
    {
        for (unsigned i = 0; i < K; i++)
        {
            lcps[i] = knownCommonLcp;
            unsigned nodeIdx = K + i;
            unsigned contenderIdx = i;

            while (nodeIdx % 2 == 1 && nodeIdx > 1)
            {
                nodeIdx >>= 1;
                contenderIdx = updateNode(nodes[nodeIdx], contenderIdx);
            }
            nodes[nodeIdx >> 1] = contenderIdx;
        }
    }

    inline void
    removeTopFromStream(unsigned streamIdx)
    {
        STREAM* stream = streams + streamIdx;

        stream->length--;
        ++(stream->elements);
        stream->isEmpty = stream->length <= 0;

        if (!stream->isEmpty)
        {
            lcps[streamIdx] = stream->elements.lcp();
        }
    }

public:
    LcpLoserTree(const LcpStringPtr& input, std::pair<size_t, size_t>* ranges, lcp_t knownCommonLcp = 0)
    {
        for (unsigned i = 0; i < K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i];
            STREAM* curr = streams + i;
            curr->elements = input + currRange.first;
            curr->length = currRange.second;
            curr->isEmpty = (curr->length <= 0);
        }

        initTree(knownCommonLcp);
    }

    inline void
    writeElementsToStream(LcpStringPtr outStream, const size_t length)
    {
        const LcpStringPtr end = outStream + length;
        unsigned contenderIdx = nodes[0];

        while (outStream < end)
        {
            outStream.set(streams[contenderIdx].elements.str(), lcps[contenderIdx]);
            ++outStream;

            removeTopFromStream(contenderIdx);

            for (unsigned nodeIdx = (K + contenderIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
            {
                contenderIdx = updateNode(nodes[nodeIdx], contenderIdx);
            }
        }

        nodes[0] = contenderIdx;
    }
};

template <unsigned K>
static inline void
lcp_mergesort_kway(string* strings, const LcpStringPtr& tmp, const LcpStringPtr& output, size_t length)
{
    if (length <= 2 * K)
    {
        memcpy(output.strings, strings, length * sizeof(string));
        return bingmann_lcp_inssort::lcp_insertion_sort(output.strings, output.lcps, length, 0);
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
        lcp_mergesort_kway<K>(strings + offset, output + offset, tmp + offset, ranges[i].second);
    }

    // K-way merge
    LcpLoserTree<K> loserTree(tmp, ranges);
    loserTree.writeElementsToStream(output, length);
}

// K must be a power of two
template <unsigned K>
static inline void
lcp_mergesort_kway(string *strings, size_t n)
{
    lcp_t* outputLcps = new lcp_t[n+1];
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n+1];

    LcpStringPtr output(strings, outputLcps);
    LcpStringPtr tmp(tmpStrings, tmpLcps);

    lcp_mergesort_kway<K>(strings, tmp, output, n);

    // check lcps
    stringtools::verify_lcp(output.strings, output.lcps, n, 0);

    delete[] outputLcps;
    delete[] tmpStrings;
    delete[] tmpLcps;
}

static inline void
lcp_mergesort_4way(string *strings, size_t n)
{
    lcp_mergesort_kway<4>(strings, n);
}

static inline void
lcp_mergesort_8way(string *strings, size_t n)
{
    lcp_mergesort_kway<8>(strings, n);
}

static inline void
lcp_mergesort_16way(string *strings, size_t n)
{
    lcp_mergesort_kway<16>(strings, n);
}

CONTESTANT_REGISTER(lcp_mergesort_4way, "bingmann/lcp_mergesort_4way",
                    "4-way LCP-Mergesort by Andreas Eberle and Timo Bingmann")
CONTESTANT_REGISTER(lcp_mergesort_8way, "bingmann/lcp_mergesort_8way",
                    "8-way LCP-Mergesort by Andreas Eberle and Timo Bingmann")
CONTESTANT_REGISTER(lcp_mergesort_16way, "bingmann/lcp_mergesort_16way",
                    "16-way LCP-Mergesort by Andreas Eberle and Timo Bingmann")

}
// namespace bingmann_lcp_mergesort

#endif // BINGMANN_LCP_MERGESORT_H_

