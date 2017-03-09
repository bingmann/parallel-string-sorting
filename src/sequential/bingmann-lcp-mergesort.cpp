/*******************************************************************************
 * src/sequential/bingmann-lcp-mergesort.cpp
 *
 * LCP aware binary and k-way mergesort, implemented to verify pseudo-code in
 * journal.  Not necessarily the fastest implementations.
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

static inline void
lcp_compare(unsigned int a, string inputA, lcp_t lcpA,
            unsigned int b, string inputB, lcp_t lcpB,
            unsigned int& outSmaller, lcp_t& outLcpSmaller,
            unsigned int& outLarger, lcp_t& outLcpAB)
{
    if (lcpA == lcpB)
    {       // CASE 1 lcps are equal => do string comparision starting at lcp+1st position
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

static inline void
lcp_merge_binary(string* input1, lcp_t* lcps1, size_t length1,
                 string* input2, lcp_t* lcps2, size_t length2,
                 string* output, lcp_t* outputLcps)
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
lcp_mergesort_binary(string* strings, const LcpStringPtr& tmp, const LcpStringPtr& out, size_t length)
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
lcp_mergesort_binary(string* strings, size_t n)
{
    // Allocate memory for LCPs and temporary string array
    lcp_t* outputLcps = new lcp_t[n];
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n];

    LcpStringPtr output(strings, outputLcps, n);
    LcpStringPtr tmp(tmpStrings, tmpLcps, n);

    // execute lcp mergesort
    lcp_mergesort_binary(strings, tmp, output, n);

    // verify result
    stringtools::verify_lcp(strings, output.lcps, n, 0);

    delete[] outputLcps;
    delete[] tmpStrings;
    delete[] tmpLcps;
}

PSS_CONTESTANT(lcp_mergesort_binary, "bingmann/lcp_mergesort_binary",
               "Binary Mergesort with LCP-merge by Andreas Eberle and Timo Bingmann")

////////////////////////////////////////////////////////////////////////////////

template <unsigned K>
class LcpLoserTree
{
    typedef LcpStringPtr Stream;

    struct Node
    {
        unsigned int idx;
        lcp_t        lcp;
    };

private:
    Stream streams[K + 1];
    Node nodes[K + 1];

    //! play one comparison edge game: contender is the node below
    //! defender. After the game, defender contains the lower index, contender
    //! the winning index, and defender.lcp = lcp(s_loser,s_winner).
    inline void
    updateNode(Node& contender, Node& defender)
    {
        const Stream& defenderStream = streams[defender.idx];

        if (defenderStream.empty())
            return;

        const Stream& contenderStream = streams[contender.idx];

        if (contenderStream.empty())
        {
            std::swap(defender, contender);
            return;
        }
#if 1
        if (defender.lcp > contender.lcp)
        {
            // CASE 2: curr->lcp > contender->lcp => curr < contender
            std::swap(defender, contender);
        }
        else if (defender.lcp == contender.lcp)
        {
            // CASE 1: compare more characters
            lcp_t lcp = defender.lcp;

            string s1 = defenderStream.firstString() + lcp;
            string s2 = contenderStream.firstString() + lcp;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != 0 && *s1 == *s2)
                s1++, s2++, lcp++;

            if (*s1 < *s2) // CASE 1.1: curr < contender
                std::swap(defender, contender);

            // update inner node with lcp(s_1,s_2)
            defender.lcp = lcp;
        }
        else {
            // CASE 3: curr->lcp < contender->lcp => contender < curr  => nothing to do
        }
#else
        lcp_compare(contender.idx, contenderStream.firstString(), contender.lcp,
                    defender.idx, defenderStream.firstString(), defender.lcp,
                    contender.idx, contender.lcp, defender.idx, defender.lcp);
#endif
        assert(scmp(streams[contender.idx].firstString(),
                    streams[defender.idx].firstString()) <= 0);

        assert(calc_lcp(streams[contender.idx].firstString(),
                        streams[defender.idx].firstString()) == defender.lcp);
    }

    inline void
    initTree(lcp_t knownCommonLcp)
    {
        //std::cout << "inittree start\n";
        for (unsigned k = 1; k <= K; k++)
        {
            Node contender;
            contender.idx = k;
            contender.lcp = knownCommonLcp;

            unsigned nodeIdx = K + k;

            //std::cout << "nodeIdx " << nodeIdx << "\n";

            while (nodeIdx % 2 == 0 && nodeIdx > 2)
            {
                nodeIdx >>= 1;
                //std::cout << "play against " << nodeIdx << "\n";
                updateNode(contender, nodes[nodeIdx]);
            }
            nodeIdx = (nodeIdx + 1) / 2;
            //std::cout << "save as " << nodeIdx << "\n";
            nodes[nodeIdx] = contender;
        }
        //std::cout << "inittree done\n";
    }

public:
    LcpLoserTree(const LcpStringPtr& input, std::pair<size_t, size_t>* ranges, lcp_t knownCommonLcp = 0)
    {
        for (unsigned i = 1; i <= K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i - 1];

            streams[i] = input.sub(currRange.first, currRange.second);
        }

        initTree(knownCommonLcp);
    }

    inline void
    writeElementsToStream(LcpStringPtr outStream, const size_t length)
    {
        const LcpStringPtr end = outStream.sub(length, 0);

        while (outStream < end)
        {
            // take winner and put into output

            unsigned winnerIdx = nodes[1].idx;
            //std::cout << "winnerIdx " << winnerIdx << std::endl;

            outStream.setFirst(streams[winnerIdx].firstString(), nodes[1].lcp);
            ++outStream;

            // advance winner stream

            Stream& stream = streams[winnerIdx];
            ++stream;

            // run new items from winner stream up the tree

            Node& contender = nodes[1];

            if (!stream.empty())
                contender.lcp = streams[winnerIdx].firstLcp();

            unsigned nodeIdx = winnerIdx + K;
            //std::cout << "nodeIdx " << nodeIdx << "\n";

            while (nodeIdx > 2) {
                nodeIdx = (nodeIdx + 1) / 2;
                //std::cout << "play against " << nodeIdx << "\n";
                updateNode(contender, nodes[nodeIdx]);
            }
            //std::cout << "play against " << nodeIdx << "\n";

            // for (unsigned nodeIdx = (K + winnerIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
            // {
            //     updateNode(contender, nodes[nodeIdx]);
            // }
        }
    }
};

template <unsigned K>
static inline void
lcp_mergesort_kway(string* strings, const LcpStringPtr& tmp,
                   const LcpStringPtr& output, size_t length)
{
    if (length <= 2 * K)
    {
        memmove(output.strings, strings, length * sizeof(string));
        return bingmann_lcp_inssort::lcp_insertion_sort(
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
        lcp_mergesort_kway<K>(strings + offset, output.sub(offset, size), tmp.sub(offset, size), size);
    }

    // K-way merge
    LcpLoserTree<K> loserTree(tmp, ranges);
    loserTree.writeElementsToStream(output, length);
}

// K must be a power of two
template <unsigned K>
static inline void
lcp_mergesort_kway(string* strings, size_t n)
{
    lcp_t* outputLcps = new lcp_t[n + 1];
    string* tmpStrings = new string[n];
    lcp_t* tmpLcps = new lcp_t[n + 1];

    LcpStringPtr output(strings, outputLcps, n);
    LcpStringPtr tmp(tmpStrings, tmpLcps, n);

    lcp_mergesort_kway<K>(strings, tmp, output, n);

    // check lcps
    stringtools::verify_lcp(output.strings, output.lcps, n, 0);

    delete[] outputLcps;
    delete[] tmpStrings;
    delete[] tmpLcps;
}

static inline void
lcp_mergesort_4way(string* strings, size_t n)
{
    lcp_mergesort_kway<4>(strings, n);
}

static inline void
lcp_mergesort_8way(string* strings, size_t n)
{
    lcp_mergesort_kway<8>(strings, n);
}

static inline void
lcp_mergesort_16way(string* strings, size_t n)
{
    lcp_mergesort_kway<16>(strings, n);
}

static inline void
lcp_mergesort_128way(string* strings, size_t n)
{
    lcp_mergesort_kway<128>(strings, n);
}

PSS_CONTESTANT(lcp_mergesort_4way, "bingmann/lcp_mergesort_4way",
               "4-way LCP-Mergesort by Andreas Eberle and Timo Bingmann")
PSS_CONTESTANT(lcp_mergesort_8way, "bingmann/lcp_mergesort_8way",
               "8-way LCP-Mergesort by Andreas Eberle and Timo Bingmann")
PSS_CONTESTANT(lcp_mergesort_16way, "bingmann/lcp_mergesort_16way",
               "16-way LCP-Mergesort by Andreas Eberle and Timo Bingmann")
PSS_CONTESTANT(lcp_mergesort_128way, "bingmann/lcp_mergesort_128way",
               "128-way LCP-Mergesort by Andreas Eberle and Timo Bingmann")

} // namespace bingmann_lcp_mergesort

#endif // BINGMANN_LCP_MERGESORT_H_

/******************************************************************************/
