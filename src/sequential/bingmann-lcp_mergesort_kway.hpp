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

#ifndef BINGMANN_LCP_MERGESORT_KWAY_H_
#define BINGMANN_LCP_MERGESORT_KWAY_H_

#include <iostream>
#include <cstring>

#include "../tools/stringtools.hpp"
#include "../tools/contest.hpp"
#include "bingmann-lcp_inssort.hpp"

namespace bingmann_lcp_mergesort {

using namespace stringtools;

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
    LcpLoserTree(const LcpStringPtr& input, std::pair<size_t, size_t>* ranges,
                 lcp_t knownCommonLcp = 0)
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
        lcp_mergesort_kway<K>(strings + offset,
                              output.sub(offset, size),
                              tmp.sub(offset, size), size);
    }

    // K-way merge
    LcpLoserTree<K> loserTree(tmp, ranges);
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
lcp_mergesort_128way(string* strings, uintptr_t* lcp, size_t n)
{
    lcp_mergesort_kway<128>(strings, lcp, n);
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

#endif // BINGMANN_LCP_MERGESORT_KWAY_H_

/******************************************************************************/
