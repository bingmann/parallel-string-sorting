/*******************************************************************************
 * src/sequential/bingmann-lcp_losertree.hpp
 *
 * Implementation of a LCP aware multiway losertree.
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

#ifndef PSS_SRC_SEQUENTIAL_BINGMANN_LCP_LOSERTREE_HEADER
#define PSS_SRC_SEQUENTIAL_BINGMANN_LCP_LOSERTREE_HEADER

#include <iostream>
#include <algorithm>
#include <utility>

#include "../tools/stringtools.hpp"
#include <tlx/define/likely.hpp>

namespace bingmann {

using namespace stringtools;

typedef unsigned char* string;

/******************************************************************************/
// LcpStringLoserTree

template <size_t K>
class LcpStringLoserTree
{
    typedef LcpStringPtr Stream;

    struct Node
    {
        size_t idx;
        lcp_t  lcp;
    };

private:
    Stream streams[K + 1];
    Node nodes[K + 1];

    //! play one comparison edge game: contender is the node below
    //! defender. After the game, defender contains the lower index, contender
    //! the winning index, and defender.lcp = lcp(s_loser,s_winner).
    void updateNode(Node& contender, Node& defender)
    {
        const Stream& defenderStream = streams[defender.idx];

        if (TLX_UNLIKELY(defenderStream.empty()))
            return;

        const Stream& contenderStream = streams[contender.idx];

        if (TLX_UNLIKELY(contenderStream.empty()))
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

    void initTree(lcp_t knownCommonLcp)
    {
        //std::cout << "inittree start\n";
        for (size_t k = 1; k <= K; k++)
        {
            Node contender;
            contender.idx = k;
            contender.lcp = knownCommonLcp;

            size_t nodeIdx = K + k;

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
    LcpStringLoserTree(const LcpStringPtr& input, std::pair<size_t, size_t>* ranges,
                       lcp_t knownCommonLcp = 0)
    {
        for (size_t i = 1; i <= K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i - 1];

            streams[i] = input.sub(currRange.first, currRange.second);
        }

        initTree(knownCommonLcp);
    }

    void writeElementsToStream(LcpStringPtr outStream, const size_t length)
    {
        const LcpStringPtr end = outStream.sub(length, 0);

        while (outStream < end)
        {
            // take winner and put into output

            size_t winnerIdx = nodes[1].idx;
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

            size_t nodeIdx = winnerIdx + K;
            //std::cout << "nodeIdx " << nodeIdx << "\n";

            while (nodeIdx > 2) {
                nodeIdx = (nodeIdx + 1) / 2;
                //std::cout << "play against " << nodeIdx << "\n";
                updateNode(contender, nodes[nodeIdx]);
            }
            //std::cout << "play against " << nodeIdx << "\n";

            // for (size_t nodeIdx = (K + winnerIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
            // {
            //     updateNode(contender, nodes[nodeIdx]);
            // }
        }
    }
};

/******************************************************************************/
// LcpCacheStringLoserTree

template <size_t K>
class LcpCacheStringLoserTree
{
    typedef LcpCacheStringPtr Stream;

public:
    //! allow public access to streams for easier initialization
    Stream streams[K];

protected:
    size_t nodes[K];
    lcp_t lcps[K];
    char_type cached[K];

    /*
     * Returns the winner of all games.
     */
    void updateNode(size_t& contenderIdx, size_t& defenderIdx)
    {
        const Stream& defenderStream = streams[defenderIdx];

        if (TLX_UNLIKELY(defenderStream.empty()))
            return;

        const Stream& contenderStream = streams[contenderIdx];

        lcp_t& contenderLcp = lcps[contenderIdx];
        lcp_t& defenderLcp = lcps[defenderIdx];

        if (TLX_UNLIKELY(contenderStream.empty()) || defenderLcp > contenderLcp)
        {
            // CASE 2: defender->lcp > contender->lcp => defender < contender
            std::swap(defenderIdx, contenderIdx);
        }
        else if (defenderLcp == contenderLcp)
        {
            // CASE 1: defender.lcp == contender.lcp
            lcp_t lcp = defenderLcp;

            char_type c1 = cached[defenderIdx];
            char_type c2 = cached[contenderIdx];

            assert(c1 == defenderStream.firstString()[lcp]);
            assert(c2 == contenderStream.firstString()[lcp]);

            // check the strings starting after lcp and calculate new lcp
            while (c1 != 0 && c1 == c2) {
                lcp++;
                c1 = defenderStream.firstString()[lcp];
                c2 = contenderStream.firstString()[lcp];
            }

            if (c1 < c2) {
                // CASE 1.1: defender < contender
                contenderLcp = lcp;
                cached[contenderIdx] = c2;
                std::swap(defenderIdx, contenderIdx);
            }
            else {
                // CASE 1.2: defender >= contender
                defenderLcp = lcp;
                cached[defenderIdx] = c1;
            }
        }
        else {
            // CASE 3: defender->lcp < contender->lcp => contender < defender
            // => nothing to do
        }
    }

    void removeTopFromStream(size_t streamIdx)
    {
        Stream& stream = streams[streamIdx];

        ++stream;

        if (TLX_LIKELY(!stream.empty()))
        {
            lcps[streamIdx] = stream.firstLcp();
            cached[streamIdx] = stream.firstCached();
        }
    }

public:
    LcpCacheStringLoserTree() { }

    LcpCacheStringLoserTree(
        const Stream& input, std::pair<size_t, size_t>* ranges,
        lcp_t knownCommonLcp = 0)
    {
        for (size_t i = 0; i < K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i];
            Stream* curr = streams + i;
            *curr = input.sub(currRange.first, currRange.second);
        }

        initTree(knownCommonLcp);
    }

    LcpCacheStringLoserTree(
        const Stream* inputs, size_t numInputs, lcp_t knownCommonLcp = 0)
    {
        assert(numInputs <= K);

        for (size_t i = 0; i < numInputs; i++)
        {
            streams[i] = inputs[i];
        }
        for (size_t i = numInputs; i < K; ++i)
        {
            streams[i] = Stream();
        }

        initTree(knownCommonLcp);
    }

    void initTree(lcp_t knownCommonLcp)
    {
        for (size_t i = 0; i < K; i++)
        {
            lcps[i] = knownCommonLcp;

            if (streams[i].size > 0)
                cached[i] = streams[i].firstString()[knownCommonLcp];

            size_t nodeIdx = K + i;
            size_t contenderIdx = i;

            while (nodeIdx % 2 == 1 && nodeIdx > 1)
            {
                nodeIdx >>= 1;
                updateNode(contenderIdx, nodes[nodeIdx]);
            }
            nodes[nodeIdx >> 1] = contenderIdx;
        }
    }

    void printTree()
    {
        size_t levelSize = 1;

        for (size_t i = 0; i < K; ++i)
        {
            if (i >= levelSize)
            {
                std::cout << "\n";
                levelSize *= 2;
            }
            std::cout << nodes[i] << " ";
        }
        std::cout << std::endl;

        for (size_t i = 0; i < K; ++i)
        {
            const Stream& stream = streams[i];
            if (stream.size > 0)
            {
                std::cout << lcps[i] << "|" << cached[i] << "|" << stream.firstString();
            }
            else
            {
                std::cout << -1;
            }
            std::cout << "(" << stream.size << ")  ";
        }

        std::cout << std::endl;
    }

    void replay(size_t& contenderIdx)
    {
#if 0
        for (size_t nodeIdx = (K + contenderIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
        {
            updateNode(contenderIdx, nodes[nodeIdx]);
        }
#else
        size_t nodeIdx = (K + contenderIdx) >> 1;

        if (K >= 256) {
            updateNode(contenderIdx, nodes[nodeIdx]);
            nodeIdx >>= 1;
        }
        if (K >= 128) {
            updateNode(contenderIdx, nodes[nodeIdx]);
            nodeIdx >>= 1;
        }
        if (K >= 64) {
            updateNode(contenderIdx, nodes[nodeIdx]);
            nodeIdx >>= 1;
        }
        if (K >= 32) {
            updateNode(contenderIdx, nodes[nodeIdx]);
            nodeIdx >>= 1;
        }
        if (K >= 16) {
            updateNode(contenderIdx, nodes[nodeIdx]);
            nodeIdx >>= 1;
        }
        if (K >= 8) {
            updateNode(contenderIdx, nodes[nodeIdx]);
            nodeIdx >>= 1;
        }
        if (K >= 4) {
            updateNode(contenderIdx, nodes[nodeIdx]);
            nodeIdx >>= 1;
        }
        if (K >= 2) {
            updateNode(contenderIdx, nodes[nodeIdx]);
        }
#endif
    }

    void writeElementsToStream(Stream outStream)
    {
        const Stream end = outStream.end();
        size_t contenderIdx = nodes[0];

        while (outStream < end)
        {
            assert(cached[contenderIdx] ==
                   streams[contenderIdx].firstString()[lcps[contenderIdx]]);

            outStream.setFirst(streams[contenderIdx].firstString(),
                               lcps[contenderIdx], cached[contenderIdx]);
            ++outStream;

            removeTopFromStream(contenderIdx);

            replay(contenderIdx);
        }

        nodes[0] = contenderIdx;
    }

    void writeElementsToStream(string* outStream, const size_t length)
    {
        const string* end = outStream + length;
        size_t contenderIdx = nodes[0];

        while (outStream < end)
        {
            *outStream = streams[contenderIdx].firstString();
            ++outStream;

            removeTopFromStream(contenderIdx);

            replay(contenderIdx);
        }

        nodes[0] = contenderIdx;
    }

    string deleteMin()
    {
        size_t contenderIdx = nodes[0];
        string min = streams[contenderIdx].firstString();
        removeTopFromStream(contenderIdx);
        replay(contenderIdx);

        return min;
    }

    void getRangesOfRemaining(
        std::pair<size_t, size_t>* ranges, const LcpStringPtr& inputBase)
    {
        for (size_t k = 0; k < K; k++)
        {
            ranges[k] = std::make_pair(streams[k] - inputBase, streams[k].length);
        }
    }

    const LcpCacheStringPtr * getRemaining()
    {
        return streams;
    }
};

} // namespace bingmann

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_LCP_LOSERTREE_HEADER

/******************************************************************************/
