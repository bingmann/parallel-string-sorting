/******************************************************************************
 * src/eberle/utils/lcp-string-losertree.h
 *
 * Implementation of a LCP aware multiway losertree.
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


#ifndef LCP_STRING_LOSERTREE_H_
#define LCP_STRING_LOSERTREE_H_

#include <iostream>
#include <algorithm>
#include <utility>

#include "../tools/stringtools.h"

namespace eberle_lcp_utils
{

using namespace stringtools;

//typedefs
typedef unsigned char* string;

//implementation follows

template<unsigned K>
class LcpStringLoserTree
{
    typedef LcpCacheStringPtr Stream;

public:
    //! allow public access to streams for easier initialization
    Stream streams[K];

protected:
    unsigned nodes[K];
    lcp_t lcps[K];

    /*
     * Returns the winner of all games.
     */
    inline void
    updateNode(unsigned &defenderIdx, unsigned &contenderIdx)
    {
        const Stream& defenderStream = streams[defenderIdx];

        if (defenderStream.empty())
            return;

        const Stream& contenderStream = streams[contenderIdx];

        lcp_t& contenderLcp = lcps[contenderIdx];
        lcp_t& defenderLcp = lcps[defenderIdx];


        if (contenderStream.empty() || defenderLcp > contenderLcp)
        { // CASE 2: defender->lcp > contender->lcp => defender < contender
            std::swap(defenderIdx, contenderIdx);
        }

        else if (defenderLcp == contenderLcp)
        { // CASE 1: defender.lcp == contender.lcp
            lcp_t lcp = defenderLcp;

            char_type c1 = defenderStream.firstCached();
            char_type c2 = contenderStream.firstCached();

            assert(c1 == defenderStream.firstString()[lcp]);
            assert(c2 == contenderStream.firstString()[lcp]);

            // check the strings starting after lcp and calculate new lcp
            while (c1 != 0 && c1 == c2) {
                lcp++;
                c1 = defenderStream.firstString()[lcp];
                c2 = contenderStream.firstString()[lcp];
            }

            if (c1 < c2)
            { 	// CASE 1.1: defender < contender
                contenderLcp = lcp;
                contenderStream.firstCached() = c2;
                std::swap(defenderIdx, contenderIdx);
            }
            else
            {	// CASE 1.2: defender >= contender
                defenderLcp = lcp;
                defenderStream.firstCached() = c1;
            }
        }
        // else
        // CASE 3: defender->lcp < contender->lcp => contender < defender  => nothing to do
    }

    inline void
    removeTopFromStream(unsigned streamIdx)
    {
        Stream& stream = streams[streamIdx];

        ++stream;

        if (!stream.empty())
        {
            lcps[streamIdx] = stream.firstLcp();
        }
    }

public:
    LcpStringLoserTree()
    {
    }

    LcpStringLoserTree(const Stream& input, std::pair<size_t, size_t>* ranges, lcp_t knownCommonLcp = 0)
    {
        for (unsigned i = 0; i < K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i];
            Stream* curr = streams + i;
            *curr = input.sub(currRange.first, currRange.second);
        }

        initTree(knownCommonLcp);
    }

    LcpStringLoserTree(const Stream* inputs, unsigned numInputs, lcp_t knownCommonLcp = 0)
    {
        assert(numInputs <= K);

        for (unsigned i = 0; i < numInputs; i++)
        {
            streams[i] = inputs[i];
        }
        for (unsigned i = numInputs; i < K; ++i)
        {
            streams[i] = Stream();
        }

        initTree(knownCommonLcp);
    }

    inline void
    initTree(lcp_t knownCommonLcp)
    {
        for (unsigned i = 0; i < K; i++)
        {
            lcps[i] = knownCommonLcp;

            if(streams[i].size > 0)
                streams[i].firstCached() = streams[i].firstString()[knownCommonLcp];

            unsigned nodeIdx = K + i;
            unsigned contenderIdx = i;

            while (nodeIdx % 2 == 1 && nodeIdx > 1)
            {
                nodeIdx >>= 1;
                updateNode(nodes[nodeIdx], contenderIdx);
            }
            nodes[nodeIdx >> 1] = contenderIdx;
        }
    }

    void
    printTree()
    {
        unsigned levelSize = 1;

        for (unsigned i = 0; i < K; ++i)
        {
            if (i >= levelSize)
            {
                std::cout << "\n";
                levelSize *= 2;
            }
            std::cout << nodes[i] << " ";
        }
        std::cout << std::endl;

        for (unsigned i = 0; i < K; ++i)
        {
            const Stream& stream = streams [i];
            if (stream.size > 0)
            {
                std::cout << lcps[i] << "|" << stream.firstString();
            }
            else
            {
                std::cout << -1;
            }
            std::cout << "(" << stream.size << ")  ";
        }

        std::cout << std::endl;
    }

    inline void
    replay(unsigned& contenderIdx)
    {
#if 0
        for (unsigned nodeIdx = (K + contenderIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
        {
            updateNode(nodes[nodeIdx], contenderIdx);
        }
#else
        unsigned nodeIdx = (K + contenderIdx) >> 1;

        if (K >= 256) {
            updateNode(nodes[nodeIdx], contenderIdx);
            nodeIdx >>= 1;
        }
        if (K >= 128) {
            updateNode(nodes[nodeIdx], contenderIdx);
            nodeIdx >>= 1;
        }
        if (K >= 64) {
            updateNode(nodes[nodeIdx], contenderIdx);
            nodeIdx >>= 1;
        }
        if (K >= 32) {
            updateNode(nodes[nodeIdx], contenderIdx);
            nodeIdx >>= 1;
        }
        if (K >= 16) {
            updateNode(nodes[nodeIdx], contenderIdx);
            nodeIdx >>= 1;
        }
        if (K >= 8) {
            updateNode(nodes[nodeIdx], contenderIdx);
            nodeIdx >>= 1;
        }
        if (K >= 4) {
            updateNode(nodes[nodeIdx], contenderIdx);
            nodeIdx >>= 1;
        }
        if (K >= 2) {
            updateNode(nodes[nodeIdx], contenderIdx);
        }
#endif
    }

    inline void
    writeElementsToStream(Stream outStream)
    {
        const Stream end = outStream.end();
        unsigned contenderIdx = nodes[0];

        while (outStream < end)
        {
            Stream currStream = streams[contenderIdx];
            outStream.setFirst(currStream.firstString(), lcps[contenderIdx], currStream.firstCached());
            ++outStream;

            removeTopFromStream(contenderIdx);

            replay(contenderIdx);
        }

        nodes[0] = contenderIdx;
    }
 
    inline void
    writeElementsToStream(string *outStream, const size_t length)
    {
        const string* end = outStream + length;
        unsigned contenderIdx = nodes[0];

        while (outStream < end)
        {
            *outStream = streams[contenderIdx].firstString();
            ++outStream;

            removeTopFromStream(contenderIdx);

            replay(contenderIdx);
        }

        nodes[0] = contenderIdx;
    }

    inline string
    deleteMin()
    {
        unsigned contenderIdx = nodes[0];
        string min = streams[contenderIdx].firstString();
        removeTopFromStream(contenderIdx);
        replay(contenderIdx);

        return min;
    }

    inline void
    getRangesOfRemaining(std::pair<size_t, size_t>* ranges, const LcpStringPtr& inputBase)
    {
        for (unsigned k = 0; k < K; k++)
        {
            ranges[k] = std::make_pair(size_t(streams[k] - inputBase), streams[k].length);
        }
    }

    inline const LcpCacheStringPtr*
    getRemaining()
    {
        return streams;
    }
}; 
 
} // namespace eberle_lcp_utils

#endif // LCP_STRING_LOSERTREE_H_
