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
    struct Stream
    {
        LcpStringPtr elements;
        unsigned length;
        bool isEmpty;
    };

private:
    Stream streams[K];
    unsigned nodes[K];
    lcp_t lcps[K];

    /*
     * Returns the winner of all games.
     */
    inline void
    updateNode(unsigned &defenderIdx, unsigned &contenderIdx)
    {
        const Stream& defenderStream = streams[defenderIdx];

        if (defenderStream.isEmpty)
            return;

        const Stream& contenderStream = streams[contenderIdx];

        lcp_t& contenderLcp = lcps[contenderIdx];
        lcp_t& defenderLcp = lcps[defenderIdx];

        if (contenderStream.isEmpty || defenderLcp > contenderLcp)
        { // CASE 2: curr->lcp > contender->lcp => curr < contender
            std::swap(defenderIdx, contenderIdx);

        }
        else if (defenderLcp == contenderLcp)
        { // CASE 1: curr.lcp == contender.lcp
            lcp_t lcp = defenderLcp;

            string s1 = defenderStream.elements.str() + lcp;
            string s2 = contenderStream.elements.str() + lcp;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != '\0' && *s1 == *s2)
                s1++, s2++, lcp++;

            if (*s1 < *s2)
            { 	// CASE 1.1: curr < contender
                contenderLcp = lcp;
                std::swap(defenderIdx, contenderIdx);
            }
            else
            {	// CASE 1.2: curr >= contender
                defenderLcp = lcp;
            }
        } // else // CASE 3: curr->lcp < contender->lcp => contender < curr  => nothing to do
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
                updateNode(nodes[nodeIdx], contenderIdx);
            }
            nodes[nodeIdx >> 1] = contenderIdx;
        }
    }

    inline void
    removeTopFromStream(unsigned streamIdx)
    {
        Stream& stream = streams[streamIdx];

        stream.length--;
        ++(stream.elements);
        stream.isEmpty = (stream.length <= 0);

        if (!stream.isEmpty)
        {
            lcps[streamIdx] = stream.elements.lcp();
        }
    }

public:
    LcpStringLoserTree(const LcpStringPtr& input, std::pair<size_t, size_t>* ranges, lcp_t knownCommonLcp = 0)
    {
        for (unsigned i = 0; i < K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i];
            Stream* curr = streams + i;
            curr->elements = input + currRange.first;
            curr->length = currRange.second;
            curr->isEmpty = (curr->length <= 0);
        }

        initTree(knownCommonLcp);
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
            Stream* stream = streams + i;
            if (stream.length > 0)
            {
                std::cout << lcps[i] << "|" << stream->elements.str();
            }
            else
            {
                std::cout << -1;
            }
            std::cout << "(" << stream->length << ")  ";
        }

        std::cout << std::endl;
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
                updateNode(nodes[nodeIdx], contenderIdx);
            }
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
            *outStream = streams[contenderIdx].elements.str();
            ++outStream;

            removeTopFromStream(contenderIdx);

            for (unsigned nodeIdx = (K + contenderIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
            {
                updateNode(nodes[nodeIdx], contenderIdx);
            }
        }

        nodes[0] = contenderIdx;
    }

    inline void
    getRangesOfRemaining(std::pair<size_t, size_t>* ranges, const LcpStringPtr& inputBase)
    {
        for (unsigned k = 0; k < K; k++)
        {
            ranges[k] = std::make_pair(size_t(streams[k].elements - inputBase), streams[k].length);
        }
    }
};

} // namespace eberle_lcp_utils

#endif // LCP_STRING_LOSERTREE_H_
