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
    typedef LcpStringPtr Stream;

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
        { // CASE 2: curr->lcp > contender->lcp => curr < contender
            std::swap(defenderIdx, contenderIdx);

        }
        else if (defenderLcp == contenderLcp)
        { // CASE 1: curr.lcp == contender.lcp
            lcp_t lcp = defenderLcp;

            string s1 = defenderStream.str() + lcp;
            string s2 = contenderStream.str() + lcp;

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
    removeTopFromStream(unsigned streamIdx)
    {
        Stream& stream = streams[streamIdx];

        ++stream;

        if (!stream.empty())
        {
            lcps[streamIdx] = stream.lcp();
        }
    }

public:
    LcpStringLoserTree()
    {
    }

    LcpStringLoserTree(const LcpStringPtr& input, std::pair<size_t, size_t>* ranges, lcp_t knownCommonLcp = 0)
    {
        for (unsigned i = 0; i < K; i++)
        {
            const std::pair<size_t, size_t> currRange = ranges[i];
            Stream* curr = streams + i;
            *curr = input.sub(currRange.first, currRange.second);
        }

        initTree(knownCommonLcp);
    }

    LcpStringLoserTree(const LcpStringPtr* inputs, unsigned numInputs, lcp_t knownCommonLcp = 0)
    {
        assert(numInputs <= K);

        for (unsigned i = 0; i < numInputs; i++)
        {
            streams[i] = inputs[i];
        }
        for (unsigned i = numInputs; i < K; ++i)
        {
            streams[i] = LcpStringPtr(NULL, NULL, 0);
        }

        initTree(knownCommonLcp);
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
            if (stream->size > 0)
            {
                std::cout << lcps[i] << "|" << stream->str();
            }
            else
            {
                std::cout << -1;
            }
            std::cout << "(" << stream->size << ")  ";
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
    writeElementsToStream(LcpStringPtr outStream)
    {
        const LcpStringPtr end = outStream.end();
        unsigned contenderIdx = nodes[0];

        while (outStream < end)
        {
            outStream.set(streams[contenderIdx].str(), lcps[contenderIdx]);
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
            *outStream = streams[contenderIdx].str();
            ++outStream;

            removeTopFromStream(contenderIdx);

            replay(contenderIdx);
        }

        nodes[0] = contenderIdx;
    }

    inline void
    getRangesOfRemaining(std::pair<size_t, size_t>* ranges, const LcpStringPtr& inputBase)
    {
        for (unsigned k = 0; k < K; k++)
        {
            ranges[k] = std::make_pair(size_t(streams[k] - inputBase), streams[k].length);
        }
    }

    inline const LcpStringPtr*
    getRemaining()
    {
        return streams;
    }
};

} // namespace eberle_lcp_utils

#endif // LCP_STRING_LOSERTREE_H_
