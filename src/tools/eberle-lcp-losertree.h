/******************************************************************************
 * src/eberle/utils/lcp-string-losertree.h
 *
 * Implementation of a LCP aware multiway losertree.
 *
 ******************************************************************************
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
 *****************************************************************************/


#ifndef LCP_STRING_LOSERTREE_H_
#define LCP_STRING_LOSERTREE_H_

#include <stdlib.h>
#include <math.h>

#include "../tools/stringtools.h"

namespace eberle_lcp_utils
{

using namespace std;
using namespace stringtools;

//typedefs
typedef unsigned char* string;

//implementation follows

template<unsigned NUMBER_OF_STREAMS>
    class LcpStringLoserTree
    {
        struct STREAM
        {
            LcpStringPtr elements;
            unsigned length;
            bool isEmpty;
        };

    private:
        STREAM streams[NUMBER_OF_STREAMS];
        unsigned nodes[NUMBER_OF_STREAMS];
        lcp_t lcps[NUMBER_OF_STREAMS];

        /*
         * Returns the winner of all games.
         */
        inline unsigned
        updateNode(unsigned &defenderIdx, unsigned contenderIdx)
        {
            const STREAM* defenderStream = streams + defenderIdx;

            if (defenderStream->isEmpty)
            {
                return contenderIdx;
            }

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
                while (*s1 != '\0' && *s1 == *s2)
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
            for (unsigned i = 0; i < NUMBER_OF_STREAMS; i++)
            {
                lcps[i] = knownCommonLcp;
                unsigned nodeIdx = NUMBER_OF_STREAMS + i;
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
        LcpStringLoserTree(const LcpStringPtr& input, pair<size_t, size_t>* ranges, lcp_t knownCommonLcp = 0)
        {
            for (unsigned i = 0; i < NUMBER_OF_STREAMS; i++)
            {
                const pair<size_t, size_t> currRange = ranges[i];
                STREAM* curr = streams + i;
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

            for (unsigned i = 0; i < NUMBER_OF_STREAMS; ++i)
            {
                if (i >= levelSize)
                {
                    cout << "\n";
                    levelSize *= 2;
                }
                cout << nodes[i] << " ";
            }
            cout << endl;

            for (unsigned i = 0; i < NUMBER_OF_STREAMS; ++i)
            {
                STREAM* stream = streams + i;
                if (stream.length > 0)
                {
                    cout << lcps[i] << "|" << stream->elements.str();
                }
                else
                {
                    cout << -1;
                }
                cout << "(" << stream->length << ")  ";
            }

            cout << endl;
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

                for (unsigned nodeIdx = (NUMBER_OF_STREAMS + contenderIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
                {
                    contenderIdx = updateNode(nodes[nodeIdx], contenderIdx);
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

                for (unsigned nodeIdx = (NUMBER_OF_STREAMS + contenderIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
                {
                    contenderIdx = updateNode(nodes[nodeIdx], contenderIdx);
                }
            }

            nodes[0] = contenderIdx;
        }

        inline void
        getRangesOfRemaining(pair<size_t, size_t>* ranges, const LcpStringPtr& inputBase)
        {
            for (unsigned k = 0; k < NUMBER_OF_STREAMS; k++)
            {
                ranges[k] = make_pair(size_t(streams[k].elements - inputBase), streams[k].length);
            }
        }
    };

} // namespace eberle_lcp_utils

#endif // LCP_STRING_LOSERTREE_H_
