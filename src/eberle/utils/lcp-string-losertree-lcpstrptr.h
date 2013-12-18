/*
 * lcp-string-losertree.h
 *
 *  Created on: Oct 30, 2013
 *      Author: aeberle
 */

#ifndef LCP_STRING_LOSERTREE_LCPPTR_H_
#define LCP_STRING_LOSERTREE_LCPPTR_H_

#include <stdlib.h>
#include <math.h>

#include "types.h"

namespace eberle_lcp_utils
{

using namespace std;
using namespace types;

//typedefs
typedef unsigned char* string;
typedef unsigned int UINT;

//implementation follows

template<unsigned NUMBER_OF_STREAMS>
    class LcpStringLcpPtrLoserTree
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
        unsigned lcps[NUMBER_OF_STREAMS];

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

            unsigned* contenderLcp = lcps + contenderIdx;
            unsigned* defenderLcp = lcps + defenderIdx;

            if (contenderStream->isEmpty || *defenderLcp > *contenderLcp)
            { // CASE 2: curr->lcp > contender->lcp => curr < contender
                std::swap(defenderIdx, contenderIdx);

            }
            else if (*defenderLcp == *contenderLcp)
            { // CASE 1: curr.lcp == contender.lcp
                unsigned lcp = *defenderLcp;
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
        initTree(unsigned knownCommonLcp)
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
        removeTopFromStream(unsigned streamIdx, const LcpStringPtr& outStream)
        {
            STREAM* stream = streams + streamIdx;
            outStream.set(stream->elements.str(), lcps[streamIdx]);

            stream->length--;
            ++(stream->elements);
            stream->isEmpty = stream->length <= 0;

            if (!stream->isEmpty)
            {
                lcps[streamIdx] = stream->elements.lcp();
            }
        }

    public:
        LcpStringLcpPtrLoserTree(const LcpStringPtr& input, pair<size_t, size_t>* ranges, unsigned knownCommonLcp = 0)
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
                removeTopFromStream(contenderIdx, outStream);
                ++outStream;

                for (unsigned nodeIdx = (NUMBER_OF_STREAMS + contenderIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
                {
                    contenderIdx = updateNode(nodes[nodeIdx], contenderIdx);
                }
            }

            nodes[0] = contenderIdx;
        }

        /*
         inline void
         writeElementsToStream(string *outStream, const size_t length)
         {
         const string* end = outStream + length;
         unsigned contenderIdx = nodes[0];

         while (outStream < end)
         {
         *outStream = removeTopFromStream(contenderIdx)->text;
         outStream++;

         for (unsigned nodeIdx = (NUMBER_OF_STREAMS + contenderIdx) >> 1; nodeIdx >= 1; nodeIdx >>= 1)
         {
         contenderIdx = updateNode(nodes[nodeIdx], contenderIdx);
         }
         }

         nodes[0] = contenderIdx;
         }
         */

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

#endif // LCP_STRING_LOSERTREE_LCPPTR_H_
