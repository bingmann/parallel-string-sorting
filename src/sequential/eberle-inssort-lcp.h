/*******************************************************************************
 * src/sequential/eberle-inssort-lcp.h
 *
 * LCP aware insertion sort.
 *
 *******************************************************************************
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
 ******************************************************************************/

#ifndef PSS_SRC_SEQUENTIAL_EBERLE_INSSORT_LCP_HEADER
#define PSS_SRC_SEQUENTIAL_EBERLE_INSSORT_LCP_HEADER

#include <iostream>
#include "../tools/stringtools.h"

//#define EBERLE_INSSORT_CHECK_LCPS

namespace eberle_inssort_lcp {

using namespace stringtools;

//typedefs
typedef unsigned char* string;

//implementations
static inline
void inssort_lcp(string* strings, const LcpStringPtr& out, size_t length)
{
    string* output = out.strings;
    lcp_t* lcps = out.lcps;

    for (size_t n = 0; n < length; n++)
    {
        const string candidateText = strings[n];
        lcp_t candidateLcp = 0;

        size_t insIdx = 0;
        for ( ; insIdx < n; insIdx++)
        {
            lcp_t currLcp = lcps[insIdx];

            if (candidateLcp == currLcp)
            {       // CASE 1 lcps are equal
                string s1 = candidateText + candidateLcp;
                string s2 = output[insIdx] + candidateLcp;

                // check the strings starting after lcp and calculate new lcp
                while (*s1 != '\0' && *s1 == *s2)
                    s1++, s2++;

                const lcp_t lcp = s1 - candidateText;

                if (*s1 <= *s2)
                {       // CASE 1.1: candidate <= curr => insert it
                    lcps[insIdx] = lcp;
                    break;
                }
                else
                {            // CASE 1.2: candidate > curr
                    candidateLcp = lcp;
                }
            }
            else if (candidateLcp > currLcp)
            {       // CASE 2: candidate < curr => insert it
                break;
            } //  CASE 3: candidate > curr => nothing to do
        }

        // move the tail one back
        for (size_t i = n; i > insIdx; i--)
        {
            output[i] = output[i - 1];
            lcps[i] = lcps[i - 1];
        }

        // insert the new element
        lcps[insIdx] = candidateLcp;
        output[insIdx] = candidateText;
    }
}

static inline
void inssort_lcp(string* strings, const LcpCacheStringPtr& out, size_t length)
{
    inssort_lcp(strings, (const LcpStringPtr&)out, length);
    out.calculateCache();
}

static inline
void eberle_lcp_inssort(string* strings, size_t n)
{
    lcp_t* lcps = new lcp_t[n];
    LcpStringPtr output(strings, lcps, n);

    inssort_lcp(strings, output, n);

    //check lcps
#ifdef EBERLE_INSSORT_CHECK_LCPS
    std::cout << "Checking LCPs" << std::endl;
    stringtools::verify_lcp(output.strings, output.lcps, n, 0);
#endif

    delete[] lcps;
}

static inline
void eberle_lcp_inssort_cache(string* strings, size_t n)
{
    lcp_t* lcps = new lcp_t[n];
    char_type* cache = new char_type[n];
    LcpCacheStringPtr output(strings, lcps, cache, n);

    inssort_lcp(strings, output, n);

    //check lcps
#ifdef EBERLE_INSSORT_CHECK_LCPS
    std::cout << "Checking LCPs" << std::endl;
    stringtools::verify_lcp_cache(output.strings, output.lcps, output.cachedChars, n, 0);
#endif

    delete[] lcps;
    delete[] cache;
}

PSS_CONTESTANT(eberle_lcp_inssort, "eberle/lcp_insertion_sort", "LCP aware inssertion sort by Andreas Eberle")
PSS_CONTESTANT(eberle_lcp_inssort_cache, "eberle/lcp_insertion_sort_cache", "LCP aware insertion sort with cached characters calculation by Andreas Eberle")

} // namespace eberle_inssort_lcp

#endif // !PSS_SRC_SEQUENTIAL_EBERLE_INSSORT_LCP_HEADER

/******************************************************************************/
