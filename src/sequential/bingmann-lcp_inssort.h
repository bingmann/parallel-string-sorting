/******************************************************************************
 * src/sequential/bingmann-lcp_inssort.h
 *
 * LCP-aware insertion sort, used in parallel variants as base case.
 *
 ******************************************************************************
 * Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
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

#include <string.h>

#include "../tools/stringtools.h"
#include "../tools/contest.h"

#ifndef BINGMANN_LCP_INSSORT_H_
#define BINGMANN_LCP_INSSORT_H_

namespace bingmann_lcp_inssort {

using namespace stringtools;

static inline
void lcp_insertion_sort(const StringPtr& strptr, size_t n, size_t depth)
{
    if (n <= 1) return;

    for (size_t j = 0; j < n - 1; ++j)
    {
        // insert strings[j] into sorted strings[0..j-1]

        string new_str = strptr.str(j);
        size_t new_lcp = depth; // start with LCP depth

        size_t i = j;
        while (i > 0)
        {
            size_t prev_lcp = new_lcp;

            string& cur_str = strptr.str(i - 1);
            size_t& cur_lcp = strptr.lcp(i);

            if (cur_lcp < new_lcp)
            {
                // CASE 1: lcp goes down -> insert string
                break;
            }
            else if (cur_lcp == new_lcp)
            {
                // CASE 2: compare more characters

                string s1 = new_str + new_lcp;
                string s2 = cur_str + new_lcp;

                while (*s1 != 0 && *s1 == *s2)
                    ++s1, ++s2, ++new_lcp;

                // if (new_str >= curr_str) -> insert string
                if (*s1 >= *s2)
                {
                    // update lcp of prev (smaller string) with inserted string
                    cur_lcp = new_lcp;
                    // lcp of inserted string with next string
                    new_lcp = prev_lcp;
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            strptr.str(i) = cur_str;
            strptr.lcp(i + 1) = cur_lcp;

            --i;
        }

        strptr.str(i) = new_str;
        strptr.lcp(i + 1) = new_lcp;
    }

    // last loop specialized with checks for out-of-bound access to lcp.
    {
        size_t j = n - 1;

        // insert strings[j] into sorted strings[0..j-1]

        string new_str = strptr.str(j);
        size_t new_lcp = depth; // start with LCP depth

        size_t i = j;
        while (i > 0)
        {
            size_t prev_lcp = new_lcp;

            string& cur_str = strptr.str(i - 1);
            size_t& cur_lcp = strptr.lcp(i);

            if (cur_lcp < new_lcp)
            {
                // CASE 1: lcp goes down -> insert string
                break;
            }
            else if (cur_lcp == new_lcp)
            {
                // CASE 2: compare more characters

                string s1 = new_str + new_lcp;
                string s2 = cur_str + new_lcp;

                while (*s1 != 0 && *s1 == *s2)
                    ++s1, ++s2, ++new_lcp;

                // if (new_str >= curr_str) -> insert string
                if (*s1 >= *s2)
                {
                    // update lcp of prev (smaller string) with inserted string
                    cur_lcp = new_lcp;
                    // lcp of inserted string with next string
                    new_lcp = prev_lcp;
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            strptr.str(i) = cur_str;

            if (i + 1 < n) // check out-of-bounds copy
                strptr.lcp(i + 1) = cur_lcp;

            --i;
        }

        strptr.str(i) = new_str;

        if (i + 1 < n) // check out-of-bounds save
            strptr.lcp(i + 1) = new_lcp;
    }
}

static inline
void do_lcp_insertion_sort(string* strings, size_t n)
{
    string* shadow = new string[n]; // allocate shadow pointer array
    StringPtr strptr(strings, shadow, false);

    strptr.lcp(0) = 42; // must keep lcp[0] unchanged

    lcp_insertion_sort(strptr, n, 0);

    std::cout << "lcp[0] " << strptr.lcp(0) << "\n";

    for (size_t j = 1; j < n; ++j)
    {
        string s1 = strptr.str(j-1), s2 = strptr.str(j);
        size_t h = 0;
        while (*s1 != 0 && *s1 == *s2)
            ++h, ++s1, ++s2;

        if (h != strptr.lcp(j)) {
            std::cout << "lcp[" << j << "] mismatch " << h << " != " << strptr.lcp(j) << std::endl;
        }
    }

    delete [] shadow;
}

CONTESTANT_REGISTER(do_lcp_insertion_sort, "bingmann/lcp_insertion_sort",
                    "bingmann/lcp_insertion_sort LCP-aware insertion sort")

} // namespace bingmann_lcp_inssort

#endif // BINGMANN_LCP_INSSORT_H_
