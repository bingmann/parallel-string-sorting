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

// Pretty tricky template magic here: we must templatize the lcp array in the
// lcp_insertion_sort() depending on whether the sorted StringPtrType contains
// a designated LCP array or not (keep temporary) one.

template <typename StringPtrType>
struct LcpArray;

template <>
struct LcpArray<StringPtr>
{
    const StringPtr& m_strptr;

    LcpArray<StringPtr>(const StringPtr& strptr)
        : m_strptr(strptr)
    { }

    uintptr_t& operator() (size_t i)
    {
        return m_strptr.lcp(i);
    }
};

template <>
struct LcpArray<StringPtrOut>
{
    const StringPtrOut& m_strptr;

    LcpArray<StringPtrOut>(const StringPtrOut& strptr)
        : m_strptr(strptr)
    { }

    uintptr_t& operator() (size_t i)
    {
        return m_strptr.lcp(i);
    }
};

template <>
struct LcpArray<StringPtrNoLcpCalc>
{
    std::vector<uintptr_t> m_lcp;

    LcpArray<StringPtrNoLcpCalc>(const StringPtrNoLcpCalc& strptr)
        : m_lcp(strptr.size())
    { }

    uintptr_t& operator() (size_t i)
    {
        assert(i < m_lcp.size());
        return m_lcp[i];
    }
};

template <>
struct LcpArray<StringPtrOutNoLcpCalc>
{
    std::vector<uintptr_t> m_lcp;

    LcpArray<StringPtrOutNoLcpCalc>(const StringPtrOutNoLcpCalc& strptr)
        : m_lcp(strptr.size())
    { }

    uintptr_t& operator() (size_t i)
    {
        assert(i < m_lcp.size());
        return m_lcp[i];
    }
};

////////////////////////////////////////////////////////////////////////////////

//! LCP insertion sort
template <typename StringPtrType>
static inline
void lcp_insertion_sort(const StringPtrType& strptr, size_t depth)
{
    size_t n = strptr.size();

    // maybe allocate temporary lcp array
    LcpArray<StringPtrType> lcp(strptr);

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

            string cur_str = strptr.out(i - 1);
            size_t cur_lcp = lcp(i);

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
                    lcp(i) = new_lcp;
                    // lcp of inserted string with next string
                    new_lcp = prev_lcp;
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            strptr.out(i) = cur_str;
            lcp(i + 1) = cur_lcp;

            --i;
        }

        strptr.out(i) = new_str;
        lcp(i + 1) = new_lcp;
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

            string cur_str = strptr.out(i - 1);
            size_t cur_lcp = lcp(i);

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
                    lcp(i) = new_lcp;
                    // lcp of inserted string with next string
                    new_lcp = prev_lcp;
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            strptr.out(i) = cur_str;

            if (i + 1 < n) // check out-of-bounds copy
                lcp(i + 1) = cur_lcp;

            --i;
        }

        strptr.out(i) = new_str;

        if (i + 1 < n) // check out-of-bounds save
            lcp(i + 1) = new_lcp;
    }
}

#if 0
static inline
void lcp_insertion_sort(const StringPtrNoLcpCalc& strptr, size_t depth)
{
    size_t n = strptr.size();
    size_t tmp_lcp[n];

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
            size_t& cur_lcp = tmp_lcp[i];

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
            tmp_lcp[i + 1] = cur_lcp;

            --i;
        }

        strptr.str(i) = new_str;
        tmp_lcp[i + 1] = new_lcp;
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
            size_t& cur_lcp = tmp_lcp[i];

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
                tmp_lcp[i + 1] = cur_lcp;

            --i;
        }

        strptr.str(i) = new_str;

        if (i + 1 < n) // check out-of-bounds save
            tmp_lcp[i + 1] = new_lcp;
    }
}
#endif

static inline
void test_lcp_insertion_sort(string* strings, size_t n)
{
    string* shadow = new string[n]; // allocate shadow pointer array
    StringPtr strptr(strings, shadow, n);

    strptr.lcp(0) = 42; // must keep lcp[0] unchanged

    lcp_insertion_sort(strptr, 0);

    std::cout << "lcp[0] " << strptr.lcp(0) << "\n";

    for (size_t j = 1; j < n; ++j)
    {
        string s1 = strptr.out(j-1), s2 = strptr.out(j);
        size_t h = calc_lcp(s1, s2);

        if (h != strptr.lcp(j)) {
            std::cout << "lcp[" << j << "] mismatch " << h << " != " << strptr.lcp(j) << std::endl;
        }
    }

    delete [] shadow;
}

static inline
void test_lcp_insertion_sort_nolcp(string* strings, size_t n)
{
    string* shadow = new string[n]; // allocate shadow pointer array
    StringPtrNoLcpCalc strptr(strings, shadow, n);

    lcp_insertion_sort(strptr, 0);

    delete [] shadow;
}

CONTESTANT_REGISTER(test_lcp_insertion_sort, "bingmann/lcp_insertion_sort",
                    "LCP-aware insertion sort")

CONTESTANT_REGISTER(test_lcp_insertion_sort_nolcp, "bingmann/lcp_insertion_sort_nolcp",
                    "LCP-aware insertion sort (without LCP output)")

static inline
void test_lcp_insertion_sort_out(string* strings, size_t n)
{
    string* shadow = new string[n]; // allocate shadow pointer array
    string* output = new string[n]; // separate output array
    StringPtrOut strptr(strings, shadow, output, n);

    strptr.lcp(0) = 42; // must keep lcp[0] unchanged

    lcp_insertion_sort(strptr, 0);

    std::cout << "lcp[0] " << strptr.lcp(0) << "\n";

    for (size_t j = 1; j < n; ++j)
    {
        string s1 = strptr.out(j-1), s2 = strptr.out(j);
        size_t h = calc_lcp(s1, s2);

        if (h != strptr.lcp(j)) {
            std::cout << "lcp[" << j << "] mismatch " << h << " != " << strptr.lcp(j) << std::endl;
        }
    }

    delete [] shadow;

    // copy output back to original array for checking
    memcpy(strings, output, n * sizeof(string));
    delete [] output;
}

static inline
void test_lcp_insertion_sort_out_nolcp(string* strings, size_t n)
{
    string* shadow = new string[n]; // allocate shadow pointer array
    string* output = new string[n]; // separate output array
    StringPtrOutNoLcpCalc strptr(strings, shadow, output, n);

    lcp_insertion_sort(strptr, 0);

    delete [] shadow;

    // copy output back to original array for checking
    memcpy(strings, output, n * sizeof(string));
    delete [] output;
}

CONTESTANT_REGISTER(test_lcp_insertion_sort_out, "bingmann/lcp_insertion_sort_out",
                    "LCP-aware insertion sort")

CONTESTANT_REGISTER(test_lcp_insertion_sort_out_nolcp, "bingmann/lcp_insertion_sort_out_nolcp",
                    "LCP-aware insertion sort (without LCP output)")

} // namespace bingmann_lcp_inssort

#endif // BINGMANN_LCP_INSSORT_H_
