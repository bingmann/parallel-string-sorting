/*******************************************************************************
 * src/sequential/bingmann-lcp_inssort.hpp
 *
 * LCP-aware insertion sort, used in parallel variants as base case.
 *
 *******************************************************************************
 * Copyright (C) 2013-2015 Timo Bingmann <tb@panthema.net>
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

#ifndef PSS_SRC_SEQUENTIAL_BINGMANN_LCP_INSSORT_HEADER
#define PSS_SRC_SEQUENTIAL_BINGMANN_LCP_INSSORT_HEADER

#include <cstring>

#include "../tools/stringtools.hpp"
#include "../tools/contest.hpp"

namespace bingmann_lcp_inssort {

using namespace stringtools;

//! LCP insertion sort
template <typename StringSet>
static inline
void lcp_insertion_sort(const StringSet& str, uintptr_t* lcp, size_t depth)
{
    typedef typename StringSet::String String;
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::CharIterator CharIterator;

    Iterator begin = str.begin();
    size_t n = str.size();

    if (n <= 1) return;

    for (size_t j = 0; j < n - 1; ++j)
    {
        // insert strings[j] into sorted strings[0..j-1]

        String new_str = std::move(str[begin + j]);
        size_t new_lcp = depth; // start with LCP depth

        size_t i = j;
        while (i > 0)
        {
            size_t prev_lcp = new_lcp;

            String cur_str = std::move(str[begin + i - 1]);
            size_t cur_lcp = lcp[i];

            if (cur_lcp < new_lcp)
            {
                // CASE 1: lcp goes down -> insert string

                // move comparison string back
                str[begin + i - 1] = std::move(cur_str);
                break;
            }
            else if (cur_lcp == new_lcp)
            {
                // CASE 2: compare more characters

                CharIterator c1 = str.get_chars(new_str, new_lcp);
                CharIterator c2 = str.get_chars(cur_str, new_lcp);

                while (str.is_equal(new_str, c1, cur_str, c2))
                    ++c1, ++c2, ++new_lcp;

                // if (new_str >= curr_str) -> insert string
                if (!str.is_less(new_str, c1, cur_str, c2))
                {
                    // update lcp of prev (smaller string) with inserted string
                    lcp[i] = new_lcp;
                    // lcp of inserted string with next string
                    new_lcp = prev_lcp;

                    // move comparison string back
                    str[begin + i - 1] = std::move(cur_str);
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            str[begin + i] = std::move(cur_str);
            lcp[i + 1] = cur_lcp;

            --i;
        }

        str[begin + i] = std::move(new_str);
        lcp[i + 1] = new_lcp;
    }

    // last loop specialized with checks for out-of-bound access to lcp.
    {
        size_t j = n - 1;

        // insert strings[j] into sorted strings[0..j-1]

        String new_str = std::move(str[begin + j]);
        size_t new_lcp = depth; // start with LCP depth

        size_t i = j;
        while (i > 0)
        {
            size_t prev_lcp = new_lcp;

            String cur_str = std::move(str[begin + i - 1]);
            size_t cur_lcp = lcp[i];

            if (cur_lcp < new_lcp)
            {
                // CASE 1: lcp goes down -> insert string

                // move comparison string back
                str[begin + i - 1] = std::move(cur_str);
                break;
            }
            else if (cur_lcp == new_lcp)
            {
                // CASE 2: compare more characters

                CharIterator c1 = str.get_chars(new_str, new_lcp);
                CharIterator c2 = str.get_chars(cur_str, new_lcp);

                while (str.is_equal(new_str, c1, cur_str, c2))
                    ++c1, ++c2, ++new_lcp;

                // if (new_str >= curr_str) -> insert string
                if (!str.is_less(new_str, c1, cur_str, c2))
                {
                    // update lcp of prev (smaller string) with inserted string
                    lcp[i] = new_lcp;
                    // lcp of inserted string with next string
                    new_lcp = prev_lcp;

                    // move comparison string back
                    str[begin + i - 1] = std::move(cur_str);
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            str[begin + i] = std::move(cur_str);

            if (i + 1 < n) // check out-of-bounds copy
                lcp[i + 1] = cur_lcp;

            --i;
        }

        str[begin + i] = std::move(new_str);

        if (i + 1 < n) // check out-of-bounds save
            lcp[i + 1] = new_lcp;
    }
}

//! LCP insertion sort, plain arguments version.
static inline
void lcp_insertion_sort(string* str, uintptr_t* lcp, size_t n, size_t depth)
{
    return lcp_insertion_sort(
        parallel_string_sorting::UCharStringSet(str, str + n), lcp, depth);
}

//! LCP insertion sort, but immediately discard the lcp
template <typename StringSet>
static inline
void lcp_insertion_sort_nolcp(const StringSet& ss, size_t depth)
{
    uintptr_t tmp_lcp[ss.size()];
    return lcp_insertion_sort(ss, tmp_lcp, depth);
}

//! LCP insertion sort, but immediately discard the lcp
template <typename StringSet>
static inline
void lcp_insertion_sort_verify(const StringSet& ss, size_t depth)
{
    std::vector<uintptr_t> tmp_lcp(ss.size());
    tmp_lcp[0] = 42;                 // must keep lcp[0] unchanged
    std::fill(tmp_lcp.begin() + 1, tmp_lcp.end(), -1);
    lcp_insertion_sort(ss, tmp_lcp.data(), depth);
    assert(stringtools::verify_lcp(ss, tmp_lcp.data(), 42));
}

//! LCP insertion sort, close to journal's pseudo-code
template <typename StringSet>
static inline
void lcp_insertion_sort_pseudocode(
    const StringSet& str, uintptr_t* lcp, size_t depth)
{
    typedef typename StringSet::String String;
    typedef typename StringSet::Iterator Iterator;

    Iterator begin = str.begin();
    size_t n = str.size();

    unsigned int cmp = 0;

    for (size_t j = 0; j < n; ++j)
    {
        // insert strings[j] into sorted strings[0..j-1]

        String snew = std::move(str[begin + j]);
        size_t h = depth; // start with LCP depth

        size_t i = j;
        while (i > 0)
        {
            if (lcp[i] < h)
            {
                // CASE 1: lcp goes down -> insert string
                break;
            }
            else if (lcp[i] == h)
            {
                // CASE 2: compare more characters

                size_t prev_lcp = h;

                const String& s2 = std::move(str[begin + i - 1]);

                while (cmp++, str.is_equal(snew, str.get_chars(snew, h),
                                           s2, str.get_chars(s2, h)))
                    ++h;

                // if (new_str >= curr_str) -> insert string
                if (!str.is_less(snew, str.get_chars(snew, h),
                                 s2, str.get_chars(s2, h)))
                {
                    // update lcp of prev (smaller string) with inserted string
                    lcp[i] = h;
                    // lcp of inserted string with next string
                    h = prev_lcp;
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            str[begin + i] = std::move(str[begin + i - 1]);

            if (i + 1 < n) // check out-of-bounds save
                lcp[i + 1] = lcp[i];

            --i;
        }

        str[begin + i] = std::move(snew);
        if (i + 1 < n) // check out-of-bounds save
            lcp[i + 1] = h;
    }

    std::cout << "lcp_inssort comparisons = " << cmp << "\n";
}

//! LCP insertion sort, but immediately discard the lcp
template <typename StringSet>
static inline
void lcp_insertion_sort_pseudocode_verify(const StringSet& ss, size_t depth)
{
    std::vector<uintptr_t> tmp_lcp(ss.size());
    tmp_lcp[0] = 42;                 // must keep lcp[0] unchanged
    std::fill(tmp_lcp.begin() + 1, tmp_lcp.end(), -1);
    lcp_insertion_sort_pseudocode(ss, tmp_lcp.data(), depth);
    assert(stringtools::verify_lcp(ss, tmp_lcp.data(), 42));
}

// *** Externally Called Functions ***

static inline
void lcp_insertion_sort(string* str, size_t n, size_t depth)
{
    uintptr_t tmp_lcp[n];
    return lcp_insertion_sort(
        parallel_string_sorting::UCharStringSet(str, str + n),
        tmp_lcp, depth);
}

template <typename StringSet>
static inline
void lcp_insertion_sort(StringShadowPtr<StringSet>& str, size_t depth)
{
    return lcp_insertion_sort_nolcp(str.active(), depth);
}

template <typename StringSet>
static inline
void lcp_insertion_sort(StringShadowLcpPtr<StringSet>& str, size_t depth)
{
    return lcp_insertion_sort(str.active(), str.lcparray(), depth);
}

// *** Registered Test Functions ***

static inline
void test_lcp_insertion_sort(string* strings, size_t n)
{
    lcp_insertion_sort(strings, n, 0);
}

PSS_CONTESTANT(test_lcp_insertion_sort,
               "bingmann/lcp_insertion_sort",
               "LCP-aware insertion sort")

static inline
void test_lcp_insertion_sort_nolcp(string* strings, size_t n)
{
    string* shadow = new string[n]; // allocate shadow pointer array
    StringShadowPtr<parallel_string_sorting::UCharStringSet> strptr(
        parallel_string_sorting::UCharStringSet(strings, strings + n),
        parallel_string_sorting::UCharStringSet(shadow, shadow + n));

    lcp_insertion_sort(strptr, 0);

    delete[] shadow;
}

PSS_CONTESTANT(test_lcp_insertion_sort_nolcp,
               "bingmann/lcp_insertion_sort_nolcp",
               "LCP-aware insertion sort (without LCP output)")

static inline
void test_lcp_insertion_sort_pseudocode(string* strings, size_t n)
{
    string* shadow = new string[n + 1]; // allocate shadow pointer array
    StringShadowPtr<parallel_string_sorting::UCharStringSet> strptr(
        parallel_string_sorting::UCharStringSet(strings, strings + n),
        parallel_string_sorting::UCharStringSet(shadow, shadow + n));

    strptr.lcp(0) = 42;                 // must keep lcp[0] unchanged

    lcp_insertion_sort_pseudocode(strptr.active(), strptr.lcparray(), 0);

    stringtools::verify_lcp(strptr.active(), strptr.lcparray(), 42);

    delete[] shadow;
}

PSS_CONTESTANT(test_lcp_insertion_sort_pseudocode,
               "bingmann/lcp_insertion_sort_pseudocode",
               "LCP-aware insertion sort close to pseudo-code, with checking")

////////////////////////////////////////////////////////////////////////////////

//! LCP insertion sort
template <typename StringSet>
static inline
void lcp_insertion_sort_cache(
    const StringSet& str, uintptr_t* lcp,
    typename StringSet::Char* cache, size_t depth)
{
    typedef typename StringSet::String String;
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::CharIterator CharIterator;
    typedef typename StringSet::Char Char;

    Iterator begin = str.begin();
    size_t n = str.size();

    if (n <= 1) return;

    for (size_t j = 0; j < n - 1; ++j)
    {
        // insert strings[j] into sorted strings[0..j-1]

        String new_str = std::move(str[begin + j]);
        size_t new_lcp = depth; // start with LCP depth
        // Char new_ch = str.get_char(new_str, new_lcp);

        size_t i = j;
        while (i > 0)
        {
            size_t prev_lcp = new_lcp;

            String cur_str = std::move(str[begin + i - 1]);
            size_t cur_lcp = lcp[i];

            if (cur_lcp < new_lcp)
            {
                // CASE 1: lcp goes down -> insert string

                // move comparison string back
                str[begin + i - 1] = std::move(cur_str);
                break;
            }
            else if (cur_lcp == new_lcp)
            {
                // CASE 2: compare more characters

                CharIterator c1 = str.get_chars(new_str, new_lcp);
                CharIterator c2 = str.get_chars(cur_str, new_lcp);

                while (str.is_equal(new_str, c1, cur_str, c2))
                {
                    ++new_lcp, ++c1, ++c2;
                }

                // if (new_str >= curr_str) -> insert string
                if (!str.is_less(new_str, c1, cur_str, c2))
                {
                    // update lcp of prev (smaller string) with inserted string
                    lcp[i] = new_lcp;
                    cache[i] = str.is_end(new_str, c1) ? 0 : *c1;

                    // lcp of inserted string with next string
                    new_lcp = prev_lcp;

                    // move comparison string back
                    str[begin + i - 1] = std::move(cur_str);
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            str[begin + i] = std::move(cur_str);
            cache[i] = cache[i - 1];

            lcp[i + 1] = cur_lcp;

            --i;
        }

        str[begin + i] = std::move(new_str);

        lcp[i + 1] = new_lcp;
        cache[i + 1] = str.get_char(str[begin + i + 1], new_lcp);
    }

    // last loop specialized with checks for out-of-bound access to lcp.
    {
        size_t j = n - 1;

        // insert strings[j] into sorted strings[0..j-1]

        String new_str = std::move(str[begin + j]);
        size_t new_lcp = depth; // start with LCP depth
        // Char new_ch = str.get_char(new_str, new_lcp);

        size_t i = j;
        while (i > 0)
        {
            size_t prev_lcp = new_lcp;

            String cur_str = std::move(str[begin + i - 1]);
            size_t cur_lcp = lcp[i];

            if (cur_lcp < new_lcp)
            {
                // CASE 1: lcp goes down -> insert string

                // move comparison string back
                str[begin + i - 1] = std::move(cur_str);
                break;
            }
            else if (cur_lcp == new_lcp)
            {
                // CASE 2: compare more characters

                CharIterator c1 = str.get_chars(new_str, new_lcp);
                CharIterator c2 = str.get_chars(cur_str, new_lcp);

                while (str.is_equal(new_str, c1, cur_str, c2))
                {
                    ++new_lcp, ++c1, ++c2;
                }

                // if (new_str >= curr_str) -> insert string
                if (!str.is_less(new_str, c1, cur_str, c2))
                {
                    // update lcp of prev (smaller string) with inserted string
                    lcp[i] = new_lcp;
                    cache[i] = str.is_end(new_str, c1) ? 0 : *c1;

                    // lcp of inserted string with next string
                    new_lcp = prev_lcp;

                    // move comparison string back
                    str[begin + i - 1] = std::move(cur_str);
                    break;
                }
            }
            // else (cur_lcp > new_lcp), CASE 3: nothing to do

            str[begin + i] = std::move(cur_str);
            cache[i] = cache[i - 1];

            if (i + 1 < n) // check out-of-bounds copy
                lcp[i + 1] = cur_lcp;

            --i;
        }

        str[begin + i] = std::move(new_str);

        if (i + 1 < n) // check out-of-bounds save
        {
            lcp[i + 1] = new_lcp;
            cache[i + 1] = str.get_char(str[begin + i + 1], new_lcp);
        }
    }
}

template <typename StringSet>
static inline
void lcp_insertion_sort_cache_verify(const StringSet& ss, size_t depth)
{
    // allocate LCP array
    std::vector<lcp_t> lcps(ss.size());
    // allocate distinguishing char cache
    std::vector<typename StringSet::Char> cache(ss.size());

    lcps[0] = -1;                 // must keep lcp[0] unchanged
    std::fill(lcps.begin() + 1, lcps.end(), -1);

    lcp_insertion_sort_cache(ss, lcps.data(), cache.data(), depth);

    assert(stringtools::verify_lcp_cache(ss, lcps.data(), cache.data(), -1));
}

static inline
void test_lcp_insertion_sort_cache(string* strings, size_t n)
{
    return lcp_insertion_sort_cache_verify(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(test_lcp_insertion_sort_cache,
               "bingmann/lcp_insertion_sort_cache",
               "LCP-aware insertion sort (with distinguishing character cache)")

} // namespace bingmann_lcp_inssort

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_LCP_INSSORT_HEADER

/******************************************************************************/
