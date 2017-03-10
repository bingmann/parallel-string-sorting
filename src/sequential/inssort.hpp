/******************************************************************************
 * src/sequential/inssort.h
 *
 * Base insertion sort from Rantala string-sorting/external/
 *
 ******************************************************************************
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

#ifndef INSSORT_H_
#define INSSORT_H_

#include "../tools/contest.hpp"
#include "../tools/stringset.hpp"
#include <tlx/define/likely.hpp>

namespace inssort {

typedef unsigned char* string;

/******************************************************************************/

static inline void
inssort(string* str, int n, int d)
{
    string *pj, s, t;

    for (string* pi = str + 1; TLX_UNLIKELY(--n > 0); pi++) {
        string tmp = *pi;

        for (pj = pi; TLX_LIKELY(pj > str); pj--) {
            for (s = *(pj-1)+d, t = tmp+d; TLX_LIKELY(*s == *t && *s != 0); ++s, ++t)
                ;
            if (TLX_UNLIKELY(*s <= *t))
                break;
            *pj = *(pj-1);
        }
        *pj = tmp;
    }
}

static inline
void insertion_sort(string* a, size_t n)
{ inssort(a, n, 0); }

PSS_CONTESTANT(insertion_sort, "insertion_sort", "String Insertion-Sort")

/******************************************************************************/

//! Generic insertion sort for objectified string sets
template <typename StringSet>
static inline void inssort_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::String String;
    typedef typename StringSet::CharIterator CharIterator;

    // this stores the begin iterator and size n, making the loops faster
    const Iterator begin = ss.begin();
    Iterator j;
    size_t n = ss.size();

    for (Iterator i = begin + 1; TLX_UNLIKELY(--n != 0); ++i)
    {
        String tmp = std::move(ss[i]);
        j = i;

        while (TLX_LIKELY(j != begin))
        {
            CharIterator s = ss.get_chars(ss[j - 1], depth);
            CharIterator t = ss.get_chars(tmp, depth);

            while (TLX_LIKELY(ss.is_equal(ss[j - 1], s, tmp, t)))
                ++s, ++t;

            if (TLX_UNLIKELY(ss.is_leq(ss[j - 1], s, tmp, t))) {
                break;
            }

            ss[j] = std::move(ss[j - 1]);
            --j;
        }

        ss[j] = std::move(tmp);
    }
}

static inline
void insertion_sort_generic(string* a, size_t n)
{
    inssort_generic(parallel_string_sorting::UCharStringSet(a, a + n), 0);
}

PSS_CONTESTANT(insertion_sort_generic, "insertion_sort_gen",
               "String Insertion-Sort generic")

/******************************************************************************/

} // namespace  inssort

#endif // INSSORT_H_
