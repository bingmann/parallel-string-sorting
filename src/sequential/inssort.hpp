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

namespace inssort {

typedef unsigned char* string;

/******************************************************************************/

static inline void
inssort(string* str, int n, int d)
{
    string *pj, s, t;

    for (string* pi = str + 1; --n > 0; pi++) {
        string tmp = *pi;

        for (pj = pi; pj > str; pj--) {
            for (s = *(pj-1)+d, t = tmp+d; *s == *t && *s != 0; ++s, ++t)
                ;
            if (*s <= *t)
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

static inline void
inssort_range(string* str_begin, string* str_end, size_t depth)
{
    for (string* i = str_begin + 1; i != str_end; ++i) {
        string* j = i;
        string tmp = *i;
        while (j > str_begin) {
            string s = *(j - 1) + depth;
            string t = tmp + depth;
            while (*s == *t && *s != 0) ++s, ++t;
            if (*s <= *t) break;
            *j = *(j - 1);
            --j;
        }
        *j = tmp;
    }
}

/******************************************************************************/

//! Generic insertion sort for objectified string sets
template <typename StringSet>
static inline void inssort_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::String String;
    typedef typename StringSet::CharIterator CharIterator;

    for (Iterator pi = ss.begin() + 1; pi != ss.end(); ++pi)
    {
        String tmp = ss[pi];
        Iterator pj = pi;

        while (pj != ss.begin())
        {
            --pj;

            CharIterator s = ss.get_chars(ss[pj], depth),
                t = ss.get_chars(tmp, depth);

            while (*s == *t && *s != 0)
                ++s, ++t;

            if (*s <= *t) {
                ++pj;
                break;
            }

            ss[pj + 1] = ss[pj];
        }

        ss[pj] = tmp;
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
