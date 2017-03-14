/*******************************************************************************
 * src/sequential/bingmann-mkqs.hpp
 *
 * Generic Multikey Quicksort
 *
 *******************************************************************************
 * Copyright (C) 2017 Timo Bingmann <tb@panthema.net>
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

#ifndef PSS_SRC_SEQUENTIAL_BINGMANN_MKQS_HEADER
#define PSS_SRC_SEQUENTIAL_BINGMANN_MKQS_HEADER

/*
  Based on Multikey quicksort, a quick sort algorithm for arrays of character
  strings by Bentley and Sedgewick.

  J. Bentley and R. Sedgewick. Fast algorithms for sorting and
  searching strings. In Proceedings of 8th Annual ACM-SIAM Symposium
  on Discrete Algorithms, 1997.

  http://www.CS.Princeton.EDU/~rs/strings/index.html
*/

#include "inssort.hpp"
#include <iostream>

namespace bingmann {

using namespace stringtools;

template <typename StringSet>
static inline void vecswap2(
    typename StringSet::Iterator a, typename StringSet::Iterator b, size_t n)
{
    while (n-- > 0)
        std::swap(*a++, *b++);
}

template <typename StringSet>
static inline typename StringSet::Iterator med3func(
    const StringSet& ss,
    typename StringSet::Iterator a, typename StringSet::Iterator b,
    typename StringSet::Iterator c, size_t depth)
{
    typename StringSet::Char va, vb, vc;
    if ((va = ss.get_char(*a, depth)) == (vb = ss.get_char(*b, depth)))
        return a;
    if ((vc = ss.get_char(*c, depth)) == va || vc == vb)
        return c;
    return va < vb
           ? (vb < vc ? b : (va < vc ? c : a))
           : (vb > vc ? b : (va < vc ? a : c));
}

template <typename StringSet>
static inline void mkqs_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::Char Char;

    const Iterator a = ss.begin();
    size_t n = ss.size();

    int r;
    Iterator pa, pb, pc, pd;

    if (ss.size() < 32) {
        return inssort::inssort_generic(ss, depth);
    }
    Iterator pl = a;
    Iterator pm = a + (n / 2);
    Iterator pn = a + (n - 1);
    if (n > 30) {
        // on big arrays: pseudomedian of 9
        size_t d = (n / 8);
        pl = med3func(ss, pl, pl + d, pl + 2 * d, depth);
        pm = med3func(ss, pm - d, pm, pm + d, depth);
        pn = med3func(ss, pn - 2 * d, pn - d, pn, depth);
    }
    pm = med3func(ss, pl, pm, pn, depth);
    std::swap(*a, *pm);
    Char pivot = ss.get_char(*a, depth);
    pa = pb = a + 1;
    pc = pd = a + n - 1;
    for ( ; ; ) {
        while (pb <= pc && (r = ss.get_char(*pb, depth) - pivot) <= 0) {
            if (r == 0) std::swap(*pa++, *pb);
            pb++;
        }
        while (pb <= pc && (r = ss.get_char(*pc, depth) - pivot) >= 0) {
            if (r == 0) std::swap(*pc, *pd--);
            pc--;
        }
        if (pb > pc) break;
        std::swap(*pb++, *pc--);
    }
    pn = a + n;

    r = std::min(pa - a, pb - pa);
    vecswap2<StringSet>(a, pb - r, r);
    r = std::min(pd - pc, pn - pd - 1);
    vecswap2<StringSet>(pb, pn - r, r);

    if ((r = pb - pa) > 1)
        mkqs_generic(ss.sub(a, a + r), depth);
    if (ss.get_char(*(a + r), depth) != 0)
        mkqs_generic(ss.sub(a + r, a + r + (pa - a) + (pn - pd - 1)), depth + 1);
    if ((r = pd - pc) > 1)
        mkqs_generic(ss.sub(a + n - r, a + n), depth);
}

static inline void mkqs_generic(unsigned char** strings, size_t n)
{
    return mkqs_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(mkqs_generic,
               "bingmann/mkqs_generic",
               "Generic Multikey-Quicksort")

} // namespace bingmann

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_MKQS_HEADER

/******************************************************************************/
