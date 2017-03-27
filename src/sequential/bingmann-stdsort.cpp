/*******************************************************************************
 * src/sequential/bingmann-stdsort.cpp
 *
 * String sorting using qsort() and std::sort with 1-, 4- and 8-bytewise
 * comparisons
 *
 *******************************************************************************
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
 ******************************************************************************/

#include <cstdlib>
#include <algorithm>

#include "../tools/debug.hpp"
#include "../tools/contest.hpp"
#include "../tools/stringtools.hpp"

namespace bingmann_qsort {

using namespace stringtools;

static const bool debug = false;

typedef unsigned char* string;

////////////////////////////////////////////////////////////////////////////////

static inline
bool stdcompare_strcmp(const string& a, const string& b)
{
    return strcmp((const char*)a, (const char*)b) < 0;
}

static inline
bool stdcompare_byte(const string& _a, const string& _b)
{
    string a = _a, b = _b;

    while (*a == *b && *a != 0)
        ++a, ++b;

    // check byte for zero -> both strings end, random tie break
    return (*a < *b);
}

template <typename key_type>
static inline
bool stdcompare_uint(const string& a, const string& b)
{
    for (size_t depth = 0; ; depth += sizeof(key_type))
    {
        key_type av = get_char<key_type>(a, depth);
        key_type bv = get_char<key_type>(b, depth);

        if (av != bv) {
            return (av > bv) ? false : true;
        }

        // check highest byte for zero -> both strings end
        const key_type mask = key_type(0xFF) << 8 * (sizeof(key_type) - 1);
        if ((av & mask) == 0)
            return false;
    }
}

void stdsort(string* strings, size_t n)
{
    std::sort(strings, strings + n, stdcompare_strcmp);
}

void stdsort1(string* strings, size_t n)
{
    std::sort(strings, strings + n, stdcompare_byte);
}

void stdsort4(string* strings, size_t n)
{
    std::sort(strings, strings + n, stdcompare_uint<uint32_t>);
}

void stdsort8(string* strings, size_t n)
{
    std::sort(strings, strings + n, stdcompare_uint<uint64_t>);
}

PSS_CONTESTANT(stdsort,
               "bingmann/stdsort",
               "Run std::sort with strcmp() string comparsions")

PSS_CONTESTANT(stdsort1,
               "bingmann/stdsort1",
               "Run std::sort with string comparsions (bytewise)")

PSS_CONTESTANT(stdsort4,
               "bingmann/stdsort4",
               "Run std::sort with string comparsions (4 bytewise)")

PSS_CONTESTANT(stdsort8,
               "bingmann/stdsort8",
               "Run std::sort with string comparsions (8 bytewise)")

} // namespace bingmann_qsort

/******************************************************************************/
