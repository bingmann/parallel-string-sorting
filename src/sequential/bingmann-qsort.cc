/*******************************************************************************
 * src/sequential/bingmann-qsort.cc
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

#include "../tools/debug.h"
#include "../tools/contest.h"
#include "../tools/stringtools.h"

namespace bingmann_qsort {

using namespace stringtools;

static const bool debug = false;

typedef unsigned char* string;

////////////////////////////////////////////////////////////////////////////////

int compare_strcmp(const void* _a, const void* _b)
{
    const char* a = *(char* const*)_a;
    const char* b = *(char* const*)_b;

    return strcmp(a, b);
}

void qsort_strcmp(string* strings, size_t n)
{
    qsort((void*)strings, n, sizeof(string), compare_strcmp);
}

PSS_CONTESTANT(qsort_strcmp,
               "bingmann/qsort_strcmp",
               "Run stdlib qsort with strcmp comparsion")

////////////////////////////////////////////////////////////////////////////////

static inline
int qcompare_byte(const void* _a, const void* _b)
{
    const char* a = *(char* const*)_a;
    const char* b = *(char* const*)_b;

    for (size_t i = 0; ; i++)
    {
        if (a[i] != b[i]) {
            return (a[i] > b[i]) ? +1 : -1;
        }

        // check byte for zero -> both strings end
        if (a[i] == 0)
            return 0;
    }
}

template <typename key_type>
static inline
int qcompare_uint(const void* _a, const void* _b)
{
    const string a = *(string const*)_a;
    const string b = *(string const*)_b;

    for (size_t depth = 0; ; depth += sizeof(key_type))
    {
        key_type av = get_char<key_type>(a, depth);
        key_type bv = get_char<key_type>(b, depth);

        if (av != bv) {
            return (av > bv) ? +1 : -1;
        }

        // check highest byte for zero -> both strings end
        const key_type mask = key_type(0xFF) << 8 * (sizeof(key_type) - 1);
        if ((av & mask) == 0)
            return 0;
    }
}

void qsort1(string* strings, size_t n)
{
    qsort((void*)strings, n, sizeof(string), qcompare_byte);
}

void qsort4(string* strings, size_t n)
{
    qsort((void*)strings, n, sizeof(string), qcompare_uint<uint32_t>);
}

void qsort8(string* strings, size_t n)
{
    qsort((void*)strings, n, sizeof(string), qcompare_uint<uint64_t>);
}

PSS_CONTESTANT(qsort1,
               "bingmann/qsort1",
               "Run stdlib qsort with string comparsions (bytewise)")

PSS_CONTESTANT(qsort4,
               "bingmann/qsort4",
               "Run stdlib qsort with string comparsions (4 bytewise)")

PSS_CONTESTANT(qsort8,
               "bingmann/qsort8",
               "Run stdlib qsort with string comparsions (8 bytewise)")

////////////////////////////////////////////////////////////////////////////////

static inline
bool stdcompare_byte(const string a, const string b)
{
    for (size_t i = 0; ; i++)
    {
        if (a[i] != b[i]) {
            return (a[i] > b[i]) ? false : true;
        }

        // check byte for zero -> both strings end
        if (a[i] == 0)
            return false;
    }
}

template <typename key_type>
static inline
bool stdcompare_uint(const string a, const string b)
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
