/*******************************************************************************
 * src/sequential/bingmann-radix_sort.hpp
 *
 * Experiments with sequential radix sort implementations.
 * Based on rantala/msd_c?.h
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

#ifndef PSS_SRC_SEQUENTIAL_BINGMANN_RADIX_SORT_HEADER
#define PSS_SRC_SEQUENTIAL_BINGMANN_RADIX_SORT_HEADER

#include "../tools/globals.hpp"
#include "../tools/stringtools.hpp"
#include "../tools/stringset.hpp"
#include "inssort.hpp"

#include <stack>

namespace bingmann_radix_sort {

static const size_t g_inssort_threshold = 32;

typedef unsigned char* string;

/******************************************************************************/

static inline void
msd_CE0(string* strings, string* sorted, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // count character occurances
    size_t bkt_size[256] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bkt_size[strings[i][depth]];

    // prefix sum
    size_t bkt_index[256];
    bkt_index[0] = 0;
    for (size_t i = 1; i < 256; ++i)
        bkt_index[i] = bkt_index[i - 1] + bkt_size[i - 1];

    // distribute
    for (size_t i = 0; i < n; ++i)
        sorted[bkt_index[strings[i][depth]]++] = strings[i];

    std::copy(sorted, sorted + n, strings);

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CE0(strings + bsum, sorted, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }
}

void bingmann_msd_CE0(string* strings, size_t n)
{
    string* sorted = new string[n];
    msd_CE0(strings, sorted, n, 0);
    delete [] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE0, "bingmann/msd_CE0",
               "bingmann/msd_CE0 (rantala CE0 baseline)")

/******************************************************************************/

template <typename StringSet>
static inline void
msd_CE0_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::String String;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    // count character occurances
    size_t bkt_size[256] = { 0 };
    for (Iterator i = ss.begin(); i != ss.end(); ++i)
        ++bkt_size[ss.get_uint8(ss[i], depth)];

    String* sorted = new String[ss.size()];

    // prefix sum
    size_t bkt_index[256];
    bkt_index[0] = 0;
    for (size_t i = 1; i < 256; ++i)
        bkt_index[i] = bkt_index[i - 1] + bkt_size[i - 1];

    // distribute
    for (Iterator i = ss.begin(); i != ss.end(); ++i)
        sorted[bkt_index[ss.get_uint8(ss[i], depth)]++] = std::move(ss[i]);

    std::move(sorted, sorted + ss.size(), ss.begin());

    delete[] sorted;

    // recursion
    Iterator bsum = ss.begin() + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        Iterator bsum_begin = bsum;
        bsum += bkt_size[i];
        if (bkt_size[i] > 1)
            msd_CE0_generic(ss.sub(bsum_begin, bsum), depth + 1);
    }
}

void bingmann_msd_CE0_generic(string* strings, size_t n)
{
    return msd_CE0_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CE0_generic, "bingmann/msd_CE0_gen",
               "bingmann/msd_CE0 generic (rantala CE baseline)")

/******************************************************************************/

static inline void
msd_CE1(string* strings, string* sorted, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // count character occurances
    size_t bkt[256 + 1] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bkt[strings[i][depth]];

    // prefix sum
    for (size_t i = 1; i <= 256; ++i)
        bkt[i] += bkt[i - 1];

    // distribute
    for (size_t i = 0; i < n; ++i)
        sorted[--bkt[strings[i][depth]]] = strings[i];

    std::copy(sorted, sorted + n, strings);

    // recursion
    for (size_t i = 1; i < 256; ++i) {
        if (bkt[i]+1 >= bkt[i+1]) continue;
        msd_CE1(strings + bkt[i], sorted, bkt[i + 1] - bkt[i], depth + 1);
    }
}

void bingmann_msd_CE1(string* strings, size_t n)
{
    string* sorted = new string[n];
    msd_CE1(strings, sorted, n, 0);
    delete [] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE1, "bingmann/msd_CE1",
               "bingmann/msd_CE1 (CE0 with reused prefix sum)")

/******************************************************************************/

template <typename StringSet>
static inline void
msd_CE1_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::String String;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    // count character occurances
    size_t bkt[256 + 1] = { 0 };
    for (Iterator i = ss.begin(); i != ss.end(); ++i)
        ++bkt[ss.get_uint8(ss[i], depth)];

    String* sorted = new String[ss.size()];

    // prefix sum
    for (size_t i = 1; i <= 256; ++i)
        bkt[i] += bkt[i - 1];

    // distribute
    for (Iterator i = ss.begin(); i != ss.end(); ++i)
        sorted[--bkt[ss.get_uint8(ss[i], depth)]] = std::move(ss[i]);

    std::move(sorted, sorted + ss.size(), ss.begin());

    delete[] sorted;

    // recursion
    for (size_t i = 1; i < 256; ++i) {
        if (bkt[i]+1 >= bkt[i+1]) continue;
        msd_CE1_generic(ss.sub(ss.begin() + bkt[i], ss.begin() + bkt[i + 1]),
                        depth + 1);
    }
}

void bingmann_msd_CE1_generic(string* strings, size_t n)
{
    return msd_CE1_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CE1_generic, "bingmann/msd_CE1_gen",
               "bingmann/msd_CE1 generic (CE with reused prefix sum)")


/******************************************************************************/

static inline void
msd_CE1_o(string* strings, string* sorted, uint8_t* charcache,
          size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // read characters and count character occurances
    size_t bkt[256 + 1] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bkt[charcache[i] = static_cast<uint8_t>(strings[i][depth])];

    // prefix sum
    for (size_t i = 1; i <= 256; ++i)
        bkt[i] += bkt[i - 1];

    // distribute
    for (size_t i = 0; i < n; ++i)
        sorted[--bkt[charcache[i]]] = strings[i];

    std::copy(sorted, sorted + n, strings);

    // recursion
    for (size_t i = 1; i < 256; ++i) {
        if (bkt[i]+1 >= bkt[i+1]) continue;
        msd_CE1_o(strings + bkt[i], sorted, charcache,
                  bkt[i + 1] - bkt[i], depth + 1);
    }
}

void bingmann_msd_CE1_o(string* strings, size_t n)
{
    string* sorted = new string[n];
    uint8_t* charcache = new uint8_t[n];

    msd_CE1_o(strings, sorted, charcache, n, 0);

    delete [] charcache;
    delete [] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE1_o, "bingmann/msd_CE1_o",
               "bingmann/msd_CE1_o (CE1 with oracle)")

/******************************************************************************/

static inline void
msd_CE1_of(string* strings, string* sorted, uint8_t* charcache,
           size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // read characters and count character occurances
    size_t bkt[256 + 1] = { 0 };
    for (size_t i = 0; i < n; ++i)
        charcache[i] = static_cast<uint8_t>(strings[i][depth]);
    for (size_t i = 0; i < n; ++i)
        ++bkt[charcache[i]];

    // prefix sum
    for (size_t i = 1; i <= 256; ++i)
        bkt[i] += bkt[i - 1];

    // distribute
    for (size_t i = 0; i < n; ++i)
        sorted[--bkt[charcache[i]]] = strings[i];

    std::copy(sorted, sorted + n, strings);

    // recursion
    for (size_t i = 1; i < 256; ++i) {
        if (bkt[i]+1 >= bkt[i+1]) continue;
        msd_CE1_of(strings + bkt[i], sorted, charcache,
                   bkt[i + 1] - bkt[i], depth + 1);
    }
}

void bingmann_msd_CE1_of(string* strings, size_t n)
{
    string* sorted = new string [n];
    uint8_t* charcache = new uint8_t[n];

    msd_CE1_of(strings, sorted, charcache, n, 0);

    delete [] charcache;
    delete [] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE1_of, "bingmann/msd_CE1_of",
               "bingmann/msd_CE1_of (CE1 with oracle, fissioned)")

/******************************************************************************/
// CI Variants
/******************************************************************************/

/******************************************************************************/
// msd_CI0

// Note: CI in-place variants cannot be done with just one prefix-sum bucket
// array, because during in-place permutation the beginning _and_ end boundaries
// of each bucket must be kept.

static inline void
msd_CI0(string* strings, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // count character occurances
    size_t bkt_size[256] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bkt_size[strings[i][depth]];

    // inclusive prefix sum
    size_t bkt_index[256];
    bkt_index[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (unsigned i = 1; i < 256; ++i) {
        bkt_index[i] = bkt_index[i - 1] + bkt_size[i];
        if (bkt_size[i]) last_bkt_size = bkt_size[i];
    }

    // premute in-place
    for (size_t i = 0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        while ((j = --bkt_index[perm[depth]]) > i)
        {
            std::swap(perm, strings[j]);
        }
        strings[i] = perm;
        i += bkt_size[perm[depth]];
    }

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CI0(strings + bsum, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }
}

void bingmann_msd_CI0(string* strings, size_t n)
{
    msd_CI0(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI0, "bingmann/msd_CI0",
               "bingmann/msd_CI0 (CE0 in-place baseline)")

/******************************************************************************/
// msd_CI0_o

static inline size_t *
msd_CI0_o_make_bkt_size(
    string * strings, uint8_t* charcache, size_t n, size_t depth)
{
    // cache and count character occurrences
    size_t* bkt_size = new size_t[256];
    memset(bkt_size, 0, 256 * sizeof(size_t));
    for (size_t i = 0; i < n; ++i)
        ++bkt_size[charcache[i] = strings[i][depth]];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (unsigned i = 1; i < 256; ++i) {
        bkt[i] = bkt[i - 1] + bkt_size[i];
        if (bkt_size[i]) last_bkt_size = bkt_size[i];
    }

    // premute in-place
    for (size_t i = 0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        uint8_t permch = charcache[i];
        while ((j = --bkt[permch]) > i)
        {
            std::swap(perm, strings[j]);
            std::swap(permch, charcache[j]);
        }
        strings[i] = perm;
        i += bkt_size[permch];
    }

    return bkt_size;
}

void
msd_CI0_o(string* strings, uint8_t* charcache, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    size_t* bkt_size = msd_CI0_o_make_bkt_size(strings, charcache, n, depth);

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CI0_o(strings + bsum, charcache, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

void bingmann_msd_CI0_o(string* strings, size_t n)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI0_o(strings, charcache, n, 0);
    delete [] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI0_o, "bingmann/msd_CI0_o",
               "bingmann/msd_CI0_o (CI4 with charcache)")

/******************************************************************************/
// msd_CI0_of

static inline size_t *
msd_CI0_of_run(
    string * strings, uint8_t* charcache, size_t n, size_t depth)
{
    // cache characters
    for (size_t i = 0; i < n; ++i)
        charcache[i] = strings[i][depth];

    // cache and count character occurrences
    size_t* bkt_size = new size_t[256];
    memset(bkt_size, 0, 256 * sizeof(size_t));
    for (size_t i = 0; i < n; ++i)
        ++bkt_size[charcache[i]];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (unsigned i = 1; i < 256; ++i) {
        bkt[i] = bkt[i - 1] + bkt_size[i];
        if (bkt_size[i]) last_bkt_size = bkt_size[i];
    }

    // premute in-place
    for (size_t i = 0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        uint8_t permch = charcache[i];
        while ((j = --bkt[permch]) > i)
        {
            std::swap(perm, strings[j]);
            std::swap(permch, charcache[j]);
        }
        strings[i] = perm;
        i += bkt_size[permch];
    }

    return bkt_size;
}

void
msd_CI0_of(string* strings, uint8_t* charcache, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    size_t* bkt_size = msd_CI0_of_run(strings, charcache, n, depth);

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CI0_of(strings + bsum, charcache, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

void bingmann_msd_CI0_of(string* strings, size_t n)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI0_of(strings, charcache, n, 0);
    delete [] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI0_of, "bingmann/msd_CI0_of",
               "bingmann/msd_CI0_of (CI4 with charcache)")

void msd_CI5(string* strings, size_t n, size_t depth)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI0_of(strings, charcache, n, depth);
    delete [] charcache;
}

/******************************************************************************/

template <typename StringSet>
static inline size_t*
msd_CI5_bktsize_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::String String;
    typedef typename StringSet::Char Char;

    // cache characters
    Char* charcache = new Char[ss.size()];
    size_t j = 0;
    for (Iterator i = ss.begin(); i != ss.end(); ++i)
        charcache[j++] = *ss.get_chars(ss[i], depth);

    // count character occurances
    size_t* bkt_size = new size_t[256];
    memset(bkt_size, 0, 256 * sizeof(size_t));
    for (size_t i = 0; i < ss.size(); ++i)
        ++bkt_size[static_cast<size_t>(charcache[i])];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (unsigned i = 1; i < 256; ++i) {
        bkt[i] = bkt[i - 1] + bkt_size[i];
        if (bkt_size[i]) last_bkt_size = bkt_size[i];
    }

    // premute in-place
    for (size_t i = 0, j; i < ss.size() - last_bkt_size; )
    {
        String perm = std::move(ss[ss.begin() + i]);
        Char permch = charcache[i];
        while ((j = --bkt[static_cast<size_t>(permch)]) > i)
        {
            std::swap(perm, ss[ss.begin() + j]);
            std::swap(permch, charcache[j]);
        }
        ss[ss.begin() + i] = std::move(perm);
        i += bkt_size[static_cast<size_t>(permch)];
    }

    delete[] charcache;

    return bkt_size;
}

template <typename StringSet>
static inline void
msd_CI5_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    size_t* bkt_size = msd_CI5_bktsize_generic(ss, depth);

    // recursion
    Iterator bsum = ss.begin() + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] == 0) continue;
        Iterator bend = bsum + bkt_size[i];
        msd_CI5_generic(ss.sub(bsum, bend), depth + 1);
        bsum = bend;
    }

    delete[] bkt_size;
}

void bingmann_msd_CI5_generic(string* strings, size_t n)
{
    msd_CI5_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CI5_generic, "bingmann/msd_CI5_gen",
               "bingmann/msd_CI0_o generic (CI4 with charcache)")

/******************************************************************************/

static inline size_t *
msd_CI0_of_16bit_run(string * strings, size_t n, size_t depth)
{
    static const size_t RADIX = 0x10000;
    using namespace stringtools;

    // cache characters
    uint16_t* charcache = new uint16_t[n];
    for (size_t i = 0; i < n; ++i)
        charcache[i] = get_char<uint16_t>(strings[i], depth);

    // count character occurances
    size_t* bkt_size = new size_t[RADIX];
    memset(bkt_size, 0, RADIX * sizeof(size_t));
    for (size_t i = 0; i < n; ++i)
        ++bkt_size[charcache[i]];

    // inclusive prefix sum
    size_t bkt[RADIX];
    bkt[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (unsigned i = 1; i < RADIX; ++i) {
        bkt[i] = bkt[i - 1] + bkt_size[i];
        if (bkt_size[i]) last_bkt_size = bkt_size[i];
    }

    // premute in-place
    for (size_t i = 0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        uint16_t permch = charcache[i];
        while ((j = --bkt[permch]) > i)
        {
            std::swap(perm, strings[j]);
            std::swap(permch, charcache[j]);
        }
        strings[i] = perm;
        i += bkt_size[permch];
    }

    delete[] charcache;

    return bkt_size;
}

static void
msd_CI0_of_16bit(string* strings, size_t n, size_t depth)
{
    if (n < 0x10000)
        return msd_CI5(strings, n, depth);

    size_t* bkt_size = msd_CI0_of_16bit_run(strings, n, depth);

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 0x10000; ++i) {
        if (bkt_size[i] == 0) continue;
        if (i & 0xFF) // not zero-terminated
            msd_CI0_of_16bit(strings + bsum, bkt_size[i], depth + 2);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

void bingmann_msd_CI0_of_16bit(string* strings, size_t n)
{
    msd_CI0_of_16bit(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI0_of_16bit, "bingmann/msd_CI0_of_16bit",
               "bingmann/msd_CI0_of_16bit (CI0_of with 16-bit radix)")

/******************************************************************************/
// Iterative Stack-Based Variants
/******************************************************************************/

struct RadixStep_CE_nr
{
    string * str;
    size_t bkt[256 + 1];
    size_t idx;

    RadixStep_CE_nr(string* strings, size_t n, size_t depth)
    {
        memset(bkt, 0, sizeof(bkt));

        for (size_t i = 0; i < n; ++i)
            ++bkt[strings[i][depth]];

        string* sorted = (string*)malloc(n * sizeof(string));

        for (unsigned i = 1; i <= 256; ++i)
            bkt[i] += bkt[i - 1];

        for (size_t i = 0; i < n; ++i)
            sorted[--bkt[strings[i][depth]]] = strings[i];

        memcpy(strings, sorted, n * sizeof(string));
        free(sorted);

        str = strings;
        idx = 0;        // will increment to 1 on first process
    }
};

static inline void
bingmann_msd_CE_nr(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    typedef RadixStep_CE_nr RadixStep;

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.push(RadixStep(strings, n, 0));

    while (radixstack.size())
    {
        RadixStep& rs = radixstack.top();

        while (++rs.idx < 256)
        {
            // process the bucket rs.idx

            if (rs.bkt[rs.idx] == rs.bkt[rs.idx + 1])
                ;
            else if (rs.bkt[rs.idx + 1] < rs.bkt[rs.idx] + g_inssort_threshold)
            {
                inssort::inssort(rs.str + rs.bkt[rs.idx],
                                 rs.bkt[rs.idx + 1] - rs.bkt[rs.idx],
                                 radixstack.size());
            }
            else
            {
                radixstack.push(RadixStep(rs.str + rs.bkt[rs.idx],
                                          rs.bkt[rs.idx + 1] - rs.bkt[rs.idx],
                                          radixstack.size()));
                break;
            }
        }

        if (radixstack.top().idx == 256) {      // rs maybe have been invalidated
            radixstack.pop();
        }
    }
}

PSS_CONTESTANT(bingmann_msd_CE_nr, "bingmann/msd_CE_nr",
               "bingmann/msd_CE_nr (CE non-recursive)")

/******************************************************************************/

struct RadixStep_CE_nr2
{
    string * str;
    size_t bkt_size[256];
    size_t idx;

    RadixStep_CE_nr2(string* strings, size_t n, size_t depth)
    {
        // count character occurances
        memset(bkt_size, 0, sizeof(bkt_size));
        for (size_t i = 0; i < n; ++i)
            ++bkt_size[strings[i][depth]];

        // prefix sum
        size_t bkt[256];
        bkt[0] = 0;
        for (size_t i = 1; i < 256; ++i)
            bkt[i] = bkt[i - 1] + bkt_size[i - 1];

        // distribute out-of-place
        string* sorted = (string*)malloc(n * sizeof(string));

        for (size_t i = 0; i < n; ++i)
            sorted[bkt[strings[i][depth]]++] = strings[i];

        memcpy(strings, sorted, n * sizeof(string));
        free(sorted);

        str = strings + bkt_size[0];
        idx = 0; // will increment to 1 on first process
    }
};

static inline void
bingmann_msd_CE_nr2(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    typedef RadixStep_CE_nr2 RadixStep;

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.push(RadixStep(strings, n, 0));

    while (radixstack.size())
    {
        while (radixstack.top().idx < 255)
        {
            RadixStep& rs = radixstack.top();
            ++rs.idx; // process the bucket rs.idx

            if (rs.bkt_size[rs.idx] == 0)
                ;
            else if (rs.bkt_size[rs.idx] < g_inssort_threshold)
            {
                inssort::inssort(rs.str, rs.bkt_size[rs.idx], radixstack.size());
                rs.str += rs.bkt_size[rs.idx];
            }
            else
            {
                rs.str += rs.bkt_size[rs.idx];
                radixstack.push(RadixStep(rs.str - rs.bkt_size[rs.idx], rs.bkt_size[rs.idx], radixstack.size()));
                // cannot add here, because rs may have invalidated
            }
        }
        radixstack.pop();
    }
}

PSS_CONTESTANT(bingmann_msd_CE_nr2, "bingmann/msd_CE_nr2",
               "bingmann/msd_CE_nr2 (CE non-recursive)")

/******************************************************************************/

struct RadixStep_CI_nr
{
    string * str;
    size_t bkt[256 + 1];
    size_t idx;

    RadixStep_CI_nr(string* strings, size_t n, size_t depth)
    {
        // count character occurances
        size_t bkt_size[256] = { 0 };
        for (size_t i = 0; i < n; ++i)
            ++bkt_size[strings[i][depth]];

        // inclusive prefix sum
        bkt[0] = bkt_size[0];
        size_t last_bkt_size = bkt_size[0], last_bkt_num = 0;
        for (unsigned i = 1; i < 256; ++i) {
            bkt[i] = bkt[i - 1] + bkt_size[i];
            if (bkt_size[i]) {
                last_bkt_size = bkt_size[i];
                last_bkt_num = i;
            }
        }
        bkt[256] = bkt[255];
        assert(bkt[256] == n);

        // premute in-place
        for (size_t i = 0, j; i < n - last_bkt_size; )
        {
            string perm = strings[i];
            if (bkt[perm[depth]] == i)
            {
                i += bkt_size[perm[depth]];
                continue;
            }
            while ((j = --bkt[perm[depth]]) > i)
            {
                std::swap(perm, strings[j]);
            }
            assert(j == i);
            strings[i] = perm;
            i += bkt_size[perm[depth]];
        }
        // update bkt boundary of last (unpermuted) bucket
        bkt[last_bkt_num] = n - last_bkt_size;

        str = strings;
        idx = 0;        // will increment to 1 on first process
    }
};

static inline void
bingmann_msd_CI_nr(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    typedef RadixStep_CI_nr RadixStep;

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.push(RadixStep(strings, n, 0));

    while (radixstack.size())
    {
        RadixStep& rs = radixstack.top();

        while (++rs.idx < 256)
        {
            // process the bucket rs.idx

            if (rs.bkt[rs.idx] == rs.bkt[rs.idx + 1])
                ;
            else if (rs.bkt[rs.idx + 1] < rs.bkt[rs.idx] + g_inssort_threshold)
            {
                inssort::inssort(rs.str + rs.bkt[rs.idx],
                                 rs.bkt[rs.idx + 1] - rs.bkt[rs.idx],
                                 radixstack.size());
            }
            else
            {
                radixstack.push(RadixStep(rs.str + rs.bkt[rs.idx],
                                          rs.bkt[rs.idx + 1] - rs.bkt[rs.idx],
                                          radixstack.size()));
                break;
            }
        }

        if (radixstack.top().idx == 256) {      // rs maybe have been invalidated
            radixstack.pop();
        }
    }
}

PSS_CONTESTANT(bingmann_msd_CI_nr, "bingmann/msd_CI_nr",
               "bingmann/msd_CI_nr (CI non-recursive)")

/******************************************************************************/

struct RadixStep_CI_nr2
{
    string * str;
    size_t idx;
    size_t bkt_size[256];

    RadixStep_CI_nr2(string* strings, size_t n, size_t depth)
    {
        // count character occurances
        memset(bkt_size, 0, sizeof(bkt_size));
        for (size_t i = 0; i < n; ++i)
            ++bkt_size[strings[i][depth]];

        // inclusive prefix sum
        size_t bkt[256];
        bkt[0] = bkt_size[0];
        size_t last_bkt_size = bkt_size[0];
        for (unsigned i = 1; i < 256; ++i) {
            bkt[i] = bkt[i - 1] + bkt_size[i];
            if (bkt_size[i]) last_bkt_size = bkt_size[i];
        }

        // premute in-place
        for (size_t i = 0, j; i < n - last_bkt_size; )
        {
            string perm = strings[i];
            while ((j = --bkt[perm[depth]]) > i)
            {
                std::swap(perm, strings[j]);
            }
            strings[i] = perm;
            i += bkt_size[perm[depth]];
        }

        str = strings + bkt_size[0];
        idx = 0; // will increment to 1 on first process, bkt 0 is not sorted further
    }
};

static inline void
bingmann_msd_CI_nr2(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    typedef RadixStep_CI_nr2 RadixStep;

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.push(RadixStep(strings, n, 0));

    while (radixstack.size())
    {
        while (radixstack.top().idx < 255)
        {
            RadixStep& rs = radixstack.top();
            ++rs.idx; // process the bucket rs.idx

            if (rs.bkt_size[rs.idx] == 0)
                ;
            else if (rs.bkt_size[rs.idx] < g_inssort_threshold)
            {
                inssort::inssort(rs.str, rs.bkt_size[rs.idx], radixstack.size());
                rs.str += rs.bkt_size[rs.idx];
            }
            else
            {
                rs.str += rs.bkt_size[rs.idx];
                radixstack.push(RadixStep(rs.str - rs.bkt_size[rs.idx], rs.bkt_size[rs.idx], radixstack.size()));
                // cannot add here, because rs may have invalidated
            }
        }
        radixstack.pop();
    }
}

PSS_CONTESTANT(bingmann_msd_CI_nr2, "bingmann/msd_CI_nr2",
               "bingmann/msd_CI_nr2 (CI non-recursive)")

/******************************************************************************/

struct RadixStep_CI_nr3
{
    string * str;
    size_t idx;
    size_t bkt_size[256];

    RadixStep_CI_nr3(string* strings, size_t n, size_t depth, uint8_t* charcache)
    {
        // cache characters
        for (size_t i = 0; i < n; ++i)
            charcache[i] = strings[i][depth];

        // count character occurances
        memset(bkt_size, 0, sizeof(bkt_size));
        for (size_t i = 0; i < n; ++i)
            ++bkt_size[charcache[i]];

        // inclusive prefix sum
        size_t bkt[256];
        bkt[0] = bkt_size[0];
        size_t last_bkt_size = bkt_size[0];
        for (unsigned i = 1; i < 256; ++i) {
            bkt[i] = bkt[i - 1] + bkt_size[i];
            if (bkt_size[i]) last_bkt_size = bkt_size[i];
        }

        // premute in-place
        for (size_t i = 0, j; i < n - last_bkt_size; )
        {
            string perm = strings[i];
            uint8_t permch = charcache[i];
            while ((j = --bkt[permch]) > i)
            {
                std::swap(perm, strings[j]);
                std::swap(permch, charcache[j]);
            }
            strings[i] = perm;
            i += bkt_size[permch];
        }

        str = strings + bkt_size[0];
        idx = 0; // will increment to 1 on first process, bkt 0 is not sorted further
    }
};

static inline void
bingmann_msd_CI_nr3(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    typedef RadixStep_CI_nr3 RadixStep;

    uint8_t* charcache = new uint8_t[n];

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.push(RadixStep(strings, n, 0, charcache));

    while (radixstack.size())
    {
        while (radixstack.top().idx < 255)
        {
            RadixStep& rs = radixstack.top();
            ++rs.idx; // process the bucket rs.idx

            if (rs.bkt_size[rs.idx] == 0)
                ;
            else if (rs.bkt_size[rs.idx] < g_inssort_threshold)
            {
                inssort::inssort(rs.str, rs.bkt_size[rs.idx], radixstack.size());
                rs.str += rs.bkt_size[rs.idx];
            }
            else
            {
                rs.str += rs.bkt_size[rs.idx];
                radixstack.push(RadixStep(rs.str - rs.bkt_size[rs.idx], rs.bkt_size[rs.idx],
                                          radixstack.size(), charcache));
                // cannot add here, because rs may have invalidated
            }
        }
        radixstack.pop();
    }

    delete[] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI_nr3, "bingmann/msd_CI_nr3",
               "bingmann/msd_CI_nr3 (CI non-recursive, charcache)")

/******************************************************************************/

} // namespace bingmann_radix_sort

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_RADIX_SORT_HEADER

/******************************************************************************/
