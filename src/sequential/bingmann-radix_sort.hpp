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
msd_CE1_c(string* strings, string* sorted, uint8_t* charcache,
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
        msd_CE1_c(strings + bkt[i], sorted, charcache,
                  bkt[i + 1] - bkt[i], depth + 1);
    }
}

void bingmann_msd_CE1_c(string* strings, size_t n)
{
    string* sorted = new string[n];
    uint8_t* charcache = new uint8_t[n];

    msd_CE1_c(strings, sorted, charcache, n, 0);

    delete [] charcache;
    delete [] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE1_c, "bingmann/msd_CE1_c",
               "bingmann/msd_CE1_c (CE1 with charcache)")

/******************************************************************************/

static inline void
msd_CE1_cf(string* strings, string* sorted, uint8_t* charcache,
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
        msd_CE1_cf(strings + bkt[i], sorted, charcache,
                   bkt[i + 1] - bkt[i], depth + 1);
    }
}

void bingmann_msd_CE1_cf(string* strings, size_t n)
{
    string* sorted = new string[n];
    uint8_t* charcache = new uint8_t[n];

    msd_CE1_cf(strings, sorted, charcache, n, 0);

    delete [] charcache;
    delete [] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE1_cf, "bingmann/msd_CE1_cf",
               "bingmann/msd_CE1_cf (CE1 with charcache, fissioned)")

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
// msd_CI0_c

static inline size_t *
msd_CI0_c_make_bkt_size(
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
msd_CI0_c(string* strings, uint8_t* charcache, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    size_t* bkt_size = msd_CI0_c_make_bkt_size(strings, charcache, n, depth);

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CI0_c(strings + bsum, charcache, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

void bingmann_msd_CI0_c(string* strings, size_t n)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI0_c(strings, charcache, n, 0);
    delete [] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI0_c, "bingmann/msd_CI0_c",
               "bingmann/msd_CI0_c (CI0 with charcache)")

/******************************************************************************/
// msd_CI0_cf

static inline size_t *
msd_CI0_cf_run(
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
msd_CI0_cf(string* strings, uint8_t* charcache, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    size_t* bkt_size = msd_CI0_cf_run(strings, charcache, n, depth);

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CI0_cf(strings + bsum, charcache, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

void bingmann_msd_CI0_cf(string* strings, size_t n)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI0_cf(strings, charcache, n, 0);
    delete [] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI0_cf, "bingmann/msd_CI0_cf",
               "bingmann/msd_CI0_cf (CI0 with charcache, fissioned)")

void msd_CI(string* strings, size_t n, size_t depth)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI0_cf(strings, charcache, n, depth);
    delete [] charcache;
}

/******************************************************************************/

template <typename StringSet>
static inline size_t*
msd_CI0_cf_bktsize_generic(
    const StringSet& ss, typename StringSet::Char* charcache, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::String String;
    typedef typename StringSet::Char Char;

    // cache characters
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

    return bkt_size;
}

template <typename StringSet>
static inline void
msd_CI0_cf_generic(
    const StringSet& ss, typename StringSet::Char* charcache, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    size_t* bkt_size = msd_CI0_cf_bktsize_generic(ss, charcache, depth);

    // recursion
    Iterator bsum = ss.begin() + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] == 0) continue;
        Iterator bend = bsum + bkt_size[i];
        msd_CI0_cf_generic(ss.sub(bsum, bend), charcache, depth + 1);
        bsum = bend;
    }

    delete[] bkt_size;
}

void bingmann_msd_CI0_cf_generic(string* strings, size_t n)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI0_cf_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), charcache, 0);
    delete[] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI0_cf_generic, "bingmann/msd_CI0_cf_gen",
               "bingmann/msd_CI0_c generic (CI0 with charcache, fissioned)")

/******************************************************************************/

static inline size_t *
msd_CI0_cf_16bit_run(string * strings, uint16_t* charcache, size_t n, size_t depth)
{
    static const size_t RADIX = 0x10000;
    using namespace stringtools;

    // cache characters
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

    return bkt_size;
}

static void
msd_CI0_cf_16bit(string* strings, uint16_t* charcache, size_t n, size_t depth)
{
    if (n < 0x10000)
        return msd_CI0_cf(
            strings, reinterpret_cast<uint8_t*>(charcache), n, depth);

    size_t* bkt_size = msd_CI0_cf_16bit_run(strings, charcache, n, depth);

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 0x10000; ++i) {
        if (bkt_size[i] == 0) continue;
        if (i & 0xFF) // not zero-terminated
            msd_CI0_cf_16bit(strings + bsum, charcache, bkt_size[i], depth + 2);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

void bingmann_msd_CI0_cf_16bit(string* strings, size_t n)
{
    uint16_t* charcache = new uint16_t[n];
    msd_CI0_cf_16bit(strings, charcache, n, 0);
    delete [] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI0_cf_16bit, "bingmann/msd_CI0_cf_16bit",
               "bingmann/msd_CI0_cf_16bit (CI0_cf with 16-bit radix)")

/******************************************************************************/
// Stack-Based Variants
/******************************************************************************/

struct RadixStep_CE0_sb
{
    string * str;
    size_t bkt_size[256];
    size_t idx;

    RadixStep_CE0_sb(string* strings, size_t n, size_t depth)
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
bingmann_msd_CE0_sb(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    typedef RadixStep_CE0_sb RadixStep;

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
                radixstack.push(
                    RadixStep(rs.str - rs.bkt_size[rs.idx],
                              rs.bkt_size[rs.idx], radixstack.size()));
                // cannot add here, because rs may have invalidated
            }
        }
        radixstack.pop();
    }
}

PSS_CONTESTANT(bingmann_msd_CE0_sb, "bingmann/msd_CE0_sb",
               "bingmann/msd_CE0_sb (CE0 stack-based)")

/******************************************************************************/

struct RadixStep_CE1_cf_sb
{
    string * str;
    size_t bkt[256 + 1];
    size_t idx;

    RadixStep_CE1_cf_sb(string* strings, size_t n, size_t depth,
                        string* sorted, uint8_t* charcache)
    {
        // cache characters
        for (size_t i = 0; i < n; ++i)
            charcache[i] = static_cast<uint8_t>(strings[i][depth]);

        // count character occurances
        memset(bkt, 0, sizeof(bkt));
        for (size_t i = 0; i < n; ++i)
            ++bkt[charcache[i]];

        // prefix sum
        for (unsigned i = 1; i <= 256; ++i)
            bkt[i] += bkt[i - 1];

        // distribute out-of-place
        for (size_t i = 0; i < n; ++i)
            sorted[--bkt[charcache[i]]] = strings[i];

        std::copy(sorted, sorted + n, strings);

        str = strings;
        idx = 0;        // will increment to 1 on first process
    }
};

static inline void
bingmann_msd_CE1_cf_sb(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    typedef RadixStep_CE1_cf_sb RadixStep;

    string* sorted = new string[n];
    uint8_t* charcache = new uint8_t[n];

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.push(RadixStep(strings, n, 0, sorted, charcache));

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
                radixstack.push(
                    RadixStep(rs.str + rs.bkt[rs.idx],
                              rs.bkt[rs.idx + 1] - rs.bkt[rs.idx],
                              radixstack.size(), sorted, charcache));
                break;
            }
        }

        // rs maybe have been invalidated
        if (radixstack.top().idx == 256)
            radixstack.pop();
    }

    delete [] charcache;
    delete [] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE1_cf_sb, "bingmann/msd_CE1_cf_sb",
               "bingmann/msd_CE1_cf_sb (CE1 stack-based, charcache, fissioned)")

/******************************************************************************/

struct RadixStep_CI0_cf_sb
{
    string * str;
    size_t idx;
    size_t bkt_size[256];

    RadixStep_CI0_cf_sb(string* strings, size_t n, size_t depth, uint8_t* charcache)
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
bingmann_msd_CI0_cf_sb(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    typedef RadixStep_CI0_cf_sb RadixStep;

    uint8_t* charcache = new uint8_t[n];

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.push(RadixStep(strings, n, 0, charcache));

    while (TLX_LIKELY(radixstack.size()))
    {
        while (TLX_LIKELY(radixstack.top().idx < 255))
        {
            RadixStep& rs = radixstack.top();
            ++rs.idx; // process the bucket rs.idx

            if (TLX_UNLIKELY(rs.bkt_size[rs.idx] == 0))
                ;
            else if (TLX_UNLIKELY(rs.bkt_size[rs.idx] < g_inssort_threshold))
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

PSS_CONTESTANT(bingmann_msd_CI0_cf_sb, "bingmann/msd_CI0_cf_sb",
               "bingmann/msd_CI0_cf_sb (CI stack-based, charcache, fissioned)")

/******************************************************************************/

} // namespace bingmann_radix_sort

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_RADIX_SORT_HEADER

/******************************************************************************/
