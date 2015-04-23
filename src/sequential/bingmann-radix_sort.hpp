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

static const size_t g_inssort_threshold = 64;

typedef unsigned char* string;

/******************************************************************************/

static inline void
msd_CE(string* strings, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // count character occurances
    size_t bktsize[256] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bktsize[strings[i][depth]];

    string* sorted = (string*)malloc(n * sizeof(string));

    // prefix sum
    size_t bktindex[256];
    bktindex[0] = 0;
    for (size_t i = 1; i < 256; ++i)
        bktindex[i] = bktindex[i - 1] + bktsize[i - 1];

    // distribute
    for (size_t i = 0; i < n; ++i)
        sorted[bktindex[strings[i][depth]]++] = strings[i];

    memcpy(strings, sorted, n * sizeof(string));
    free(sorted);

    // recursion
    size_t bsum = bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        msd_CE(strings + bsum, bktsize[i], depth + 1);
        bsum += bktsize[i];
    }
}

void bingmann_msd_CE(string* strings, size_t n)
{
    return msd_CE(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CE, "bingmann/msd_CE",
               "bingmann/msd_CE (rantala CE original)")

/******************************************************************************/

template <typename StringSet>
static inline void
msd_CE_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::String String;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    // count character occurances
    size_t bktsize[256] = { 0 };
    for (Iterator i = ss.begin(); i != ss.end(); ++i)
        ++bktsize[ss.get_uint8(ss[i], depth)];

    String* sorted = new String[ss.size()];

    // prefix sum
    size_t bktindex[256];
    bktindex[0] = 0;
    for (size_t i = 1; i < 256; ++i)
        bktindex[i] = bktindex[i - 1] + bktsize[i - 1];

    // distribute
    for (Iterator i = ss.begin(); i != ss.end(); ++i)
        sorted[bktindex[ss.get_uint8(ss[i], depth)]++] = ss[i];

    size_t i = 0;
    for (Iterator o = ss.begin(); o != ss.end(); ++o, ++i)
        ss[o] = sorted[i];
    delete[] sorted;

    // recursion
    Iterator bsum = ss.begin() + bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        Iterator bend = bsum + bktsize[i];
        msd_CE_generic(ss.sub(bsum, bend), depth + 1);
        bsum = bend;
    }
}

void bingmann_msd_CE_generic(string* strings, size_t n)
{
    return msd_CE_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CE_generic, "bingmann/msd_CE_gen",
               "bingmann/msd_CE generic (rantala CE original)")

/******************************************************************************/

static inline void
msd_CE2(string* strings, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // count character occurances
    size_t bkt[256 + 1] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bkt[strings[i][depth]];

    string* sorted = (string*)malloc(n * sizeof(string));

    // prefix sum
    for (size_t i = 1; i <= 256; ++i)
        bkt[i] += bkt[i - 1];

    // distribute
    for (size_t i = 0; i < n; ++i)
        sorted[--bkt[strings[i][depth]]] = strings[i];

    memcpy(strings, sorted, n * sizeof(string));
    free(sorted);

    // recursion
    for (size_t i = 1; i < 256; ++i) {
        if (bkt[i] == bkt[i + 1]) continue;
        //if (bkt[i]+1 >= bkt[i+1]) continue;
        msd_CE2(strings + bkt[i], bkt[i + 1] - bkt[i], depth + 1);
    }
}

void bingmann_msd_CE2(string* strings, size_t n)
{
    return msd_CE2(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CE2, "bingmann/msd_CE2",
               "bingmann/msd_CE2 (CE with reused prefix sum)")

/******************************************************************************/

template <typename StringSet>
static inline void
msd_CE2_generic(const StringSet& ss, size_t depth)
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
        sorted[--bkt[ss.get_uint8(ss[i], depth)]] = ss[i];

    size_t i = 0;
    for (Iterator o = ss.begin(); o != ss.end(); ++o, ++i)
        ss[o] = sorted[i];
    delete[] sorted;

    // recursion
    for (size_t i = 1; i < 256; ++i) {
        if (bkt[i] == bkt[i + 1]) continue;
        //if (bkt[i]+1 >= bkt[i+1]) continue;
        msd_CE2_generic(ss.sub(ss.begin() + bkt[i], ss.begin() + bkt[i + 1]),
                        depth + 1);
    }
}

void bingmann_msd_CE2_generic(string* strings, size_t n)
{
    return msd_CE2_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CE2_generic, "bingmann/msd_CE2_gen",
               "bingmann/msd_CE2 generic (CE with reused prefix sum)")

/******************************************************************************/

static inline void
msd_CE3(string* str_begin, string* str_end, size_t depth)
{
    if (str_begin + g_inssort_threshold > str_end)
        return inssort::inssort_range(str_begin, str_end, depth);

    // count character occurances
    size_t bkt[256 + 1] = { 0 };
    for (string* str = str_begin; str != str_end; ++str)
        ++bkt[(*str)[depth]];

    string* sorted = (string*)malloc((str_end - str_begin) * sizeof(string));

    // prefix sum
    for (size_t i = 1; i <= 256; ++i)
        bkt[i] += bkt[i - 1];

    // distribute
    for (string* str = str_begin; str != str_end; ++str)
        sorted[--bkt[(*str)[depth]]] = *str;

    memcpy(str_begin, sorted, (str_end - str_begin) * sizeof(string));
    free(sorted);

    // recursion
    for (size_t i = 1; i < 256; ++i) {
        if (bkt[i] == bkt[i + 1]) continue;
        //if (bkt[i]+1 >= bkt[i+1]) continue;
        msd_CE3(str_begin + bkt[i], str_begin + bkt[i + 1], depth + 1);
    }
}

void bingmann_msd_CE3(string* strings, size_t n)
{
    return msd_CE3(strings, strings + n, 0);
}

PSS_CONTESTANT(bingmann_msd_CE3, "bingmann/msd_CE3",
               "bingmann/msd_CE3 (CE2 with iterators)")

/******************************************************************************/

template <typename BucketType>
struct distblock {
    string     ptr;
    BucketType bkt;
};

template <typename BucketsizeType>
static void msd_CI(string* strings, size_t n, size_t depth)
{
    if (n < g_inssort_threshold) {
        inssort::inssort(strings, n, depth);
        return;
    }
    BucketsizeType bktsize[256] = { 0 };
    string restrict oracle = (string)malloc(n);
    for (size_t i = 0; i < n; ++i)
        oracle[i] = strings[i][depth];
    for (size_t i = 0; i < n; ++i)
        ++bktsize[oracle[i]];
    static size_t bktindex[256];
    bktindex[0] = bktsize[0];
    BucketsizeType last_bkt_size = bktsize[0];
    for (unsigned i = 1; i < 256; ++i) {
        bktindex[i] = bktindex[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }
    for (size_t i = 0; i < n - last_bkt_size; ) {
        distblock<uint8_t> tmp = { strings[i], oracle[i] };
        while (1) {
            // Continue until the current bucket is completely in
            // place
            if (--bktindex[tmp.bkt] <= i)
                break;
            // backup all information of the position we are about
            // to overwrite
            size_t backup_idx = bktindex[tmp.bkt];
            distblock<uint8_t> tmp2 = { strings[backup_idx], oracle[backup_idx] };
            // overwrite everything, ie. move the string to correct
            // position
            strings[backup_idx] = tmp.ptr;
            oracle[backup_idx] = tmp.bkt;
            tmp = tmp2;
        }
        // Commit last pointer to place. We don't need to copy the
        // oracle entry, it's not read after this.
        strings[i] = tmp.ptr;
        i += bktsize[tmp.bkt];
    }
    free(oracle);
    size_t bsum = bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        msd_CI<BucketsizeType>(strings + bsum, bktsize[i], depth + 1);
        bsum += bktsize[i];
    }
}

void bingmann_msd_CI(string* strings, size_t n)
{
    msd_CI<size_t>(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI, "bingmann/msd_CI",
               "bingmann/msd_CI (rantala CI original with oracle)")

/******************************************************************************/

template <typename BucketsizeType>
static inline void msd_CI2(string* strings, size_t n, size_t depth)
{
    if (n < g_inssort_threshold) {
        inssort::inssort(strings, n, depth);
        return;
    }
    BucketsizeType bktsize[256] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bktsize[strings[i][depth]];
    size_t bktindex[256];
    bktindex[0] = bktsize[0];
    BucketsizeType last_bkt_size = bktsize[0];
    for (unsigned i = 1; i < 256; ++i) {
        bktindex[i] = bktindex[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }
    for (size_t i = 0; i < n - last_bkt_size; ) {
        distblock<uint8_t> tmp = { strings[i], strings[i][depth] };
        while (1) {
            // Continue until the current bkt is completely in
            // place
            if (--bktindex[tmp.bkt] <= i)
                break;
            // backup all information of the position we are about
            // to overwrite
            size_t backup_idx = bktindex[tmp.bkt];
            distblock<uint8_t> tmp2 = { strings[backup_idx], strings[backup_idx][depth] };
            // overwrite everything, ie. move the string to correct
            // position
            strings[backup_idx] = tmp.ptr;
            tmp = tmp2;
        }
        // Commit last pointer to place. We don't need to copy the
        // oracle entry, it's not read after this.
        strings[i] = tmp.ptr;
        i += bktsize[tmp.bkt];
    }
    size_t bsum = bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        msd_CI2<BucketsizeType>(strings + bsum, bktsize[i], depth + 1);
        bsum += bktsize[i];
    }
}

void bingmann_msd_CI2(string* strings, size_t n)
{
    msd_CI2<size_t>(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI2, "bingmann/msd_CI2",
               "bingmann/msd_CI2 (CI without oracle)")

/******************************************************************************/

static inline void
msd_CI3(string* strings, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // count character occurances
    size_t bktsize[256] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bktsize[strings[i][depth]];

    // prefix sum
    size_t bktindex[256];
    bktindex[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned i = 1; i < 256; ++i) {
        bktindex[i] = bktindex[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }

    // premute in-place
    for (size_t i = 0, j; i < n - last_bkt_size; )
    {
        while ((j = --bktindex[strings[i][depth]]) > i)
        {
            std::swap(strings[i], strings[j]);
        }
        i += bktsize[strings[i][depth]];
    }

    // recursion
    size_t bsum = bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        msd_CI3(strings + bsum, bktsize[i], depth + 1);
        bsum += bktsize[i];
    }
}

void bingmann_msd_CI3(string* strings, size_t n)
{
    msd_CI3(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI3, "bingmann/msd_CI3",
               "bingmann/msd_CI3 (CI2 with swap operations)")

/******************************************************************************/

// Note: CI in-place variant cannot be done with just one prefix-sum bucket
// array, because during in-place permutation the beginning _and_ end
// boundaries of each bucket must be kept.

static inline void
msd_CI4(string* strings, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // count character occurances
    size_t bktsize[256] = { 0 };
    for (size_t i = 0; i < n; ++i)
        ++bktsize[strings[i][depth]];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned i = 1; i < 256; ++i) {
        bkt[i] = bkt[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
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
        i += bktsize[perm[depth]];
    }

    // recursion
    size_t bsum = bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        msd_CI4(strings + bsum, bktsize[i], depth + 1);
        bsum += bktsize[i];
    }
}

void bingmann_msd_CI4(string* strings, size_t n)
{
    msd_CI4(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI4, "bingmann/msd_CI4",
               "bingmann/msd_CI4 (CI3 with swap cache)")

/******************************************************************************/

static inline size_t *
msd_CI5_bktsize(string * strings, size_t n, size_t depth)
{
    // cache characters
    uint8_t* charcache = new uint8_t[n];
    for (size_t i = 0; i < n; ++i)
        charcache[i] = strings[i][depth];

    // count character occurances
    size_t* bktsize = new size_t[256];
    memset(bktsize, 0, 256 * sizeof(size_t));
    for (size_t i = 0; i < n; ++i)
        ++bktsize[charcache[i]];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned i = 1; i < 256; ++i) {
        bkt[i] = bkt[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
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
        i += bktsize[permch];
    }

    delete[] charcache;

    return bktsize;
}

void
msd_CI5(string* strings, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    size_t* bktsize = msd_CI5_bktsize(strings, n, depth);

    // recursion
    size_t bsum = bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        msd_CI5(strings + bsum, bktsize[i], depth + 1);
        bsum += bktsize[i];
    }

    delete[] bktsize;
}

void bingmann_msd_CI5(string* strings, size_t n)
{
    msd_CI5(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI5, "bingmann/msd_CI5",
               "bingmann/msd_CI5 (CI4 with charcache)")

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
    size_t* bktsize = new size_t[256];
    memset(bktsize, 0, 256 * sizeof(size_t));
    for (size_t i = 0; i < ss.size(); ++i)
        ++bktsize[static_cast<size_t>(charcache[i])];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned i = 1; i < 256; ++i) {
        bkt[i] = bkt[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }

    // premute in-place
    for (size_t i = 0, j; i < ss.size() - last_bkt_size; )
    {
        String perm = ss[ss.begin() + i];
        Char permch = charcache[i];
        while ((j = --bkt[static_cast<size_t>(permch)]) > i)
        {
            std::swap(perm, ss[ss.begin() + j]);
            std::swap(permch, charcache[j]);
        }
        ss[ss.begin() + i] = perm;
        i += bktsize[static_cast<size_t>(permch)];
    }

    delete[] charcache;

    return bktsize;
}

template <typename StringSet>
static inline void
msd_CI5_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    size_t* bktsize = msd_CI5_bktsize_generic(ss, depth);

    // recursion
    Iterator bsum = ss.begin() + bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        Iterator bend = bsum + bktsize[i];
        msd_CI5_generic(ss.sub(bsum, bend), depth + 1);
        bsum = bend;
    }

    delete[] bktsize;
}

void bingmann_msd_CI5_generic(string* strings, size_t n)
{
    msd_CI5_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CI5_generic, "bingmann/msd_CI5_gen",
               "bingmann/msd_CI5 generic (CI4 with charcache)")

/******************************************************************************/

static inline size_t *
msd_CI5_16bit_bktsize(string * strings, size_t n, size_t depth)
{
    static const size_t RADIX = 0x10000;
    using namespace stringtools;

    // cache characters
    uint16_t* charcache = new uint16_t[n];
    for (size_t i = 0; i < n; ++i)
        charcache[i] = get_char<uint16_t>(strings[i], depth);

    // count character occurances
    size_t* bktsize = new size_t[RADIX];
    memset(bktsize, 0, RADIX * sizeof(size_t));
    for (size_t i = 0; i < n; ++i)
        ++bktsize[charcache[i]];

    // inclusive prefix sum
    size_t bkt[RADIX];
    bkt[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned i = 1; i < RADIX; ++i) {
        bkt[i] = bkt[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
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
        i += bktsize[permch];
    }

    delete[] charcache;

    return bktsize;
}

static void
msd_CI5_16bit(string* strings, size_t n, size_t depth)
{
    if (n < 0x10000)
        return msd_CI5(strings, n, depth);

    size_t* bktsize = msd_CI5_16bit_bktsize(strings, n, depth);

    // recursion
    size_t bsum = bktsize[0];
    for (size_t i = 1; i < 0x10000; ++i) {
        if (bktsize[i] == 0) continue;
        if (i & 0xFF) // not zero-terminated
            msd_CI5_16bit(strings + bsum, bktsize[i], depth + 2);
        bsum += bktsize[i];
    }

    delete[] bktsize;
}

void bingmann_msd_CI5_16bit(string* strings, size_t n)
{
    msd_CI5_16bit(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI5_16bit, "bingmann/msd_CI5_16bit",
               "bingmann/msd_CI5_16bit (CI5 with 16-bit radix)")

/******************************************************************************/

template <typename StringSet>
static inline size_t*
msd_CI6_bktsize_generic(const StringSet& ss, size_t depth)
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
    size_t* bktsize = new size_t[256];
    memset(bktsize, 0, 256 * sizeof(size_t));
    for (size_t i = 0; i < ss.size(); ++i)
        ++bktsize[static_cast<size_t>(charcache[i])];

    // inclusive prefix sum
    Iterator bkt[256];
    bkt[0] = ss.begin() + bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned i = 1; i < 256; ++i) {
        bkt[i] = bkt[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }

    // premute in-place
    for (Iterator i = ss.begin(), j; i != ss.end() - last_bkt_size; )
    {
        String perm = ss[i];
        Char permch = charcache[i - ss.begin()];
        while ((j = --bkt[static_cast<size_t>(permch)]) > i)
        {
            std::swap(perm, ss[j]);
            std::swap(permch, charcache[j - ss.begin()]);
        }
        ss[i] = perm;
        i += bktsize[static_cast<size_t>(permch)];
    }

    delete[] charcache;

    return bktsize;
}

template <typename StringSet>
static inline void
msd_CI6_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    size_t* bktsize = msd_CI6_bktsize_generic(ss, depth);

    // recursion
    Iterator bsum = ss.begin() + bktsize[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bktsize[i] == 0) continue;
        Iterator bend = bsum + bktsize[i];
        msd_CI6_generic(ss.sub(bsum, bend), depth + 1);
        bsum = bend;
    }

    delete[] bktsize;
}

void bingmann_msd_CI6_generic(string* strings, size_t n)
{
    msd_CI6_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CI6_generic, "bingmann/msd_CI6_gen",
               "bingmann/msd_CI6 generic (CI5 with ptr bkts)")

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
    size_t bktsize[256];
    size_t idx;

    RadixStep_CE_nr2(string* strings, size_t n, size_t depth)
    {
        // count character occurances
        memset(bktsize, 0, sizeof(bktsize));
        for (size_t i = 0; i < n; ++i)
            ++bktsize[strings[i][depth]];

        // prefix sum
        size_t bkt[256];
        bkt[0] = 0;
        for (size_t i = 1; i < 256; ++i)
            bkt[i] = bkt[i - 1] + bktsize[i - 1];

        // distribute out-of-place
        string* sorted = (string*)malloc(n * sizeof(string));

        for (size_t i = 0; i < n; ++i)
            sorted[bkt[strings[i][depth]]++] = strings[i];

        memcpy(strings, sorted, n * sizeof(string));
        free(sorted);

        str = strings + bktsize[0];
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

            if (rs.bktsize[rs.idx] == 0)
                ;
            else if (rs.bktsize[rs.idx] < g_inssort_threshold)
            {
                inssort::inssort(rs.str, rs.bktsize[rs.idx], radixstack.size());
                rs.str += rs.bktsize[rs.idx];
            }
            else
            {
                rs.str += rs.bktsize[rs.idx];
                radixstack.push(RadixStep(rs.str - rs.bktsize[rs.idx], rs.bktsize[rs.idx], radixstack.size()));
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
        size_t bktsize[256] = { 0 };
        for (size_t i = 0; i < n; ++i)
            ++bktsize[strings[i][depth]];

        // inclusive prefix sum
        bkt[0] = bktsize[0];
        size_t last_bkt_size = bktsize[0], last_bkt_num = 0;
        for (unsigned i = 1; i < 256; ++i) {
            bkt[i] = bkt[i - 1] + bktsize[i];
            if (bktsize[i]) {
                last_bkt_size = bktsize[i];
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
                i += bktsize[perm[depth]];
                continue;
            }
            while ((j = --bkt[perm[depth]]) > i)
            {
                std::swap(perm, strings[j]);
            }
            assert(j == i);
            strings[i] = perm;
            i += bktsize[perm[depth]];
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
    size_t bktsize[256];

    RadixStep_CI_nr2(string* strings, size_t n, size_t depth)
    {
        // count character occurances
        memset(bktsize, 0, sizeof(bktsize));
        for (size_t i = 0; i < n; ++i)
            ++bktsize[strings[i][depth]];

        // inclusive prefix sum
        size_t bkt[256];
        bkt[0] = bktsize[0];
        size_t last_bkt_size = bktsize[0];
        for (unsigned i = 1; i < 256; ++i) {
            bkt[i] = bkt[i - 1] + bktsize[i];
            if (bktsize[i]) last_bkt_size = bktsize[i];
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
            i += bktsize[perm[depth]];
        }

        str = strings + bktsize[0];
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

            if (rs.bktsize[rs.idx] == 0)
                ;
            else if (rs.bktsize[rs.idx] < g_inssort_threshold)
            {
                inssort::inssort(rs.str, rs.bktsize[rs.idx], radixstack.size());
                rs.str += rs.bktsize[rs.idx];
            }
            else
            {
                rs.str += rs.bktsize[rs.idx];
                radixstack.push(RadixStep(rs.str - rs.bktsize[rs.idx], rs.bktsize[rs.idx], radixstack.size()));
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
    size_t bktsize[256];

    RadixStep_CI_nr3(string* strings, size_t n, size_t depth, uint8_t* charcache)
    {
        // cache characters
        for (size_t i = 0; i < n; ++i)
            charcache[i] = strings[i][depth];

        // count character occurances
        memset(bktsize, 0, sizeof(bktsize));
        for (size_t i = 0; i < n; ++i)
            ++bktsize[charcache[i]];

        // inclusive prefix sum
        size_t bkt[256];
        bkt[0] = bktsize[0];
        size_t last_bkt_size = bktsize[0];
        for (unsigned i = 1; i < 256; ++i) {
            bkt[i] = bkt[i - 1] + bktsize[i];
            if (bktsize[i]) last_bkt_size = bktsize[i];
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
            i += bktsize[permch];
        }

        str = strings + bktsize[0];
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

            if (rs.bktsize[rs.idx] == 0)
                ;
            else if (rs.bktsize[rs.idx] < g_inssort_threshold)
            {
                inssort::inssort(rs.str, rs.bktsize[rs.idx], radixstack.size());
                rs.str += rs.bktsize[rs.idx];
            }
            else
            {
                rs.str += rs.bktsize[rs.idx];
                radixstack.push(RadixStep(rs.str - rs.bktsize[rs.idx], rs.bktsize[rs.idx],
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
