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

namespace bingmann {

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
    for (string* s = strings; s != strings + n; ++s)
        ++bkt_size[(*s)[depth]];

    {
        // prefix sum
        string* bkt[256];
        bkt[0] = sorted;
        for (size_t i = 1; i < 256; ++i)
            bkt[i] = bkt[i - 1] + bkt_size[i - 1];

        // distribute
        for (string* s = strings; s != strings + n; ++s)
            *(bkt[(*s)[depth]]++) = *s;

        std::copy(sorted, sorted + n, strings);
    }

    // recursion
    string* bsum = strings + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CE0(bsum, sorted, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }
}

static inline void bingmann_msd_CE0(string* strings, size_t n)
{
    string* sorted = new string[n];
    msd_CE0(strings, sorted, n, 0);
    delete[] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE0, "bingmann/msd_CE0",
               "bingmann/msd_CE0 (rantala CE0 baseline)")

/******************************************************************************/

template <typename StringSet>
static inline void
msd_CE0_generic(const StringSet& ss, const StringSet& sorted, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    // count character occurances
    size_t bkt_size[256] = { 0 };
    for (Iterator i = ss.begin(); i != ss.end(); ++i)
        ++bkt_size[ss.get_uint8(ss[i], depth)];

    {
        // prefix sum
        Iterator bkt_index[256];
        bkt_index[0] = sorted.begin();
        for (size_t i = 1; i < 256; ++i)
            bkt_index[i] = bkt_index[i - 1] + bkt_size[i - 1];

        // distribute
        for (Iterator i = ss.begin(); i != ss.end(); ++i)
            *(bkt_index[ss.get_uint8(ss[i], depth)]++) = std::move(ss[i]);

        std::move(sorted.begin(), sorted.begin() + ss.size(), ss.begin());
    }

    // recursion
    Iterator bsum = ss.begin() + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        Iterator bsum_begin = bsum;
        bsum += bkt_size[i];
        if (bkt_size[i] > 1)
            msd_CE0_generic(ss.sub(bsum_begin, bsum), sorted, depth + 1);
    }
}

template <typename StringSet>
static inline void
msd_CE0_generic(const StringSet& ss, size_t depth)
{
    typename StringSet::Container sorted = ss.allocate(ss.size());
    msd_CE0_generic(ss, StringSet(sorted), depth);
    StringSet::deallocate(sorted);
}

static inline void bingmann_msd_CE0_generic(string* strings, size_t n)
{
    msd_CE0_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CE0_generic, "bingmann/msd_CE0_gen",
               "bingmann/msd_CE0 generic (rantala CE baseline)")

/******************************************************************************/

static inline void
msd_CE1(string* strings, string* sorted, uint8_t* charcache,
        size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // read characters and count character occurances
    size_t bkt_size[256] = { 0 };
    uint8_t* cc = charcache;
    for (string* s = strings; s != strings + n; ++s, ++cc)
        ++bkt_size[*cc = (*s)[depth]];

    {
        // prefix sum
        string* bkt[256];
        bkt[0] = sorted;
        for (size_t i = 1; i < 256; ++i)
            bkt[i] = bkt[i - 1] + bkt_size[i - 1];

        // distribute
        uint8_t* cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *(bkt[*cc]++) = *s;

        std::copy(sorted, sorted + n, strings);
    }

    // recursion
    string* bsum = strings + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CE1(bsum, sorted, charcache, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }
}

static inline void bingmann_msd_CE1(string* strings, size_t n)
{
    string* sorted = new string[n];
    uint8_t* charcache = new uint8_t[n];

    msd_CE1(strings, sorted, charcache, n, 0);

    delete[] charcache;
    delete[] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE1, "bingmann/msd_CE1",
               "bingmann/msd_CE1 (with charcache, fused loop)")

/******************************************************************************/

static inline void
msd_CE2(string* strings, string* sorted, uint8_t* charcache,
        size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    // read characters and count character occurances
    size_t bkt_size[256] = { 0 };
    uint8_t* cc = charcache;
    for (string* s = strings; s != strings + n; ++s, ++cc)
        *cc = (*s)[depth];
    for (cc = charcache; cc != charcache + n; ++cc)
        ++bkt_size[*cc];

    {
        // prefix sum
        string* bkt[256];
        bkt[0] = sorted;
        for (size_t i = 1; i < 256; ++i)
            bkt[i] = bkt[i - 1] + bkt_size[i - 1];

        // distribute out-of-place
        uint8_t* cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *(bkt[*cc]++) = *s;

        std::copy(sorted, sorted + n, strings);
    }

    // recursion
    string* bsum = strings + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CE2(bsum, sorted, charcache, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }
}

static inline void bingmann_msd_CE2(string* strings, size_t n)
{
    string* sorted = new string[n];
    uint8_t* charcache = new uint8_t[n];

    msd_CE2(strings, sorted, charcache, n, 0);

    delete[] charcache;
    delete[] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE2, "bingmann/msd_CE2",
               "bingmann/msd_CE2 (with charcache, fissioned loop)")

/******************************************************************************/

static inline void
msd_CE3(string* strings, string* sorted, uint16_t* charcache,
        size_t n, size_t depth)
{
    static const size_t RADIX = 0x10000;
    using namespace stringtools;

    if (n < 0x10000)
        return msd_CE2(strings, sorted,
                       reinterpret_cast<uint8_t*>(charcache), n, depth);

    // read characters and count character occurances
    size_t bkt_size[RADIX] = { 0 };
    uint16_t* cc = charcache;
    for (string* s = strings; s != strings + n; ++s, ++cc)
        *cc = get_char<uint16_t>(*s, depth);
    for (cc = charcache; cc != charcache + n; ++cc)
        ++bkt_size[*cc];

    {
        // prefix sum
        string* bkt[RADIX];
        bkt[0] = sorted;
        for (size_t i = 1; i < RADIX; ++i)
            bkt[i] = bkt[i - 1] + bkt_size[i - 1];

        // distribute
        uint16_t* cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *(bkt[*cc]++) = *s;

        std::copy(sorted, sorted + n, strings);
    }

    // recursion
    string* bsum = strings + bkt_size[0];
    for (size_t i = 1; i < RADIX; ++i) {
        if ((i & 0xFF) != 0 && bkt_size[i] > 1) // not zero-terminated
            msd_CE3(bsum, sorted, charcache, bkt_size[i], depth + 2);
        bsum += bkt_size[i];
    }
}

static inline void bingmann_msd_CE3(string* strings, size_t n)
{
    string* sorted = new string[n];
    uint16_t* charcache = new uint16_t[n];

    msd_CE3(strings, sorted, charcache, n, 0);

    delete[] charcache;
    delete[] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE3, "bingmann/msd_CE3",
               "bingmann/msd_CE3 (with charcache, fissioned loop, 16-bit adaptive)")

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
    for (string* s = strings; s != strings + n; ++s)
        ++bkt_size[(*s)[depth]];

    {
        // inclusive prefix sum
        string* bkt[256];
        bkt[0] = strings + bkt_size[0];
        size_t last_bkt_size = bkt_size[0];
        for (size_t i = 1; i < 256; ++i) {
            bkt[i] = bkt[i - 1] + bkt_size[i];
            if (bkt_size[i]) last_bkt_size = bkt_size[i];
        }

        // premute in-place
        for (string* i = strings, * j; i < strings + n - last_bkt_size; )
        {
            string perm = *i;
            while ((j = --bkt[perm[depth]]) > i)
            {
                std::swap(perm, *j);
            }
            *i = perm;
            i += bkt_size[perm[depth]];
        }
    }

    // recursion
    string* bsum = strings + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CI0(bsum, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }
}

static inline void bingmann_msd_CI0(string* strings, size_t n)
{
    msd_CI0(strings, n, 0);
}

PSS_CONTESTANT(bingmann_msd_CI0, "bingmann/msd_CI0",
               "bingmann/msd_CI0 (in-place baseline)")

/******************************************************************************/
// msd_CI1

static inline size_t *
msd_CI1_make_bkt_size(
    string * strings, uint8_t * charcache, size_t n, size_t depth)
{
    // cache and count character occurrences
    size_t* bkt_size = new size_t[256]();
    uint8_t* cc = charcache;
    for (string* s = strings; s != strings + n; ++s, ++cc)
        ++bkt_size[*cc = (*s)[depth]];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
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

static inline
void msd_CI1(string* strings, uint8_t* charcache, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    size_t* bkt_size = msd_CI1_make_bkt_size(strings, charcache, n, depth);

    // recursion
    string* bsum = strings + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CI1(bsum, charcache, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

static inline void bingmann_msd_CI1(string* strings, size_t n)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI1(strings, charcache, n, 0);
    delete[] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI1, "bingmann/msd_CI1",
               "bingmann/msd_CI1 (with charcache, fused loop)")

/******************************************************************************/
// msd_CI2

static inline size_t *
msd_CI2_run(string * strings, uint8_t * charcache, size_t n, size_t depth)
{
    // cache characters
    uint8_t* cc = charcache;
    for (string* s = strings; s != strings + n; ++s, ++cc)
        *cc = (*s)[depth];

    // cache and count character occurrences
    size_t* bkt_size = new size_t[256]();
    for (cc = charcache; cc != charcache + n; ++cc)
        ++bkt_size[*cc];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
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

static inline void msd_CI2(
    string* strings, uint8_t* charcache, size_t n, size_t depth)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, depth);

    size_t* bkt_size = msd_CI2_run(strings, charcache, n, depth);

    // recursion
    string* bsum = strings + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        if (bkt_size[i] > 1)
            msd_CI2(bsum, charcache, bkt_size[i], depth + 1);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

static inline void bingmann_msd_CI2(string* strings, size_t n)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI2(strings, charcache, n, 0);
    delete[] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI2, "bingmann/msd_CI2",
               "bingmann/msd_CI2 (with charcache, fissioned loop)")

static inline void msd_CI(string* strings, size_t n, size_t depth)
{
    uint8_t* charcache = new uint8_t[n];
    msd_CI2(strings, charcache, n, depth);
    delete[] charcache;
}

/******************************************************************************/

template <typename StringSet>
static inline size_t*
msd_CI2_bktsize_generic(
    const StringSet& ss, typename StringSet::Char* charcache, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::String String;
    typedef typename StringSet::Char Char;

    // cache characters
    Char* cc = charcache;
    for (Iterator i = ss.begin(); i != ss.end(); ++i, ++cc)
        *cc = *ss.get_chars(ss[i], depth);

    // count character occurances
    size_t* bkt_size = new size_t[256]();
    for (cc = charcache; cc != charcache + ss.size(); ++cc)
        ++bkt_size[static_cast<size_t>(*cc)];

    // inclusive prefix sum
    size_t bkt[256];
    bkt[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
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
msd_CI2_generic(
    const StringSet& ss, typename StringSet::Char* charcache, size_t depth)
{
    typedef typename StringSet::Iterator Iterator;

    if (ss.size() < g_inssort_threshold)
        return inssort::inssort_generic(ss, depth);

    size_t* bkt_size = msd_CI2_bktsize_generic(ss, charcache, depth);

    // recursion
    Iterator bsum = ss.begin() + bkt_size[0];
    for (size_t i = 1; i < 256; ++i) {
        Iterator bend = bsum + bkt_size[i];
        if (bkt_size[i] > 1)
            msd_CI2_generic(ss.sub(bsum, bend), charcache, depth + 1);
        bsum = bend;
    }

    delete[] bkt_size;
}

template <typename StringSet>
static inline void
msd_CI2_generic(const StringSet& ss, size_t depth)
{
    typedef typename StringSet::Char Char;
    Char* charcache = new Char[ss.size()];
    msd_CI2_generic(ss, charcache, depth);
    delete[] charcache;
}

static inline void bingmann_msd_CI2_generic(string* strings, size_t n)
{
    msd_CI2_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_msd_CI2_generic, "bingmann/msd_CI2_gen",
               "bingmann/msd_CI2 generic (with charcache, fissioned loop)")

/******************************************************************************/

static inline
size_t * msd_CI3_run(
    string * strings, uint16_t * charcache, size_t n, size_t depth)
{
    static const size_t RADIX = 0x10000;

    // cache characters
    uint16_t* cc = charcache;
    for (string* s = strings; s != strings + n; ++s, ++cc)
        *cc = stringtools::get_char<uint16_t>(*s, depth);

    // count character occurances
    size_t* bkt_size = new size_t[RADIX]();
    for (cc = charcache; cc != charcache + n; ++cc)
        ++bkt_size[*cc];

    // inclusive prefix sum
    size_t bkt[RADIX];
    bkt[0] = bkt_size[0];
    size_t last_bkt_size = bkt_size[0];
    for (size_t i = 1; i < RADIX; ++i) {
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

static inline
void msd_CI3(string* strings, uint16_t* charcache, size_t n, size_t depth)
{
    if (n < 0x10000)
        return msd_CI2(
            strings, reinterpret_cast<uint8_t*>(charcache), n, depth);

    size_t* bkt_size = msd_CI3_run(strings, charcache, n, depth);

    // recursion
    size_t bsum = bkt_size[0];
    for (size_t i = 1; i < 0x10000; ++i) {
        if ((i & 0xFF) != 0 && bkt_size[i] > 1) // not zero-terminated
            msd_CI3(strings + bsum, charcache, bkt_size[i], depth + 2);
        bsum += bkt_size[i];
    }

    delete[] bkt_size;
}

static inline void bingmann_msd_CI3(string* strings, size_t n)
{
    uint16_t* charcache = new uint16_t[n];
    msd_CI3(strings, charcache, n, 0);
    delete[] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI3, "bingmann/msd_CI3",
               "bingmann/msd_CI3 (with charcache, fissioned loop, 16-bit adaptive)")

/******************************************************************************/
// Stack-Based Variants
/******************************************************************************/

struct RadixStep_CE0_sb
{
    string * str;
    size_t idx;
    size_t bkt_size[256];

    RadixStep_CE0_sb(string* strings, string* sorted, size_t n, size_t depth)
    {
        // count character occurances
        memset(bkt_size, 0, sizeof(bkt_size));
        for (string* s = strings; s != strings + n; ++s)
            ++bkt_size[(*s)[depth]];

        // prefix sum
        string* bkt[256];
        bkt[0] = sorted;
        for (size_t i = 1; i < 256; ++i)
            bkt[i] = bkt[i - 1] + bkt_size[i - 1];

        // distribute out-of-place

        for (string* s = strings; s != strings + n; ++s)
            *(bkt[(*s)[depth]]++) = *s;

        std::copy(sorted, sorted + n, strings);

        str = strings + bkt_size[0];
        idx = 0; // will increment to 1 on first process
    }
};

static inline void
bingmann_msd_CE0_sb(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    string* sorted = new string[n];

    typedef RadixStep_CE0_sb RadixStep;

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.emplace(strings, sorted, n, 0);

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
                radixstack.emplace(
                    rs.str - rs.bkt_size[rs.idx], sorted,
                    rs.bkt_size[rs.idx], radixstack.size());
                // cannot add here, because rs may have invalidated
            }
        }
        radixstack.pop();
    }

    delete[] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE0_sb, "bingmann/msd_CE0_sb",
               "bingmann/msd_CE0_sb (CE0 stack-based)")

/******************************************************************************/

struct RadixStep_CE2_sb
{
    string * str;
    size_t idx;
    size_t bkt_size[256];

    RadixStep_CE2_sb(string* strings, size_t n, size_t depth,
                     string* sorted, uint8_t* charcache)
    {
        // read characters and count character occurances
        memset(bkt_size, 0, sizeof(bkt_size));
        uint8_t* cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *cc = (*s)[depth];
        for (cc = charcache; cc != charcache + n; ++cc)
            ++bkt_size[*cc];

        // exclusive pointer prefix sum
        string* bkt[256];
        bkt[0] = sorted;
        for (size_t i = 1; i < 256; ++i) {
            bkt[i] = bkt[i - 1] + bkt_size[i - 1];
        }

        // distribute out-of-place
        cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *(bkt[*cc]++) = *s;

        std::copy(sorted, sorted + n, strings);

        str = strings + bkt_size[0];
        idx = 0;        // will increment to 1 on first process
    }
};

static inline void
bingmann_msd_CE2_sb(
    string* strings, size_t n, string* sorted, uint8_t* charcache, size_t depth)
{
    typedef RadixStep_CE2_sb RadixStep;

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.emplace(strings, n, depth, sorted, charcache);

    while (TLX_LIKELY(!radixstack.empty()))
    {
        while (TLX_LIKELY(radixstack.top().idx < 255))
        {
            RadixStep& rs = radixstack.top();
            ++rs.idx; // process the bucket rs.idx

            if (TLX_UNLIKELY(rs.bkt_size[rs.idx] == 0))
                ;
            else if (TLX_UNLIKELY(rs.bkt_size[rs.idx] < g_inssort_threshold))
            {
                inssort::inssort(rs.str, rs.bkt_size[rs.idx],
                                 depth + radixstack.size());
                rs.str += rs.bkt_size[rs.idx];
            }
            else
            {
                // have to increment first, as rs may be invalidated
                rs.str += rs.bkt_size[rs.idx];
                radixstack.emplace(
                    rs.str - rs.bkt_size[rs.idx], rs.bkt_size[rs.idx],
                    depth + radixstack.size(), sorted, charcache);
            }
        }
        radixstack.pop();
    }
}

static inline void
bingmann_msd_CE2_sb(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, /* depth */ 0);

    string* sorted = new string[n];
    uint8_t* charcache = new uint8_t[n];

    bingmann_msd_CE2_sb(strings, n, sorted, charcache, /* depth */ 0);

    delete[] charcache;
    delete[] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE2_sb, "bingmann/msd_CE2_sb",
               "bingmann/msd_CE2_sb (CE2 stack-based, charcache, fissioned loop)")

/******************************************************************************/

struct RadixStep_CE3_sb
{
    static const size_t RADIX = 0x10000;

    string              * str;
    size_t              idx;
    size_t              bkt_size[RADIX];

    RadixStep_CE3_sb(string* strings, size_t n, size_t depth,
                     string* sorted, uint16_t* charcache)
    {
        // read characters and count character occurrences
        memset(bkt_size, 0, sizeof(bkt_size));
        uint16_t* cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *cc = stringtools::get_char<uint16_t>(*s, depth);
        for (cc = charcache; cc != charcache + n; ++cc)
            ++bkt_size[*cc];

        // exclusive pointer prefix sum
        string* bkt[RADIX];
        bkt[0] = sorted;
        for (size_t i = 1; i < RADIX; ++i) {
            bkt[i] = bkt[i - 1] + bkt_size[i - 1];
        }

        // distribute out-of-place
        cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *(bkt[*cc]++) = *s;

        std::copy(sorted, sorted + n, strings);

        str = strings + bkt_size[0];
        idx = 0;        // will increment to 1 on first process
    }
};

static inline void
bingmann_msd_CE3_sb(string* strings, size_t n)
{
    static const size_t RADIX = 0x10000;

    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    if (n < RADIX)
        return bingmann_msd_CE2_sb(strings, n);

    typedef RadixStep_CE3_sb RadixStep;

    string* sorted = new string[n];
    uint16_t* charcache = new uint16_t[n];

    typedef std::stack<RadixStep, std::vector<RadixStep> > radixstack_type;
    radixstack_type radixstack;
    radixstack.emplace(strings, n, 0, sorted, charcache);

    while (TLX_LIKELY(!radixstack.empty()))
    {
        while (TLX_LIKELY(radixstack.top().idx < RADIX))
        {
            RadixStep& rs = radixstack.top();
            ++rs.idx;                                      // process the bucket rs.idx

            if (TLX_UNLIKELY(rs.bkt_size[rs.idx] == 0))
                ;
            else if (TLX_UNLIKELY((rs.idx & 0xFF) == 0)) { // zero-termination
                rs.str += rs.bkt_size[rs.idx];
            }
            else if (TLX_UNLIKELY(rs.bkt_size[rs.idx] < g_inssort_threshold))
            {
                inssort::inssort(
                    rs.str, rs.bkt_size[rs.idx], 2* radixstack.size());
                rs.str += rs.bkt_size[rs.idx];
            }
            else if (rs.bkt_size[rs.idx] < RADIX)
            {
                bingmann_msd_CE2_sb(rs.str, rs.bkt_size[rs.idx], sorted,
                                    reinterpret_cast<uint8_t*>(charcache),
                                    /* depth */ 2 * radixstack.size());
                rs.str += rs.bkt_size[rs.idx];
            }
            else
            {
                // have to increment first, as rs may be invalidated
                rs.str += rs.bkt_size[rs.idx];
                radixstack.emplace(
                    rs.str - rs.bkt_size[rs.idx], rs.bkt_size[rs.idx],
                    2 * radixstack.size(), sorted, charcache);
            }
        }
        radixstack.pop();
    }

    delete[] charcache;
    delete[] sorted;
}

PSS_CONTESTANT(bingmann_msd_CE3_sb, "bingmann/msd_CE3_sb",
               "bingmann/msd_CE3_sb (CE3 stack-based, charcache, fissioned loop, 16-bit adaptive)")

/******************************************************************************/

struct RadixStep_CI2_sb
{
    string * str;
    size_t idx;
    size_t bkt_size[256];

    RadixStep_CI2_sb(string* strings, size_t n, size_t depth, uint8_t* charcache)
    {
        // cache characters
        uint8_t* cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *cc = (*s)[depth];

        // count character occurances
        memset(bkt_size, 0, sizeof(bkt_size));
        for (cc = charcache; cc != charcache + n; ++cc)
            ++bkt_size[*cc];

        // inclusive prefix sum
        size_t bkt[256];
        bkt[0] = bkt_size[0];
        size_t last_bkt_size = bkt_size[0];
        for (size_t i = 1; i < 256; ++i) {
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
bingmann_msd_CI2_sb(string* strings, size_t n, size_t depth, uint8_t* charcache)
{
    typedef RadixStep_CI2_sb RadixStep;

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.emplace(strings, n, depth, charcache);

    while (TLX_LIKELY(!radixstack.empty()))
    {
        while (TLX_LIKELY(radixstack.top().idx < 255))
        {
            RadixStep& rs = radixstack.top();
            ++rs.idx; // process the bucket rs.idx

            if (TLX_UNLIKELY(rs.bkt_size[rs.idx] == 0))
                ;
            else if (TLX_UNLIKELY(rs.bkt_size[rs.idx] < g_inssort_threshold))
            {
                inssort::inssort(rs.str, rs.bkt_size[rs.idx],
                                 depth + radixstack.size());
                rs.str += rs.bkt_size[rs.idx];
            }
            else
            {
                // have to increment first, as rs may be invalidated
                rs.str += rs.bkt_size[rs.idx];
                radixstack.emplace(
                    rs.str - rs.bkt_size[rs.idx], rs.bkt_size[rs.idx],
                    depth + radixstack.size(), charcache);
            }
        }
        radixstack.pop();
    }
}

static inline void
bingmann_msd_CI2_sb(string* strings, size_t n)
{
    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    uint8_t* charcache = new uint8_t[n];

    bingmann_msd_CI2_sb(strings, n, /* depth */ 0, charcache);

    delete[] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI2_sb, "bingmann/msd_CI2_sb",
               "bingmann/msd_CI2_sb (CI stack-based, charcache, fissioned)")

/******************************************************************************/

struct RadixStep_CI3_sb
{
    static const size_t RADIX = 0x10000;

    string              * str;
    size_t              idx;
    size_t              bkt_size[RADIX];

    RadixStep_CI3_sb(string* strings, size_t n, size_t depth, uint16_t* charcache)
    {
        // read characters and count character occurrences
        memset(bkt_size, 0, sizeof(bkt_size));
        uint16_t* cc = charcache;
        for (string* s = strings; s != strings + n; ++s, ++cc)
            *cc = stringtools::get_char<uint16_t>(*s, depth);
        for (cc = charcache; cc != charcache + n; ++cc)
            ++bkt_size[*cc];

        // inclusive prefix sum
        size_t bkt[RADIX];
        bkt[0] = bkt_size[0];
        size_t last_bkt_size = bkt_size[0];
        for (size_t i = 1; i < RADIX; ++i) {
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

        str = strings + bkt_size[0];
        idx = 0; // will increment to 1 on first process, bkt 0 is not sorted further
    }
};

static inline void
bingmann_msd_CI3_sb(string* strings, size_t n)
{
    static const size_t RADIX = 0x10000;

    if (n < g_inssort_threshold)
        return inssort::inssort(strings, n, 0);

    if (n < RADIX)
        return bingmann_msd_CI2_sb(strings, n);

    typedef RadixStep_CI3_sb RadixStep;

    uint16_t* charcache = new uint16_t[n];

    std::stack<RadixStep, std::vector<RadixStep> > radixstack;
    radixstack.emplace(strings, n, 0, charcache);

    while (TLX_LIKELY(!radixstack.empty()))
    {
        while (TLX_LIKELY(radixstack.top().idx < RADIX))
        {
            RadixStep& rs = radixstack.top();
            ++rs.idx;                                      // process the bucket rs.idx

            if (TLX_UNLIKELY(rs.bkt_size[rs.idx] == 0))
                ;
            else if (TLX_UNLIKELY((rs.idx & 0xFF) == 0)) { // zero-termination
                rs.str += rs.bkt_size[rs.idx];
            }
            else if (TLX_UNLIKELY(rs.bkt_size[rs.idx] < g_inssort_threshold))
            {
                inssort::inssort(rs.str, rs.bkt_size[rs.idx], 2* radixstack.size());
                rs.str += rs.bkt_size[rs.idx];
            }
            else if (rs.bkt_size[rs.idx] < RADIX)
            {
                bingmann_msd_CI2_sb(rs.str, rs.bkt_size[rs.idx],
                                    /* depth */ 2 * radixstack.size(),
                                    reinterpret_cast<uint8_t*>(charcache));
                rs.str += rs.bkt_size[rs.idx];
            }
            else
            {
                rs.str += rs.bkt_size[rs.idx];
                radixstack.emplace(
                    rs.str - rs.bkt_size[rs.idx], rs.bkt_size[rs.idx],
                    2 * radixstack.size(), charcache);
                // cannot add here, because rs may have invalidated
            }
        }
        radixstack.pop();
    }

    delete[] charcache;
}

PSS_CONTESTANT(bingmann_msd_CI3_sb, "bingmann/msd_CI3_sb",
               "bingmann/msd_CI3_sb (CI stack-based, charcache, fissioned, 16-bit adaptive)")

/******************************************************************************/

} // namespace bingmann

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_RADIX_SORT_HEADER

/******************************************************************************/
