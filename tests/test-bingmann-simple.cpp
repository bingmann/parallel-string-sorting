/*******************************************************************************
 * tests/test-bingmann-simple.cpp
 *
 * Parallel string sorting test program
 *
 *******************************************************************************
 * Copyright (C) 2015 Timo Bingmann <tb@panthema.net>
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

#include <sequential/inssort.hpp>
#include <sequential/bingmann-lcp_inssort.hpp>
#include <sequential/bingmann-radix_sort.hpp>
#include <parallel/bingmann-parallel_mkqs.hpp>
#include <tools/stringset.hpp>
#include <tools/lcgrandom.hpp>

using namespace parallel_string_sorting;

template <typename Iterator>
void fill_random(LCGRandom& rng, const std::string& letters,
                 Iterator begin, Iterator end)
{
    for (Iterator i = begin; i != end; ++i)
        *i = letters[(rng() / 100) % letters.size()];
}

template <void(* Algo)(const UCharStringSet& ss, size_t depth)>
void TestUCharString(const size_t nstrings, const size_t nchars,
                     const std::string& letters)
{
    typedef unsigned char* string;

    LCGRandom rng(1234567);

    // array of string pointers
    string* cstrings = new string[nstrings];

    // generate random strings of length nchars
    for (size_t i = 0; i < nstrings; ++i)
    {
        size_t slen = nchars + (rng() >> 8) % (nchars / 4);

        cstrings[i] = new unsigned char[slen + 1];
        fill_random(rng, letters, cstrings[i], cstrings[i] + slen);
        cstrings[i][slen] = 0;
    }

    // run sorting algorithm
    UCharStringSet ss(cstrings, cstrings + nstrings);
    Algo(ss, 0);
    if (0) ss.print();

    // check result
    if (!ss.check_order())
        abort();

    // free memory.
    for (size_t i = 0; i < nstrings; ++i)
        delete[] cstrings[i];

    delete[] cstrings;
}

template <void(* Algo)(const VectorStringSet& ss, size_t depth)>
void TestVectorString(const size_t nstrings, const size_t nchars,
                      const std::string& letters)
{
    LCGRandom rng(1234567);

    // vector of std::string objects
    std::vector<std::string> strings(nstrings);

    // generate random strings of length nchars
    for (size_t i = 0; i < nstrings; ++i)
    {
        size_t slen = nchars + (rng() >> 8) % (nchars / 4);

        strings[i].resize(slen);
        fill_random(rng, letters, strings[i].begin(), strings[i].end());
    }

    // run sorting algorithm
    VectorStringSet ss(strings.begin(), strings.end());
    Algo(ss, 0);
    //if (0) ss.print();

    // check result
    if (!ss.check_order())
        abort();
}

template <void(* Algo)(const VectorPtrStringSet& ss, size_t depth)>
void TestVectorPtrString(const size_t nstrings, const size_t nchars,
                         const std::string& letters)
{
    LCGRandom rng(1234567);

    // vector of pointers to std::string objects
    typedef std::unique_ptr<std::string> unique_ptr_string;
    std::vector<unique_ptr_string> strings(nstrings);

    // generate random strings of length nchars
    for (size_t i = 0; i < nstrings; ++i)
    {
        size_t slen = nchars + (rng() >> 8) % (nchars / 4);

        strings[i] = unique_ptr_string(new std::string(slen, 0));
        fill_random(rng, letters, strings[i]->begin(), strings[i]->end());
    }

    // run sorting algorithm
    VectorPtrStringSet ss(strings.begin(), strings.end());
    Algo(ss, 0);
    //ss.print();

    // check result
    if (!ss.check_order())
        abort();
}

template <void(* Algo)(const StringSuffixSet& ss, size_t depth)>
void TestSuffixString(const size_t nchars, const std::string& letters)
{
    LCGRandom rng(1234567);

    // std::string text object
    std::string text(nchars, 0);
    fill_random(rng, letters, text.begin(), text.end());

    std::vector<size_t> suffixarray;
    StringSuffixSet ss = StringSuffixSet::Initialize(text, suffixarray);

    // run sorting algorithm
    Algo(ss, 0);
    if (0) ss.print();

    // check result
    if (!ss.check_order())
        abort();
}

const std::string letters_alnum
    = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

// use macro because one cannot pass template functions as template parameters:
#define run_tests(func)                                  \
    TestUCharString<func>(nstrings, 16, letters_alnum);  \
    TestVectorString<func>(nstrings, 16, letters_alnum); \
    TestSuffixString<func>(nstrings, letters_alnum);     \
    TestVectorPtrString<func>(nstrings, 16, letters_alnum);

void test_all(const size_t nstrings)
{
    if (nstrings <= 1024) {
        run_tests(inssort::inssort_generic);
        run_tests(bingmann_lcp_inssort::lcp_insertion_sort_verify);
    }
    run_tests(bingmann_radix_sort::msd_CE_generic);
    run_tests(bingmann_radix_sort::msd_CE2_generic);
    run_tests(bingmann_radix_sort::msd_CI5_generic);
    run_tests(bingmann_parallel_mkqs::bingmann_sequential_mkqs_cache8);
    run_tests(bingmann_parallel_mkqs::bingmann_parallel_mkqs);
}

int main()
{
    test_all(16);
    test_all(256);
    test_all(65550);

    return 0;
}

/******************************************************************************/
