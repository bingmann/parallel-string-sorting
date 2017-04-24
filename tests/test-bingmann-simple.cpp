/*******************************************************************************
 * tests/test-bingmann-simple.cpp
 *
 * String sorting test program
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

#include <sequential/bingmann-sample_sort.hpp>
#include <tools/stringset.hpp>
#include <tools/lcgrandom.hpp>

template <typename Iterator>
void fill_random(LCGRandom& rng, const std::string& letters,
                 Iterator begin, Iterator end)
{
    for (Iterator i = begin; i != end; ++i)
        *i = letters[(rng() / 100) % letters.size()];
}

void TestUCharString(const char* name,
                     void (* algo)(unsigned char** ss, size_t size),
                     const size_t nstrings, const size_t nchars,
                     const std::string& letters)
{
    typedef unsigned char* string;

    LCGRandom rng(1234567);

    std::cout << "Running " << name
              << " on " << nstrings << " uchar* strings" << std::endl;

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
    algo(cstrings, nstrings);

    // check result
    if (!parallel_string_sorting::UCharStringSet(
            cstrings, cstrings + nstrings).check_order()) {
        std::cout << "Result is not sorted!" << std::endl;
        abort();
    }

    // free memory.
    for (size_t i = 0; i < nstrings; ++i)
        delete[] cstrings[i];

    delete[] cstrings;
}

static const char* letters_alnum
    = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

// use macro because one cannot pass template functions as template parameters:
#define run_tests(func) \
    TestUCharString(#func, func, nstrings, 16, letters_alnum);

void test_all(const size_t nstrings)
{
    run_tests(bingmann_sample_sort::bingmann_sample_sortBSC);

    run_tests(bingmann_sample_sort::bingmann_sample_sortBTC);
    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCA);
    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCU);
    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCU1);
    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCU2);
    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCU4);

    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCT);
    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCTU);

    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCE);
    run_tests(bingmann_sample_sort::bingmann_sample_sortBTCEA);
}

int main()
{
    test_all(16);
    test_all(256);
    test_all(65550);
    test_all(1024 * 1024);
    test_all(16 * 1024 * 1024);

    return 0;
}

/******************************************************************************/
