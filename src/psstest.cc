/******************************************************************************
 * src/psstest.cc
 *
 * Parallel string sorting test program
 *
 ******************************************************************************
 * Copyright (C) 2012 Timo Bingmann <tb@panthema.net>
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
 *****************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <string.h>
#include <errno.h>

#include <string>
#include <bitset>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <omp.h>
#include <getopt.h>

#include <boost/static_assert.hpp>
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/integral_constant.hpp>

#include "zio.h"

// *** Global Input Data Structures ***

size_t          gopt_inputsizelimit = 0;
const char*     gopt_algorithm = NULL;

const char*     g_stringdata = NULL;
size_t          g_stringdatasize = 0;
std::vector<size_t> g_stringoffsets;

// *** Tools and Algorithms

#include "tools/contest.h"
#include "tools/input.h"
#include "tools/checker.h"

#include "sequential/inssort.h"
#include "sequential/mkqs.h"
#include "sequential/burstsortA.h"
#include "sequential/burstsortL.h"
#include "sequential/cradix.h"
#include "sequential/cradix-rantala.h"

#include "rantala/debug.h"
#include "rantala/get_char.h"
#include "rantala/insertion_sort.h"
#include "rantala/msd_ce.h"
#include "rantala/msd_ci.h"
#include "rantala/vector_malloc.h"
#include "rantala/vector_realloc.h"
#include "rantala/vector_block.h"
#include "rantala/vector_bagwell.h"
#include "rantala/vector_brodnik.h"
#include "rantala/burstsort.h"
#include "rantala/burstsort2.h"
#include "rantala/mergesort.h"
#include "rantala/mergesort_unstable.h"
#include "rantala/losertree.h"
#include "rantala/mergesort_losertree.h"
//#include "rantala/mergesort_lcp.h"
#include "rantala/funnelsort.h"
#undef debug

void Contest::run_contest(const char* path)
{
    // read input datafile
    readfile_lines(path);

    // iterate over all contestants
    for (list_type::iterator c = m_list.begin(); c != m_list.end(); ++c)
    {
        if (gopt_algorithm)
        {
            if ((*c)->m_funcname.find(gopt_algorithm) == std::string::npos)
                continue;
        }

        (*c)->run();
    }
}

void Contestant_UCArray::run()
{
    // create unsigned char* array from offsets
    std::vector<unsigned char*> stringptr;

    for (size_t i = 0; i < g_stringoffsets.size(); ++i)
    {
        stringptr.push_back( (unsigned char*)g_stringdata + g_stringoffsets[i] );
    }

    // save permutation check evaluation result
    PermutationCheck pc(stringptr);

    std::cout << "Running " << m_funcname << " - " << m_description << "\n";

    double ts1, ts2;

    ts1 = omp_get_wtime();
    m_func(stringptr.data(), stringptr.size());
    ts2 = omp_get_wtime();

    (std::cout << ts2-ts1 << " checking ").flush();

    if (check_sorted_order(stringptr, pc))
        std::cout << "ok" << std::endl;
}

void print_usage(const char* prog)
{
    std::cerr << "Usage: " << prog << " [-s <input size limit>] filename" << std::endl;
}

int main(int argc, char* argv[])
{
    static const struct option longopts[] = {
        { "help",    no_argument,        0, 'h' },
        { "size",    required_argument,  0, 's' },
        { "algo",    required_argument,  0, 'a' },
        { 0,0,0,0 },
    };

    while (1)
    {
        int index;
        int argi = getopt_long(argc, argv, "hs:a:", longopts, &index);

        if (argi < 0) break;

        switch (argi)
        {
        case 'h':
            print_usage(argv[0]);
            break;

        case 'a':
            gopt_algorithm = optarg;
            std::cerr << "Selecting algorithms containing " << optarg << std::endl;
            break;

        case 's':
            if (!parse_filesize(optarg, gopt_inputsizelimit)) {
                std::cerr << "Invalid size parameter: " << optarg << std::endl;
                exit(EXIT_FAILURE);                
            }
            std::cerr << "Limiting input size to " << gopt_inputsizelimit << std::endl;
            break;

        default:
            std::cerr << "Invalid parameter: " << argi << std::endl;
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    while (optind < argc)
    {
        getContestSingleton()->run_contest(argv[optind++]);
    }

    return 0;
}
