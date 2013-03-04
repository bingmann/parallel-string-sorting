/******************************************************************************
 * src/psstest.cc
 *
 * Parallel string sorting test program
 *
 ******************************************************************************
 * Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
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
#include <math.h>

#include <string>
#include <bitset>
#include <vector>
#include <list>
#include <deque>
#include <set>
#include <map>
#include <stack>
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

#include "src/config.h"
#include "zio.h"

#include "tools/statsfile.h"

#ifdef MALLOC_COUNT
#include "tools/malloc_count.h"
#include "tools/stack_count.h"
#include "tools/memprofile.h"
#endif

// *** Global Input Data Structures ***

size_t          gopt_inputsizelimit = 0;
std::vector<const char*> gopt_algorithm;

const char*     g_stringdata = NULL;
size_t          g_stringdatasize = 0;
std::vector<size_t> g_stringoffsets;

static StatsCache g_statscache;

// file name of statistics output
static const char* statsfile = "pss-runs1.txt";

// *** Tools and Algorithms

#include "tools/debug.h"
#include "tools/lcgrandom.h"
#include "tools/contest.h"
#include "tools/input.h"
#include "tools/checker.h"
#include "tools/stringtools.h"
#include "tools/logfloor.h"

#include "sequential/inssort.h"
#include "sequential/mbm-radix.h"
#include "sequential/mkqs.h"
#include "sequential/burstsortA.h"
#include "sequential/burstsortL.h"
#include "sequential/cradix.h"
#include "sequential/cradix-rantala.h"
#include "sequential/bingmann-radix_sort.h"
#include "sequential/bingmann-sample_sort.h"

#include "rantala/tools/debug.h"
#include "rantala/tools/get_char.h"
#include "rantala/tools/median.h"
#include "rantala/tools/vector_malloc.h"
#include "rantala/tools/vector_realloc.h"
#include "rantala/tools/vector_block.h"
#include "rantala/tools/vector_bagwell.h"
#include "rantala/tools/vector_brodnik.h"
#include "rantala/tools/insertion_sort.h"

#include "rantala/multikey_block.h"
#include "rantala/multikey_cache.h"
#include "rantala/multikey_dynamic.h"
//#include "rantala/multikey_multipivot.h"
//#include "rantala/multikey_simd.h"

#include "rantala/msd_a.h"
#include "rantala/msd_a2.h"
#include "rantala/msd_ce.h"
#include "rantala/msd_ci.h"
#include "rantala/msd_dyn_block.h"
#include "rantala/msd_dyn_vector.h"

#include "rantala/burstsort.h"
#include "rantala/burstsort2.h"
#include "rantala/burstsort_mkq.h"
#include "rantala/mergesort.h"
#include "rantala/mergesort_unstable.h"
#include "rantala/tools/losertree.h"
#include "rantala/mergesort_losertree.h"
//#include "rantala/mergesort_lcp.h"

#undef debug

int pss_num_threads = 0;

Contest* getContestSingleton()
{
    static Contest* c = NULL;
    if (!c) c = new Contest;
    return c;
}

void Contest::run_contest(const char* path)
{
    // read input datafile
    if (!input::load(path))
        return;

    std::cerr << "Sorting " << g_stringoffsets.size() << " strings composed of " << g_stringdatasize << " bytes." << std::endl;

    std::sort(m_list.begin(), m_list.end(), sort_contestants);

    // iterate over all contestants
    for (list_type::iterator c = m_list.begin(); c != m_list.end(); ++c)
    {
        if (!gopt_algorithm.size()) {
            (*c)->run();
        }
        else {
            // iterate over gopt_algorithm list as a filter
            for (size_t ai = 0; ai < gopt_algorithm.size(); ++ai) {
                if (strstr((*c)->m_funcname, gopt_algorithm[ai]) != NULL) {
                    (*c)->run();
                    break;
                }
            }
        }
    }
}

void Contestant_UCArray::run()
{
    // create unsigned char* array from offsets
    std::vector<unsigned char*> stringptr;
    stringptr.reserve( g_stringoffsets.size() );

    for (size_t i = 0; i < g_stringoffsets.size(); ++i)
    {
        stringptr.push_back( (unsigned char*)g_stringdata + g_stringoffsets[i] );
    }

    // save permutation check evaluation result
    PermutationCheck pc(stringptr);

    std::cerr << "Running " << m_funcname << " - " << m_description << "\n";

    g_statscache >> "algo" << m_funcname
                 >> "char_count" << g_stringdatasize
                 >> "string_count" << stringptr.size();

    //(std::cerr << m_funcname << "\t").flush();

#ifdef MALLOC_COUNT
    MemProfile memprofile( m_funcname );
    size_t memuse = malloc_count_current();
    //void* stack = stack_count_clear();
#endif

    double ts1, ts2;

    ts1 = omp_get_wtime();
    m_func(stringptr.data(), stringptr.size());
    ts2 = omp_get_wtime();

#ifdef MALLOC_COUNT
    //stack_count_usage(stack);
    memprofile.finish();

    if (memuse < malloc_count_current())
    {
        std::cerr << "MEMORY LEAKED: " << (malloc_count_current() - memuse) << " B\n";
    }
#endif

    g_statscache >> "time" << (ts2-ts1);
    (std::cerr << ts2-ts1 << "\tchecking ").flush();

    if (check_sorted_order(stringptr, pc)) {
        std::cerr << "ok" << std::endl;
        g_statscache >> "status" << "ok";
    }
    else {
        g_statscache >> "status" << "failed";
    }

    // print timing data out to results file
    //StatsWriter(statsfile).append_statsmap(g_statscache);
    g_statscache.clear();
}

void Contestant_UCArray_Parallel::run()
{
    for (int p = 1; p <= omp_get_num_procs(); p *= 2)
    {
        pss_num_threads = p;
        std::cerr << "threads=" << p << " ";
        omp_set_num_threads(p);
        g_statscache >> "threads" << p;

        Contestant_UCArray::run();
    }
}

void print_usage(const char* prog)
{
    std::cerr << "Usage: " << prog << " [-s <input size limit>] [-a <algorithm match>] filename" << std::endl;
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
            gopt_algorithm.push_back(optarg);
            std::cerr << "Selecting algorithms containing " << optarg << std::endl;
            break;

        case 's':
            if (!input::parse_filesize(optarg, gopt_inputsizelimit)) {
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

    if (optind == argc) { // no input data parameter given
        print_usage(argv[0]);
    }

    while (optind < argc)
    {
        getContestSingleton()->run_contest(argv[optind++]);
    }

    if (g_stringdata)
        free((void*)g_stringdata);

    return 0;
}
