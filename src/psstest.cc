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

static const char* memprofile_path = "memprofile.txt";
#endif

// *** Global Input Data Structures ***

size_t          gopt_inputsize = 0;
size_t          gopt_inputsize_minlimit = 0;
size_t          gopt_inputsize_maxlimit = 0;
size_t          gopt_repeats = 0;
std::vector<const char*> gopt_algorithm;

const char*     g_dataname = NULL;
const char*     g_string_data = NULL;
size_t          g_string_datasize = 0;
std::vector<size_t> g_string_offsets;
size_t          g_string_dprefix = 0;

static StatsCache g_statscache;

// file name of statistics output
static const char* statsfile = "pss-runs1.txt";

size_t g_smallsort = 64;

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
#include "rantala/multikey_multipivot.h"
#include "rantala/multikey_simd.h"

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

static inline bool gopt_algorithm_match(const char* funcname)
{
    if (!gopt_algorithm.size()) return true;

    // iterate over gopt_algorithm list as a filter
    for (size_t ai = 0; ai < gopt_algorithm.size(); ++ai) {
        if (strstr(funcname, gopt_algorithm[ai]) != NULL)
            return true;
    }

    return false;
}

void Contest::run_contest(const char* path)
{
    // read input datafile
    if (!input::load(path))
        return;

    g_dataname = path;
    g_string_dprefix = 0;

    std::cerr << "Sorting " << g_string_offsets.size() << " strings composed of " << g_string_datasize << " bytes." << std::endl;

    std::sort(m_list.begin(), m_list.end(), sort_contestants);

    // iterate over all contestants
    for (list_type::iterator c = m_list.begin(); c != m_list.end(); ++c)
    {
        if (gopt_algorithm_match( (*c)->m_funcname )) {
            (*c)->run();
        }
    }
}

void Contest::list_contentants()
{
    std::cout << "Available string sorting algorithms:" << std::endl;

    std::sort(m_list.begin(), m_list.end(), sort_contestants);

    size_t w_funcname = 0;
    for (list_type::iterator c = m_list.begin(); c != m_list.end(); ++c)
    {
        if (!gopt_algorithm_match( (*c)->m_funcname )) continue;
        w_funcname = std::max(w_funcname, strlen( (*c)->m_funcname ));
    }

    // iterate over all contestants
    for (list_type::iterator c = m_list.begin(); c != m_list.end(); ++c)
    {
        if (!gopt_algorithm_match( (*c)->m_funcname )) continue;
        std::cout << std::left << std::setw(w_funcname) << (*c)->m_funcname << "  " << (*c)->m_description << std::endl;
    }
}

void Contestant_UCArray::run()
{
    size_t repeats = gopt_repeats ? gopt_repeats : 1;

    // create unsigned char* array from offsets
    std::vector<unsigned char*> stringptr;
    stringptr.reserve( g_string_offsets.size() );

    for (size_t i = 0; i < g_string_offsets.size(); ++i)
    {
        stringptr.push_back( (unsigned char*)g_string_data + g_string_offsets[i] );
    }

    // save permutation check evaluation result
    PermutationCheck pc(stringptr);

    std::cerr << "Running " << m_funcname << " - " << m_description << "\n";

    g_statscache >> "algo" << m_funcname
                 >> "data" << g_dataname
                 >> "char_count" << g_string_datasize
                 >> "string_count" << stringptr.size();

    g_statscache >> "smallsort" << g_smallsort;

    if (repeats > 1)
        g_statscache >> "repeats" << repeats;

    //(std::cerr << m_funcname << "\t").flush();

#ifdef MALLOC_COUNT
    //MemProfile memprofile( m_funcname, memprofile_path );
    size_t memuse = malloc_count_current();
    void* stack = stack_count_clear();
    malloc_count_reset_peak();
#endif

    double ts1, ts2;

    ts1 = omp_get_wtime();
    do
    {
        m_func(stringptr.data(), stringptr.size());

        if (repeats > 1) { // refill stringptr array for next repeat
            for (size_t i = 0; i < g_string_offsets.size(); ++i) {
                stringptr[i] = (unsigned char*)g_string_data + g_string_offsets[i];
            }
        }
    } while (--repeats);
    ts2 = omp_get_wtime();

#ifdef MALLOC_COUNT
    std::cerr << "Max stack usage: " << stack_count_usage(stack) << "\n";
    std::cerr << "Max heap usage: " << malloc_count_peak() - memuse << "\n";
    //memprofile.finish();

    g_statscache >> "heapuse" << (malloc_count_peak() - memuse)
                 >> "stackuse" << stack_count_usage(stack);

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

        if (!g_string_dprefix)
            g_string_dprefix = calc_distinguishing_prefix(stringptr, g_string_datasize);

        g_statscache >> "dprefix" << g_string_dprefix;
    }
    else {
        g_statscache >> "status" << "failed";
    }

    // print timing data out to results file
    StatsWriter(statsfile).append_statsmap(g_statscache);
    g_statscache.clear();
}

void Contestant_UCArray_Parallel::run()
{
    int p = 1;

    while (1)
    {
        pss_num_threads = p;
        std::cerr << "threads=" << p << " ";
        omp_set_num_threads(p);
        g_statscache >> "threads" << p;

        Contestant_UCArray::run();

        if (p == omp_get_num_procs()) break;
        p = std::min( omp_get_num_procs(), 2 * p );
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
        { "maxsize", required_argument,  0, 'S' },
        { "algo",    required_argument,  0, 'a' },
        { "repeat",  required_argument,  0, 'r' },
        { 0,0,0,0 },
    };

#ifdef MALLOC_COUNT
    if (truncate(memprofile_path, 0)) {
        perror("Cannot truncate memprofile datafile");
    }
#endif

    while (1)
    {
        int index;
        int argi = getopt_long(argc, argv, "hs:S:a:r:", longopts, &index);

        if (argi < 0) break;

        switch (argi)
        {
        case 'h':
            print_usage(argv[0]);
            break;

        case 'a':
            if (strcmp(optarg,"list") == 0)
            {
                getContestSingleton()->list_contentants();
                return 0;
            }
            gopt_algorithm.push_back(optarg);
            std::cerr << "Selecting algorithms containing " << optarg << std::endl;
            break;

        case 's':
            if (!input::parse_filesize(optarg, gopt_inputsize_minlimit)) {
                std::cerr << "Invalid size parameter: " << optarg << std::endl;
                exit(EXIT_FAILURE);
            }
            std::cerr << "Limiting input size to " << gopt_inputsize_minlimit << std::endl;
            break;

        case 'S':
            if (!input::parse_filesize(optarg, gopt_inputsize_maxlimit)) {
                std::cerr << "Invalid maxsize parameter: " << optarg << std::endl;
                exit(EXIT_FAILURE);
            }
            std::cerr << "Limiting maximum input size to " << gopt_inputsize_maxlimit << std::endl;
            break;

        case 'r':
            gopt_repeats = atoi(optarg);
            std::cerr << "Repeating string sorting algorithms " << gopt_repeats << std::endl;
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

    if (gopt_inputsize_maxlimit < gopt_inputsize_minlimit)
        gopt_inputsize_maxlimit = gopt_inputsize_minlimit;

    for (; optind < argc; ++optind)
    {
        // iterate over input size range
        for (gopt_inputsize = gopt_inputsize_minlimit; gopt_inputsize <= gopt_inputsize_maxlimit;
             gopt_inputsize += 16)
        {
            // iterate over small sort size
            //for (g_smallsort = 64; g_smallsort < 512; g_smallsort += 8)
            {
                getContestSingleton()->run_contest(argv[optind]);
            }
        }
    }

    if (g_string_data)
        free((void*)g_string_data);

    return 0;
}
