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
#include <sys/wait.h>

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
#include <gmpxx.h>

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
int             gopt_timeout = 0;

const char*     g_dataname = NULL;
const char*     g_string_data = NULL;
size_t          g_string_datasize = 0;
char*           g_string_databuff = NULL;
std::vector<size_t> g_string_offsets;
size_t          g_string_dprefix = 0;

const char*     gopt_output = NULL;

bool            gopt_suffixsort = false;

static StatsCache g_statscache;

// file name of statistics output
static const char* statsfile = "pss-runs1.txt";

//size_t g_smallsort = 64;

static const bool use_forkrun = true;
static const bool use_forkdataload = false;

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
#include "sequential/ng-cradix.h"
#include "sequential/ng-cradix-rantala.h"
#include "sequential/ng-lcpmergesort.h"
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
    g_dataname = path;

    if (!use_forkdataload)
    {
        // read input datafile
        if (!input::load(g_dataname))
            return;

        g_string_dprefix = 0;

        std::cerr << "Sorting " << g_string_offsets.size() << " strings composed of " << g_string_datasize << " bytes." << std::endl;
    }

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
    if (!use_forkrun) return real_run();

    pid_t p = fork();
    if (p == 0)
    {
        if (use_forkdataload)
        {
            // read input datafile
            if (!input::load(g_dataname))
                return;

            g_string_dprefix = 0;

            std::cerr << "Sorting " << g_string_offsets.size() << " strings composed of " << g_string_datasize << " bytes." << std::endl;
        }

        if (gopt_timeout) alarm(gopt_timeout); // terminate child program after use_timeout seconds
        real_run();

        if (g_string_data) free(g_string_databuff);
        exit(0);
    }

    int status;
    wait(&status);

    if (WIFEXITED(status)) {
        //std::cout << "Child normally with return code " << WEXITSTATUS(status) << std::endl;
    }
    else if (WIFSIGNALED(status)) {
        std::cout << "Child terminated abnormally with signal " << WTERMSIG(status) << std::endl;

        // write out exit status information to results file

        g_statscache >> "algo" << m_funcname
                     >> "data" << g_dataname
                     >> "char_count" << g_string_datasize
                     >> "string_count" << g_string_offsets.size();

        if (gopt_repeats > 1)
            g_statscache >> "repeats" << gopt_repeats;

        if (WTERMSIG(status) == SIGALRM)
        {
            g_statscache >> "status" << "timeout";
        }
        else if (WTERMSIG(status) == SIGSEGV)
        {
            g_statscache >> "status" << "segfault";
        }
        else if (WTERMSIG(status) == SIGABRT)
        {
            g_statscache >> "status" << "aborted";
        }
        else
        {
            g_statscache >> "status" << "SIG" << WTERMSIG(status);
        }

        StatsWriter(statsfile).append_statsmap(g_statscache);
    }
    else {
        std::cout << "Child wait returned with status " << status << std::endl;

        g_statscache >> "algo" << m_funcname
                     >> "data" << g_dataname
                     >> "char_count" << g_string_datasize
                     >> "string_count" << g_string_offsets.size()
                     >> "status" << "weird";

        StatsWriter(statsfile).append_statsmap(g_statscache);
    }

    if (gopt_output) exit(0);
    g_statscache.clear();
}

void Contestant_UCArray::real_run()
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

    //g_statscache >> "smallsort" << g_smallsort;

    if (repeats > 1)
        g_statscache >> "repeats" << repeats;

    omp_set_num_threads(pss_num_threads);

#ifdef MALLOC_COUNT
    //MemProfile memprofile( m_funcname, memprofile_path );
    size_t memuse = malloc_count_current();
    void* stack = stack_count_clear();
    malloc_count_reset_peak();
#endif

    if (m_prepare_func)
        m_prepare_func(stringptr.data(), stringptr.size());

    double ts1, ts2;

    ts1 = omp_get_wtime();
    do
    {
        m_run_func(stringptr.data(), stringptr.size());

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
    }
    else {
        g_statscache >> "status" << "failed";
    }
    
    if (!g_string_dprefix)
        g_string_dprefix = calc_distinguishing_prefix(stringptr, g_string_datasize);

    g_statscache >> "dprefix" << g_string_dprefix;

    // print timing data out to results file
    StatsWriter(statsfile).append_statsmap(g_statscache);
    g_statscache.clear();

    if (gopt_output)
    {
        std::cout << "Writing sorted output to " << gopt_output << std::endl;
        std::ofstream f(gopt_output);
        for (size_t i = 0; i < stringptr.size(); ++i)
            f << stringptr[i] << "\n";
        f.close();
        exit(0);
    }
}

void Contestant_UCArray_Parallel::run()
{
    int p = 1;

    while (1)
    {
        pss_num_threads = p;
        std::cerr << "threads=" << p << " ";
        g_statscache >> "threads" << p;

        Contestant_UCArray::run();

        if (p == omp_get_num_procs()) break;
        p = std::min( omp_get_num_procs(), 2 * p );
    }
}

void print_usage(const char* prog)
{
    std::cerr << "Usage: " << prog << " [options] filename" << std::endl
              << "Options:" << std::endl
              << "  -a, --algo <match>          Run only algorithms containing this substring, can be used multile times. Try \"list\"." << std::endl
              << "  -o, --output <path>         Write sorted strings to output file, terminate after first algorithm run." << std::endl
              << "  -r, --repeat <num>          Repeat experiment a number of times and divide by repetition count." << std::endl
              << "  -s, --size <size>           Limit the input size to this number of characters." << std::endl
              << "  -S, --maxsize <size>        Run through powers of two for input size limit." << std::endl
              << "      --suffix                Suffix sort the input file." << std::endl
              << "  -T, --timeout <sec>         Abort algorithms after this timeout (default: disabled)." << std::endl
        ;
}

int main(int argc, char* argv[])
{
    static const struct option longopts[] = {
        { "algo",    required_argument,  0, 'a' },
        { "help",    no_argument,        0, 'h' },
        { "maxsize", required_argument,  0, 'S' },
        { "output",  required_argument,  0, 'o' },
        { "repeat",  required_argument,  0, 'r' },
        { "size",    required_argument,  0, 's' },
        { "timeout", required_argument,  0, 'T' },
        { "suffix",  no_argument,        0, 1 },
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
        int argi = getopt_long(argc, argv, "hs:S:a:r:o:T:", longopts, &index);

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

        case 'o':
            gopt_output = optarg;
            std::cerr << "Will write output strings to " << gopt_output << std::endl;
            break;

        case 'T':
            gopt_timeout = atoi(optarg);
            std::cerr << "Aborting algorithms after " << gopt_timeout << " seconds timeout." << std::endl;
            break;

        case 1:
            gopt_suffixsort = true;
            std::cerr << "Running as suffix sorter on input file." << std::endl;
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
             gopt_inputsize *= 2)
        {
            // iterate over small sort size
            //for (g_smallsort = 64; g_smallsort < 512; g_smallsort += 8)
            {
                getContestSingleton()->run_contest(argv[optind]);
            }

            if (gopt_inputsize == 0) break;
        }
    }

    if (g_string_data) free(g_string_databuff);

    return 0;
}
