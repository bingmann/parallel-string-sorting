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
#include <sys/mman.h>

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

const char*     g_datapath = NULL;      // path of input
std::string     g_dataname;             // stripped name of input
const char*     g_string_data = NULL;   // pointer to first string
char*           g_string_databuff = NULL; // pointer to data buffer (with padding)
size_t          g_string_datasize = 0;  // total size of string data (without padding)
size_t          g_string_count = 0;     // number of strings.
size_t          g_string_dprefix = 0;   // calculated distinguishing prefix

const char*     gopt_inputwrite = NULL; // argument -i, --input
const char*     gopt_output = NULL; // argument -o, --output

bool            gopt_suffixsort = false; // argument --suffix
bool            gopt_threads = false; // argument --threads
bool            gopt_allthreads = false; // argument --allthreads
bool            gopt_nocheck = false;   // argument --nocheck

StatsCache      g_statscache;

// file name of statistics output
static const char* statsfile = "pss-runs1.txt";

size_t          g_smallsort = 0;

bool            gopt_forkrun = false;
bool            gopt_forkdataload = false;

bool            gopt_sequential_only = false; // argument --sequential

const size_t    g_stacklimit = 64*1024*1024; // increase from 8 MiB

// *** Tools and Algorithms

#include "tools/debug.h"
#include "tools/membuffer.h"
#include "tools/lcgrandom.h"
#include "tools/contest.h"
#include "tools/input.h"
#include "tools/checker.h"
#include "tools/stringtools.h"
#include "tools/logfloor.h"

#include "sequential/inssort.h"
#include "sequential/mbm-radix.h"
#include "sequential/bs-mkqs.h"
#include "sequential/sinha-burstsortA.h"
#include "sequential/sinha-burstsortL.h"
#include "sequential/ng-cradix.h"
#include "sequential/ng-cradix-rantala.h"
#include "sequential/ng-lcpmergesort.h"

#include "sequential/bingmann-radix_sort.h"
#include "sequential/bingmann-sample_sort.h"
#include "sequential/bingmann-sample_sortBS.h"
#include "sequential/bingmann-sample_sortBSC.h"
#include "sequential/bingmann-sample_sortBT.h"
#include "sequential/bingmann-sample_sortBTC.h"
#include "sequential/bingmann-sample_sortBTCE.h"
#include "sequential/bingmann-sample_sortRBTCE.h"

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
#include "rantala/mergesort_lcp.h"

#undef debug

int pss_num_threads = 0;

Contest* getContestSingleton()
{
    static Contest* c = NULL;
    if (!c) c = new Contest;
    return c;
}

static inline bool gopt_algorithm_match(const char* algoname)
{
    if (!gopt_algorithm.size()) return true;

    // iterate over gopt_algorithm list as a filter
    for (size_t ai = 0; ai < gopt_algorithm.size(); ++ai) {
        if (strstr(algoname, gopt_algorithm[ai]) != NULL)
            return true;
    }

    return false;
}

static inline void maybe_inputwrite()
{
    if (!gopt_inputwrite) return;

    std::cout << "Writing unsorted input to " << gopt_inputwrite << std::endl;
    std::ofstream f(gopt_inputwrite);

    if (!gopt_suffixsort)
    {
        const char* sbegin = g_string_data;
        for (size_t i = 0; i < g_string_datasize; ++i)
        {
            if (g_string_data[i] == 0) {
                f << sbegin << "\n";
                sbegin = g_string_data + i+1;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < g_string_datasize; ++i)
            f << g_string_data + i << "\n";
    }
}

void Contest::run_contest(const char* path)
{
    g_datapath = path;

    if (!gopt_forkdataload)
    {
        // read input datafile
        if (!input::load(g_datapath))
            return;

        g_string_dprefix = 0;

        maybe_inputwrite();
        std::cout << "Sorting " << g_string_count << " strings composed of " << g_string_datasize << " bytes." << std::endl;
    }

    std::sort(m_list.begin(), m_list.end(), sort_contestants);

    // iterate over all contestants
    for (list_type::iterator c = m_list.begin(); c != m_list.end(); ++c)
    {
        if (gopt_algorithm_match( (*c)->m_algoname )) {
                (*c)->run();
        }
    }
}

void Contest::list_contentants()
{
    std::cout << "Available string sorting algorithms:" << std::endl;

    std::sort(m_list.begin(), m_list.end(), sort_contestants);

    size_t w_algoname = 0;
    for (list_type::iterator c = m_list.begin(); c != m_list.end(); ++c)
    {
        if (!gopt_algorithm_match( (*c)->m_algoname )) continue;
        w_algoname = std::max(w_algoname, strlen( (*c)->m_algoname ));
    }

    // iterate over all contestants
    for (list_type::iterator c = m_list.begin(); c != m_list.end(); ++c)
    {
        if (!gopt_algorithm_match( (*c)->m_algoname )) continue;
        std::cout << std::left << std::setw(w_algoname) << (*c)->m_algoname << "  " << (*c)->m_description << std::endl;
    }
}

void Contestant_UCArray::run()
{
    if (!gopt_forkrun) return real_run();

    pid_t p = fork();
    if (p == 0)
    {
        std::cout << "fork() ------------------------------------------------------------" << std::endl;

        if (gopt_forkdataload)
        {
            // read input datafile
            if (!input::load(g_datapath))
                return;

            g_string_dprefix = 0;

            maybe_inputwrite();
            std::cout << "Sorting " << g_string_count << " strings composed of " << g_string_datasize << " bytes." << std::endl;
        }

        if (gopt_timeout) alarm(gopt_timeout); // terminate child program after gopt_timeout seconds
        real_run();

        if (g_string_databuff) free(g_string_databuff);
        exit(0);
    }

    int status = 0;
    wait(&status);

    if (WIFEXITED(status)) {
        //std::cout << "Child normally with return code " << WEXITSTATUS(status) << std::endl;
    }
    else if (WIFSIGNALED(status)) {
        std::cout << "Child terminated abnormally with signal " << WTERMSIG(status) << std::endl;

        // write out exit status information to results file

        g_statscache >> "algo" << m_algoname
                     >> "data" << g_dataname
                     >> "char_count" << gopt_inputsize;

        if (gopt_repeats > 1)
            g_statscache >> "repeats" << gopt_repeats;

        if (WTERMSIG(status) == SIGALRM)
        {
            g_statscache >> "status" << "timeout"
                         >> "timeout" << gopt_timeout;
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

        StatsWriter(statsfile).append_stats(g_statscache);
    }
    else {
        std::cout << "Child wait returned with status " << status << std::endl;

        g_statscache >> "algo" << m_algoname
                     >> "data" << g_dataname
                     >> "char_count" << g_string_datasize
                     >> "string_count" << g_string_count
                     >> "status" << "weird";

        StatsWriter(statsfile).append_stats(g_statscache);
    }

    if (gopt_output) exit(0);
    g_statscache.clear();
}

void Contestant_UCArray::real_run()
{
    typedef unsigned char* string;

    size_t repeats = gopt_repeats ? gopt_repeats : 1;

    // lock process into memory (on Linux)
    if (mlockall(MCL_CURRENT | MCL_FUTURE) != 0) {
        std::cout << "Error locking process into memory: " << strerror(errno) << std::endl;
    }
    else {
        std::cout << "Successfully locked process into memory." << std::endl;
    }

    // create unsigned char* array from offsets
    membuffer<string> stringptr( g_string_count );

    if (!gopt_suffixsort)
    {
        size_t j = 0;
        stringptr[j++] = (string)g_string_data;
        for (size_t i = 0; i < g_string_datasize; ++i)
        {
            if (g_string_data[i] == 0 && i+1 < g_string_datasize) {
                assert(j < stringptr.size());
                stringptr[j++] = (string)g_string_data + i+1;
            }
        }
        assert(j == g_string_count);
    }
    else
    {
        assert(g_string_count == g_string_datasize);
        for (size_t i = 0; i < g_string_datasize; ++i)
            stringptr[i] = (string)g_string_data + i;
    }

    membuffer<string> stringptr_copy;
    if (repeats > 1) stringptr_copy.copy(stringptr);

    // save permutation check evaluation result
    PermutationCheck pc;
    if (!gopt_nocheck) pc = PermutationCheck(stringptr);

    std::cout << "Running " << m_algoname << " - " << m_description << std::endl;

    g_statscache >> "algo" << m_algoname
                 >> "data" << g_dataname
                 >> "char_count" << g_string_datasize
                 >> "string_count" << stringptr.size();

    if (g_smallsort)
        g_statscache >> "smallsort" << g_smallsort;

    if (repeats > 1)
        g_statscache >> "repeats" << repeats;

    omp_set_num_threads(pss_num_threads);

#ifdef MALLOC_COUNT
    //MemProfile memprofile( m_algoname, memprofile_path );
    size_t memuse = malloc_count_current();
    void* stack = stack_count_clear();
    malloc_count_reset_peak();
#endif

    if (m_prepare_func)
        m_prepare_func(stringptr.data(), stringptr.size());

    MeasureTime timer;

    timer.start();
    do
    {
        m_run_func(stringptr.data(), stringptr.size());

        if (repeats > 1) // refill stringptr array for next repeat
            stringptr.copy(stringptr_copy);

    } while (--repeats);
    timer.stop();

#ifdef MALLOC_COUNT
    std::cout << "Max stack usage: " << stack_count_usage(stack) << std::endl;
    std::cout << "Max heap usage: " << malloc_count_peak() - memuse << std::endl;
    //memprofile.finish();

    g_statscache >> "heapuse" << (malloc_count_peak() - memuse)
                 >> "stackuse" << stack_count_usage(stack);

    if (memuse < malloc_count_current())
    {
        std::cout << "MEMORY LEAKED: " << (malloc_count_current() - memuse) << " B" << std::endl;
    }
#endif

    g_statscache >> "time" << timer.delta();
    (std::cout << timer.delta() << "\tchecking ").flush();

    if (!gopt_nocheck)
    {
        if (check_sorted_order(stringptr, pc)) {
            std::cout << "ok" << std::endl;
            g_statscache >> "status" << "ok";
        }
        else {
            g_statscache >> "status" << "failed";
        }

        if (!g_string_dprefix)
            g_string_dprefix = calc_distinguishing_prefix(stringptr, g_string_datasize);

        g_statscache >> "dprefix" << g_string_dprefix;
    }
    else
    {
        std::cout << "skipped" << std::endl;
    }

    // print timing data out to results file
    StatsWriter(statsfile).append_stats(g_statscache);
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
    if (gopt_sequential_only) return;

    int p = (gopt_threads ? 1 : omp_get_num_procs());

    while (1)
    {
        pss_num_threads = p;
        std::cout << "threads=" << p << std::endl;
        g_statscache >> "threads" << p;

        Contestant_UCArray::run();

        if (p == omp_get_num_procs()) break;

        if (!gopt_allthreads)
            p = std::min( omp_get_num_procs(), 2 * p );
        else
            p = std::min( omp_get_num_procs(), p+1 );
    }
}

static inline void increase_stacklimit(size_t stacklimit)
{
    struct rlimit rl;

    if (getrlimit(RLIMIT_STACK, &rl)) {
        std::cout << "Error getrlimit(RLIMIT_STACK): " << strerror(errno) << std::endl;
    }
    else if (rl.rlim_cur < stacklimit)
    {
        rl.rlim_cur = stacklimit;
        if (setrlimit(RLIMIT_STACK, &rl)) {
            std::cout << "Error increasing stack limit with setrlimit(RLIMIT_STACK): " << strerror(errno) << std::endl;
        }
        else {
            std::cout << "Successfully increased stack limit to " << stacklimit << std::endl;
        }
    }
}

void print_usage(const char* prog)
{
    std::cout << "Usage: " << prog << " [options] filename" << std::endl
              << "Options:" << std::endl
              << "  -a, --algo <match>          Run only algorithms containing this substring, can be used multile times. Try \"list\"." << std::endl
              << "  --allthreads                Run linear thread increase test from 1 to max_processors." << std::endl
              << "  -D, --datafork              Fork before running algorithm and load data within fork!" << std::endl
              << "  -i, --input <path>          Write unsorted input strings to file, usually for checking." << std::endl
              << "      --nocheck               Skip checking of sorted order and distinguishing prefix calculation." << std::endl
              << "  -o, --output <path>         Write sorted strings to output file, terminate after first algorithm run." << std::endl
              << "  -r, --repeat <num>          Repeat experiment a number of times and divide by repetition count." << std::endl
              << "  -s, --size <size>           Limit the input size to this number of characters." << std::endl
              << "  -S, --maxsize <size>        Run through powers of two for input size limit." << std::endl
              << "      --sequential            Run only sequential algorithms." << std::endl
              << "      --suffix                Suffix sort the input file." << std::endl
              << "  -T, --timeout <sec>         Abort algorithms after this timeout (default: disabled)." << std::endl
              << "  --threads                   Run tests with doubling number of threads from 1 to max_processors." << std::endl
        ;
}

int main(int argc, char* argv[])
{
    static const struct option longopts[] = {
        { "help",    no_argument,        0, 'h' },
        { "algo",    required_argument,  0, 'a' },
        { "datafork", no_argument,       0, 'D' },
        { "input",   required_argument,  0, 'i' },
        { "nocheck", no_argument,        0, 'N' },
        { "output",  required_argument,  0, 'o' },
        { "repeat",  required_argument,  0, 'r' },
        { "size",    required_argument,  0, 's' },
        { "maxsize", required_argument,  0, 'S' },
        { "timeout", required_argument,  0, 'T' },
        { "suffix",  no_argument,        0, 1 },
        { "sequential", no_argument,     0, 2 },
        { "threads", no_argument,        0, 3 },
        { "allthreads", no_argument,     0, 4 },
        { 0,0,0,0 },
    };

    {
        char hostname[128];
        gethostname(hostname, sizeof(hostname));
        std::cout << "Running parallel-string-sorting test on " << hostname << std::endl;

        std::cout << "Called as";
        for (int i = 0; i < argc; ++i)
            std::cout << " " << argv[i];
        std::cout << std::endl;
    }

#ifdef MALLOC_COUNT
    if (truncate(memprofile_path, 0)) {
        perror("Cannot truncate memprofile datafile");
    }
#endif

    while (1)
    {
        int index;
        int argi = getopt_long(argc, argv, "hs:S:a:r:o:i:T:DN", longopts, &index);

        if (argi < 0) break;

        switch (argi)
        {
        case 'h':
            print_usage(argv[0]);
            return 0;

        case 'a':
            if (strcmp(optarg,"list") == 0)
            {
                getContestSingleton()->list_contentants();
                return 0;
            }
            gopt_algorithm.push_back(optarg);
            std::cout << "Option -a: selecting algorithms containing " << optarg << std::endl;
            break;

        case 'D':
            gopt_forkrun = gopt_forkdataload = true;
            std::cout << "Option -D: forking before each algorithm run and loading data after fork." << std::endl;
            break;

        case 'i':
            gopt_inputwrite = optarg;
            std::cout << "Option -i: will write input strings to \"" << gopt_inputwrite << "\"" << std::endl;
            break;

        case 'N':
            gopt_nocheck = true;
            std::cout << "Option --nocheck: skipping checking of sorted order and distinguishing prefix calculation." << std::endl;
            break;

        case 'o':
            gopt_output = optarg;
            std::cout << "Option -o: will write output strings to \"" << gopt_output << "\"" << std::endl;
            break;

        case 'r':
            gopt_repeats = atoi(optarg);
            std::cout << "Option -r: repeating string sorting algorithms " << gopt_repeats << std::endl;
            break;

        case 's':
            if (!input::parse_filesize(optarg, gopt_inputsize_minlimit)) {
                std::cout << "Option -s: invalid size parameter: " << optarg << std::endl;
                exit(EXIT_FAILURE);
            }
            std::cout << "Option -s: limiting input size to " << gopt_inputsize_minlimit << std::endl;
            break;

        case 'S':
            if (!input::parse_filesize(optarg, gopt_inputsize_maxlimit)) {
                std::cout << "Option -S: invalid maxsize parameter: " << optarg << std::endl;
                exit(EXIT_FAILURE);
            }
            std::cout << "Option -S: limiting maximum input size to " << gopt_inputsize_maxlimit << std::endl;
            break;

        case 'T':
            gopt_timeout = atoi(optarg);
            std::cout << "Option -T: aborting algorithms after " << gopt_timeout << " seconds timeout." << std::endl;
            break;

        case 1: // --suffix
            gopt_suffixsort = true;
            std::cout << "Option --suffix: running as suffix sorter on input file." << std::endl;
            break;

        case 2: // --sequential
            gopt_sequential_only = true;
            std::cout << "Option --sequential: running only sequential algorithms." << std::endl;
            break;

        case 3: // --threads
            gopt_threads = true;
            std::cout << "Option --threads: running test with exponentially increasing thread count." << std::endl;
            break;

        case 4: // --allthreads
            gopt_allthreads = true;
            std::cout << "Option --allthreads: running test with linear increasing thread count." << std::endl;
            break;

        default:
            std::cout << "Invalid parameter: " << argi << std::endl;
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (optind == argc) { // no input data parameter given
        print_usage(argv[0]);
        return 0;
    }

    increase_stacklimit(g_stacklimit);

    std::cout << "Using CLOCK_MONOTONIC with resolution: " << MeasureTime().resolution() << std::endl;

    if (gopt_inputsize_maxlimit < gopt_inputsize_minlimit)
        gopt_inputsize_maxlimit = gopt_inputsize_minlimit;

    for (; optind < argc; ++optind)
    {
        // iterate over input size range
        for (gopt_inputsize = gopt_inputsize_minlimit; gopt_inputsize <= gopt_inputsize_maxlimit;
             gopt_inputsize *= 2)
        {
            // iterate over small sort size
            //for (g_smallsort = 1*1024*1024; g_smallsort <= 1*1024*1024; g_smallsort *= 2)
            {
                getContestSingleton()->run_contest(argv[optind]);
            }

            if (gopt_inputsize == 0) break;
        }
    }

    if (g_string_databuff) free(g_string_databuff);

    return 0;
}
