/*******************************************************************************
 * src/parallel/eberle-parallel-lcp-merge.hpp
 *
 * Parallel LCP aware merge implementation.
 *
 *******************************************************************************
 * Copyright (C) 2013-2014 Andreas Eberle <email@andreas-eberle.com>
 * Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
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

#ifndef PSS_SRC_PARALLEL_EBERLE_PARALLEL_LCP_MERGE_HEADER
#define PSS_SRC_PARALLEL_EBERLE_PARALLEL_LCP_MERGE_HEADER

#include <limits>
#include <utility>

#include "../tools/eberle-lcp-losertree.hpp"
#include "../sequential/bingmann-lcp_mergesort.hpp"

#include "../tools/jobqueue.hpp"
#include "../tools/stringtools.hpp"

#include "../tools/debug.hpp"
#undef DBGX
#define DBGX DBGX_OMP

namespace eberle_parallel_lcp_merge {

using std::numeric_limits;

using namespace eberle_lcp_utils;

using namespace jobqueue;
using namespace stringtools;

// due to ambigous symbol
using stringtools::string;

// debugging constants
static const bool debug_jobtype_on_creation = false;
static const bool debug_merge_start_message = true;

// global variables
static string* g_outputBase;  // used for debugging as base pointer to create delta output
static size_t g_lengthOfLongestJob = 0;

static unsigned g_splittingsExecuted;
static unsigned g_mergeJobsCreated;
static double g_splittingTime;

// constants
static const bool USE_WORK_SHARING = true;
static const size_t MERGE_BULK_SIZE = 4 * 1024;
static const size_t SHARE_WORK_THRESHOLD = 4 * MERGE_BULK_SIZE;

//structs defining the jobs
struct CopyDataJob : public Job
{
    LcpCacheStringPtr input;
    string            * output;

    CopyDataJob(const LcpCacheStringPtr& input, string* output)
        : input(input), output(output)
    {
        g_mergeJobsCreated++;
        DBG(debug_jobtype_on_creation,
            "CopyDataJob (output: " << (output - g_outputBase) << ", length: " << input.size << ")");
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        (void)jobQueue;

        input.copyStringsTo(output, input.size);

        return true;
    }
};

struct BinaryMergeJob : public Job
{
    LcpCacheStringPtr input1;
    LcpCacheStringPtr input2;
    lcp_t             firstLcp;
    string            * output;

    BinaryMergeJob(const LcpCacheStringPtr& input1, const LcpCacheStringPtr& input2, lcp_t firstLcp, string* output)
        : input1(input1), input2(input2), firstLcp(firstLcp), output(output)
    {
        g_mergeJobsCreated++;
        DBG(debug_jobtype_on_creation,
            "BinaryMergeJob (length1: " << input1.size << ", length2: " << input2.size << ", output: " << (output - g_outputBase) << ")");
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        (void)jobQueue;
        assert(!input1.empty() && !input2.empty());

        input1.firstLcp() = firstLcp;
        input2.firstLcp() = firstLcp;

        input1.firstCached() = input1.firstString()[firstLcp];
        input2.firstCached() = input2.firstString()[firstLcp];

        bingmann_lcp_mergesort::lcp_merge_binary_cache_opt(input1, input2, output);

        return true;
    }
};

} // namespace eberle_parallel_lcp_merge

#endif // !PSS_SRC_PARALLEL_EBERLE_PARALLEL_LCP_MERGE_HEADER

/******************************************************************************/
