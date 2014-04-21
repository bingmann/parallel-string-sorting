/******************************************************************************
 * src/parallel/eberle-parallel-lcp-mergesort.h
 *
 * Parallel LCP aware merge sort.
 *
 ******************************************************************************
 * Copyright (C) 2014 Andreas Eberle <email@andreas-eberle.com>
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

#ifndef EBERLE_PARALLEL_LCP_MERGESORT_H_
#define EBERLE_PARALLEL_LCP_MERGESORT_H_

#include <iostream>
#include <utility>

#include "eberle-parallel-lcp-merge-lcp-splitting.h"
#include "eberle-parallel-lcp-merge-standard-splitting.h"
#include "eberle-parallel-lcp-merge-binary-splitting.h"
#include "../sequential/eberle-mergesort-lcp-losertree.h"

#include "../tools/eberle-utilities.h"

#include "../tools/jobqueue.h"
#include "../tools/stringtools.h"

#include "../tools/debug.h"
#undef DBGX
#define DBGX DBGX_OMP


namespace eberle_parallel_lcp_mergesort
{
using std::numeric_limits;

using namespace eberle_utils;

using namespace jobqueue;
using namespace stringtools;

//due to ambigous symbol
using stringtools::string;

//debugging constants
static const bool debug_toplevel_merge_duration = true;

//constants
static const unsigned MERGESORT_BRANCHES = 64;

//method definitions

static inline void
eberle_parallel_lcp_mergesort(string *strings, size_t n, void (*parallelMerge)(const LcpCacheStringPtr*, unsigned, string*, size_t))
{
    int realNumaNodes = numa_num_configured_nodes();
    if (realNumaNodes < 1) realNumaNodes = 1;
    const unsigned topLevelBranches = omp_get_max_threads();

    unsigned threadsPerNode = ceil(float(omp_get_max_threads()) / realNumaNodes);
    if (threadsPerNode < 1) threadsPerNode = 1;

    // calculate ranges
    std::pair<size_t, size_t> ranges[topLevelBranches];
    calculateRanges(ranges, topLevelBranches, n);

    LcpCacheStringPtr stringPtr[topLevelBranches];

#pragma omp parallel for
    for(unsigned k = 0; k < topLevelBranches; k++)
    {
        const size_t length = ranges[k].second;
        const unsigned numaNode = k / threadsPerNode;

        // tie thread to a NUMA node
        numa_run_on_node(numaNode);
        numa_set_preferred(numaNode);

        // allocate memory on numa node
        stringPtr[k].allocateNumaMemory(numaNode, length);

        // execute mergesort
        eberle_mergesort::eberle_mergesort_losertree_lcp_kway<MERGESORT_BRANCHES>(strings + ranges[k].first, stringPtr[k]);
    }

    // do top level merge
    parallelMerge(stringPtr, topLevelBranches, strings, n);

    // free memory
    for(unsigned k = 0; k < topLevelBranches; k++)
    {
        stringPtr[k].freeNumaMemory();
    }
}


void eberle_parallel_lcp_mergesort_lcp_splitting(string* strings, size_t n)
{
    eberle_parallel_lcp_mergesort(strings, n, eberle_parallel_lcp_merge::parallelLcpMerge);
}

void eberle_parallel_lcp_mergesort_standard_splitting(string* strings, size_t n)
{
    eberle_parallel_lcp_mergesort(strings, n, eberle_parallel_lcp_merge::parallelLcpMergeStandardSplitting);
}

void eberle_parallel_lcp_mergesort_binary_splitting(string* strings, size_t n)
{
    eberle_parallel_lcp_mergesort(strings, n, eberle_parallel_lcp_merge::parallelLcpMergeBinarySplitting);
}

CONTESTANT_REGISTER_PARALLEL(eberle_parallel_lcp_mergesort_lcp_splitting,
    "eberle/parallel-lcp-mergesort-lcp-splitting",
    "parallel LCP aware mergesort by Andreas Eberle")

CONTESTANT_REGISTER_PARALLEL(eberle_parallel_lcp_mergesort_standard_splitting,
    "eberle/parallel-lcp-mergesort-standard-splitting",
    "parallel LCP aware mergesort by Andreas Eberle")

CONTESTANT_REGISTER_PARALLEL(eberle_parallel_lcp_mergesort_binary_splitting,
    "eberle/parallel-lcp-mergesort-binary-splitting",
    "parallel LCP aware mergesort by Andreas Eberle")


} // namespace eberle_parallel_lcp_mergesort

#endif // EBERLE_PARALLEL_LCP_MERGESORT_H_
