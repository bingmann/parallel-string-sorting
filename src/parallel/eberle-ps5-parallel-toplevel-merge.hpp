/*******************************************************************************
 * src/parallel/eberle-ps5-parallel-toplevel-merge.hpp
 *
 * NUMA aware parallel sorting algorithm running separate pS5s on every NUMA
 * node and then combining the result via parallel LCP aware merge.
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

#ifndef PSS_SRC_PARALLEL_EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_HEADER
#define PSS_SRC_PARALLEL_EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_HEADER

#include <utility>

#include "eberle-parallel-lcp-merge-lcp-splitting.hpp"
#include "eberle-parallel-lcp-merge-standard-splitting.hpp"
#include "eberle-parallel-lcp-merge-binary-splitting.hpp"
#include "bingmann-parallel_sample_sort.hpp"

#include "../tools/jobqueue.hpp"
#include "../tools/stringtools.hpp"

namespace eberle_ps5_parallel_toplevel_merge {

using std::numeric_limits;

using namespace jobqueue;
using namespace stringtools;

using namespace bingmann_parallel_sample_sort;

//due to ambigous symbol
using stringtools::string;

//debugging constants
static const bool debug_toplevel_merge_duration = true;
static const bool debug_verify_ps5_lcp_cache = false;

//method definitions
static inline
void eberle_ps5_parallel_toplevel_merge(
    string* strings, size_t n,
    void (* parallelMerge)(const LcpCacheStringPtr*, unsigned, string*, size_t))
{
    int realNumaNodes = numa_num_configured_nodes();
    if (realNumaNodes < 1) realNumaNodes = 1;

    if (realNumaNodes == 1) {
        LOG1 << "No or just one NUMA nodes detected on the system.";
        LOG1 << "pS5-LCP-Mergesort is designed for NUMA systems.";
        LOG1 << "Continuing anyway, at your own peril!";
    }

    g_stats >> "num_numa_nodes" << realNumaNodes;

    // maybe raise number of used NUMA nodes for developement
    int numNumaNodes = g_numa_nodes;
    if (numNumaNodes < 1) numNumaNodes = 1;

    if (numNumaNodes != realNumaNodes || g_numa_nodes == 0)
    {
        LOG1 << "!!! WARNING !!! emulating " << numNumaNodes << " NUMA nodes! "
             << "Remove --numa-nodes for REAL EXPERIMENTS.";
    }

    // create outputs
    LcpCacheStringPtr outputs[numNumaNodes];

    // distribute threads amoung NUMA nodes
    int numThreadsPerNode = omp_get_max_threads() / numNumaNodes;
    int remainThreads = omp_get_max_threads() % numNumaNodes;

    ClockTimer timer;

    if (numThreadsPerNode == 0)
    {
        LOG1 << "Fewer threads than NUMA nodes detected.";
        LOG1 << "Strange things will happen now.";
        LOG1 << "Continuing anyway, at your own peril!";

        // set num threads = 1 for nodes, and limit number of active nodes via
        // normal num_threads() in next loop

#pragma omp parallel for schedule(dynamic)
        for (int k = 0; k < numNumaNodes; k++)
        {
            size_t start = g_numa_strings[k];
            size_t length = g_numa_string_count[k];

            int nodeThreads = 1;
            int numaNode = k % realNumaNodes;

            ClockTimer timer;

            outputs[k].allocateNumaMemory(numaNode, length);

            parallel_sample_sort_numa(strings + start, length,
                                      numaNode, nodeThreads,
                                      outputs[k]);

            LOGC(debug_toplevel_merge_duration)
                << "node[" << k << "] took : " << timer.elapsed() << " s";
        }
    }
    else
    {
#pragma omp parallel for num_threads(numNumaNodes) schedule(static)
        for (int k = 0; k < numNumaNodes; k++)
        {
            size_t start = g_numa_strings[k];
            size_t length = g_numa_string_count[k];

            int nodeThreads = numThreadsPerNode;
            int numaNode = k % realNumaNodes;
            if (k < remainThreads) nodeThreads++; // distribute extra threads

            LOG1 << "node[" << numaNode << "] gets " << nodeThreads << " threads";

            ClockTimer timer;

            outputs[k].allocateNumaMemory(numaNode, length);

            parallel_sample_sort_numa(strings + start, length,
                                      numaNode, nodeThreads,
                                      outputs[k]);

            LOGC(debug_toplevel_merge_duration)
                << "node[" << k << "] took : " << timer.elapsed() << " s";
        }
    }

    LOGC(debug_toplevel_merge_duration)
        << "all nodes took : " << timer.elapsed() << " s";

    if (debug_verify_ps5_lcp_cache)
    {
        for (int k = 0; k < numNumaNodes; k++)
        {
            LcpCacheStringPtr& out = outputs[k];

            stringtools::verify_lcp_cache(
                out.strings, out.lcps, out.cachedChars,
                out.size, 0);
        }
    }

    // do top level merge

    timer.start();
    //eberle_parallel_lcp_merge::sequentialLcpMerge(outputs, numNumaNodes, strings, n);
    parallelMerge(outputs, numNumaNodes, strings, n);

    LOGC(debug_toplevel_merge_duration)
        << "top level merge needed: " << timer.elapsed() << " s";

    for (int k = 0; k < numNumaNodes; k++)
    {
        outputs[k].freeNumaMemory();
    }
}

static inline
void eberle_ps5_parallel_toplevel_merge_lcp_splitting(string* strings, size_t n)
{
    eberle_ps5_parallel_toplevel_merge(
        strings, n, eberle_parallel_lcp_merge::parallelLcpMerge);
}

static inline
void eberle_ps5_parallel_toplevel_merge_standard_splitting(string* strings, size_t n)
{
    eberle_ps5_parallel_toplevel_merge(
        strings, n, eberle_parallel_lcp_merge::parallelLcpMergeStandardSplitting);
}

static inline
void eberle_ps5_parallel_toplevel_merge_binary_splitting(string* strings, size_t n)
{
    eberle_ps5_parallel_toplevel_merge(
        strings, n, eberle_parallel_lcp_merge::parallelLcpMergeBinarySplitting);
}

PSS_CONTESTANT_PARALLEL(
    eberle_ps5_parallel_toplevel_merge_lcp_splitting,
    "eberle/ps5-parallel-toplevel-merge-lcp-splitting",
    "NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle")

PSS_CONTESTANT_PARALLEL(
    eberle_ps5_parallel_toplevel_merge_standard_splitting,
    "eberle/ps5-parallel-toplevel-merge-standard-splitting",
    "NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle")

PSS_CONTESTANT_PARALLEL(
    eberle_ps5_parallel_toplevel_merge_binary_splitting,
    "eberle/ps5-parallel-toplevel-merge-binary-splitting",
    "NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle")

static inline
void eberle_ps5_parallel_toplevel_merge_assisting(
    string* strings, size_t n,
    void (* parallelMerge)(const LcpCacheStringPtr*, unsigned, string*, size_t))
{
    int realNumaNodes = numa_num_configured_nodes();
    if (realNumaNodes < 1) realNumaNodes = 1;

    if (realNumaNodes == 1) {
        LOG1 << "No or just one NUMA nodes detected on the system.";
        LOG1 << "pS5-LCP-Mergesort is designed for NUMA systems.";
        LOG1 << "Continuing anyway, at your own peril!";
    }

    g_stats >> "num_numa_nodes" << realNumaNodes;

    // maybe raise number of used NUMA nodes for developement
    int numNumaNodes = g_numa_nodes;
    if (numNumaNodes < 1) numNumaNodes = 1;

    if (numNumaNodes != realNumaNodes || g_numa_nodes == 0)
    {
        LOG1 << "!!! WARNING !!! emulating " << numNumaNodes << " NUMA nodes! "
             << "Remove --numa-nodes for REAL EXPERIMENTS.";
    }

    ClockTimer timer;

    // create outputs
    LcpCacheStringPtr outputs[numNumaNodes];

    // create in/output StringPtrs
    typedef UCharStringSet StringSet;
    std::vector<UCharStringShadowLcpCacheOutPtr> inputs;

    for (int k = 0; k < numNumaNodes; k++)
    {
        size_t start = g_numa_strings[k];
        size_t length = g_numa_string_count[k];

        outputs[k].allocateNumaMemory(k % realNumaNodes, length);

        inputs.emplace_back(
            StringSet(strings + start, strings + start + length),
            StringSet((string*)outputs[k].lcps, (string*)outputs[k].lcps + length),
            StringSet(outputs[k].strings, outputs[k].strings + length),
            outputs[k].lcps, outputs[k].cachedChars);
    }

    LOGC(debug_toplevel_merge_duration)
        << "allocation needed: " << timer.elapsed() << " s";

    parallel_sample_sort_numa2(inputs.data(), numNumaNodes);

    if (debug_verify_ps5_lcp_cache)
    {
        for (int k = 0; k < numNumaNodes; k++)
        {
            LcpCacheStringPtr& out = outputs[k];

            stringtools::verify_lcp_cache(
                out.strings, out.lcps, out.cachedChars,
                out.size, 0);
        }
    }

    // do top level merge

    timer.start();
    //eberle_parallel_lcp_merge::sequentialLcpMerge(outputs, numNumaNodes, strings, n);
    parallelMerge(outputs, numNumaNodes, strings, n);

    LOGC(debug_toplevel_merge_duration)
        << "top level merge needed: " << timer.elapsed() << " s";

    for (int k = 0; k < numNumaNodes; k++)
    {
        outputs[k].freeNumaMemory();
    }
}

static inline
void eberle_ps5_parallel_toplevel_merge_assisting_lcp_splitting(
    string* strings, size_t n)
{
    eberle_ps5_parallel_toplevel_merge_assisting(
        strings, n, eberle_parallel_lcp_merge::parallelLcpMerge);
}

static inline
void eberle_ps5_parallel_toplevel_merge_assisting_standard_splitting(
    string* strings, size_t n)
{
    eberle_ps5_parallel_toplevel_merge_assisting(
        strings, n, eberle_parallel_lcp_merge::parallelLcpMergeStandardSplitting);
}

static inline
void eberle_ps5_parallel_toplevel_merge_assisting_binary_splitting(
    string* strings, size_t n)
{
    eberle_ps5_parallel_toplevel_merge_assisting(
        strings, n, eberle_parallel_lcp_merge::parallelLcpMergeBinarySplitting);
}

PSS_CONTESTANT_PARALLEL(
    eberle_ps5_parallel_toplevel_merge_assisting_lcp_splitting,
    "eberle/ps5-parallel-toplevel-merge-assisting-lcp-splitting",
    "pS5-LCP-Merge with JobQueue assisting each other by Andreas Eberle and Timo Bingmann")

PSS_CONTESTANT_PARALLEL(
    eberle_ps5_parallel_toplevel_merge_assisting_standard_splitting,
    "eberle/ps5-parallel-toplevel-merge-assisting-standard-splitting",
    "pS5-LCP-Merge with JobQueue assisting each other by Andreas Eberle and Timo Bingmann")

PSS_CONTESTANT_PARALLEL(
    eberle_ps5_parallel_toplevel_merge_assisting_binary_splitting,
    "eberle/ps5-parallel-toplevel-merge-assisting-binary-splitting",
    "pS5-LCP-Merge with JobQueue assisting each other by Andreas Eberle and Timo Bingmann")

} // namespace eberle_ps5_parallel_toplevel_merge

#endif // !PSS_SRC_PARALLEL_EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_HEADER

/******************************************************************************/
