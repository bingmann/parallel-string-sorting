/******************************************************************************
 * src/parallel/eberle-ps5-parallel-toplevel-merge.h
 *
 * NUMA aware parallel sorting algorithm running separate pS5s on every NUMA
 * node and then combining the result via parallel LCP aware merge.
 *
 ******************************************************************************
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
 *****************************************************************************/

#ifndef EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_H_
#define EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_H_

#include <iostream>
#include <utility>

#include "eberle-parallel-lcp-merge.h"
#include "bingmann-parallel_sample_sort.h"

#include "../tools/eberle-utilities.h"

#include "../tools/jobqueue.h"
#include "../tools/stringtools.h"

#include "../tools/debug.h"
#undef DBGX
#define DBGX DBGX_OMP


namespace eberle_ps5_parallel_toplevel_merge
{
using std::numeric_limits;

using namespace eberle_utils;

using namespace jobqueue;
using namespace stringtools;

using namespace bingmann_parallel_sample_sort_lcp;

//due to ambigous symbol
using stringtools::string;

//debugging constants
static const bool debug_toplevel_merge_duration = true;

//method definitions
void
eberle_ps5_parallel_toplevel_merge(string *strings, size_t n)
{
    int realNumaNodes = numa_num_configured_nodes();
    if (realNumaNodes < 1) realNumaNodes = 1;

    if (realNumaNodes == 1) {
        DBG(1, "No or just one NUMA nodes detected on the system.");
        DBG(1, "pS5-LCP-Mergesort is designed for NUMA systems.");
        DBG(1, "Continuing anyway, at your own peril!");
    }

    g_statscache >> "num_numa_nodes" << realNumaNodes;

    // this max ensures a parallel merge on developer machine
    int numNumaNodes = realNumaNodes;
    numNumaNodes = std::max(4, realNumaNodes);

    // calculate ranges
    std::pair<size_t, size_t> ranges[numNumaNodes];
    calculateRanges(ranges, numNumaNodes, n);

    LcpCacheStringPtr outputs[numNumaNodes];

    // distribute threads amoung NUMA nodes
    int numThreadsPerNode = omp_get_max_threads() / numNumaNodes;
    int remainThreads = omp_get_max_threads() % numNumaNodes;

    if (numThreadsPerNode == 0)
    {
        DBG(1, "Fewer threads than NUMA nodes detected.");
        DBG(1, "Strange things will happen now.");
        DBG(1, "Continuing anyway, at your own peril!");

        // set num threads = 1 for nodes, and limit number of active nodes via
        // normal num_threads() in next loop
        numThreadsPerNode = 1;

#pragma omp parallel for schedule(dynamic)
        for (int k = 0; k < numNumaNodes; k++)
        {
            size_t start = ranges[k].first;
            size_t length = ranges[k].second;

            parallel_sample_sort_numa(strings + start, length,
                                      k % realNumaNodes, numThreadsPerNode,
                                      outputs[k]);

            outputs[k].size = length;
        }
    }
    else
    {
        ClockTimer all_timer;

#pragma omp parallel for num_threads(numNumaNodes) schedule(static)
        for (int k = 0; k < numNumaNodes; k++)
        {
            size_t start = ranges[k].first;
            size_t length = ranges[k].second;

            int nodeThreads = numThreadsPerNode;
            if (k < remainThreads) nodeThreads++; // distribute extra threads

            DBG(1, "node[" << k << "] gets " << nodeThreads << " threads");

            ClockTimer timer;

            parallel_sample_sort_numa(strings + start, length,
                                      k % realNumaNodes, numThreadsPerNode,
                                      outputs[k]);

            outputs[k].size = length;

            DBG(debug_toplevel_merge_duration, "node[" << k << "] took : " << timer.elapsed() << " s");
        }

        DBG(debug_toplevel_merge_duration, "all nodes took : " << all_timer.elapsed() << " s");
    }

    // calculate cache characters
#pragma omp parallel for num_threads(numNumaNodes) schedule(static)
    for (int k = 0; k < numNumaNodes; k++)
    {
        LcpCacheStringPtr& outputPtr = outputs[k];
        for(unsigned i = 0; i < outputPtr.size; i++)
        {
            string s = outputPtr.strings[i];
            lcp_t lcp = outputPtr.lcps[i];
            char c = s[lcp];
            outputPtr.cachedChars[i] = c;
        }

//        verify_lcp_cache(outputPtr.strings, outputPtr.lcps, outputPtr.cachedChars, outputPtr.size, 0);
    }

    // do top level merge

    ClockTimer timer;
   // eberle_parallel_lcp_merge::sequentialLcpMerge(outputs, numNumaNodes, strings, n);
    eberle_parallel_lcp_merge::parallelLcpMerge(outputs, numNumaNodes, strings, n);

    DBG(debug_toplevel_merge_duration, std::endl << "top level merge needed: " << timer.elapsed() << " s" << std::endl);

    for (int k = 0; k < numNumaNodes; k++)
    {
        numa_free(outputs[k].strings, ranges[k].second * sizeof(string));
        numa_free(outputs[k].lcps, ranges[k].second * sizeof(string));
        numa_free(outputs[k].cachedChars, ranges[k].second * sizeof(char));
    }
}

CONTESTANT_REGISTER_PARALLEL(eberle_ps5_parallel_toplevel_merge,
    "eberle/ps5-parallel-toplevel-merge",
    "NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle")

} // namespace eberle_ps5_parallel_toplevel_merge

#endif // EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_H_
