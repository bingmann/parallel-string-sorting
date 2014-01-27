/******************************************************************************
 * src/eberle/parallel/eberle-ps5-parallel-toplevel-merge.h
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
#include <algorithm>
#include <utility>
#include <limits>

#include "../tools/eberle-utilities.h"
#include "../tools/eberle-lcp-losertree.h"
#include "../sequential/eberle-mergesort-lcp.h"

#include "../tools/jobqueue.h"
#include "../tools/stringtools.h"

#include "../tools/debug.h"
#undef DBGX
#define DBGX DBGX_OMP

#include "../parallel/bingmann-parallel_sample_sort.h"

namespace eberle_ps5_parallel_toplevel_merge
{
using std::numeric_limits;

using namespace eberle_lcp_utils;
using namespace eberle_utils;
using namespace eberle_mergesort_lcp;

using namespace jobqueue;
using namespace stringtools;

using namespace bingmann_parallel_sample_sort_lcp;

//debugging constants
static const bool debug_jobtype_on_creation = false;
static const bool debug_job_details = false;
static const bool debug_job_creation = false;

static const bool debug_created_jobs_count = true;
static const bool debug_toplevel_merge_duration = true;

//constants
static const bool USE_WORK_SHARING = true;
static const size_t MERGE_BULK_SIZE = 3000;
static const size_t SHARE_WORK_THRESHOLD = 3 * MERGE_BULK_SIZE;

//method definitions

static inline void
createJobs(JobQueue &jobQueue, const LcpStringPtr* input, unsigned numInputs,
           string* output, size_t numberOfElements,
           lcp_t baseLcp);

// debug variables (for delta calculations)
string * g_outputBase;

// variable definitions
typedef uint64_t CHAR_TYPE;
size_t g_lengthOfLongestJob = 0;

//structs defining the jobs

struct CopyDataJob : public Job
{
    LcpStringPtr input;
    string* output;
    size_t length;

    CopyDataJob(const LcpStringPtr& input, string* output, size_t length)
        : input(input), output(output), length(length)
    {
        DBG(debug_jobtype_on_creation,
                "CopyDataJob (output: " << (output - g_outputBase) << ", length: " << length << ")");
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        (void) jobQueue;

        input.copyStringsTo(output, length);

        return true;
    }
};

struct BinaryMergeJob : public Job
{
    LcpStringPtr input1;
    size_t length1;
    LcpStringPtr input2;
    size_t length2;
    string* output;

    BinaryMergeJob(const LcpStringPtr& input1, size_t length1, const LcpStringPtr& input2, size_t length2, string* output) :
            input1(input1), length1(length1), input2(input2), length2(length2), output(output)
    {
        DBG(debug_jobtype_on_creation,
                "BinaryMergeJob (length1: " << length1 << ", length2: " << length2 << ", output: " << (output - g_outputBase) << ")");
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        (void) jobQueue;
        input1.setLcp(0, 0);
        input2.setLcp(0, 0);
        eberle_lcp_merge(input1, length1, input2, length2, output);

        return true;
    }
};

template <unsigned K>
struct MergeJob : public Job
{
    LcpStringLoserTree<K> loserTree;

    string* output;
    size_t length;

    lcp_t baseLcp;
    lcp_t nextBaseLcp;

    MergeJob(string* output, size_t length, lcp_t baseLcp, lcp_t nextBaseLcp)
        : output(output), length(length), baseLcp(baseLcp), nextBaseLcp(nextBaseLcp)
    {
        DBG(debug_jobtype_on_creation, "MergeJob<" << K << "> (output: " << (output - g_outputBase) << ", baseLcp: " << baseLcp << ", nextBaseLcp: " << nextBaseLcp << ", length: " << length << ")");
    }

    /*
     * returns true if all elements have been written to output
     * false if the merge has been stopped to free work.
     */
    inline bool
    mergeToOutput(JobQueue& jobQueue)
    {
        for (size_t lastLength = length; length >= MERGE_BULK_SIZE; length -= MERGE_BULK_SIZE, output += MERGE_BULK_SIZE)
        {
            if (g_lengthOfLongestJob == lastLength)
                g_lengthOfLongestJob = length;

            if (g_lengthOfLongestJob < length)
                g_lengthOfLongestJob = length; // else if to prevent work sharing when we just increased g_lengthOfLongestJob
            else if (USE_WORK_SHARING &&
                     jobQueue.has_idle() &&
                     length > SHARE_WORK_THRESHOLD &&
                     g_lengthOfLongestJob == length)
                return false;

            loserTree.writeElementsToStream(output, MERGE_BULK_SIZE);
            lastLength = length;
        }

        loserTree.writeElementsToStream(output, length);

        return true;
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        loserTree.initTree(baseLcp);

        // merge

        if (!mergeToOutput(jobQueue))
        {
            // share work by splitting remaining streams

            createJobs(jobQueue, loserTree.getRemaining(), K,
                       output, length, nextBaseLcp);

            if (g_lengthOfLongestJob == length)
                g_lengthOfLongestJob = 0;
        }

        return true;
    }
};

struct InitialSplitJob : public Job
{
    const LcpStringPtr* input;
    unsigned numInputs;
    string* output;
    size_t length;

    InitialSplitJob(const LcpStringPtr* input, unsigned numInputs, string* output, size_t length)
        : input(input), numInputs(numInputs), output(output), length(length)
    {
        g_lengthOfLongestJob = length; // prevents that the first MergeJob immediately starts splitting itself

        g_outputBase = output;
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        createJobs(jobQueue, input, numInputs, output, length, lcp_t(0));
        g_lengthOfLongestJob = 0;
        return true;
    }
};

static inline size_t
findNextSplitter(LcpStringPtr& inputStream,
                 lcp_t baseLcp, lcp_t maxAllowedLcp, CHAR_TYPE& lastCharacter, CHAR_TYPE keyMask)
{
    LcpStringPtr end = inputStream.end();

    size_t length = 1;
    ++inputStream;

    for (; inputStream < end; ++inputStream, ++length)
    {
        lcp_t lcp = inputStream.lcp();

        if (lcp <= maxAllowedLcp)
        {
            CHAR_TYPE character = get_char<CHAR_TYPE>(inputStream.str(), baseLcp);
            if ((character & keyMask) != (lastCharacter & keyMask))
            {
                lastCharacter = character;
                return length;
            }
        }
    }

    lastCharacter = numeric_limits<CHAR_TYPE>::max();
    return length;
}

static inline void
createJobs(JobQueue &jobQueue, const LcpStringPtr* inputStreams, unsigned numInputs,
           string* output, size_t numberOfElements, lcp_t baseLcp)
{
    DBG(debug_job_creation, std::endl << "CREATING JOBS at baseLcp: " << baseLcp << ", numberOfElements: " << numberOfElements);

    LcpStringPtr inputs[numInputs];
    CHAR_TYPE splitterCharacter[numInputs];

    for (unsigned k = 0; k < numInputs; ++k)
    {
        inputs[k] = inputStreams[k];

        if (inputStreams[k].size > 0)
        {
            splitterCharacter[k] = get_char<CHAR_TYPE>(inputs[k].str(), baseLcp);
        }
        else
        {
            splitterCharacter[k] = numeric_limits<CHAR_TYPE>::max();
        }
    }

    const unsigned overProvFactor = 500;
    const size_t expectedJobLength = std::max(MERGE_BULK_SIZE, numberOfElements / (overProvFactor * numa_num_configured_cpus()));

    DBG(debug_job_creation, "Expected job length: " << expectedJobLength);

    unsigned keyWidth = 8;
    unsigned createdJobsCtr = 0;
    size_t elementsProcessed = 0;

    unsigned indexesOfFound[numInputs];

    while (true)
    {
        lcp_t maxAllowedLcp = baseLcp + keyWidth - 1;
        CHAR_TYPE keyMask = numeric_limits<CHAR_TYPE>::max() << ((key_traits<CHAR_TYPE>::add_depth - keyWidth) * 8);

        CHAR_TYPE currBucket = numeric_limits<CHAR_TYPE>::max();
        unsigned numberOfFoundBuckets = 0;

        for (unsigned k = 0; k < numInputs; ++k)
        {
            CHAR_TYPE splitter = splitterCharacter[k] & keyMask;
            if (splitter <= currBucket)
            {
                if (splitter < currBucket)
                {
                    currBucket = splitter;

                    indexesOfFound[0] = k;
                    numberOfFoundBuckets = 1;
                }
                else
                {
                    indexesOfFound[numberOfFoundBuckets] = k;
                    numberOfFoundBuckets++;
                }
            }
        }

        if (currBucket == (numeric_limits<CHAR_TYPE>::max() & keyMask))
        {
            break;
        }

        size_t length = 0;

        switch (numberOfFoundBuckets)
        {
        case 1:
        {
            const unsigned idx = indexesOfFound[0];
            LcpStringPtr start = inputs[idx];
            length += findNextSplitter(inputs[idx],
                                       baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
            jobQueue.enqueue(new CopyDataJob(start, output, length));
            break;
        }
        case 2:
        {
            const unsigned idx1 = indexesOfFound[0];
            const LcpStringPtr start1 = inputs[idx1];
            size_t length1 = findNextSplitter(inputs[idx1],
                                              baseLcp, maxAllowedLcp, splitterCharacter[idx1],
                                              keyMask);

            const unsigned idx2 = indexesOfFound[1];
            const LcpStringPtr start2 = inputs[idx2];
            size_t length2 = findNextSplitter(inputs[idx2],
                                              baseLcp, maxAllowedLcp, splitterCharacter[idx2],
                                              keyMask);

            jobQueue.enqueue(new BinaryMergeJob(start1, length1, start2, length2, output));
            length = length1 + length2;
            break;
        }
        case 3: case 4:
        {
            static const unsigned K = 4;

            MergeJob<K>* job = new MergeJob<K>(output, length, baseLcp, maxAllowedLcp + 1);

            unsigned k = 0;
            for (; k < numberOfFoundBuckets; ++k)
            {
                const unsigned idx = indexesOfFound[k];
                const LcpStringPtr start = inputs[idx];
                size_t currLength =
                    findNextSplitter(inputs[idx],
                                     baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
                job->loserTree.streams[k] = start.sub(0, currLength);
                length += currLength;
            }
            for (; k < K; k++)
            {
                // this stream is not used
                job->loserTree.streams[k] = LcpStringPtr(NULL, NULL, 0);
            }

            job->length = length;
            jobQueue.enqueue(job);
            break;
        }
        case 5: case 6: case 7: case 8:
        {
            static const unsigned K = 8;

            MergeJob<K>* job = new MergeJob<K>(output, length, baseLcp, maxAllowedLcp + 1);

            unsigned k = 0;
            for (; k < numberOfFoundBuckets; ++k)
            {
                const unsigned idx = indexesOfFound[k];
                const LcpStringPtr start = inputs[idx];
                size_t currLength =
                    findNextSplitter(inputs[idx],
                                     baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
                job->loserTree.streams[k] = start.sub(0, currLength);
                length += currLength;
            }
            for (; k < K; k++)
            {
                // this stream is not used
                job->loserTree.streams[k] = LcpStringPtr(NULL, NULL, 0);
            }

            job->length = length;
            jobQueue.enqueue(job);
            break;
        }
        default:
            DBG(1, "Found " << numberOfFoundBuckets << ", which is more NUMA nodes than expected, ADD MORE CASES IN SWITCH.");
            abort();
        }

        output += length;
        elementsProcessed += length;
        createdJobsCtr++;

        const unsigned expectedCreatedJobs = elementsProcessed / expectedJobLength;
        const int diffExpectedReal = int(expectedCreatedJobs) - int(createdJobsCtr);

        const int tollerance = expectedCreatedJobs / 30 + 5;

        if (diffExpectedReal <= -tollerance)
        {
            keyWidth = std::max(unsigned(1), keyWidth - 1);

            DBG(debug_job_creation, "decreased key to " << keyWidth << "  diff: " << diffExpectedReal);
        }
        else if (diffExpectedReal >= tollerance)
        {
            keyWidth = std::min(unsigned(8), keyWidth + 1);

            DBG(debug_job_creation, "increased key to " << keyWidth << "  diff: " << diffExpectedReal);
        }
    }

    DBG(debug_created_jobs_count, "Created " << createdJobsCtr << " Jobs!");
}

static inline
void
parallelMerge(const LcpStringPtr* input, unsigned numInputs, string* output, size_t length)
{
    JobQueue jobQueue;
    DBG(debug_toplevel_merge_duration, "doing parallel merge for " << numInputs << " input streams using " << omp_get_max_threads() << " threads");
    jobQueue.enqueue(new InitialSplitJob(input, numInputs, output, length));
    jobQueue.numaLoop(-1, omp_get_max_threads());
}

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

    LcpStringPtr output[numNumaNodes];

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
                                      &output[k].strings, &output[k].lcps);

            output[k].size = length;
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
                                      &output[k].strings, &output[k].lcps);

            output[k].size = length;

            DBG(debug_toplevel_merge_duration, "node[" << k << "] took : " << timer.elapsed() << " s");
        }

        DBG(debug_toplevel_merge_duration, "all nodes took : " << all_timer.elapsed() << " s");
    }

    // do top level merge

    ClockTimer timer;
    parallelMerge(output, numNumaNodes, strings, n);

    DBG(debug_toplevel_merge_duration, std::endl << "top level merge needed: " << timer.elapsed() << " s" << std::endl);

    for (int k = 0; k < numNumaNodes; k++)
    {
        delete[] output[k].strings;
        delete[] output[k].lcps;
    }
}

CONTESTANT_REGISTER_PARALLEL(eberle_ps5_parallel_toplevel_merge,
    "eberle/ps5-parallel-toplevel-merge",
    "NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle")

} // namespace eberle_ps5_parallel_toplevel_merge

#endif // EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_H_
