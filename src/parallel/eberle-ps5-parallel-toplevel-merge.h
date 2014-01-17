/******************************************************************************
 * src/eberle/parallel/eberle-ps5-parallel-toplevel-merge.h
 *
 * NUMA aware parallel sorting algorithm running separate pS5s on every NUMA
 * node and then combining the result via parallel LCP aware merge.
 *
 ******************************************************************************
 * Copyright (C) 2013-2014 Andreas Eberle <email@andreas-eberle.com>
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
using std::pair;
using std::numeric_limits;

using namespace eberle_lcp_utils;
using namespace eberle_utils;
using namespace eberle_mergesort_lcp;

using namespace jobqueue;
using namespace stringtools;

using namespace bingmann_parallel_sample_sort_lcp;

//typedefs
typedef unsigned char* string;

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
createJobs(JobQueue &jobQueue, const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams, size_t numberOfElements,
        lcp_t baseLcp);

// debug variables (for delta calculations)
string * g_outputBase;
LcpStringPtr g_inputBase;

// variable definitions
typedef uint64_t CHAR_TYPE;
size_t g_lengthOfLongestJob(0);

//structs defining the jobs

struct CopyDataJob : public Job
{
    LcpStringPtr input;
    string* output;
    size_t length;

    CopyDataJob(const LcpStringPtr& input, string* output, size_t length) :
            input(input), output(output), length(length)
    {
        DBG(debug_jobtype_on_creation,
                "CopyDataJob (input: " << (input - g_inputBase) << ", output: " << (output - g_outputBase) << ", length: " << length << ")");
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
                "BinaryMergeJob (input1: " << (input1 - g_inputBase) << ", length1: " << length1 << ", input2: " << (input2 - g_inputBase) << ", length2: " << length2 << ", output: " << (output - g_outputBase) << ")");
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

template<unsigned K>
struct MergeJob : public Job
{
    LcpStringPtr input;
    string* output;
    pair<size_t, size_t>* ranges;
    size_t length;
    lcp_t baseLcp;
    lcp_t nextBaseLcp;

    MergeJob(const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, size_t length, lcp_t baseLcp, lcp_t nextBaseLcp) :
            input(input), output(output), ranges(ranges), length(length), baseLcp(baseLcp), nextBaseLcp(nextBaseLcp)
    {
        DBG(debug_jobtype_on_creation,
                "MergeJob<" << K << "> (output: " << (output - g_outputBase) << ", baseLcp: " << baseLcp << ", nextBaseLcp: " << nextBaseLcp << ", length: " << length << ")");

        for (unsigned k = 0; k < K; ++k)
        {
            DBG(debug_job_details, "" << k << ": " << ranges[k].first << " length: " << ranges[k].second);
        }
    }

    /*
     * returns true if all elements have been written to output
     * false if the merge has been stopped to free work.
     */
    inline bool
    mergeToOutput(JobQueue& jobQueue, LcpStringLoserTree<K> & loserTree)
    {
        for (size_t lastLength = length; length >= MERGE_BULK_SIZE; length -= MERGE_BULK_SIZE, output += MERGE_BULK_SIZE)
        {
            if (g_lengthOfLongestJob == lastLength)
                g_lengthOfLongestJob = length;

            if (g_lengthOfLongestJob < length)
                g_lengthOfLongestJob = length; // else if to prevent work sharing when we just increased g_lengthOfLongestJob
            else if (USE_WORK_SHARING && jobQueue.has_idle() && length > SHARE_WORK_THRESHOLD && g_lengthOfLongestJob == length)
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
        for (unsigned k = 0; k < K; k++)
        { // this ensures that the streams are all equally compared at the start
            input.setLcp(ranges[k].first, baseLcp);
        }

        //merge
        LcpStringLoserTree<K> loserTree(input, ranges);

        if (!mergeToOutput(jobQueue, loserTree))
        {
            // share work
            pair<size_t, size_t> newRanges[K];
            loserTree.getRangesOfRemaining(newRanges, input);

            createJobs(jobQueue, input, output, newRanges, K, length, nextBaseLcp);

            if (g_lengthOfLongestJob == length)
                g_lengthOfLongestJob = 0;
        }

        return true;
    }

    ~MergeJob()
    {
        delete [] ranges;
    }
};

struct InitialSplitJob : public Job
{
    const LcpStringPtr& input;
    string* output;
    pair<size_t, size_t>* ranges;
    size_t length;
    unsigned numStreams;

    InitialSplitJob(const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, size_t length, unsigned numStreams) :
            input(input), output(output), ranges(ranges), length(length), numStreams(numStreams)
    {
        g_lengthOfLongestJob = length; // prevents that the first MergeJob immediately starts splitting itself

        g_inputBase = input;
        g_outputBase = output;
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        createJobs(jobQueue, input, output, ranges, numStreams, length, lcp_t(0));
        g_lengthOfLongestJob = 0;
        return true;
    }

    ~InitialSplitJob()
    {
        delete [] ranges;
    }
};

// Implementations follow.

static inline void
enqueueMergeJob(JobQueue& jobQueue, const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, size_t length, unsigned numStreams,
        lcp_t baseLcp, lcp_t nextBaseLcp)
{
    switch (numStreams)
    {
    case 2:
        jobQueue.enqueue(new BinaryMergeJob(input + ranges[0].first, ranges[0].second, input + ranges[1].first, ranges[1].second, output));
        break;
    case 4:
        jobQueue.enqueue(new MergeJob<4>(input, output, ranges, length, baseLcp, nextBaseLcp));
        break;
    case 8:
        jobQueue.enqueue(new MergeJob<8>(input, output, ranges, length, baseLcp, nextBaseLcp));
        break;
    case 16:
        jobQueue.enqueue(new MergeJob<16>(input, output, ranges, length, baseLcp, nextBaseLcp));
        break;
    case 32:
        jobQueue.enqueue(new MergeJob<32>(input, output, ranges, length, baseLcp, nextBaseLcp));
        break;
    case 64:
        jobQueue.enqueue(new MergeJob<64>(input, output, ranges, length, baseLcp, nextBaseLcp));
        break;
    default:
        DBG(1, "CANNOT MERGE! TO MANY STREAMS: " << numStreams);
        abort();
    }
}

static inline size_t
findNextSplitter(LcpStringPtr& inputStream, const LcpStringPtr& end, lcp_t baseLcp, lcp_t maxAllowedLcp, CHAR_TYPE &lastCharacter, CHAR_TYPE keyMask)
{
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
createJobs(JobQueue &jobQueue, const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams, size_t numberOfElements, lcp_t baseLcp)
{
    DBG(debug_job_creation, std::endl << "CREATING JOBS at baseLcp: " << baseLcp << ", numberOfElements: " << numberOfElements);

    LcpStringPtr inputStreams[numStreams];
    LcpStringPtr ends[numStreams];
    CHAR_TYPE splitterCharacter[numStreams];

    for (unsigned k = 0; k < numStreams; ++k)
    {
        if (ranges[k].second > 0)
        {
            inputStreams[k] = input + ranges[k].first;
            ends[k] = inputStreams[k] + ranges[k].second;
            splitterCharacter[k] = get_char<CHAR_TYPE>(inputStreams[k].str(), baseLcp);
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

    unsigned indexesOfFound[numStreams];

    while (true)
    {
        lcp_t maxAllowedLcp(baseLcp + keyWidth - 1);
        CHAR_TYPE keyMask = numeric_limits < CHAR_TYPE > ::max() << ((key_traits<CHAR_TYPE>::add_depth - keyWidth) * 8);

        CHAR_TYPE currBucket = numeric_limits < CHAR_TYPE > ::max();
        unsigned numberOfFoundBuckets = 0;

        for (unsigned k = 0; k < numStreams; ++k)
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
            const unsigned streamIdx = indexesOfFound[0];
            LcpStringPtr inputStart = inputStreams[streamIdx];
            length += findNextSplitter(inputStreams[streamIdx], ends[streamIdx], baseLcp, maxAllowedLcp, splitterCharacter[streamIdx], keyMask);
            jobQueue.enqueue(new CopyDataJob(inputStart, output, length));
            break;
        }
        case 2:
        {
            const unsigned streamIdx1 = indexesOfFound[0];
            const LcpStringPtr inputStart1 = inputStreams[streamIdx1];
            size_t length1 = findNextSplitter(inputStreams[streamIdx1], ends[streamIdx1], baseLcp, maxAllowedLcp, splitterCharacter[streamIdx1],
                    keyMask);

            const unsigned streamIdx2 = indexesOfFound[1];
            const LcpStringPtr inputStart2 = inputStreams[streamIdx2];
            size_t length2 = findNextSplitter(inputStreams[streamIdx2], ends[streamIdx2], baseLcp, maxAllowedLcp, splitterCharacter[streamIdx2],
                    keyMask);

            jobQueue.enqueue(new BinaryMergeJob(inputStart1, length1, inputStart2, length2, output));
            length = length1 + length2;
            break;
        }

        default:
        {
            unsigned numNewStreams = 1 << ilog2_ceil(numberOfFoundBuckets);
            pair < size_t, size_t > *newRange = new pair<size_t, size_t> [numNewStreams];

            unsigned k = 0;
            for (; k < numberOfFoundBuckets; ++k)
            {
                const unsigned idx = indexesOfFound[k];
                const LcpStringPtr inputStart = inputStreams[idx];
                size_t currLength = findNextSplitter(inputStreams[idx], ends[idx], baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
                newRange[k] = std::make_pair(inputStart - input, currLength);
                length += currLength;
            }
            for (; k < numNewStreams; k++)
            {
                newRange[k] = std::make_pair(0, 0); // this stream is not used
            }

            enqueueMergeJob(jobQueue, input, output, newRange, length, numNewStreams, baseLcp, maxAllowedLcp + 1);
            break;
        }
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
parallelMerge(const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, size_t length, unsigned numStreams)
{
    JobQueue jobQueue;
    DBG(debug_toplevel_merge_duration, "doing parallel merge for " << numStreams << " streams");
    jobQueue.enqueue(new InitialSplitJob(input, output, ranges, length, numStreams));
    jobQueue.loop();
}

void
eberle_ps5_parallel_toplevel_merge(string *strings, size_t n)
{
    int realNumaNodes = numa_num_configured_nodes();
    if (realNumaNodes < 1) realNumaNodes = 1;

    if (realNumaNodes == 1) {
        DBG(1, "No or just one NUMA nodes detected on the system.");
        DBG(1, "pS5-LCP-Mergesort is designed for NUMA systems.");
    }

    unsigned numNumaNodes = std::max(unsigned(4), unsigned(realNumaNodes)); // this max ensures a parallel merge on developer machine
    int numThreadsPerNode = numa_num_configured_cpus() / numNumaNodes;
    if (numThreadsPerNode < 1) numThreadsPerNode = 1;

    //allocate memory for lcps and temporary strings
    string* shadow = new string[n];
    string* tmp = new string[n];

    pair < size_t, size_t > *ranges = new pair<size_t, size_t> [numNumaNodes];
    calculateRanges(ranges, numNumaNodes, n);

    // enable nested parallel regions
    omp_set_nested(true);

    for (unsigned k = 0; k < numNumaNodes; k++)
    {
        size_t start = ranges[k].first;
        size_t length = ranges[k].second;

        const StringPtrOut strptr(strings + start, shadow + start, tmp + start, length);
        parallel_sample_sort_numa(strptr, k % realNumaNodes, numThreadsPerNode);
    }

    // do top level merge
    LcpStringPtr lcpStringPtr(tmp, (lcp_t*) shadow);

    MeasureTime<0> timer;
    timer.start();
    parallelMerge(lcpStringPtr, strings, ranges, n, numNumaNodes);
    timer.stop();

    DBG(debug_toplevel_merge_duration, std::endl << "top level merge needed: " << timer.delta() << " s" << std::endl);

    delete[] shadow;
    delete[] tmp;
}

CONTESTANT_REGISTER_PARALLEL(eberle_ps5_parallel_toplevel_merge, "eberle/ps5-parallel-toplevel-merge",
        "NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle")

}
// namespace eberle_parallel_mergesort_lcp_loosertree

#endif // EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_H_

