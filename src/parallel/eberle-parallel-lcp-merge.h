/******************************************************************************
 * src/parallel/eberle-parallel-lcp-merge.h
 *
 * Parallel LCP aware merge implementation.
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

#ifndef EBERLE_PARALLEL_LCP_MERGE_H_
#define EBERLE_PARALLEL_LCP_MERGE_H_

#include <iostream>
#include <limits>
#include <utility>

#include "../tools/eberle-lcp-losertree.h"
#include "../sequential/eberle-mergesort-lcp.h"

#include "../tools/jobqueue.h"
#include "../tools/stringtools.h"

#include "../tools/debug.h"
#undef DBGX
#define DBGX DBGX_OMP

namespace eberle_parallel_lcp_merge
{
using std::numeric_limits;

using namespace eberle_lcp_utils;
using namespace eberle_mergesort_lcp;

using namespace jobqueue;
using namespace stringtools;

//due to ambigous symbol
using stringtools::string;

//debugging constants
static const bool debug_jobtype_on_creation = false;
static const bool debug_job_details = false;
static const bool debug_job_creation = false;

static const bool debug_created_jobs_count = true;
static const bool debug_merge_start_message = true;

//constants
static const bool USE_WORK_SHARING = true;
static const size_t MERGE_BULK_SIZE = 3000;
static const size_t SHARE_WORK_THRESHOLD = 3 * MERGE_BULK_SIZE;

//method definitions

static inline void
createJobs(JobQueue &jobQueue, const LcpCacheStringPtr* input, unsigned numInputs,
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

    CopyDataJob(const LcpStringPtr& input, string* output)
        : input(input), output(output)
    {
        DBG(debug_jobtype_on_creation,
                "CopyDataJob (output: " << (output - g_outputBase) << ", length: " << input.size << ")");
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        (void) jobQueue;

        input.copyStringsTo(output, input.size);

        return true;
    }
};

struct BinaryMergeJob : public Job
{
    LcpStringPtr input1;
    LcpStringPtr input2;
    string* output;

    BinaryMergeJob(const LcpStringPtr& input1, const LcpStringPtr& input2, string* output) :
            input1(input1), input2(input2), output(output)
    {
        DBG(debug_jobtype_on_creation,
                "BinaryMergeJob (length1: " << input1.size << ", length2: " << input2.size << ", output: " << (output - g_outputBase) << ")");
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        (void) jobQueue;
        input1.lcp() = 0;
        input2.lcp() = 0;
        eberle_lcp_merge(input1, input2, output);

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
    const LcpCacheStringPtr* input;
    unsigned numInputs;
    string* output;
    size_t length;

    InitialSplitJob(const LcpCacheStringPtr* input, unsigned numInputs, string* output, size_t length)
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
findNextSplitter(LcpCacheStringPtr& inputStream,
                 lcp_t baseLcp, lcp_t maxAllowedLcp, CHAR_TYPE& lastCharacter, CHAR_TYPE keyMask)
{
    LcpCacheStringPtr end = inputStream.end();

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
createJobs(JobQueue &jobQueue, const LcpCacheStringPtr* inputStreams, unsigned numInputs,
           string* output, size_t numberOfElements, lcp_t baseLcp)
{
    DBG(debug_job_creation, std::endl << "CREATING JOBS at baseLcp: " << baseLcp << ", numberOfElements: " << numberOfElements);

    LcpCacheStringPtr inputs[numInputs];
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
            const LcpCacheStringPtr start = inputs[idx];
            length += findNextSplitter(inputs[idx],
                                       baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
            jobQueue.enqueue(new CopyDataJob(LcpStringPtr(start, length), output));
            break;
        }
        case 2:
        {
            const unsigned idx1 = indexesOfFound[0];
            const LcpCacheStringPtr start1 = inputs[idx1];
            size_t length1 = findNextSplitter(inputs[idx1],
                                              baseLcp, maxAllowedLcp, splitterCharacter[idx1],
                                              keyMask);

            const unsigned idx2 = indexesOfFound[1];
            const LcpCacheStringPtr start2 = inputs[idx2];
            size_t length2 = findNextSplitter(inputs[idx2],
                                              baseLcp, maxAllowedLcp, splitterCharacter[idx2],
                                              keyMask);

            jobQueue.enqueue(new BinaryMergeJob(LcpStringPtr(start1, length1), LcpStringPtr(start2, length2), output));
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
                const LcpCacheStringPtr start = inputs[idx];
                size_t currLength =
                    findNextSplitter(inputs[idx],
                                     baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
                job->loserTree.streams[k] = start.sub(0, currLength);
                length += currLength;
            }
            for (; k < K; k++)
            {
                // this stream is not used
                job->loserTree.streams[k] = LcpCacheStringPtr(NULL, NULL, NULL, 0);
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
                const LcpCacheStringPtr start = inputs[idx];
                size_t currLength =
                    findNextSplitter(inputs[idx],
                                     baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
                job->loserTree.streams[k] = start.sub(0, currLength);
                length += currLength;
            }
            for (; k < K; k++)
            {
                // this stream is not used
                job->loserTree.streams[k] = LcpCacheStringPtr(NULL, NULL, NULL, 0);
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
sequentialLcpMerge(const LcpCacheStringPtr* input, unsigned numInputs, string* output, size_t length)
{
    DBG(debug_merge_start_message, "doing sequential lcp merge for " << numInputs << " input streams");

    JobQueue jobQueue;

    switch (numInputs)
    {
    case 1:
    {
        jobQueue.enqueue(new CopyDataJob(LcpStringPtr(input[0]), output));
        break;
    }
    case 2:
    {
        jobQueue.enqueue(new BinaryMergeJob(LcpStringPtr(input[0]), LcpStringPtr(input[1]), output));
        break;
    }
    case 3: case 4:
    {
        static const unsigned K = 4;

        MergeJob<K>* job = new MergeJob<K>(output, length, 0, 1);
        job->length = 0;

        unsigned k = 0;
        for (; k < numInputs; ++k)
        {
            job->loserTree.streams[k] = input[k];
            job->length += input[k].size;
        }
        for (; k < K; k++)
        {
            // this stream is not used
            job->loserTree.streams[k] = LcpCacheStringPtr(NULL, NULL, NULL, 0);
        }

        jobQueue.enqueue(job);
        break;
    }
    default:
        abort();
    }

    jobQueue.numaLoop(-1, 1); // use just one thread
}

static inline
void
parallelLcpMerge(const LcpCacheStringPtr* input, unsigned numInputs, string* output, size_t length)
{
    JobQueue jobQueue;
    DBG(debug_merge_start_message, "doing parallel lcp merge for " << numInputs << " input streams using " << omp_get_max_threads() << " threads");
    jobQueue.enqueue(new InitialSplitJob(input, numInputs, output, length));
    jobQueue.numaLoop(-1, omp_get_max_threads());
}


} // namespace eberle_parallel_lcp_merge

#endif // EBERLE_PARALLEL_LCP_MERGE_H_
