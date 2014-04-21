/******************************************************************************
 * src/parallel/eberle-parallel-lcp-merge-lcp-splitting.h
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

#ifndef EBERLE_PARALLEL_LCP_MERGE_LCP_SPLITTING_H_
#define EBERLE_PARALLEL_LCP_MERGE_LCP_SPLITTING_H_

#include "eberle-parallel-lcp-merge.h"

#include "../tools/debug.h"
#undef DBGX
#define DBGX DBGX_OMP

namespace eberle_parallel_lcp_merge
{

// debugging constants
static const bool debug_job_details = false;
static const bool debug_job_creation = false;
static const bool debug_splitter_key_mask = false;

static const bool debug_created_jobs_count = true;

// variable definitions
typedef uint64_t CHAR_TYPE;

// method definitions
static inline void
createJobsLcpSplitting(JobQueue &jobQueue, const LcpCacheStringPtr* input, unsigned numInputs,
                        string* output, size_t numberOfElements, lcp_t baseLcp);

//structs defining the jobs
template <unsigned K>
struct MergeJobLcpSplitting : public Job
{
    LcpStringLoserTree<K> loserTree;

    string* output;
    size_t length;

    lcp_t baseLcp;
    lcp_t nextBaseLcp;

    MergeJobLcpSplitting(string* output, size_t length, lcp_t baseLcp, lcp_t nextBaseLcp)
        : output(output), length(length), baseLcp(baseLcp), nextBaseLcp(nextBaseLcp)
    {
        g_mergeJobsCreated++;
        DBG(debug_jobtype_on_creation, "MergeJobLcpSplitting<" << K << "> (output: " << (output - g_outputBase) << ", baseLcp: " << baseLcp << ", nextBaseLcp: " << nextBaseLcp << ", length: " << length << ")");
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

            createJobsLcpSplitting(jobQueue, loserTree.getRemaining(), K,
                       output, length, nextBaseLcp);

            if (g_lengthOfLongestJob == length)
                g_lengthOfLongestJob = 0;
        }

        return true;
    }
};

struct InitialJobLcpSplitting : public Job
{
    const LcpCacheStringPtr* input;
    unsigned numInputs;
    string* output;
    size_t length;

    InitialJobLcpSplitting(const LcpCacheStringPtr* input, unsigned numInputs, string* output, size_t length)
        : input(input), numInputs(numInputs), output(output), length(length)
    {
        g_lengthOfLongestJob = length; // prevents that the first MergeJob immediately starts splitting itself

        g_outputBase = output;
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        createJobsLcpSplitting(jobQueue, input, numInputs, output, length, lcp_t(0));
        g_lengthOfLongestJob = 0;
        return true;
    }
};

// Search forward in input streams, skip all strings with larger LCP, test
// equal positions if they still match the splitter.
static inline size_t
findNextSplitter(LcpCacheStringPtr& inputStream,
                 lcp_t baseLcp, lcp_t maxAllowedLcp, CHAR_TYPE& lastCharacter, CHAR_TYPE keyMask)
{
    LcpCacheStringPtr end = inputStream.end();

    size_t length = 1;
    ++inputStream;

    for (; inputStream < end; ++inputStream, ++length)
    {
        lcp_t lcp = inputStream.firstLcp();

        if (lcp <= maxAllowedLcp)
        {
            CHAR_TYPE character = get_char<CHAR_TYPE>(inputStream.firstString(), baseLcp);
            if ((character & keyMask) != (lastCharacter & keyMask))
            {
                lastCharacter = character;
                return length;
            }
        }
    }

    return length;
}

// Create K-way MergeJob for the selected input streams.
template <unsigned K>
static inline unsigned
createMergeJob(JobQueue &jobQueue, string* output, LcpCacheStringPtr* inputs, unsigned* indexesOfFound,
        unsigned numberOfFoundBuckets, unsigned baseLcp, unsigned maxAllowedLcp, CHAR_TYPE* splitterCharacter, CHAR_TYPE keyMask)
{
    unsigned length = 0;
    MergeJobLcpSplitting<K>* job = new MergeJobLcpSplitting<K>(output, 0, baseLcp, maxAllowedLcp + 1);

    unsigned k = 0;
    for (; k < numberOfFoundBuckets; ++k)
    {
        const unsigned idx = indexesOfFound[k];
        const LcpCacheStringPtr start = inputs[idx];
        size_t currLength = findNextSplitter(inputs[idx], baseLcp,
                             maxAllowedLcp, splitterCharacter[idx], keyMask);
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

    return length;
}

static inline void
createJobsLcpSplitting(JobQueue &jobQueue, const LcpCacheStringPtr* inputStreams, unsigned numInputs,
           string* output, size_t numberOfElements, lcp_t baseLcp)
{
    DBG(debug_job_creation, "CREATING JOBS at baseLcp: " << baseLcp << ", numberOfElements: " << numberOfElements);
    g_splittingsExecuted++;
    ClockTimer splittingTimer;
    splittingTimer.start();

    LcpCacheStringPtr inputs[numInputs];
    CHAR_TYPE splitterCharacter[numInputs];

    for (unsigned k = 0; k < numInputs; ++k)
    {
        inputs[k] = inputStreams[k];

        if (!inputs[k].empty())
        {
            splitterCharacter[k] = get_char<CHAR_TYPE>(inputs[k].firstString(), baseLcp);
        }
    }

    const unsigned overProvFactor = 100;
    const size_t expectedJobLength = std::max(MERGE_BULK_SIZE, numberOfElements / (overProvFactor * numa_num_configured_cpus()));

    DBG(debug_job_creation, "Expected job length: " << expectedJobLength);

    unsigned keyWidth = 8;
    unsigned createdJobsCtr = 0;
    size_t elementsProcessed = 0;

    unsigned indexesOfFound[numInputs];

    while (true)
    {
        lcp_t maxAllowedLcp = baseLcp + keyWidth - 1;

        // construct current splitter mask
        CHAR_TYPE keyMask = numeric_limits<CHAR_TYPE>::max() << ((key_traits<CHAR_TYPE>::add_depth - keyWidth) * 8);

        // *** Select input streams with smallest current splitter ***

        CHAR_TYPE currBucket = numeric_limits<CHAR_TYPE>::max();
        unsigned numberOfFoundBuckets = 0;

        DBG(debug_splitter_key_mask, std::hex << keyMask);

        for (unsigned k = 0; k < numInputs; ++k)
        {
            if (inputs[k].empty()) continue;

            CHAR_TYPE splitter = splitterCharacter[k] & keyMask;

            if (splitter < currBucket) // smaller splitter found, reset.
            {
                currBucket = splitter;

                indexesOfFound[0] = k;
                numberOfFoundBuckets = 1;
            }
            else if (splitter == currBucket) // additional splitter found
            {
                indexesOfFound[numberOfFoundBuckets] = k;
                numberOfFoundBuckets++;
            }
        }

        if (numberOfFoundBuckets == 0)
            break;

        // *** Create Job to Merge Buckets, Depending on Input Streams ***

        size_t length = 0;

        switch (numberOfFoundBuckets)
        {
        case 1:
        {
            const unsigned idx = indexesOfFound[0];
            const LcpCacheStringPtr start = inputs[idx];
            length += findNextSplitter(inputs[idx],
                                       baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
            jobQueue.enqueue(new CopyDataJob(start.sub(0, length), output));
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

            jobQueue.enqueue(new BinaryMergeJob(start1.sub(0, length1), start2.sub(0, length2),
                                                baseLcp, output));
            length = length1 + length2;
            break;
        }
        case 3: case 4:
        {
            length = createMergeJob<4>(jobQueue, output, inputs, indexesOfFound, numberOfFoundBuckets,
                                       baseLcp, maxAllowedLcp, splitterCharacter, keyMask);
            break;
        }
        case 5: case 6: case 7: case 8:
        {
            length = createMergeJob<8>(jobQueue, output, inputs, indexesOfFound, numberOfFoundBuckets,
                                       baseLcp, maxAllowedLcp, splitterCharacter, keyMask);
            break;
        }
        case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16:
        {
            length = createMergeJob<16>(jobQueue, output, inputs, indexesOfFound, numberOfFoundBuckets,
                                        baseLcp, maxAllowedLcp, splitterCharacter, keyMask);
            break;
        }
        default:
            if (numberOfFoundBuckets <= 32) {
                length = createMergeJob<32>(jobQueue, output, inputs, indexesOfFound, numberOfFoundBuckets,
                                            baseLcp, maxAllowedLcp, splitterCharacter, keyMask);
            }
            else if (numberOfFoundBuckets <= 64) {
                length = createMergeJob<64>(jobQueue, output, inputs, indexesOfFound, numberOfFoundBuckets,
                                            baseLcp, maxAllowedLcp, splitterCharacter, keyMask);
            }
            else {
                DBG(1, "Found " << numberOfFoundBuckets << ", which is more input streams than expected, ADD MORE CASES IN SWITCH.");
                abort();
            }
        }

        output += length;
        elementsProcessed += length;
        createdJobsCtr++;

        // *** Heuristic to Adapt Size of Splitter Keys ***

        const unsigned expectedCreatedJobs = elementsProcessed / expectedJobLength;
        const int diffExpectedReal = int(expectedCreatedJobs) - int(createdJobsCtr);

        const int tolerance = expectedCreatedJobs / 30 + 5;

        if (diffExpectedReal <= -tolerance)
        {
            keyWidth = std::max(unsigned(1), keyWidth - 1);

            DBG(debug_job_creation, "decreased key to " << keyWidth << "  diff: " << diffExpectedReal);
        }
        else if (diffExpectedReal >= tolerance)
        {
            keyWidth = std::min(unsigned(8), keyWidth + 1);

            DBG(debug_job_creation, "increased key to " << keyWidth << "  diff: " << diffExpectedReal);
        }
    }

    DBG(debug_created_jobs_count, "Created " << createdJobsCtr << " jobs!");
    g_splittingTime += splittingTimer.elapsed();
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
        jobQueue.enqueue(new CopyDataJob(input[0].sub(0, length), output));
        break;
    }
    case 2:
    {
        jobQueue.enqueue(new BinaryMergeJob(input[0], input[1], 0, output));
        break;
    }
    case 3: case 4:
    {
        static const unsigned K = 4;

        MergeJobLcpSplitting<K>* job = new MergeJobLcpSplitting<K>(output, length, 0, 1);
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
    case 5: case 6: case 7: case 8:
    {
        static const unsigned K = 8;

        MergeJobLcpSplitting<K>* job = new MergeJobLcpSplitting<K>(output, length, 0, 1);
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
    g_outputBase = output;
    g_splittingsExecuted = 0;
    g_mergeJobsCreated = 0;
    g_splittingTime = 0;

	ClockTimer timer;
	timer.start();

    JobQueue jobQueue;
    DBG(debug_merge_start_message, "doing parallel lcp merge for " << numInputs << " input streams using " << omp_get_max_threads() << " threads");
    jobQueue.enqueue(new InitialJobLcpSplitting(input, numInputs, output, length));
    jobQueue.loop();
	
	
    g_stats >> "toplevelmerge_time" << timer.elapsed();
    g_stats >> "splittings_executed" << g_splittingsExecuted;
    g_stats >> "mergejobs_created" << g_mergeJobsCreated;
    g_stats >> "splitting_time" << g_splittingTime;
}


} // namespace eberle_parallel_lcp_merge

#endif // EBERLE_PARALLEL_LCP_MERGE_LCP_SPLITTING_H_
