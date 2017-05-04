/*******************************************************************************
 * src/parallel/eberle-parallel-lcp-merge-binary-splitting.hpp
 *
 * Parallel LCP aware merge implementation with binary splitting algorithm.
 *
 *******************************************************************************
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
 ******************************************************************************/

#ifndef PSS_SRC_PARALLEL_EBERLE_PARALLEL_LCP_MERGE_BINARY_SPLITTING_HEADER
#define PSS_SRC_PARALLEL_EBERLE_PARALLEL_LCP_MERGE_BINARY_SPLITTING_HEADER

#include "eberle-parallel-lcp-merge.hpp"

namespace eberle_parallel_lcp_merge {

// debugging constants
static const bool debug_binary_splitting = false;
static const bool debug_binary_splitting_splits_count = true;

// method definitions
static inline void
createJobsBinarySplitting(JobQueue& jobQueue, const LcpCacheStringPtr* inputStreams, unsigned numInputs, string* output, size_t numberOfElements);

//structs defining the jobs
template <unsigned K>
struct MergeJobBinarySplitting : public Job
{
    LcpCacheStringLoserTree<K> loserTree;

    string                     * output;
    size_t                     length;
    bool                       splittable;

    MergeJobBinarySplitting(const LcpCacheStringPtr* inputs, unsigned numInputs, string* output, size_t length, bool splittable)
        : loserTree(inputs, numInputs), output(output), length(length), splittable(splittable)
    {
        g_mergeJobsCreated++;
        LOGC(debug_jobtype_on_creation)
            << "MergeJobStandardSplitting<" << K << "> (output: "
            << (output - g_outputBase) << ", length: " << length << ")";
    }

    inline bool shouldShareWork(JobQueue& jobQueue)
    {
        return USE_WORK_SHARING && splittable && jobQueue.has_idle() && length > SHARE_WORK_THRESHOLD;
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
            else if (shouldShareWork(jobQueue) && g_lengthOfLongestJob == length)
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
        if (shouldShareWork(jobQueue))
        {
            createJobsBinarySplitting(jobQueue, loserTree.getRemaining(), K, output, length);
        }
        else {
            loserTree.initTree(0);

            // merge

            if (!mergeToOutput(jobQueue))
            {
                // share work by splitting remaining streams

                createJobsBinarySplitting(jobQueue, loserTree.getRemaining(), K, output, length);

                if (g_lengthOfLongestJob == length)
                    g_lengthOfLongestJob = 0;
            }
        }

        return true;
    }
};

static inline void
enqueueBinarySplittingJob(JobQueue& jobQueue, const LcpCacheStringPtr* inputs, unsigned numInputs, string* output, size_t jobLength, bool splittable)
{
    if (numInputs == 1)
        jobQueue.enqueue(new CopyDataJob(inputs[0], output));

    else if (numInputs <= 2)
        jobQueue.enqueue(new BinaryMergeJob(inputs[0], inputs[1], 0, output));

    else if (numInputs <= 4)
        jobQueue.enqueue(new MergeJobBinarySplitting<4>(inputs, numInputs, output, jobLength, splittable));

    else if (numInputs <= 8)
        jobQueue.enqueue(new MergeJobBinarySplitting<8>(inputs, numInputs, output, jobLength, splittable));

    else if (numInputs <= 16)
        jobQueue.enqueue(new MergeJobBinarySplitting<16>(inputs, numInputs, output, jobLength, splittable));

    else if (numInputs <= 32)
        jobQueue.enqueue(new MergeJobBinarySplitting<32>(inputs, numInputs, output, jobLength, splittable));

    else if (numInputs <= 64)
        jobQueue.enqueue(new MergeJobBinarySplitting<64>(inputs, numInputs, output, jobLength, splittable));

    else
    {
        LOG1 << "Can't create job with that many streams. Add more cases.";
        abort();
    }
}

static inline void
createJobsBinarySplitting(JobQueue& jobQueue, const LcpCacheStringPtr* inputStreams, unsigned numInputs, string* output, size_t numberOfElements)
{
    LOGC(debug_binary_splitting)
        << "CREATING JOBS for numberOfElements: " << numberOfElements;
    g_splittingsExecuted++;
    ClockTimer splittingTimer;
    splittingTimer.start();

    const unsigned numSplitters = numInputs;

    string splitters[numSplitters];
    LcpCacheStringPtr streams[numSplitters];
    unsigned nonEmptyStreams = 0;

    for (unsigned i = 0; i < numInputs; i++)
    {
        streams[i] = inputStreams[i];

        if (!streams[i].empty())
        {
            splitters[nonEmptyStreams] = streams[i].strings[streams[i].size / 2];
            ++nonEmptyStreams;
        }
    }

    string splitterString = splitters[unsigned(rand() % nonEmptyStreams)];

    LOGC(debug_binary_splitting)
        << "SplitterString: " << splitterString;

    LcpCacheStringPtr jobStreams[2][numInputs];
    unsigned nonEmptyCtr[2] = { 0, 0 };
    unsigned jobLength[2] = { 0, 0 };

    for (unsigned i = 0; i < numInputs; i++)
    {
        LcpCacheStringPtr stream = streams[i];

        if (!stream.empty())
        {
            const size_t idx = stream.binarySearch(splitterString);

            jobStreams[0][nonEmptyCtr[0]] = stream.sub(0, idx);
            nonEmptyCtr[0]++;
            jobLength[0] += idx;

            LOGC(debug_binary_splitting) << "Found at [" << idx << "]: ";

            const size_t restLength = stream.size - idx;
            if (restLength > 0)
            {
                jobStreams[1][nonEmptyCtr[1]] = stream.sub(idx, restLength);
                nonEmptyCtr[1]++;
                jobLength[1] += restLength;
            }
        }
    }

    if (jobLength[0] > 0) {
        // only splittable, if we already split this job (think of long equal string sequences)
        enqueueBinarySplittingJob(jobQueue, jobStreams[0], nonEmptyCtr[0], output, jobLength[0], jobLength[1] > 0);
        output += jobLength[0];
    }
    if (jobLength[1] > 0) {
        // only splittable, if we already split this job (think of long equal string sequences)
        enqueueBinarySplittingJob(jobQueue, jobStreams[1], nonEmptyCtr[1], output, jobLength[1], jobLength[0] > 0);
    }

    g_splittingTime += splittingTimer.elapsed();
}

static inline void
parallelLcpMergeBinarySplitting(const LcpCacheStringPtr* inputs, unsigned numInputs, string* output, size_t length)
{
    g_outputBase = output;
    g_splittingsExecuted = 0;
    g_mergeJobsCreated = 0;
    g_splittingTime = 0;

    ClockTimer timer;
    timer.start();

    JobQueue jobQueue;
    LOGC(debug_merge_start_message)
        << "doing parallel lcp merge for " << numInputs << " input streams using " << omp_get_max_threads() << " threads with binary splitting";
    enqueueBinarySplittingJob(jobQueue, inputs, numInputs, output, length, true);
    jobQueue.numaLoop(-1, omp_get_max_threads());

    LOGC(debug_binary_splitting_splits_count)
        << "Binary Splitting executed " << g_splittingsExecuted << " splittings; created " << g_mergeJobsCreated << " jobs";

    g_stats >> "toplevelmerge_time" << timer.elapsed();
    g_stats >> "splittings_executed" << g_splittingsExecuted;
    g_stats >> "mergejobs_created" << g_mergeJobsCreated;
    g_stats >> "splitting_time" << g_splittingTime;
}

} // namespace eberle_parallel_lcp_merge

#endif // !PSS_SRC_PARALLEL_EBERLE_PARALLEL_LCP_MERGE_BINARY_SPLITTING_HEADER

/******************************************************************************/
