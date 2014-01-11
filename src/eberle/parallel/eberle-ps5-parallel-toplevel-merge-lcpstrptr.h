#ifndef EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_LCPSTRPTR_H_
#define EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_LCPSTRPTR_H_

#include <iostream>

#include "../utils/types.h"
#include "../utils/utility-functions.h"
#include "../utils/lcp-string-losertree-lcpstrptr.h"
#include "../sequential/eberle-mergesort-lcp.h"

#include "../../tools/jobqueue.h"
#include "../../tools/stringtools.h"

#include "../../parallel/bingmann-parallel_sample_sort.h"

//#define PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
//#define PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS_DETAILED
#define PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION
#define PARALLEL_LCP_MERGE_DEBUG_TOP_LEVEL_MERGE_DURATION

namespace eberle_parallel_mergesort_lcp_loosertree
{

using namespace std;

using namespace types;
using namespace eberle_lcp_utils;
using namespace eberle_utils;
using namespace eberle_mergesort_lcp;

using namespace jobqueue;
using namespace stringtools;

using namespace bingmann_parallel_sample_sort_lcp;

//typedefs
typedef unsigned char* string;

//constants
static const bool USE_WORK_SHARING = true;
static const size_t MERGE_BULK_SIZE = 3000;
static const size_t SHARE_WORK_THRESHOLD = 3 * MERGE_BULK_SIZE;

//method definitions

static inline void
createJobs(JobQueue &jobQueue, const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams, size_t numberOfElements,
        lcp_t baseLcp);

// debug variables
string * outputBase;
LcpStringPtr inputBase;

// variable definitions
typedef uint64_t CHAR_TYPE;
size_t lengthOfLongestJob(0);

//structs defining the jobs

struct CopyDataJob : public Job
{
    LcpStringPtr input;
    string* output;
    size_t length;

    CopyDataJob(const LcpStringPtr& input, string* output, size_t length) :
            input(input), output(output), length(length)
    {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
        {
            cout << "CopyDataJob (input: " << (input - inputBase) << ", output: " << (output - outputBase) << ", length: " << length << ")" << endl;
        }
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
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
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
        {
            cout << "BinaryMergeJob (input1: " << (input1 - inputBase) << ", length1: " << length1 << ", input2: " << (input2 - inputBase)
            << ", length2: " << length2 << ", output: " << (output - outputBase) << ")" << endl;
        }
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
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
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
            {
                cout << "MergeJob<" << K << "> (output: " << (output - outputBase) << ", baseLcp: " << baseLcp << ", nextBaseLcp: " << nextBaseLcp
                << ", length: " << length << ")" << endl;
#ifdef PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS_DETAILED
                for (unsigned k = 0; k < K; ++k)
                {
                    cout << k << ": " << ranges[k].first << " length: " << ranges[k].second << endl;
                }
                cout << endl;
#endif // PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS_DETAILED
            }
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
        }

        /*
         * returns true if all elements have been written to output
         * false if the merge has been stopped to free work.
         */
        inline bool
        mergeToOutput(JobQueue& jobQueue, LcpStringLcpPtrLoserTree<K> * loserTree)
        {
            for (size_t lastLength = length; length >= MERGE_BULK_SIZE; length -= MERGE_BULK_SIZE, output += MERGE_BULK_SIZE)
            {
                if (lengthOfLongestJob == lastLength)
                    lengthOfLongestJob = length;

                if (lengthOfLongestJob < length)
                    lengthOfLongestJob = length; // else if to prevent work sharing when we just increased lengthOfLongestJob
                else if (USE_WORK_SHARING && jobQueue.has_idle() && length > SHARE_WORK_THRESHOLD && lengthOfLongestJob == length)
                    return false;

                loserTree->writeElementsToStream(output, MERGE_BULK_SIZE);
                lastLength = length;
            }

            loserTree->writeElementsToStream(output, length);

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
            LcpStringLcpPtrLoserTree<K> *loserTree = new LcpStringLcpPtrLoserTree<K>(input, ranges);

            if (!mergeToOutput(jobQueue, loserTree))
            {
                // share work
                pair < size_t, size_t > newRanges[K];
                loserTree->getRangesOfRemaining(newRanges, input);

                createJobs(jobQueue, input, output, newRanges, K, length, nextBaseLcp);

                if (lengthOfLongestJob == length)
                    lengthOfLongestJob = 0;
            }

            delete loserTree;

            return true;
        }

        ~MergeJob()
        {
            delete ranges;
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
        lengthOfLongestJob = length; // prevents that the first MergeJob immediately starts splitting itself

        inputBase = input;
        outputBase = output;
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        createJobs(jobQueue, input, output, ranges, numStreams, length, lcp_t(0));
        lengthOfLongestJob = 0;
        return true;
    }

    ~InitialSplitJob()
    {
        delete ranges;
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
        cerr << "CANNOT MERGE! TO MANY STREAMS: " << numStreams;
        break;
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

    lastCharacter = numeric_limits < CHAR_TYPE > ::max();
    return length;
}

static inline void
createJobs(JobQueue &jobQueue, const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams, size_t numberOfElements,
        lcp_t baseLcp)
{
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION
    cout << endl << "CREATING JOBS at baseLcp: " << baseLcp << ", numberOfElements: " << numberOfElements << endl;
#else
    (void) numberOfElements;
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION

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
            splitterCharacter[k] = numeric_limits < CHAR_TYPE > ::max();
        }
    }

    const unsigned overProvFactor = 100;
    const size_t expectedJobLength = max(MERGE_BULK_SIZE, numberOfElements / (overProvFactor * numa_num_configured_cpus()));
    cout << "Expected job length: " << expectedJobLength << endl;

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

        if (currBucket == (numeric_limits < CHAR_TYPE > ::max() & keyMask))
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
            unsigned numNewStreams = getNextHigherPowerOfTwo(numberOfFoundBuckets);
            pair < size_t, size_t > *newRange = new pair<size_t, size_t> [numNewStreams];

            unsigned k = 0;
            for (; k < numberOfFoundBuckets; ++k)
            {
                const unsigned idx = indexesOfFound[k];
                const LcpStringPtr inputStart = inputStreams[idx];
                size_t currLength = findNextSplitter(inputStreams[idx], ends[idx], baseLcp, maxAllowedLcp, splitterCharacter[idx], keyMask);
                newRange[k] = make_pair(inputStart - input, currLength);
                length += currLength;
            }
            for (; k < numNewStreams; k++)
            {
                newRange[k] = make_pair(0, 0); // this stream is not used
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
            keyWidth = max(unsigned(1), keyWidth - 1);
            cout << "decreased key to " << keyWidth << "  diff: " << diffExpectedReal << endl;
        }
        else if (diffExpectedReal >= tollerance)
        {
            keyWidth = min(unsigned(8), keyWidth + 1);
            cout << "increased key to " << keyWidth << "  diff: " << diffExpectedReal << endl;
        }
    }

#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION
    cout << "Created " << createdJobsCtr << " Jobs!" << endl;
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION
}

static inline
void
parallelMerge(const LcpStringPtr& input, string* output, pair<size_t, size_t>* ranges, size_t length, unsigned numStreams)
{
    JobQueue jobQueue;
    cout << "doing parallel merge for " << numStreams << " streams" << endl;
    jobQueue.enqueue(new InitialSplitJob(input, output, ranges, length, numStreams));
    jobQueue.loop();
}

void
eberle_parallel_mergesort_lcp_loosertree(string *strings, size_t n)
{
    int realNumaNodes = numa_num_configured_nodes();
    unsigned numNumaNodes = max(unsigned(4), unsigned(realNumaNodes)); // this max ensures a parallel merge on developer machine
    int numThreadsPerPart = numa_num_configured_cpus() / numNumaNodes;

//allocate memory for lcps and temporary strings
    string* shadow = new string[n]; // allocate shadow pointer array
    string* tmp = new string[n];

    pair < size_t, size_t > *ranges = new pair<size_t, size_t> [numNumaNodes];
    calculateRanges(ranges, numNumaNodes, n);

// enable nested parallel regions
    omp_set_nested(true);

#pragma omp parallel for
    for (unsigned k = 0; k < numNumaNodes; k++)
    {
        size_t start = ranges[k].first;
        size_t length = ranges[k].second;

        const StringPtrOut strptr(strings + start, shadow + start, tmp + start, length);
        parallel_sample_sort_numa(strptr, k % realNumaNodes, numThreadsPerPart);
    }

    LcpStringPtr lcpStringPtr(tmp, (lcp_t*) shadow);

#ifdef PARALLEL_LCP_MERGE_DEBUG_TOP_LEVEL_MERGE_DURATION
    MeasureTime<0> timer;
    timer.start();
    parallelMerge(lcpStringPtr, strings, ranges, n, numNumaNodes);
    timer.stop();
    cout << endl << "top level merge needed: " << timer.delta() << " s" << endl << endl;
#else
    parallelMerge(tmp, output, ranges, n, numNumaNodes);
#endif

    delete[] shadow;
    delete[] tmp;
}

CONTESTANT_REGISTER_PARALLEL(eberle_parallel_mergesort_lcp_loosertree, "eberle/ps5-parallel-toplevel-merge-lcpstrptr",
        "NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle")

}
// namespace eberle_parallel_mergesort_lcp_loosertree

#endif // EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_LCPSTRPTR_H_

