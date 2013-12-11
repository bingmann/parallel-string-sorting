#ifndef EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_H_
#define EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_H_

#include <iostream>
#include <vector>

#include "../utils/types.h"
#include "../utils/utility-functions.h"
#include "../utils/lcp-string-losertree.h"
#include "../sequential/eberle-lcp-mergesort.h"

#include "../../tools/jobqueue.h"
#include "../../tools/stringtools.h"

#include "../../parallel/bingmann-parallel_sample_sort.h"

//#define PARALLEL_LCP_MERGE_DEBUG_SPLITTER_DETECTION
//#define PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
//#define PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS_DETAILED
//#define PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION
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
typedef unsigned int UINT;

//constants
static const bool USE_WORK_SHARING = true;
static const size_t MERGE_BULK_SIZE = 3000;
static const size_t SHARE_WORK_THRESHOLD = 3 * MERGE_BULK_SIZE;

//method definitions

static inline void
createJobs(JobQueue& jobQueue, AS* input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams);

static inline void
createJobs2(JobQueue &jobQueue, AS* input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams, size_t numberOfElements,
        unsigned baseLcp);

//definitions

typedef uint32_t CHAR_TYPE;

struct Splitter
{
    size_t pos;
    CHAR_TYPE character;

    Splitter(size_t pos, CHAR_TYPE character) :
            pos(pos), character(character)
    {
    }

    Splitter()
    {
    }
};

//structs defining the jobs

struct CopyDataJob : public Job
{
    AS* input;
    string* output;
    size_t length;

    CopyDataJob(AS* input, string* output, size_t length) :
            input(input), output(output), length(length)
    {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
        {
            cout << "CopyDataJob (output: " << output << ", length: " << length << ")" << endl;
        }
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        (void) jobQueue;
        //memcpy(output, input, length * sizeof(AS));

        for (string* end = output + length; output < end; output++, input++)
        {
            *output = input->text;
        }

        return true;
    }
};

struct BinaryMergeJob : public Job
{
    AS* input1;
    size_t length1;
    AS* input2;
    size_t length2;
    string* output;

    BinaryMergeJob(AS* input1, size_t length1, AS* input2, size_t length2, string* output) :
            input1(input1), length1(length1), input2(input2), length2(length2), output(output)
    {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
        {
            cout << "BinaryMergeJob (output: " << output << ", length1: " << length1 << ", length2: " << length2 << ")" << endl;
        }
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        (void) jobQueue;
        input1->lcp = 0;
        input2->lcp = 0;
        eberle_lcp_merge(input1, length1, input2, length2, output);

        return true;
    }
};
size_t lengthOfLongestJob(0);

template<unsigned K>
    struct MergeJob : public Job
    {
        AS* input;
        string* output;
        pair<size_t, size_t>* ranges;
        size_t length;
        unsigned baseLcp;

        MergeJob(AS* input, string* output, pair<size_t, size_t>* ranges, size_t length, unsigned baseLcp) :
                input(input), output(output), ranges(ranges), length(length), baseLcp(baseLcp)
        {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
            {
                cout << "MergeJob<" << K << "> (output: " << output << ", baseLcp: " << baseLcp << ", length: " << length << ")" << endl;
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
        mergeToOutput(JobQueue& jobQueue, LcpStringLoserTree<K> * loserTree)
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
                input[ranges[k].first].lcp = baseLcp;
            }

            //merge
            LcpStringLoserTree<K> *loserTree = new LcpStringLoserTree<K>(input, ranges);

            if (!mergeToOutput(jobQueue, loserTree))
            {
                // share work
                pair < size_t, size_t > newRanges[K];
                loserTree->getRangesOfRemaining(newRanges, input);

                //createJobs(jobQueue, input, output, newRanges, K);
                createJobs2(jobQueue, input, output, newRanges, K, length, baseLcp + key_traits<CHAR_TYPE>::add_depth);

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
    AS* input;
    string* output;
    pair<size_t, size_t>* ranges;
    size_t length;
    unsigned numStreams;
    unsigned splitStartLcp;

    InitialSplitJob(AS* input, string* output, pair<size_t, size_t>* ranges, size_t length, unsigned numStreams, unsigned splitStartLcp) :
            input(input), output(output), ranges(ranges), length(length), numStreams(numStreams), splitStartLcp(splitStartLcp)
    {
        lengthOfLongestJob = length; // prevents that the first MergeJob immediately starts splitting itself
    }

    virtual bool
    run(JobQueue& jobQueue)
    {
        createJobs2(jobQueue, input, output, ranges, numStreams, length, splitStartLcp);
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
enqueueMergeJob(JobQueue& jobQueue, AS* input, string* output, pair<size_t, size_t>* ranges, size_t length, unsigned numStreams, unsigned baseLcp)
{
    switch (numStreams)
    {
    case 2:
        jobQueue.enqueue(new BinaryMergeJob(input + ranges[0].first, ranges[0].second, input + ranges[1].first, ranges[1].second, output));
        break;
    case 4:
        jobQueue.enqueue(new MergeJob<4>(input, output, ranges, length, baseLcp));
        break;
    case 8:
        jobQueue.enqueue(new MergeJob<8>(input, output, ranges, length, baseLcp));
        break;
    case 16:
        jobQueue.enqueue(new MergeJob<16>(input, output, ranges, length, baseLcp));
        break;
    case 32:
        jobQueue.enqueue(new MergeJob<32>(input, output, ranges, length, baseLcp));
        break;
    case 64:
        jobQueue.enqueue(new MergeJob<64>(input, output, ranges, length, baseLcp));
        break;
    default:
        cerr << "CANNOT MERGE! TO MANY STREAMS: " << numStreams;
        break;
    }
}

static inline unsigned
findNextSplitter(AS* &inputStream, AS* end, unsigned baseLcp, unsigned maxAllowedLcp, CHAR_TYPE &lastCharacter)
{
    AS* streamStart = inputStream;
    inputStream++;

    for (; inputStream < end; ++inputStream)
    {
        unsigned lcp = inputStream->lcp;

        if (lcp <= maxAllowedLcp)
        {
            CHAR_TYPE character = get_char<CHAR_TYPE>(inputStream->text, baseLcp);
            if (character != lastCharacter)
            {
                lastCharacter = character;
                return inputStream - streamStart;
            }
        }
    }

    lastCharacter = key_traits<CHAR_TYPE>::maxValue;
    return inputStream - streamStart;
}

static inline void
createJobs2(JobQueue &jobQueue, AS* input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams, size_t numberOfElements,
        unsigned baseLcp)
{
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION
    cout << endl << "CREATING JOBS at baseLcp: " << baseLcp << ", numberOfElements: " << numberOfElements << endl;
#else
    (void) numberOfElements;
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION

    AS* inputStreams[numStreams];
    AS* ends[numStreams];
    CHAR_TYPE splitterCharacter[numStreams];

    for (unsigned k = 0; k < numStreams; ++k)
    {
        if (ranges[k].second > 0)
        {
            inputStreams[k] = input + ranges[k].first;
            ends[k] = inputStreams[k] + ranges[k].second;
            splitterCharacter[k] = get_char<CHAR_TYPE>(inputStreams[k][0].text, baseLcp);
        }
        else
        {
            splitterCharacter[k] = key_traits<CHAR_TYPE>::maxValue;
        }
    }

    unsigned maxAllowedLcp(baseLcp + key_traits<CHAR_TYPE>::add_depth - 1);
    unsigned indexesOfFound[numStreams];

    unsigned createdJobsCtr = 0;

    while (true)
    {
        CHAR_TYPE currBucket = key_traits<CHAR_TYPE>::maxValue;
        unsigned numberOfFoundBuckets = 0;

        for (unsigned k = 0; k < numStreams; ++k)
        {
            if (splitterCharacter[k] <= currBucket)
            {
                if (splitterCharacter[k] < currBucket)
                {
                    currBucket = splitterCharacter[k];

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

        if (currBucket == key_traits<CHAR_TYPE>::maxValue)
        {
            break;
        }

        size_t length = 0;

        switch (numberOfFoundBuckets)
        {
        case 1:
        {
            unsigned streamIdx = indexesOfFound[0];
            AS* inputStart = inputStreams[streamIdx];
            length += findNextSplitter(inputStreams[streamIdx], ends[streamIdx], baseLcp, maxAllowedLcp, splitterCharacter[streamIdx]);
            jobQueue.enqueue(new CopyDataJob(inputStart, output, length));
            break;
        }
        case 2:
        {
            unsigned streamIdx1 = indexesOfFound[0];
            AS* inputStart1 = inputStreams[streamIdx1];
            size_t length1 = findNextSplitter(inputStreams[streamIdx1], ends[streamIdx1], baseLcp, maxAllowedLcp, splitterCharacter[streamIdx1]);

            unsigned streamIdx2 = indexesOfFound[1];
            AS* inputStart2 = inputStreams[streamIdx2];
            size_t length2 = findNextSplitter(inputStreams[streamIdx2], ends[streamIdx2], baseLcp, maxAllowedLcp, splitterCharacter[streamIdx2]);

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
                unsigned idx = indexesOfFound[k];
                AS* inputStart = inputStreams[idx];
                size_t currLength = findNextSplitter(inputStreams[idx], ends[idx], baseLcp, maxAllowedLcp, splitterCharacter[idx]);
                newRange[k] = make_pair(inputStart - input, currLength);
                length += currLength;
            }
            for (; k < numNewStreams; k++)
            {
                newRange[k] = make_pair(0, 0); // this stream is not used
            }

            enqueueMergeJob(jobQueue, input, output, newRange, length, numNewStreams, baseLcp);
            break;
        }
        }

        output += length;
        createdJobsCtr++;
    }

#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION
    cout << "Created " << createdJobsCtr << " Jobs!" << endl;
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_CREATION
}

static inline vector<Splitter*>*
findSplitters(AS* input, pair<size_t, size_t>* ranges, unsigned numStreams)
{
    unsigned smallestLcp = UINT_MAX;

//calculate the minimum lcp between the leading elements
    string lastText = NULL;
    for (unsigned k = 0; k < numStreams; ++k)
    {
        if (ranges[k].second > 0)
        {
            string text = input[ranges[k].first].text;

            if (lastText != NULL)
            {
                unsigned lcp = stringtools::calc_lcp(lastText, text);
                smallestLcp = min(smallestLcp, lcp);
            }

            lastText = text;
        }
    }

//create list to store the position and lcps of the minimum lcps in the streams
    vector < pair<size_t, unsigned> > minimas[numStreams];

// find minimas in each stream
    for (unsigned k = 0; k < numStreams; ++k)
    {
        size_t currLength = ranges[k].second;

        if (currLength <= 1)
        { // the stream has no or just one element => no splitters can be found => skip it.
            continue;
        }

        vector<pair<size_t, unsigned>>* currList = minimas + k;
        AS* currInput = input + ranges[k].first;

        smallestLcp = min(smallestLcp, currInput[1].lcp);
        char lastCharacter = currInput[0].text[smallestLcp];
        unsigned lastLcp = smallestLcp;

        for (size_t i = 1; i < currLength; ++i)
        {
            unsigned lcp = currInput[i].lcp;

            if (lcp <= smallestLcp + key_traits<CHAR_TYPE>::add_depth - 1)
            {
                smallestLcp = min(smallestLcp, lcp);

                char character = currInput[i].text[lcp];

                if (lastLcp != lcp || character != lastCharacter)
                { // prevents multiple insertions for equal strings at start
                    currList->push_back(make_pair(i, lcp));
                    lastCharacter = character;
                    lastLcp = lcp;
                }
            }
        }
    }

// Positions with lcps with at max this value are allowed as splitter points
    unsigned maxAllowedLcp = smallestLcp + key_traits<CHAR_TYPE>::add_depth - 1;

// calculate result and get distinguishing characters
    vector<Splitter*>* splitters = new vector<Splitter*> [numStreams];

    for (unsigned k = 0; k < numStreams; ++k)
    {
        if (ranges[k].second > 0)
        {
            vector < pair<size_t, unsigned> > *currList = minimas + k;
            AS* currInput = input + ranges[k].first;

            splitters[k].push_back(new Splitter(0, get_char<CHAR_TYPE>(currInput[0].text, smallestLcp)));

            for (vector<pair<size_t, unsigned> >::iterator it = currList->begin(); it != currList->end(); ++it)
            {
                if (it->second <= maxAllowedLcp)
                {
                    splitters[k].push_back(new Splitter(it->first, get_char<CHAR_TYPE>(currInput[it->first].text, smallestLcp)));
                }
            }
        }

        splitters[k].push_back(new Splitter(ranges[k].second, key_traits<CHAR_TYPE>::maxValue));
    }

#ifdef PARALLEL_LCP_MERGE_DEBUG_SPLITTER_DETECTION
#pragma omp critical (OUTPUT)
    {
        // now we have all minimas of minimaHeight.
        cout << endl << "minLcp: " << smallestLcp << " at positions:" << endl;
        for (unsigned k = 0; k < numStreams; k++)
        {
            AS* currInput = input + ranges[k].first;
            vector < pair<size_t, unsigned> > *v = minimas + k;

            for (vector<pair<size_t, unsigned> >::iterator it = v->begin(); it != v->end(); ++it)
            {
                cout << it->first << ":" << (it->second - smallestLcp) << "(";
                string text = currInput[it->first].text;
                char lastCharacter = text[smallestLcp];
                for (unsigned lcp = smallestLcp + 1; lcp <= maxAllowedLcp && lastCharacter != 0; lcp++)
                {
                    cout << lastCharacter;
                    lastCharacter = text[lcp];
                }
                cout << lastCharacter << "), ";
            }
            cout << endl;
        }
        cout << endl;
    }
#endif // PARALLEL_LCP_MERGE_DEBUG_SPLITTER_DETECTION

    return splitters;
}

static inline void
createJobs(JobQueue& jobQueue, AS* input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams)
{
    vector<Splitter*>* splitters = findSplitters(input, ranges, numStreams);

    vector<Splitter*>::iterator iterators[numStreams]; // get iterators of minimas of each streams
    for (unsigned k = 0; k < numStreams; ++k)
    {
        iterators[k] = splitters[k].begin();

        if ((*iterators[k])->character == 0)
        {
            size_t start = (*iterators[k])->pos;
            iterators[k]++;
            size_t length = (*iterators[k])->pos - start;

            AS* newInput = input + ranges[k].first + start;
            jobQueue.enqueue(new CopyDataJob(newInput, output, length));
            output += length;
        }
    }

//Go through all found minimas, find the matching buckets and create jobs for them
    while (true)
    {
        // find all matching buckets with the smallest character
        CHAR_TYPE currBucket = key_traits<CHAR_TYPE>::maxValue;

        unsigned numberOfFoundBuckets = 0;
        unsigned indexesOfFound[numStreams];

        for (unsigned k = 0; k < numStreams; ++k)
        {
            if ((*iterators[k])->character <= currBucket)
            {
                if ((*iterators[k])->character < currBucket)
                {
                    currBucket = (*iterators[k])->character;

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

        if (currBucket == key_traits<CHAR_TYPE>::maxValue)
        { // no bucket < maxValue found => we reached the end
            break;
        }

        size_t length = 0; // length of the current job (needed to increase output)

        if (numberOfFoundBuckets == 1)
        { // only one stream => copy instead of merge
            unsigned index = indexesOfFound[0];
            size_t start = (*iterators[index])->pos;
            iterators[index]++;
            length = (*iterators[index])->pos - start;

            AS* newInput = input + ranges[index].first + start;
            jobQueue.enqueue(new CopyDataJob(newInput, output, length));

        }
        else if (numberOfFoundBuckets == 2)
        { // only two steams => binary merge
            unsigned index1 = indexesOfFound[0];
            size_t start1 = (*iterators[index1])->pos;
            iterators[index1]++;
            size_t length1 = (*iterators[index1])->pos - start1;
            AS* newInput1 = input + ranges[index1].first + start1;

            unsigned index2 = indexesOfFound[1];
            size_t start2 = (*iterators[index2])->pos;
            iterators[index2]++;
            size_t length2 = (*iterators[index2])->pos - start2;
            AS* newInput2 = input + ranges[index2].first + start2;

            jobQueue.enqueue(new BinaryMergeJob(newInput1, length1, newInput2, length2, output));

            length = length1 + length2;

        }
        else
        { // default case: use a loosertree to merge the streams.
            unsigned numStreams = getNextHigherPowerOfTwo(numberOfFoundBuckets);
            pair < size_t, size_t > *newRange = new pair<size_t, size_t> [numStreams];

            unsigned k = 0;
            for (; k < numberOfFoundBuckets; ++k)
            {
                unsigned idx = indexesOfFound[k];
                size_t start = (*iterators[idx])->pos;
                iterators[idx]++;
                size_t currLength = (*iterators[idx])->pos - start;
                newRange[k] = make_pair(ranges[idx].first + start, currLength);
                length += currLength;
            }
            for (; k < numStreams; k++)
            {
                newRange[k] = make_pair(0, 0); // this stream is not used
            }

            enqueueMergeJob(jobQueue, input, output, newRange, length, numStreams, 0);

        }

        output += length;
    }

    delete[] splitters;
}

static inline
void
parallelMerge(AS* input, string* output, pair<size_t, size_t>* ranges, size_t length, unsigned numStreams)
{
    JobQueue jobQueue;
    cout << "doing parallel merge for " << numStreams << " streams" << endl;
    jobQueue.enqueue(new InitialSplitJob(input, output, ranges, length, numStreams, unsigned(0)));
    jobQueue.loop();
}

void
eberle_parallel_mergesort_lcp_loosertree(string *strings, size_t n)
{
    int realNumaNodes = numa_num_configured_nodes();
    unsigned numNumaNodes = max(unsigned(4), unsigned(realNumaNodes)); // this max ensures a parallel merge on developer machine
    int numThreadsPerPart = numa_num_configured_cpus() / numNumaNodes;

//allocate memory for annotated strings
    AS *tmp = static_cast<AS *>(malloc(n * sizeof(AS)));
    string* shadow = new string[n]; // allocate shadow pointer array

    pair < size_t, size_t > *ranges = new pair<size_t, size_t> [numNumaNodes];
    calculateRanges(ranges, numNumaNodes, n);

// enable nested parallel regions
    omp_set_nested(true);

#pragma omp parallel for
    for (unsigned k = 0; k < numNumaNodes; k++)
    {
        size_t start = ranges[k].first;
        size_t length = ranges[k].second;

        StringPtr strptr(strings + start, shadow + start, length);
        parallel_sample_sort_numa(strptr, k % realNumaNodes, numThreadsPerPart);

        //create AS* array
        MeasureTime < 0 > timer;
        timer.start();

        for (size_t pos = 0; pos < length; pos++)
        {
            tmp[start + pos].text = strptr.str(pos);
            tmp[start + pos].lcp = strptr.lcp(pos);
        }

        timer.stop();
        cout << endl << "Creating AS* needed: " << timer.delta() << " s" << endl << endl;
    }

    delete[] shadow;

#ifdef PARALLEL_LCP_MERGE_DEBUG_TOP_LEVEL_MERGE_DURATION
    MeasureTime < 0 > timer;
    timer.start();
    parallelMerge(tmp, strings, ranges, n, numNumaNodes);
    timer.stop();
    cout << endl << "top level merge needed: " << timer.delta() << " s" << endl << endl;
#else
    parallelMerge(tmp, output, ranges, n, numNumaNodes);
#endif

    free(tmp);
}

CONTESTANT_REGISTER_PARALLEL(eberle_parallel_mergesort_lcp_loosertree, "eberle/ps5-parallel-toplevel-merge",
        "NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle")

}
// namespace eberle_parallel_mergesort_lcp_loosertree

#endif // EBERLE_PS5_PARALLEL_TOPLEVEL_MERGE_H_

