#ifndef PARALLEL_LCP_LOOSERTREE_MERGESORT_H_
#define PARALLEL_LCP_LOOSERTREE_MERGESORT_H_

#include <iostream>
#include <list>
#include <vector>
#include <climits>

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
static const int SHARE_WORK_THRESHOLD = 3 * MERGE_BULK_SIZE;

//method definitions

static inline
void
createJobs(JobQueue& jobQueue, AS* input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams);

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
}
;

template<unsigned K>
    struct MergeJob : public Job
    {
        AS* input;
        string* output;
        pair<size_t, size_t>* ranges;
        size_t length;

        MergeJob(AS* input, string* output, pair<size_t, size_t>* ranges, size_t length) :
                input(input), output(output), ranges(ranges), length(length)
        {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
            {
                cout << "MergeJob<" << K << "> (output: " << output << ", length: " << length << ")" << endl;
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

            const string* end = output + length - MERGE_BULK_SIZE;
            for (; output < end; output += MERGE_BULK_SIZE)
            {
                if (USE_WORK_SHARING && jobQueue.has_idle() && (end - output) > SHARE_WORK_THRESHOLD)
                {
                    return false;
                }

                loserTree->writeElementsToStream(output, MERGE_BULK_SIZE);
            }

            loserTree->writeElementsToStream(output, end - output + MERGE_BULK_SIZE);

            return true;
        }

        virtual bool
        run(JobQueue& jobQueue)
        {
            for (unsigned k = 0; k < K; k++)
            { // this is a temporary fix, because the losertree modifies the lcps.
                input[ranges[k].first].lcp = 0;
            }

            //merge
            LcpStringLoserTree<K> *loserTree = new LcpStringLoserTree<K>(input, ranges);

            if (!mergeToOutput(jobQueue, loserTree))
            {
                // share work
                pair < size_t, size_t > newRanges[K];
                loserTree->getRangesOfRemaining(newRanges, input);
                createJobs(jobQueue, input, output, newRanges, K);
            }

            delete loserTree;

            return true;
        }

        ~MergeJob()
        {
            delete ranges;
        }
    };

typedef uint32_t CHAR_TYPE;
const CHAR_TYPE CHAR_TYPE_MAX = CHAR_TYPE(-1);
const unsigned CHARS_PER_KEY = 4;

//implementations follow
struct Splitter
{
    size_t pos;
    CHAR_TYPE character;

    Splitter(size_t pos, CHAR_TYPE character) :
            pos(pos), character(character)
    {
    }
};

static inline list<Splitter*>*
findSplitters(AS* input, pair<size_t, size_t>* ranges, unsigned numStreams)
{
    //create list to store the minimas in the streams
    list < pair<size_t, unsigned> > minimas[numStreams];
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
                if (lcp < smallestLcp)
                {
                    smallestLcp = lcp;
                }
            }

            lastText = text;
        }
    }

    // find minimas in each stream
    for (unsigned k = 0; k < numStreams; ++k)
    {
        size_t currLength = ranges[k].second;

        if (currLength <= 1)
        { // the stream has no or just one element => no splitters can be found => skip it.
            continue;
        }

        list<pair<size_t, unsigned>>* currList = minimas + k;

        AS* currInput = input + ranges[k].first;
        if (currInput[1].lcp < smallestLcp)
        {
            smallestLcp = currInput[1].lcp;
        }
        char lastCharacter = currInput[0].text[smallestLcp];
        unsigned lastLcp = smallestLcp;

        for (size_t i = 1; i < currLength; ++i)
        {
            unsigned lcp = currInput[i].lcp;

            if (lcp <= smallestLcp + CHARS_PER_KEY - 1)
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

    unsigned maxAllowedLcp = smallestLcp + CHARS_PER_KEY - 1;

#ifdef PARALLEL_LCP_MERGE_DEBUG_SPLITTER_DETECTION
#pragma omp critical (OUTPUT)
    {
        // now we have all minimas of minimaHeight.
        cout << endl << "minLcp: " << smallestLcp << " at positions:" << endl;
        for (unsigned k = 0; k < numStreams; k++)
        {
            AS* currInput = input + ranges[k].first;
            list < pair<size_t, unsigned> > *v = minimas + k;

            for (list<pair<size_t, unsigned> >::iterator it = v->begin(); it != v->end(); ++it)
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

    // calculate result and get distinguishing characters
    list<Splitter*>* splitters = new list<Splitter*> [numStreams];

    for (unsigned k = 0; k < numStreams; ++k)
    {
        if (ranges[k].second > 0)
        {
            list < pair<size_t, unsigned> > *currList = minimas + k;
            AS* currInput = input + ranges[k].first;

            currList->push_front(make_pair(0, smallestLcp));

            for (list<pair<size_t, unsigned> >::iterator it = currList->begin(); it != currList->end(); ++it)
            {
                if (it->second <= maxAllowedLcp)
                {
                    splitters[k].push_back(new Splitter(it->first, get_char<CHAR_TYPE>(currInput[it->first].text, smallestLcp)));
                }
            }
        }

        splitters[k].push_back(new Splitter(ranges[k].second, CHAR_TYPE_MAX));
    }

    return splitters;
}

static inline
void
createJobs(JobQueue& jobQueue, AS* input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams)
{
    list<Splitter*>* splitters = findSplitters(input, ranges, numStreams);

    list<Splitter*>::iterator iterators[numStreams]; // get iterators of minimas of each streams
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
        CHAR_TYPE currBucket = CHAR_TYPE_MAX;

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

        if (currBucket == CHAR_TYPE_MAX)
        { // no bucket < CHAR_MAX found => we reached the end
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

            switch (numStreams)
            {
            case 4:
                jobQueue.enqueue(new MergeJob<4>(input, output, newRange, length));
                break;
            case 8:
                jobQueue.enqueue(new MergeJob<8>(input, output, newRange, length));
                break;
            case 16:
                jobQueue.enqueue(new MergeJob<16>(input, output, newRange, length));
                break;
            case 32:
                jobQueue.enqueue(new MergeJob<32>(input, output, newRange, length));
                break;
            }

        }

        output += length;
    }

    delete[] splitters;
}

static inline
void
parallelMerge(AS* input, string* output, pair<size_t, size_t>* ranges, unsigned numStreams)
{
    JobQueue jobQueue;
    cout << "doing parallel merge for " << numStreams << " streams" << endl;
    createJobs(jobQueue, input, output, ranges, numStreams);
    jobQueue.loop();
}

void
eberle_parallel_mergesort_lcp_loosertree(string *strings, size_t n)
{
    int realNumaNodes = numa_num_configured_nodes();
    unsigned numNumaNodes = max(unsigned(8), unsigned(realNumaNodes)); // this max ensures a parallel merge on developer machine
    int numThreadsPerPart = numa_num_configured_cpus() / numNumaNodes;

//allocate memory for annotated strings
    AS *tmp = static_cast<AS *>(malloc(n * sizeof(AS)));
    string* shadow = new string[n]; // allocate shadow pointer array

    std::pair < size_t, size_t > ranges[numNumaNodes];
    calculateRanges(ranges, numNumaNodes, n);

#pragma omp parallel for
    for (unsigned k = 0; k < numNumaNodes; k++)
    {
        size_t start = ranges[k].first;
        size_t length = ranges[k].second;

        StringPtr strptr(strings + start, shadow + start, length);
        parallel_sample_sort_numa(strptr, k % realNumaNodes, numThreadsPerPart);

        //calculate lcps
        for (size_t pos = 0; pos < length; pos++)
        {
            tmp[start + pos].text = strptr.str(pos);
            tmp[start + pos].lcp = strptr.lcp(pos);
        }
    }

    delete[] shadow;

#ifdef PARALLEL_LCP_MERGE_DEBUG_TOP_LEVEL_MERGE_DURATION
    MeasureTime < 0 > timer;
    timer.start();
    parallelMerge(tmp, strings, ranges, numNumaNodes);
    timer.stop();
    cout << endl << "top level merge needed: " << timer.delta() << " s" << endl << endl;
#else
    parallelMerge(tmp, output, ranges, numNumaNodes);
#endif

    free(tmp);
}

CONTESTANT_REGISTER_PARALLEL(eberle_parallel_mergesort_lcp_loosertree, "eberle/parallel_mergesort_lcp_losertree",
        "Parallel Mergesort with LCP-Losertree-usage by Andreas Eberle")

}
// namespace eberle_parallel_mergesort_lcp_loosertree

#endif // PARALLEL_LCP_LOOSERTREE_MERGESORT_H_

