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

#include "../../parallel/bingmann-parallel_sample_sort.h"

//#define PARALLEL_LCP_MERGE_DEBUG_MINIMA_DETECTION
//#define PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS
//#define PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#define PARALLEL_LCP_MERGE_DEBUG_TOP_LEVEL_MERGE_DURATION

namespace eberle_parallel_mergesort_lcp_loosertree {

using namespace std;

using namespace types;
using namespace eberle_lcp_utils;
using namespace eberle_utils;
using namespace eberle_mergesort_lcp;

using namespace jobqueue;
using namespace stringtools;

using namespace bingmann_parallel_sample_sort;

//typedefs
typedef unsigned char* string;
typedef unsigned int UINT;

//constants
static const bool USE_WORK_SHARING = true;
static const size_t MERGE_BULK_SIZE = 50;
static const int SHARE_WORK_THRESHOLD = 5 * MERGE_BULK_SIZE;

//method definitions

template<unsigned K>
static inline
void createJobs(JobQueue& jobQueue, AS* input, AS* output,
		pair<size_t, size_t>* ranges);

//structs defining the jobs

struct CopyDataJob: public Job {
	AS* input;
	AS* output;
	size_t length;

	CopyDataJob(AS* input, AS* output, size_t length) :
			input(input), output(output), length(length) {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
		{
			cout << "CopyDataJob (output: " << output
			<< ", length: " << length << ")" << endl;
		}
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
	}

	virtual bool run(JobQueue& jobQueue) {
		(void) jobQueue;
		memcpy(output, input, length * sizeof(AS));
                return true;
	}
};

struct BinaryMergeJob: public Job {
	AS* input1;
	size_t length1;
	AS* input2;
	size_t length2;
	AS* output;

	BinaryMergeJob(AS* input1, size_t length1, AS* input2, size_t length2,
			AS* output) :
			input1(input1), length1(length1), input2(input2), length2(length2), output(
					output) {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#pragma omp critical (OUTPUT)
		{
			cout << "BinaryMergeJob (output: " << output
			<< ", length1: " << length1 << ", length2: " << length2
			<< ")" << endl;
		}
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
	}

	virtual bool run(JobQueue& jobQueue) {
		(void) jobQueue;
		input1->lcp = 0;
		input2->lcp = 0;
		eberle_lcp_merge(input1, length1, input2, length2, output);
                return true;
	}
}
;

template<unsigned K>
struct MergeJob: public Job {
	AS* input;
	AS* output;
	pair<size_t, size_t>* ranges;
	size_t length;

	MergeJob(AS* input, AS* output, pair<size_t, size_t>* ranges, size_t length) :
			input(input), output(output), ranges(ranges), length(length) {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
#ifndef PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS
#pragma omp critical (OUTPUT)
		{
			cout << "MergeJob<" << K << "> (output: " << output << ", length: " << length << ")" << endl;
		}
#endif // PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS
#endif // PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION

#ifdef PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS
#pragma omp critical (OUTPUT)
		{
			cout << "MergeJob<" << K << "> (output: " << output << ",  length: " << length << endl;
			for (unsigned k = 0; k < K; ++k) {
				cout << k << ": " << ranges[k].first << " length: "
				<< ranges[k].second << endl;
			}
			cout << endl;
		}
#endif // PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS
	}

	/*
	 * returns true if all elements have been written to output
	 * false if the merge has been stopped to free work.
	 */
	inline bool mergeToOutput(JobQueue& jobQueue,
			LcpStringLoserTree<K> * loserTree) {

		const AS* end = output + length - MERGE_BULK_SIZE;
		for (; output < end; output += MERGE_BULK_SIZE) {
			if (USE_WORK_SHARING && jobQueue.has_idle()
					&& (end - output) > SHARE_WORK_THRESHOLD) {
				return false;
			}

			loserTree->writeElementsToStream(output, MERGE_BULK_SIZE);
		}

		loserTree->writeElementsToStream(output,
				end - output + MERGE_BULK_SIZE);

		return true;
	}

	virtual bool run(JobQueue& jobQueue) {
		for (unsigned k = 0; k < K; k++) { // this is a temporary fix, because the losertree modifies the lcps.
			input[ranges[k].first].lcp = 0;
		}

		//merge
		LcpStringLoserTree<K> *loserTree = new LcpStringLoserTree<K>(input,
				ranges);

		if (!mergeToOutput(jobQueue, loserTree)) {
			// share work
			pair < size_t, size_t > newRanges[K];
			loserTree->getRangesOfRemaining(newRanges, input);
			createJobs<K>(jobQueue, input, output, newRanges);
		}

		delete loserTree;

                return true;
	}

	~MergeJob() {
		delete ranges;
	}
};

//implementations follow

template<unsigned K>
static inline list<pair<size_t, char> > ** findMinimas(AS* input,
		pair<size_t, size_t>* ranges) {
//create list to store the minimas in the streams
	list < pair<size_t, char> > **minimas = new list<pair<size_t, char> >*[4];
	unsigned minimaHeights[K];
	unsigned minLcp = UINT_MAX;

// find minimas in each stream
	for (unsigned k = 0; k < K; ++k) {
		list<pair<size_t, char>> * currList = new list<pair<size_t, char>>();
		size_t currLength = ranges[k].second;
		AS* currInput = input + ranges[k].first;

		if (currLength <= 1) { // the stream has no or just one element => no splitters can be found => skip it.
			minimas[k] = currList;
			continue;
		}

		if (currInput[1].lcp < minLcp) {
			minLcp = currInput[1].lcp;
		}
		char lastCharacter = currInput[0].text[minLcp];

		for (size_t i = 1; i < currLength; ++i) {
			unsigned lcp = currInput[i].lcp;

			if (lcp <= minLcp) {
				char character = currInput[i].text[lcp];

				if (lcp < minLcp) {
					currList->clear();
					minLcp = lcp;

					currList->push_back(make_pair(i, character)); // minLcp has been reduced => insert split point
					lastCharacter = character;

				} else if (character != lastCharacter) { // prevents multiple insertions for equal strings at start

					currList->push_back(make_pair(i, character));
					lastCharacter = character;
				}
			}
		}

		minimaHeights[k] = minLcp;
		minimas[k] = currList;
	}

	string lastText = NULL;
	for (unsigned k = 0; k < K; ++k) {
		if (ranges[k].second > 0) {
			string text = input[ranges[k].first].text;

			if (lastText != NULL) {
				unsigned lcp = calc_lcp(lastText, text);
				if (lcp < minLcp) {
					minLcp = lcp;
				}
			}

			lastText = text;
		}
	}

// ensure that all minimas of the different streams are of the same lcp-height and add start and end
	for (unsigned k = 0; k < K; ++k) {
		if (ranges[k].second <= 0) {
			continue;
		}

		if (minimaHeights[k] > minLcp) {
			minimas[k]->clear();
		}

		minimas[k]->push_front(
				make_pair(0, input[ranges[k].first].text[minLcp]));
		minimas[k]->push_back(make_pair(ranges[k].second, CHAR_MAX));
	}

#ifdef PARALLEL_LCP_MERGE_DEBUG_MINIMA_DETECTION

#pragma omp critical (OUTPUT)
	{
		// now we have all minimas of minimaHeight.
		cout << endl << "minLcp: " << minLcp << " at positions:" << endl;
		for (unsigned k = 0; k < K; k++) {
			list < pair<size_t, char> > *v = minimas[k];
			for (list<pair<size_t, char> >::iterator it = v->begin();
					it != v->end(); ++it) {
				cout << it->first << ":" << it->second << ", ";
			}
			cout << endl;
		}
		cout << endl;
	}
#endif // PARALLEL_LCP_MERGE_DEBUG_MINIMA_DETECTION

	return minimas;
}

template<unsigned K>
static inline
void createJobs(JobQueue& jobQueue, AS* input, AS* output,
		pair<size_t, size_t>* ranges) {
	list < pair<size_t, char> > **minimas = findMinimas<K>(input, ranges);

	list<pair<size_t, char> >::iterator iterators[K]; // get iterators of minimas of each streams
	for (unsigned k = 0; k < K; ++k) {
		iterators[k] = minimas[k]->begin();

		if (iterators[k]->second == 0) {
			size_t start = iterators[k]->first;
			iterators[k]++;
			size_t length = iterators[k]->first - start;

			AS* newInput = input + ranges[k].first + start;
			jobQueue.enqueue(new CopyDataJob(newInput, output, length));
			output += length;
		}
	}

//Go through all found minimas, find the matching buckets and create jobs for them
	while (true) {
		// find all matching buckets with the smallest character
		char currBucket = CHAR_MAX;
		unsigned numberOfFoundBuckets = 0;
		unsigned indexesOfFound[K];

		for (unsigned k = 0; k < K; ++k) {
			if (iterators[k] != minimas[k]->end()
					&& iterators[k]->second <= currBucket) {
				if (iterators[k]->second < currBucket) {
					currBucket = iterators[k]->second;

					indexesOfFound[0] = k;
					numberOfFoundBuckets = 1;
				} else {
					indexesOfFound[numberOfFoundBuckets] = k;
					numberOfFoundBuckets++;
				}
			}
		}

		if (currBucket == CHAR_MAX) { // no bucket < CHAR_MAX found => we reached the end
			break;
		}

		size_t length = 0; // length of the current job (needed to increase output)

		if (numberOfFoundBuckets == 1) { // only one stream => copy instead of merge
			unsigned index = indexesOfFound[0];
			size_t start = iterators[index]->first;
			iterators[index]++;
			length = iterators[index]->first - start;

			AS* newInput = input + ranges[index].first + start;
			jobQueue.enqueue(new CopyDataJob(newInput, output, length));

		} else if (numberOfFoundBuckets == 2) { // only two steams => binary merge
			unsigned index1 = indexesOfFound[0];
			size_t start1 = iterators[index1]->first;
			iterators[index1]++;
			size_t length1 = iterators[index1]->first - start1;
			AS* newInput1 = input + ranges[index1].first + start1;

			unsigned index2 = indexesOfFound[1];
			size_t start2 = iterators[index2]->first;
			iterators[index2]++;
			size_t length2 = iterators[index2]->first - start2;
			AS* newInput2 = input + ranges[index2].first + start2;

			jobQueue.enqueue(
					new BinaryMergeJob(newInput1, length1, newInput2, length2,
							output));

			length = length1 + length2;

		} else { // default case: use a loosertree to merge the streams.
			pair < size_t, size_t > *newRange = new pair<size_t, size_t> [K];

			for (unsigned k = 0; k < K; ++k) {
				if (iterators[k] != minimas[k]->end()
						&& iterators[k]->second == currBucket) {
					size_t start = iterators[k]->first;
					iterators[k]++;
					size_t currLength = iterators[k]->first - start;
					newRange[k] = make_pair(ranges[k].first + start,
							currLength);
					length += currLength;
				} else {
					newRange[k] = make_pair(0, 0); // this stream is not used
				}
			}

			jobQueue.enqueue(new MergeJob<K>(input, output, newRange, length));
		}

		output += length;
	}

	delete[] minimas;
}

template<unsigned K>
static inline
void parallelMerge(AS* input, AS* output, pair<size_t, size_t>* ranges) {
	JobQueue jobQueue;
	createJobs<K>(jobQueue, input, output, ranges);
	jobQueue.loop();
}

void eberle_parallel_mergesort_lcp_loosertree(string *strings, size_t n) {
	const unsigned K = 4;

	int numNumaNodes = numa_num_configured_nodes();
	int numThreadsPerPart = numa_num_configured_cpus() / K;

//allocate memory for annotated strings
	AS *tmp = static_cast<AS *>(malloc(n * sizeof(AS)));
	AS *output = static_cast<AS *>(malloc(n * sizeof(AS)));

	std::pair < size_t, size_t > ranges[K];
	calculateRanges(ranges, K, n);

#pragma omp parallel for
	for (unsigned k = 0; k < K; k++) {
		size_t start = ranges[k].first;
		size_t length = ranges[k].second;

		parallel_sample_sort_numa(strings + start, length, k % numNumaNodes,
				numThreadsPerPart);

		//calculate lcps
		size_t end = start + length;
		tmp[start].text = strings[start];
		tmp[start].lcp = 0;
		for (size_t pos = start + 1; pos < end; pos++) {
			tmp[pos].text = strings[pos];
			tmp[pos].lcp = calc_lcp(strings[pos - 1], strings[pos]);
		}
	}

#ifdef PARALLEL_LCP_MERGE_DEBUG_TOP_LEVEL_MERGE_DURATION
	MeasureTime < 0 > timer;
	timer.start();
	parallelMerge<K>(tmp, output, ranges);
	timer.stop();
	cout << endl << "top level merge needed: " << timer.delta() << " s" << endl
			<< endl;
#else
	parallelMerge<K>(tmp, output, ranges);
#endif

	for (size_t i = 0; i < n; i++) {
		strings[i] = output[i].text;
	}

	free(tmp);
	free(output);
}

CONTESTANT_REGISTER_PARALLEL(eberle_parallel_mergesort_lcp_loosertree,
		"eberle/parallel_mergesort_lcp_loosertree",
		"Parallel Mergesort with LCP-usage by Andreas Eberle")

}
// namespace eberle_parallel_mergesort_lcp_loosertree

#endif // PARALLEL_LCP_LOOSERTREE_MERGESORT_H_

