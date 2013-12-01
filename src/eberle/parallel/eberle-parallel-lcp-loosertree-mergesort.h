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

#define PARALLEL_LCP_MERGE_DEBUG_MINIMA_DETECTION
//#define PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS
#define PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION

namespace eberle_parallel_mergesort_lcp_loosertree {

using namespace std;

using namespace types;
using namespace eberle_lcp_utils;
using namespace eberle_mergesort_lcp;

using namespace jobqueue;

using namespace bingmann_parallel_sample_sort;

//typedefs
typedef unsigned char* string;
typedef unsigned int UINT;

//constants
static const size_t MERGE_BULK_SIZE = 50;
static const bool USE_WORK_SHARING = true;
static const int SHARE_WORK_THRESHOLD = 3 * MERGE_BULK_SIZE;

//method definitionstemplate<unsigned K>
template<unsigned K>
static inline
void parallelMerge(AS* input, AS* output,
		boost::array<std::pair<size_t, size_t>, K> &ranges);

struct CopyDataJob: public Job {
	AS* input;
	AS* output;
	size_t length;

	CopyDataJob(AS* input, AS* output, size_t length) :
			input(input), output(output), length(length) {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
		cout << "CopyDataJob (length: " << length << ")" << endl;
#endif
	}

	virtual void run(JobQueue& jobQueue) {
		(void) jobQueue;
		memcpy(output, input, length * sizeof(AS));
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
		cout << "BinaryMergeJob (length1: " << length1 << ", length2: "
				<< length2 << ")" << endl;
#endif
	}

	virtual void run(JobQueue& jobQueue) {
		(void) jobQueue;
		eberle_lcp_merge(input1, length1, input2, length2, output);
	}
};

template<unsigned K>
struct MergeJob: public Job {
	AS* input;
	AS* output;
	boost::array<pair<size_t, size_t>, K>* ranges;
	size_t length;

	MergeJob(AS* input, AS* output,
			boost::array<pair<size_t, size_t>, K>* ranges, size_t length) :
			input(input), output(output), ranges(ranges), length(length) {
#ifdef PARALLEL_LCP_MERGE_DEBUG_JOB_TYPE_ON_CREATION
		cout << "MergeJob<" << K << "> (length: " << length << ")" << endl;
#endif

#ifdef PARALLEL_LCP_MERGE_DEBUG_MERGE_JOBS
		cout << "MergeJob: in: " << input << "  out: " << output << "  length: "
		<< length << endl;
		for (unsigned k = 0; k < K; ++k) {
			cout << k << ": " << (*ranges)[k].first << " length: "
			<< (*ranges)[k].second << endl;
		}
		cout << endl;
#endif
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

	virtual void run(JobQueue& jobQueue) {
		for (unsigned k = 0; k < K; k++) {
			input[(*ranges)[k].first].lcp = 0;
		}

		//merge
		LcpStringLoserTree<K> *loserTree = new LcpStringLoserTree<K>(input,
				ranges);

		if (!mergeToOutput(jobQueue, loserTree)) {
			// share work
			boost::array<pair<size_t, size_t>, K> newRanges;
			loserTree->getRangesOfRemaining(newRanges, input);
			parallelMerge<K>(input, output, newRanges);
		}

		delete loserTree;
	}

	~MergeJob() {
		delete ranges;
	}
};

//implementations follow

template<unsigned K>
static inline list<pair<size_t, char> > ** findMinimas(AS* input,
		boost::array<std::pair<size_t, size_t>, K> &ranges) {
	//find minimas in first stream

	list < pair<size_t, char> > **minimas = new list<pair<size_t, char> >*[4];
	unsigned minimaHeights[K];
	unsigned currMinimaHeight = input[ranges[0].first + 1].lcp;

	// find minimas in each stream
	for (unsigned k = 0; k < K; k++) {
		list < pair<size_t, char> > *currList = new list<pair<size_t, char> >();
		AS* currInput = input + ranges[k].first;
		const size_t currLength = ranges[k].second;

		for (size_t i = 1; i < currLength; i++) {
			unsigned lcp = currInput[i].lcp;
			if (lcp <= currMinimaHeight) {
				if (lcp < currMinimaHeight) {
					currList->clear();
					currMinimaHeight = lcp;
				}

				currList->push_back(make_pair(i, currInput[i].text[lcp]));
			}
		}

		minimaHeights[k] = currMinimaHeight;
		minimas[k] = currList;
	}

	// ensure that all minimas of the different streams are of the same lcp-height and add start and end
	for (unsigned k = 0; k < K; k++) {
		if (minimaHeights[k] > currMinimaHeight) { // remove higher minimas
			minimas[k]->clear();
		}
		//add start and end to the lists
		minimas[k]->push_front(
				make_pair(0, input[ranges[k].first].text[currMinimaHeight]));
		minimas[k]->push_back(make_pair(ranges[k].second, CHAR_MAX));
	}

#ifdef PARALLEL_LCP_MERGE_DEBUG_MINIMA_DETECTION
	// now we have all minimas of minimaHeight.
	cout << endl << "minimaHeight: " << currMinimaHeight << " at positions:"
			<< endl;
	for (unsigned k = 0; k < K; k++) {
		list < pair<size_t, char> > *v = minimas[k];
		for (list<pair<size_t, char> >::iterator it = v->begin();
				it != v->end(); ++it) {
			cout << it->first << ":" << it->second << "(" << "), ";
		}
		cout << endl;
	}
	cout << endl;
#endif

	return minimas;
}

template<unsigned K>
static inline
void createJobs(JobQueue& jobQueue, AS* input, AS* output,
		boost::array<pair<size_t, size_t>, K> &ranges) {
	list < pair<size_t, char> > **minimas = findMinimas<K>(input, ranges);

	list<pair<size_t, char> >::iterator iterators[K]; // get iterators of minimas of each streams
	for (unsigned k = 0; k < K; ++k) {
		iterators[k] = minimas[k]->begin();
	}

	//Go through all found minimas, find the matching buckets and create jobs for them
	while (true) {
		// find all matching buckets with the smallest character
		char currBucket = CHAR_MAX;
		unsigned numberOfFoundBuckets = 0;
		boost::array<unsigned, K> indexesOfFound;

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
			boost::array<pair<size_t, size_t>, K> *newRange = new boost::array<
					std::pair<size_t, size_t>, K>;

			for (unsigned k = 0; k < K; ++k) {
				if (iterators[k] != minimas[k]->end()
						&& iterators[k]->second == currBucket) {
					size_t start = iterators[k]->first;
					iterators[k]++;
					size_t currLength = iterators[k]->first - start;
					(*newRange)[k] = make_pair(ranges[k].first + start,
							currLength);
					length += currLength;
				} else {
					(*newRange)[k] = make_pair(0, 0); // this stream is not used
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
void parallelMerge(AS* input, AS* output,
		boost::array<std::pair<size_t, size_t>, K> &ranges) {
	JobQueue jobQueue;
	createJobs<K>(jobQueue, input, output, ranges);

	jobQueue.loop();
}

static inline unsigned calculateLcp(string s1, string s2) {
	unsigned lcp = 0;
	while (*s1 != '\0' && *s1 == *s2)
		s1++, s2++, lcp++;

	return lcp;
}

void eberle_parallel_mergesort_lcp_loosertree(string *strings, size_t n) {
	const unsigned K = 4;

	int numNumaNodes = numa_num_configured_nodes();
	int numThreadsPerPart = numa_num_configured_cpus() / K;

	//allocate memory for annotated strings
	AS *tmp = static_cast<AS *>(malloc(n * sizeof(AS)));
	AS *output = static_cast<AS *>(malloc(n * sizeof(AS)));

	boost::array<std::pair<size_t, size_t>, K> ranges;
	eberle_utils::calculateRanges<K>(ranges, n);

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
			tmp[pos].lcp = calculateLcp(strings[pos - 1], strings[pos]);
		}
	}

	MeasureTime < 0 > timer;
	timer.start();

	parallelMerge<K>(tmp, output, ranges);

	timer.stop();
	cout << "top level merge needed: " << timer.delta() << "s" << endl;

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

