/*
 * string-losertree.h
 *
 *  Created on: Oct 30, 2013
 *      Author: aeberle
 */

#ifndef STRING_LOSERTREE_H_
#define STRING_LOSERTREE_H_

#include <stdlib.h>
#include "types.h"

namespace eberle_utils {

using namespace std;

//typedefs
typedef unsigned char* string;
typedef unsigned int UINT;

//implementation follows

static inline
int scmp(string s1, string s2) {
	while (*s1 != '\0' && *s1 == *s2)
		s1++, s2++;
	return (*s1 - *s2);
}

template<unsigned NUMBER_OF_STREAMS>
class StringLoserTree {
	struct STREAM {
		string* stream;
		unsigned length;
		bool isEmpty;
	};

private:
	STREAM streams[NUMBER_OF_STREAMS];
	unsigned nodes[NUMBER_OF_STREAMS];
	unsigned remainingElements;

	inline void initTree() {
		for (unsigned i = 0; i < NUMBER_OF_STREAMS; i++) {
			unsigned parentIdx = NUMBER_OF_STREAMS + i;
			unsigned streamIdx = i;

			while (parentIdx % 2 == 1 && parentIdx > 1) {
				parentIdx = parentIdx / 2;

				const unsigned streamIdxAtParent = nodes[parentIdx];

				if (scmp(streams[streamIdxAtParent].stream[0],
						streams[streamIdx].stream[0]) < 0) {
					nodes[parentIdx] = streamIdx;
					streamIdx = streamIdxAtParent;
				}
			}
			nodes[parentIdx / 2] = streamIdx;
		}
	}

	inline string removeTopFromStream(unsigned streamIdx) {
		STREAM* stream = &(streams[streamIdx]);
		const string top = *(stream->stream);
		stream->length--;
		stream->stream++;
		stream->isEmpty = stream->length <= 0;
		remainingElements--;

		return top;
	}

public:
	StringLoserTree(string * strings, pair<size_t, size_t>* ranges) {

		remainingElements = 0;

		for (unsigned i = 0; i < NUMBER_OF_STREAMS; i++) {
			pair < size_t, size_t > currRange = ranges[i];
			STREAM* curr = &(this->streams[i]);
			curr->stream = strings + currRange.first;
			curr->length = currRange.second;
			curr->isEmpty = (curr->length <= 0);
			remainingElements += currRange.second;
		}

		initTree();
	}

	void printTree() {
		unsigned levelSize = 1;

		for (unsigned i = 0; i < NUMBER_OF_STREAMS; ++i) {
			if (i >= levelSize) {
				cout << "\n";
				levelSize *= 2;
			}
			cout << nodes[i] << " ";
		}
		cout << endl;

		for (unsigned i = 0; i < NUMBER_OF_STREAMS; ++i) {
			STREAM stream = streams[i];
			if (stream.length > 0) {
				cout << stream.stream[0];
			} else {
				cout << -1;
			}
			cout << "(" << stream.length << ")  ";
		}

		cout << endl << endl;
	}

	inline string getMin() {
		return streams[nodes[0]].stream[0];
	}

	inline string deleteMin() {
		unsigned winnerIdx = nodes[0];
		const string min = removeTopFromStream(winnerIdx);

		unsigned nodeIdx = NUMBER_OF_STREAMS + winnerIdx;
		STREAM winnerStream = streams[winnerIdx];

		while (nodeIdx > 1) {
			nodeIdx = nodeIdx / 2;
			unsigned lastLoserIdx = nodes[nodeIdx];
			STREAM lastLoserStream = streams[lastLoserIdx];

			if (winnerStream.isEmpty
					|| (!lastLoserStream.isEmpty
							&& scmp(lastLoserStream.stream[0],
									winnerStream.stream[0]) < 0)) {
				// last loser won this time
				nodes[nodeIdx] = winnerIdx;
				winnerIdx = lastLoserIdx;
				winnerStream = streams[winnerIdx];

			} else {
				// last loser lost again, nothing to do
			}
		}

		nodes[0] = winnerIdx;

		return min;
	}

	inline
	bool isEmpty() {
		return remainingElements <= 0;
	}

	bool isValidTree() {
		string top = streams[nodes[0]].stream[0];

		for (unsigned i = 0; i < NUMBER_OF_STREAMS; i++) {
			string curr = streams[i].stream[0];
			if (scmp(top, curr) > 0) {
				return false;
			}
		}

		return true;
	}

	void writeElementsToStream(string *outStream, const size_t length) {
		const string* end = outStream + length;
		while (outStream < end) {
			*outStream = deleteMin();
			outStream++;
		}
	}
};

}

#endif // STRING_LOSERTREE_H_
