/*
 * Copyright 2008 by Tommi Rantala <tt.rantala@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/*
 * A variant of the Multi-Key-Quicksort using dynamic arrays to store the three
 * buckets.
 */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include <vector>

#include <boost/array.hpp>

#include "tools/debug.h"
#include "tools/get_char.h"
#include "tools/median.h"
#include "tools/vector_malloc.h"
#include "tools/vector_realloc.h"
#include "tools/vector_block.h"
#include "tools/vector_bagwell.h"
#include "tools/vector_brodnik.h"

#include "../tools/contest.h"
#include "../sequential/bs-mkqs.h"

namespace rantala {

static inline void
copy(const std::vector<unsigned char*>& bucket, unsigned char** dst)
{
	std::copy(bucket.begin(), bucket.end(), dst);
}

template <typename CharT>
static inline unsigned
get_bucket(CharT c, CharT pivot)
{
        return ((c > pivot) << 1) | (c == pivot);
}

template <typename BucketT>
static inline void
clear_bucket(BucketT& bucket)
{ bucket.clear(); }

static inline void
clear_bucket(std::vector<unsigned char*>& bucket)
{ bucket.clear(); std::vector<unsigned char*>().swap(bucket); }

template <typename BucketT, typename CharT>
static void
multikey_dynamic(unsigned char** strings, size_t N, size_t depth)
{
	if (N < 10000) {
            mkqsort(strings, N, depth);
            return;
	}
	boost::array<BucketT, 3> buckets;
	CharT partval = pseudo_median<CharT>(strings, N, depth);
	// Use a small cache to reduce memory stalls.
	size_t i=0;
	for (; i < N-N%32; i+=32) {
		boost::array<CharT, 32> cache;
		for (unsigned j=0; j<32; ++j) {
			cache[j] = get_char<CharT>(strings[i+j], depth);
		}
		for (unsigned j=0; j<32; ++j) {
			const unsigned b = get_bucket(cache[j], partval);
			buckets[b].push_back(strings[i+j]);
		}
	}
	for (; i < N; ++i) {
		const CharT c = get_char<CharT>(strings[i], depth);
		const unsigned b = get_bucket(c, partval);
		buckets[b].push_back(strings[i]);
	}
	const size_t bucketsize0 = buckets[0].size();
	const size_t bucketsize1 = buckets[1].size();
	const size_t bucketsize2 = buckets[2].size();
	assert(bucketsize0 + bucketsize1 + bucketsize2 == N);
	if (bucketsize0) copy(buckets[0], strings);
	if (bucketsize1) copy(buckets[1], strings+bucketsize0);
	if (bucketsize2) copy(buckets[2], strings+bucketsize0+bucketsize1);
	clear_bucket(buckets[0]);
	clear_bucket(buckets[1]);
	clear_bucket(buckets[2]);
	multikey_dynamic<BucketT, CharT>(strings, bucketsize0, depth);
	if (not is_end(partval))
		multikey_dynamic<BucketT, CharT>(strings+bucketsize0,
				bucketsize1, depth+sizeof(CharT));
	multikey_dynamic<BucketT, CharT>(strings+bucketsize0+bucketsize1,
			bucketsize2, depth);
}

void multikey_dynamic_vector1(unsigned char** strings, size_t n)
{
	typedef std::vector<unsigned char*> BucketT;
	typedef unsigned char CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}
void multikey_dynamic_vector2(unsigned char** strings, size_t n)
{
	typedef std::vector<unsigned char*> BucketT;
	typedef uint16_t CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}
void multikey_dynamic_vector4(unsigned char** strings, size_t n)
{
	typedef std::vector<unsigned char*> BucketT;
	typedef uint32_t CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}

void multikey_dynamic_brodnik1(unsigned char** strings, size_t n)
{
	typedef vector_brodnik<unsigned char*> BucketT;
	typedef unsigned char CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}
void multikey_dynamic_brodnik2(unsigned char** strings, size_t n)
{
	typedef vector_brodnik<unsigned char*> BucketT;
	typedef uint16_t CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}
void multikey_dynamic_brodnik4(unsigned char** strings, size_t n)
{
	typedef vector_brodnik<unsigned char*> BucketT;
	typedef uint32_t CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}

void multikey_dynamic_bagwell1(unsigned char** strings, size_t n)
{
	typedef vector_bagwell<unsigned char*> BucketT;
	typedef unsigned char CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}
void multikey_dynamic_bagwell2(unsigned char** strings, size_t n)
{
	typedef vector_bagwell<unsigned char*> BucketT;
	typedef uint16_t CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}
void multikey_dynamic_bagwell4(unsigned char** strings, size_t n)
{
	typedef vector_bagwell<unsigned char*> BucketT;
	typedef uint32_t CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}

void multikey_dynamic_vector_block1(unsigned char** strings, size_t n)
{
	typedef vector_block<unsigned char*> BucketT;
	typedef unsigned char CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}
void multikey_dynamic_vector_block2(unsigned char** strings, size_t n)
{
	typedef vector_block<unsigned char*> BucketT;
	typedef uint16_t CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}
void multikey_dynamic_vector_block4(unsigned char** strings, size_t n)
{
	typedef vector_block<unsigned char*> BucketT;
	typedef uint32_t CharT;
	multikey_dynamic<BucketT, CharT>(strings, n, 0);
}

CONTESTANT_REGISTER(multikey_dynamic_vector1,
                    "rantala/multikey_dynamic_vector1",
                    "multikey_dynamic with std::vector bucket type and 1byte alphabet")
CONTESTANT_REGISTER(multikey_dynamic_vector2,
                    "rantala/multikey_dynamic_vector2",
                    "multikey_dynamic with std::vector bucket type and 2byte alphabet")
CONTESTANT_REGISTER(multikey_dynamic_vector4,
                    "rantala/multikey_dynamic_vector4",
                    "multikey_dynamic with std::vector bucket type and 4byte alphabet")

CONTESTANT_REGISTER(multikey_dynamic_brodnik1,
                    "rantala/multikey_dynamic_brodnik1",
                    "multikey_dynamic with vector_brodnik bucket type and 1byte alphabet")
CONTESTANT_REGISTER(multikey_dynamic_brodnik2,
                    "rantala/multikey_dynamic_brodnik2",
                    "multikey_dynamic with vector_brodnik bucket type and 2byte alphabet")
CONTESTANT_REGISTER(multikey_dynamic_brodnik4,
                    "rantala/multikey_dynamic_brodnik4",
                    "multikey_dynamic with vector_brodnik bucket type and 4byte alphabet")

CONTESTANT_REGISTER(multikey_dynamic_bagwell1,
                    "rantala/multikey_dynamic_bagwell1",
                    "multikey_dynamic with vector_bagwell bucket type and 1byte alphabet")
CONTESTANT_REGISTER(multikey_dynamic_bagwell2,
                    "rantala/multikey_dynamic_bagwell2",
                    "multikey_dynamic with vector_bagwell bucket type and 2byte alphabet")
CONTESTANT_REGISTER(multikey_dynamic_bagwell4,
                    "rantala/multikey_dynamic_bagwell4",
                    "multikey_dynamic with vector_bagwell bucket type and 4byte alphabet")

CONTESTANT_REGISTER(multikey_dynamic_vector_block1,
                    "rantala/multikey_dynamic_vector_block1",
                    "multikey_dynamic with vector_block bucket type and 1byte alphabet")
CONTESTANT_REGISTER(multikey_dynamic_vector_block2,
                    "rantala/multikey_dynamic_vector_block2",
                    "multikey_dynamic with vector_block bucket type and 2byte alphabet")
CONTESTANT_REGISTER(multikey_dynamic_vector_block4,
                    "rantala/multikey_dynamic_vector_block4",
                    "multikey_dynamic with vector_block bucket type and 4byte alphabet")

} // namespace rantala
