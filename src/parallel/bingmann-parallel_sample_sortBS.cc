/******************************************************************************
 * src/parallel/bingmann-parallel_sample_sortBS.h
 *
 * Parallel Super Scalar String Sample-Sort with work-balancing, variant 1.
 *
 ******************************************************************************
 * Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
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

#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <vector>
#include <algorithm>

#include "../tools/debug.h"
#include "../tools/lcgrandom.h"
#include "../tools/contest.h"
#include "../tools/stringtools.h"
#include "../tools/jobqueue.h"

namespace bingmann_parallel_radix_sort3 {

using namespace stringtools;
using namespace jobqueue;

// prototype to call radix sort for small buckets
extern void EnqueueSmall(JobQueue& jobqueue, string* strings, size_t n, size_t depth);

}

namespace bingmann_parallel_sample_sortBS {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = false;

static const bool debug_splitter = false;
static const bool debug_bucketsize = false;
static const bool debug_recursion = false;

static const bool use_work_stealing = true;

/// Prototype called to schedule deeper sorts
void Enqueue(JobQueue& jobqueue, string* strings, size_t n, size_t depth);

// ****************************************************************************
// *** SampleSortStep out-of-place parallel sample sort with separate Jobs

typedef uint64_t key_type;

struct SampleSortStep
{
#if 0
    static const size_t leaves = 2*16;       // TODO: calculate to match L2 cache
#else
    static const size_t l2cache = 256*1024;

// bounding equations:
// splitters            + bktsize
// n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    static const size_t leaves = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));

//std::cout << "leaves: " << leaves << "\n";

#endif

    static const size_t bktnum = 2*leaves+1;

    string*             strings;
    size_t              n;
    size_t              depth;
    
    unsigned int        parts;
    size_t              psize;
    std::atomic<unsigned int> pwork;
    key_type            splitter[leaves];
    unsigned char       splitter_lcp[leaves];
    size_t*             bkt;

    string*             sorted;

    void sample(JobQueue& jobqueue);

    void count(unsigned int p, JobQueue& jobqueue);
    void count_finished(JobQueue& jobqueue);

    void distribute(unsigned int p, JobQueue& jobqueue);
    void distribute_finished(JobQueue& jobqueue);

    void copyback(unsigned int p, JobQueue& jobqueue);
    void copyback_finished(JobQueue& jobqueue);
};

struct SampleJob : public Job
{
    SampleSortStep*   step;

    SampleJob(SampleSortStep* _step)
        : step(_step) { }

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process SampleJob " << step);
        step->sample(jobqueue);
    }
};

struct CountJob : public Job
{
    SampleSortStep*   step;
    unsigned int        p;

    CountJob(SampleSortStep* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process CountJob " << p << " @ " << step);
        step->count(p, jobqueue);
    }
};

struct DistributeJob : public Job
{
    SampleSortStep*   step;
    unsigned int        p;

    DistributeJob(SampleSortStep* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process DistributeJob " << p << " @ " << step);
        step->distribute(p, jobqueue);
    }
};

struct CopybackJob : public Job
{
    SampleSortStep*        step;
    unsigned int        p;

    CopybackJob(SampleSortStep* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process CopybackJob " << p << " @ " << step);
        step->copyback(p, jobqueue);
    }
};

inline unsigned int find_bkt(const key_type& key, const key_type* splitter, size_t leaves)
{
    unsigned int lo = 0, hi = leaves;

    while ( lo < hi )
    {
        unsigned int mid = (lo + hi) >> 1;
        assert(mid < leaves);

        if (key <= splitter[mid])
            hi = mid;
        else // (key > splitter[mid])
            lo = mid + 1;
    }

    size_t b = lo * 2;                              // < bucket
    if (lo < leaves && splitter[lo] == key) b += 1; // equal bucket

    return b;
}

void SampleSortStep::sample(JobQueue& jobqueue)
{
    const size_t oversample_factor = 4;
    size_t samplesize = oversample_factor * leaves;

    key_type samples[ samplesize ];

    LCGRandom rng(9384234);

    for (unsigned int i = 0; i < samplesize; ++i)
    {
        samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
    }

    std::sort(samples, samples + samplesize);

    DBG(debug_splitter, "splitter:");
    splitter_lcp[0] = 0; // sentinel for first < everything bucket
    for (size_t i = 0, j = oversample_factor/2; i < leaves; ++i)
    {
        splitter[i] = samples[j];
        DBG(debug_splitter, "key " << toHex(splitter[i]));

        if (i != 0) {
            key_type andSplit = splitter[i-1] ^ splitter[i];

            DBG1(debug_splitter, "    XOR -> " << toHex(andSplit) << " - ");

            DBG3(debug_splitter, count_high_zero_bits(andSplit) << " bits = "
                << count_high_zero_bits(andSplit) / 8 << " chars lcp");

            splitter_lcp[i] = count_high_zero_bits(andSplit) / 8;
        }

        j += oversample_factor;
    }

    bkt = new size_t[bktnum * parts + 1];

    // create new jobs
    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new CountJob(this, p) );
}

void SampleSortStep::count(unsigned int p, JobQueue& jobqueue)
{
    string* strB = strings + p * psize;
    string* strE = strings + std::min((p+1) * psize, n);
    if (strE < strB) strE = strB;

    // TODO: check if processor-local stack + copy is faster
#if 1
    size_t* mybkt = bkt + p * bktnum;

    memset(mybkt, 0, bktnum * sizeof(size_t));
    for (string* str = strB; str != strE; ++str)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(*str, depth);

        unsigned int b = find_bkt(key, splitter, leaves);
        assert(b < bktnum);
        ++mybkt[ b ];
    }
#else
    size_t mybkt[bktnum] = { 0 };

    for (string* str = strB; str != strE; ++str)
        ++mybkt[ (*str)[depth] ];

    memcpy(bkt + p * bktnum, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        count_finished(jobqueue);
}

void SampleSortStep::count_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing CountJob " << this << " with prefixsum");

    // inclusive prefix sum over bkt
    size_t sum = 0;
    for (unsigned int i = 0; i < bktnum; ++i)
    {
        for (unsigned int p = 0; p < parts; ++p)
        {
            bkt[p * bktnum + i] = (sum += bkt[p * bktnum + i]);
        }
    }
    assert(sum == n);

    // allocate out-of-place memory
    sorted = new string[n];

    // create new jobs
    pwork = parts;

    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new DistributeJob(this, p) );
}

void SampleSortStep::distribute(unsigned int p, JobQueue& jobqueue)
{
    string* strB = strings + p * psize;
    string* strE = strings + std::min((p+1) * psize, n);
    if (strE < strB) strE = strB;

    // TODO: check if processor-local stack + copy is faster
#if 1
    size_t* mybkt = bkt + p * bktnum;

    for (string* str = strB; str != strE; ++str)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(*str, depth);

        unsigned int b = find_bkt(key, splitter, leaves);
        assert(b < bktnum);
        sorted[ --mybkt[b] ] = *str;
    }
#else
    size_t mybkt[bktnum];
    memcpy(mybkt, bkt + p * bktnum, sizeof(mybkt));

    for (string* str = strB; str != strE; ++str)
        sorted[ --mybkt[(*str)[depth]] ] = *str;

    if (p == 0) // these are needed for recursion into bkts
        memcpy(bkt, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        distribute_finished(jobqueue);
}

void SampleSortStep::distribute_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing DistributeJob " << this << " with copy to original");

    // create new jobs
    pwork = parts;

    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new CopybackJob(this, p) );
}

void SampleSortStep::copyback(unsigned int p, JobQueue& jobqueue)
{
    size_t strB = p * psize;
    size_t strE = std::min((p+1) * psize, n);
    if (strE < strB) strE = strB;

    memcpy(strings + strB, sorted + strB, (strE - strB) * sizeof(string));

    if (--pwork == 0)
        copyback_finished(jobqueue);
}

void SampleSortStep::copyback_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing CopybackJob " << this << " with copy to original");

    delete [] sorted;

    // first processor's bkt pointers are boundaries between bkts, just add sentinel:
    bkt[bktnum] = n;

    for (size_t i=0; i < bktnum-1; ++i) {
        // i is even -> bkt[i] is less-than bucket
        size_t bktsize = bkt[i+1] - bkt[i];
        if (bktsize > 1)
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bkt[i] << " size " << bktsize << " lcp " << int(splitter_lcp[i/2]));
            Enqueue(jobqueue, strings+bkt[i], bktsize, depth + splitter_lcp[i/2]);
        }
        ++i;
        bktsize = bkt[i+1] - bkt[i];
        // i is odd -> bkt[i] is equal bucket
        if (bktsize > 1)
        {
            if ( (splitter[i/2] & 0xFF) == 0 ) { // equal-bucket has NULL-terminated key
                // done.
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " is done!");
            }
            else {
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " lcp keydepth!");
                Enqueue(jobqueue, strings+bkt[i], bktsize, depth + sizeof(key_type));
            }
        }
    }

    {
        size_t bktsize = n - bkt[bktnum-1];

        if (bktsize > 0)
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bkt[bktnum-1] << " size " << bktsize << " no lcp");
            Enqueue(jobqueue, strings+bkt[bktnum-1], bktsize, depth);
        }
    }

    delete [] bkt;
    delete this;
}

void EnqueueBig(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
{
    unsigned int parts = omp_get_max_threads();
    size_t psize = (n + parts-1) / parts;

    SampleSortStep* step = new SampleSortStep;
    step->strings = strings;
    step->n = n;
    step->depth = depth;

    step->parts = parts;
    step->psize = psize;
    step->pwork = parts;

    jobqueue.enqueue( new SampleJob(step) );
}

void Enqueue(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
{
    // TODO: tune parameter
    if (n > 1024*1024)
        return EnqueueBig(jobqueue, strings, n, depth);
    else
        return bingmann_parallel_radix_sort3::EnqueueSmall(jobqueue, strings, n, depth);
}

void bingmann_parallel_sample_sortBS(string* strings, size_t n)
{
    JobQueue jobqueue;
    Enqueue(jobqueue, strings,n,0);
    jobqueue.loop();
}

CONTESTANT_REGISTER_UCARRAY_PARALLEL(bingmann_parallel_sample_sortBS,
                                     "bingmann/parallel_sample_sortBS: binary search, no cache")

} // namespace bingmann_parallel_sample_sortBS
