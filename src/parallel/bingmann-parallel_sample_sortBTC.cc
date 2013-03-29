/******************************************************************************
 * src/parallel/bingmann-parallel_sample_sortBTC.h
 *
 * Parallel Super Scalar String Sample-Sort with work-balancing, variant BT:
 * with binary splitter tree and no cache.
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
#include <math.h>

#include <iostream>
#include <vector>
#include <algorithm>

#include "../tools/debug.h"
#include "../tools/lcgrandom.h"
#include "../tools/contest.h"
#include "../tools/stringtools.h"
#include "../tools/jobqueue.h"
#include "../tools/logfloor.h"

extern size_t g_smallsort;

namespace bingmann_parallel_radix_sort3 {

using namespace stringtools;
using namespace jobqueue;

// prototype to call radix sort for small buckets
extern void EnqueueSmall(JobQueue& jobqueue, string* strings, size_t n, size_t depth);

}

namespace bingmann_parallel_sample_sortBTC {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = false;

static const bool debug_splitter = false;
static const bool debug_bucketsize = false;
static const bool debug_recursion = false;
static const bool debug_splitter_tree = false;

static const size_t MAXPROCS = 128;

/// Prototype called to schedule deeper sorts
void Enqueue(JobQueue& jobqueue, string* strings, size_t n, size_t depth);

// ****************************************************************************
// *** SampleSortStep out-of-place parallel sample sort with separate Jobs

typedef uint64_t key_type;

struct SampleSortStep
{
#if 0
    static const size_t numsplitters = 2*16;       // TODO: calculate to match L2 cache
#else
    static const size_t l2cache = 256*1024;

// bounding equations:
// splitters            + bktsize
// n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));
    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / ( sizeof(key_type) );

    static const size_t numsplitters = (1 << logfloor_<numsplitters2>::value) - 1;

#endif

    static const size_t bktnum = 2*numsplitters+1;

    string*             strings;
    size_t              n;
    size_t              depth;
    
    unsigned int        parts;
    size_t              psize;
    std::atomic<unsigned int> pwork;
    key_type            splitter[numsplitters];
    key_type            splitter_tree[numsplitters];
    unsigned char       splitter_lcp[numsplitters];
    size_t*             bkt;

    string*             sorted;

    uint16_t*           bktcache[MAXPROCS];

    SampleSortStep() {}

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

/// binary search on splitter array for bucket number
static inline unsigned int
find_bkt_splittertree(const key_type& key, const key_type* splitter, const key_type* splitter_tree0, size_t numsplitters)
{
    // binary tree traversal without left branch

    const key_type* splitter_tree = splitter_tree0 - 1;
    unsigned int i = 1;

    while ( i <= numsplitters )
    {
        if (key <= splitter_tree[i])
            i = 2*i + 0;
        else // (key > splitter_tree[i])
            i = 2*i + 1;
    }

    i -= numsplitters+1;

    size_t b = i * 2;                                   // < bucket
    if (i < numsplitters && splitter[i] == key) b += 1; // equal bucket

    return b;
}

void SampleSortStep::sample(JobQueue& jobqueue)
{
    const size_t oversample_factor = 4;
    size_t samplesize = oversample_factor * numsplitters;

    key_type samples[ samplesize ];

    LCGRandom rng(&samples);

    for (unsigned int i = 0; i < samplesize; ++i)
    {
        samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
    }

    std::sort(samples, samples + samplesize);

    DBG(debug_splitter, "splitter:");
    splitter_lcp[0] = 0; // sentinel for first < everything bucket
    for (size_t i = 0, j = oversample_factor/2; i < numsplitters; ++i)
    {
        splitter[i] = samples[j];
        DBG(debug_splitter, "key " << toHex(splitter[i]));

        if (i != 0) {
            key_type xorSplit = splitter[i-1] ^ splitter[i];

            DBG1(debug_splitter, "    XOR -> " << toHex(xorSplit) << " - ");

            DBG3(debug_splitter, count_high_zero_bits(xorSplit) << " bits = "
                << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            splitter_lcp[i] = count_high_zero_bits(xorSplit) / 8;
        }

        j += oversample_factor;
    }

    // step 2.1: construct splitter tree to perform binary search

    {
        size_t t = 0;
        size_t highbit = (numsplitters+1) / 2;

        while ( highbit > 0 )
        {
            DBG(debug_splitter_tree, "highbit = " << highbit);

            size_t p = highbit-1;
            size_t inc = highbit << 1;

            while ( p <= numsplitters )
            {
                DBG(debug_splitter_tree, "p = " << p);

                splitter_tree[t++] = splitter[p];

                p += inc;
            }

            highbit >>= 1;
        }
    }

    if (debug_splitter_tree)
    {
        DBG1(1, "splitter_tree: ");
        for (size_t i = 0; i < numsplitters; ++i)
        {
            DBG2(1, splitter_tree[i] << " ");
        }
        DBG3(1, "");
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
    size_t strN = strE-strB;

    // TODO: check if processor-local stack + copy is faster
#if 1
    size_t* mybkt = bkt + p * bktnum;
    uint16_t* mybktcache = bktcache[p] = new uint16_t[strN];

    memset(mybkt, 0, bktnum * sizeof(size_t));
    for (size_t i = 0; i < strN; ++i)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(strB[i], depth);

        unsigned int b = find_bkt_splittertree(key, splitter, splitter_tree, numsplitters);
        assert(b < bktnum);
        mybktcache[ i ] = b;
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
    size_t strN = strE-strB;

    // TODO: check if processor-local stack + copy is faster
#if 1
    size_t* mybkt = bkt + p * bktnum;
    uint16_t* mybktcache = bktcache[p];

    for (size_t i = 0; i < strN; ++i)
    {
        // binary search in splitter with equal check
        unsigned int b = mybktcache[i];
        assert(b < bktnum);
        sorted[ --mybkt[b] ] = strB[i];
    }

    delete bktcache[p];
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
    if (n > 128*1024)
        return EnqueueBig(jobqueue, strings, n, depth);
    else
        return bingmann_parallel_radix_sort3::EnqueueSmall(jobqueue, strings, n, depth);
}

void parallel_sample_sortBTC(string* strings, size_t n)
{
    JobQueue jobqueue;
    Enqueue(jobqueue, strings,n,0);
    jobqueue.loop();
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTC,
                             "bingmann/parallel_sample_sortBTC",
                             "bingmann/parallel_sample_sortBTC: binary tree, bktcache")

} // namespace bingmann_parallel_sample_sortBTC
