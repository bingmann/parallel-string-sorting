/******************************************************************************
 * src/parallel/bingmann-parallel_sample_sortBTC.h
 *
 * Parallel Super Scalar String Sample-Sort with work-balancing, variant BTC:
 * with binary splitter tree and cache.
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

namespace bingmann_parallel_radix_sort {

using namespace stringtools;
using namespace jobqueue;

// prototype to call radix sort for small buckets
extern void EnqueueSmall(JobQueue& jobqueue, const StringPtr& strptr, size_t n, size_t depth);

}

namespace bingmann_parallel_sample_sortBTC {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = false;

static const bool debug_splitter = false;
static const bool debug_bucketsize = false;
static const bool debug_recursion = false;
static const bool debug_splitter_tree = false;

static const size_t MAXPROCS = 64;

static const size_t g_inssort_threshold = 64;

size_t g_totalsize;             // total size of input
size_t g_sequential_threshold;  // calculated threshold for sequential sorting
size_t g_threadnum;             // number of threads overall

/// Prototype called to schedule deeper sorts
void Enqueue(JobQueue& jobqueue, const StringPtr& strptr, size_t n, size_t depth);

// ****************************************************************************
// *** SampleSortStep out-of-place parallel sample sort with separate Jobs

typedef uint64_t key_type;

struct SampleSortStep
{
#if 0
    static const size_t numsplitters = 2*16;       // TODO: calculate to match L2 cache
#else
    static const size_t l2cache = 256*1024;

    // bounding equation: 2*K+1 buckets when counting bkt occurances.
    static const size_t numsplitters2 = (l2cache - sizeof(size_t)) / (2 * sizeof(size_t));

    static const size_t treebits = logfloor_<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;
#endif

    static const size_t bktnum = 2*numsplitters+1;

    StringPtr           strptr;
    size_t              n, depth;
    
    unsigned int        parts;
    size_t              psize;
    std::atomic<unsigned int> pwork;
    key_type            splitter[numsplitters];
    key_type            splitter_tree[numsplitters];
    unsigned char       splitter_lcp[numsplitters];
    size_t*             bkt;

    uint16_t*           bktcache[MAXPROCS];

    SampleSortStep(JobQueue& jobqueue, const StringPtr& _strptr, size_t _n, size_t _depth);

    void sample(unsigned int p, JobQueue& jobqueue);

    void count(unsigned int p, JobQueue& jobqueue);
    void count_finished(JobQueue& jobqueue);

    void distribute(unsigned int p, JobQueue& jobqueue);
    void distribute_finished(JobQueue& jobqueue);

    static inline void put_stats()
    {
        g_statscache >> "l2cache" << size_t(l2cache)
                     >> "splitter_treebits" << size_t(treebits)
                     >> "numsplitters" << size_t(numsplitters);
    }
};

template <void (SampleSortStep::*func)(unsigned int p, JobQueue& jobqueue)>
struct GenericJob : public Job
{
    SampleSortStep*     step;
    unsigned int        p;

    GenericJob(SampleSortStep* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        (step->*func)(p, jobqueue);
    }
};

typedef GenericJob<&SampleSortStep::sample> SampleJob;
typedef GenericJob<&SampleSortStep::count> CountJob;
typedef GenericJob<&SampleSortStep::distribute> DistributeJob;

SampleSortStep::SampleSortStep(JobQueue& jobqueue, const StringPtr& _strptr, size_t _n, size_t _depth)
    : strptr(_strptr), n(_n), depth(_depth)
{
    parts = g_threadnum * n / g_totalsize;
    if (parts == 0) parts = 1;

    psize = (n + parts-1) / parts;

    jobqueue.enqueue( new SampleJob(this, 0) );
}

/// binary search on splitter array for bucket number
static inline unsigned int
find_bkt_splittertree(const key_type& key, const key_type* splitter, const key_type* splitter_tree0, size_t numsplitters)
{
    // binary tree traversal without left branch

    const key_type* splitter_tree = splitter_tree0 - 1;
    unsigned int i = 1;

    while ( i <= numsplitters )
    {
        i = 2*i + (key <= splitter_tree[i] ? 0 : 1);
    }

    i -= numsplitters+1;

    size_t b = i * 2;                                   // < bucket
    if (i < numsplitters && splitter[i] == key) b += 1; // equal bucket

    return b;
}

void SampleSortStep::sample(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process SampleJob " << p << " @ " << this);

    const size_t oversample_factor = 1;
    size_t samplesize = oversample_factor * numsplitters;

    string* strings = strptr.active();
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
    pwork = parts;
    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new CountJob(this, p) );
}

void SampleSortStep::count(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process CountJob " << p << " @ " << this);

    string* strB = strptr.active() + p * psize;
    string* strE = strptr.active() + std::min((p+1) * psize, n);
    if (strE < strB) strE = strB;
    size_t strN = strE-strB;

    uint16_t* mybktcache = bktcache[p] = new uint16_t[strN];

    for (size_t i = 0; i < strN; ++i)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(strB[i], depth);

        unsigned int b = find_bkt_splittertree(key, splitter, splitter_tree, numsplitters);
        assert(b < bktnum);
        mybktcache[ i ] = b;
    }

    size_t* mybkt = bkt + p * bktnum;
    memset(mybkt, 0, bktnum * sizeof(size_t));

    for (uint16_t* bc = mybktcache; bc != mybktcache + strN; ++bc)
        ++mybkt[ *bc ];

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

    // create new jobs
    pwork = parts;
    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new DistributeJob(this, p) );
}

void SampleSortStep::distribute(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process DistributeJob " << p << " @ " << this);

    string* strB = strptr.active() + p * psize;
    string* strE = strptr.active() + std::min((p+1) * psize, n);
    if (strE < strB) strE = strB;

    string* sorted = strptr.shadow(); // get alternative shadow pointer array

    size_t* mybkt = bkt + p * bktnum;
    uint16_t* mybktcache = bktcache[p];

    for (string* str = strB; str != strE; ++str, ++mybktcache)
        sorted[ --mybkt[ *mybktcache ] ] = *str;

    delete bktcache[p];

    if (--pwork == 0)
        distribute_finished(jobqueue);
}

void SampleSortStep::distribute_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing DistributeJob " << this << " with enqueuing subjobs");

    // first processor's bkt pointers are boundaries between bkts, just add sentinel:
    assert(bkt[0] == 0);
    bkt[bktnum] = n;

    size_t i = 0;
    while (i < bktnum-1)
    {
        // i is even -> bkt[i] is less-than bucket
        size_t bktsize = bkt[i+1] - bkt[i];
        if (bktsize == 0)
            ;
        else if (bktsize == 1) // just one string pointer, copyback
            strptr.flip_ptr(bkt[i]).to_original(1);
        else
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bkt[i] << " size " << bktsize << " lcp " << int(splitter_lcp[i/2]));
            Enqueue(jobqueue, strptr.flip_ptr(bkt[i]), bktsize, depth + splitter_lcp[i/2]);
        }
        ++i;
        bktsize = bkt[i+1] - bkt[i];
        // i is odd -> bkt[i] is equal bucket
        if (bktsize == 0)
            ;
        else if (bktsize == 1) // just one string pointer, copyback
            strptr.flip_ptr(bkt[i]).to_original(1);
        else
        {
            if ( (splitter[i/2] & 0xFF) == 0 ) { // equal-bucket has NULL-terminated key
                // done.
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " is done!");
                strptr.flip_ptr(bkt[i]).to_original(bktsize);
            }
            else {
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " lcp keydepth!");
                Enqueue(jobqueue, strptr.flip_ptr(bkt[i]), bktsize, depth + sizeof(key_type));
            }
        }
        ++i;
    }

    size_t bktsize = bkt[i+1] - bkt[i];

    if (bktsize == 0)
        ;
    else if (bktsize == 1) // just one string pointer, copyback
        strptr.flip_ptr(bkt[i]).to_original(1);
    else
    {
        DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bkt[i] << " size " << bktsize << " no lcp");
        Enqueue(jobqueue, strptr.flip_ptr(bkt[i]), bktsize, depth);
    }

    delete [] bkt;
    delete this;
}

void Enqueue(JobQueue& jobqueue, const StringPtr& strptr, size_t n, size_t depth)
{
    if (n > g_sequential_threshold)
        new SampleSortStep(jobqueue, strptr, n, depth);
    else
        bingmann_parallel_radix_sort::EnqueueSmall(jobqueue, strptr, n, depth);
}

void parallel_sample_sortBTC(string* strings, size_t n)
{
    g_totalsize = n;
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);

    SampleSortStep::put_stats();

    string* shadow = new string[n]; // allocate shadow pointer array

    JobQueue jobqueue;
    Enqueue(jobqueue, StringPtr(strings, shadow), n, 0);
    jobqueue.loop();

    delete [] shadow;
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTC,
                             "bingmann/parallel_sample_sortBTC",
                             "bingmann/parallel_sample_sortBTC: binary tree, bktcache")

} // namespace bingmann_parallel_sample_sortBTC
