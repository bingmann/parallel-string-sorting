/******************************************************************************
 * src/parallel/bingmann-parallel_mkqs.cc
 *
 * Parallel multikey-quicksort.
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

#include "../tools/debug.h"
#include "../tools/contest.h"
#include "../tools/stringtools.h"
#include "../tools/jobqueue.h"

#include "../sequential/inssort.h"

extern size_t g_smallsort;

namespace bingmann_parallel_mkqs {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = true;
static const bool debug_blocks = true;
static const bool debug_cmp = false;

static const bool use_work_stealing = true;

typedef unsigned char* string;

static const size_t g_inssort_threshold = 64;

/// Prototype called to schedule deeper sorts
void Enqueue(JobQueue& jobqueue, string* strings, size_t n, size_t depth);

/// Multi-character cache block for each string

struct StrCache
{
    uint64_t cache;
    string str;

    friend std::ostream& operator<< (std::ostream& os, const StrCache& sc)
    {
        return os << toHex(sc.cache);
    }
};

template <typename CharT>
CharT med3char(CharT a, CharT b, CharT c)
{
    if (a == b)           return a;
    if (c == a || c == b) return c;
    if (a < b) {
        if (b < c) return b;
        if (a < c) return c;
        return a;
    }
    if (b > c) return b;
    if (a < c) return a;
    return c;
}

// TODO: make faster
static inline int cmp(const uint64_t& a, const uint64_t& b)
{
    if (a > b) return 1;
    if (a < b) return -1;
    return 0;
}

// ****************************************************************************
// *** ParallelMKQS - split ternary with 8-byte superalphabet

struct ParallelMKQS
{
    StrCache*           strcache;
    size_t              n;
    size_t              depth;

    unsigned int        procs;
    std::atomic<unsigned int> pwork;
    size_t              partsize;

    uint64_t            pivot;

    // Variables for blocking scheme on strcache. The strcache array is split
    // into blocks of block_size.

    static const unsigned int   block_size = 1024;
    unsigned int                block_count;
    std::atomic<uint_fast64_t>  block_lt_gt;

    // The atomic variable block_lt_gt holds both pointers, to the beginning of
    // the first unused less-than block and the end of the last unused
    // greater-than block. A lt- or gt-block is reserved by atomic exchange of
    // the block_lt_gt variable.

    // function to retreive blocks from an end
    StrCache*   get_block_lt(unsigned int& pos, unsigned int& size);
    StrCache*   get_block_gt(unsigned int& pos, unsigned int& size);

    ParallelMKQS() {}

    void fillcache(unsigned int p, JobQueue& jobqueue);
    void fillcache_finished(JobQueue& jobqueue);

    void partition(unsigned int p, JobQueue& jobqueue);
    void partition_finished(JobQueue& jobqueue);
};

static inline uint64_t pack64(const uint32_t& lo, const uint32_t& hi)
{ return (((uint64_t)hi) << 32) | lo; }

static inline uint32_t unpack64_low(const uint64_t& x)
{ return (uint32_t)x; }

static inline uint32_t unpack64_high(const uint64_t& x)
{ return (uint32_t)(x >> 32); }

StrCache* ParallelMKQS::get_block_lt(unsigned int& pos, unsigned int& size)
{
    while (1)
    {
        uint64_t blk = block_lt_gt;

        uint64_t blk_lt = unpack64_low(blk);
        uint64_t blk_gt = unpack64_high(blk);
        DBG(debug_blocks, "current blocks: lt " << blk_lt << " gt " << blk_gt);

        // if block pointers overlap -> no block left.
        if ( blk_lt >= blk_gt ) return NULL;

        uint64_t nblk = pack64( blk_lt+1, blk_gt );

        // atomic swap in nblk, which reserves one lt-block
        if (block_lt_gt.compare_exchange_weak(blk, nblk))
        {
            // calculate size of reserved block, last one is smaller
            pos = 0;
            size = std::min<size_t>(block_size, n - (blk_lt + block_count-blk_gt) * block_size);
            StrCache* sc = strcache + blk_lt * block_size;
            DBG(debug_blocks, "reserved lt-block " << blk_lt << " @ " << (sc - strcache) << " size " << size);
            return sc;
        }
    }
}
 
StrCache* ParallelMKQS::get_block_gt(unsigned int& pos, unsigned int& size)
{
    while (1)
    {
        uint64_t blk = block_lt_gt;

        uint64_t blk_lt = unpack64_low(blk);
        uint64_t blk_gt = unpack64_high(blk);
        DBG(debug_blocks, "current blocks: lt " << blk_lt << " gt " << blk_gt);

        // if block pointers overlap -> no block left.
        if ( blk_lt >= blk_gt ) return NULL;

        uint64_t nblk = pack64( blk_lt, blk_gt-1 );

        // atomic swap in nblk, which reserves one gt-block
        if (block_lt_gt.compare_exchange_weak(blk, nblk))
        {
            // calculate size of reserved block, last one is smaller
            pos = 0;
            size = std::min<size_t>(block_size, n - (blk_lt + block_count-blk_gt) * block_size);
            StrCache* sc = strcache + n - (block_count-blk_gt) * block_size - size;
            DBG(debug_blocks, "reserved gt-block " << blk_gt << " @ " << (sc - strcache) << " size " << size);
            return sc;
        }
    }
}

template < void (ParallelMKQS::*func)(unsigned int, JobQueue&) >
struct GenericJob : public Job
{
    ParallelMKQS*       step;
    unsigned int        p;

    GenericJob(ParallelMKQS* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        (step->*func)(p, jobqueue);
    }
};

typedef GenericJob<&ParallelMKQS::fillcache> FillCacheJob;
typedef GenericJob<&ParallelMKQS::partition> PartitionJob;

void ParallelMKQS::fillcache(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process FillCacheJob " << p << " @ " << this);

    StrCache* strB = strcache + p * partsize;
    StrCache* strE = strcache + std::min((p+1) * partsize, n);
    if (strE < strB) strE = strB;

    // fill cache for processor's part
    for (StrCache* s = strB; s != strE; ++s)
        s->cache = get_char<uint64_t>(s->str, depth);

    if (--pwork == 0)
        fillcache_finished(jobqueue);
}

void ParallelMKQS::fillcache_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "finished FillCacheJobs @ " << this);

    // select pivot from median of 9, no swapping
    pivot = med3char(
        med3char(strcache[0].cache,       strcache[n/8].cache,     strcache[n/4].cache),
        med3char(strcache[n/2-n/8].cache, strcache[n/2].cache,     strcache[n/2+n/8].cache),
        med3char(strcache[n-1-n/4].cache, strcache[n-1-n/8].cache, strcache[n-3].cache)
        );

    DBG(debug_jobs, "selected pivot " << toHex(pivot));

    // create new jobs
    pwork = procs;
    for (unsigned int p = 0; p < procs; ++p)
        jobqueue.enqueue( new PartitionJob(this, p) );
}

void ParallelMKQS::partition(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "process PartitionJob " << p << " @ " << this);

    unsigned int pos_lt = 0, pos_eq = 0, pos_gt = 0;
    unsigned int size_lt = 0, size_eq = 0, size_gt = 0;
    StrCache *blk_lt = NULL, *blk_eq = NULL, *blk_gt = NULL;

    while (1)
    {
        while ( ( pos_lt < size_lt || (blk_lt = get_block_lt(pos_lt, size_lt)) ) &&
                ( pos_eq < size_eq || (blk_eq = get_block_lt(pos_eq, size_eq)) ) )
        {
            int res = cmp(blk_lt[pos_lt].cache, pivot);
            if (res > 0) { // > than pivot
                DBG(debug_cmp, "blk_lt[" << pos_lt << "] = " << blk_lt[pos_lt] << " > pivot " << toHex(pivot) << ", break.");
                break;
            }
            else if (res == 0) { // = pivot
                DBG(debug_cmp, "blk_lt[" << pos_lt << "] = " << blk_lt[pos_lt] << " = pivot " << toHex(pivot) << ", swap to blk_eq");
                std::swap(blk_lt[pos_lt], blk_eq[pos_eq++]);
            }
            else { // < than pivot
                DBG(debug_cmp, "blk_lt[" << pos_lt << "] = " << blk_lt[pos_lt] << " < pivot " << toHex(pivot));
                pos_lt++;
            }
        }

        if (!blk_lt || !blk_eq) break;

        while ( ( pos_gt < size_gt || (blk_gt = get_block_gt(pos_gt, size_gt)) ) &&
                ( pos_eq < size_eq || (blk_eq = get_block_lt(pos_eq, size_eq)) ) )
        {
            int res = cmp(blk_gt[pos_gt].cache, pivot);
            if (res < 0) { // < than pivot
                DBG(debug_cmp, "blk_gt[" << pos_gt << "] = " << blk_gt[pos_gt] << " < pivot " << toHex(pivot) << ", break.");
                break;
            }
            else if (res == 0) { // = pivot
                DBG(debug_cmp, "blk_gt[" << pos_gt << "] = " << blk_gt[pos_gt] << " = pivot " << toHex(pivot) << ", swap to blk_eq");
                std::swap(blk_gt[pos_gt], blk_eq[pos_eq++]);
            }
            else { // > than pivot
                DBG(debug_cmp, "blk_gt[" << pos_gt << "] = " << blk_gt[pos_gt] << " > pivot " << toHex(pivot));
                pos_gt++;
            }
        }

        if (!blk_gt || !blk_eq) break;

        DBG(debug_cmp, "swap blk_lt[" << pos_lt << "] = " << blk_lt[pos_lt]
            << " and blk_gt[" << pos_gt << "] = " << blk_gt[pos_gt]);

        assert( blk_lt[pos_lt].cache > pivot && blk_gt[pos_gt].cache < pivot );
        std::swap(blk_lt[pos_lt++], blk_gt[pos_gt++]);
    }


    abort();




    if (--pwork == 0)
        partition_finished(jobqueue);
}

void ParallelMKQS::partition_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "finished PartitionJobs @ " << this);

    abort();
}

size_t totalsize;

void EnqueueParallelMKQS(JobQueue& jobqueue, StrCache* str, size_t n, size_t depth)
{
    size_t procs = omp_get_max_threads() * n / totalsize;
    if (procs == 0) procs = 1;

    procs = 1;

    ParallelMKQS* step = new ParallelMKQS;
    step->strcache = str;
    step->n = n;
    step->depth = depth;

    step->procs = procs;
    step->pwork = procs;
    step->partsize = (n + procs-1) / procs; // size for each processor

    // need to pack block_count into a uint32_t
    uint64_t block_count = (n + step->block_size-1) / step->block_size;

    DBG(debug_blocks, "input has block_count " << block_count);
    if (block_count >= (((uint64_t)1) << 32)) abort();

    step->block_count = block_count;
    step->block_lt_gt = pack64( 0, block_count );

    for (unsigned int p = 0; p < procs; ++p)
        jobqueue.enqueue( new FillCacheJob(step, p) );
}

void parallel_mkqs(string* strings, size_t n)
{
    totalsize = n;

    StrCache* cache = new StrCache[n];

    // copy pointers to charcache
//#pragma omp parallel for schedule(static)
    for (size_t i=0; i < n; ++i)
        cache[i].str = strings[i];
    
    JobQueue jobqueue;
    EnqueueParallelMKQS(jobqueue, cache, n, 0);
    jobqueue.loop();

    // copy pointers from charcache
//#pragma omp parallel for schedule(static)
    for (size_t i=0; i < n; ++i)
        strings[i] = cache[i].str;

    free(cache);
}

CONTESTANT_REGISTER_PARALLEL(parallel_mkqs,
                             "bingmann/parallel_mkqs",
                             "Parallel MKQS with Blocks and Cache8")

} // namespace bingmann_parallel_mkqs
