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

#undef DBGX
#define DBGX DBGX_OMP

#include <tr1/memory>
#include <tbb/concurrent_queue.h>

namespace bingmann_parallel_mkqs {

using namespace stringtools;
using namespace jobqueue;

static const bool debug = false;
static const bool debug_parajobs = false;
static const bool debug_blocks = false;
static const bool debug_cmp1 = false;
static const bool debug_cmp2 = false;
static const bool debug_seqjobs = false;

static const bool debug_selfcheck = false;

static const bool use_work_sharing = true;

typedef unsigned char* string;
typedef uint64_t key_type;

static const size_t g_inssort_threshold = 64;

size_t g_totalsize;             // total size of input
size_t g_sequential_threshold;  // calculated threshold for sequential sorting
size_t g_threadnum;             // number of threads overall

static const unsigned int block_size = 128*1024;

/// output string ranges for debugging

string* g_strings;

static inline std::string srange(string* str, size_t n)
{
    std::ostringstream oss;
    oss << "[" << (str - g_strings) << "," << (str + n - g_strings) << ")=" << n;
    return oss.str();
}

/// Multi-character cache block for each string

struct StrCache
{
    uint64_t key;
    string str;

    friend std::ostream& operator<< (std::ostream& os, const StrCache& sc)
    {
        return os << toHex(sc.key);
    }
};

struct StrCacheBlock
{
    unsigned int fill;
    StrCache    cache[block_size];

    inline string& str(unsigned int i)
    {
        assert(i < block_size);
        return cache[i].str;
    }

    inline uint64_t& key(unsigned int i)
    {
        assert(i < block_size);
        return cache[i].key;
    }
};

typedef tbb::concurrent_queue<StrCacheBlock*> BlockQueueType;
typedef tbb::concurrent_queue<key_type> PivotKeyQueueType;
typedef tbb::concurrent_queue<string> PivotStrQueueType;

template <typename CharT>
static inline const CharT&
med3char(const CharT& a, const CharT& b, const CharT& c)
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

template <typename CharT>
static inline CharT&
med3charC(CharT& a, CharT& b, CharT& c)
{
    if (cmp(a, b) == 0)                   return a;
    if (cmp(c, a) == 0 || cmp(c, b) == 0) return c;
    if (cmp(a, b) < 0) {
        if (cmp(b, c) < 0) return b;
        if (cmp(a, c) < 0) return c;
        return a;
    }
    if (cmp(b, c) > 0) return b;
    if (cmp(a, c) < 0) return a;
    return c;
}

// TODO: make faster
static inline int cmp(const uint64_t& a, const uint64_t& b)
{
    if (a > b) return 1;
    if (a < b) return -1;
    return 0;
}

static inline int cmp(const StrCache& a, const StrCache& b)
{
    return cmp(a.key, b.key);
}

// ****************************************************************************
// *** Rantala's Multikey-Quicksort with cached 8-byte superalphabet

// Insertion sort, ignores any cached characters.
static inline void
insertion_sort_nocache(StrCache* cache, int n, size_t depth)
{
    StrCache *pi, *pj;
    unsigned char *s, *t;
    for (pi = cache + 1; --n > 0; ++pi) {
        unsigned char* tmp = pi->str;
        for (pj = pi; pj > cache; --pj) {
            for (s=(pj-1)->str+depth, t=tmp+depth; *s==*t && *s!=0;
                 ++s, ++t)
                ;
            if (*s <= *t)
                break;
            pj->str = (pj-1)->str;
        }
        pj->str = tmp;
    }
}

// Insertion sorts the strings only based on the cached characters.
static inline void
insertion_sort_cache_block(StrCache* cache, int n)
{
    StrCache *pi, *pj;
    for (pi = cache + 1; --n > 0; ++pi) {
        StrCache tmp = *pi;
        for (pj = pi; pj > cache; --pj) {
            if (cmp((pj-1)->key, tmp.key) <= 0)
                break;
            *pj = *(pj-1);
        }
        *pj = tmp;
    }
}

template <bool CacheDirty>
static inline void
insertion_sort(StrCache* cache, int n, size_t depth)
{
    if (n == 0) return;
    if (CacheDirty) return insertion_sort_nocache(cache, n, depth);

    insertion_sort_cache_block(cache, n);

    size_t start=0, cnt=1;
    for (int i=0; i < n-1; ++i) {
        if (cmp(cache[i], cache[i+1]) == 0) {
            ++cnt;
            continue;
        }
        if (cnt > 1 && cache[start].key & 0xFF)
            insertion_sort_nocache(cache + start, cnt,
                                   depth + sizeof(key_type));
        cnt = 1;
        start = i+1;
    }
    if (cnt > 1 && cache[start].key & 0xFF)
        insertion_sort_nocache(cache+start, cnt, depth + sizeof(key_type));
}

// ****************************************************************************
// *** SequentialMKQS - split ternary with 8-byte superalphabet

struct MKQSStep
{
    StrCache* cache;
    size_t num_lt, num_eq, num_gt, n, depth;
    size_t idx;
    bool eq_recurse;

    MKQSStep(StrCache* cache, size_t n, size_t depth, bool CacheDirty)
    {
        if (CacheDirty)
        {
            for (size_t i=0; i < n; ++i) {
                cache[i].key = get_char<key_type>(cache[i].str, depth);
            }
        }
        // Move pivot to first position to avoid wrapping the unsigned values
        // we are using in the main loop from zero to max.
        std::swap(cache[0], med3charC(
                      med3charC(cache[0],       cache[n/8],     cache[n/4]),
                      med3charC(cache[n/2-n/8], cache[n/2],     cache[n/2+n/8]),
                      med3charC(cache[n-1-n/4], cache[n-1-n/8], cache[n-3])
                      ));
        StrCache pivot = cache[0];
        size_t first   = 1;
        size_t last    = n-1;
        size_t beg_ins = 1;
        size_t end_ins = n-1;
        while (true) {
            while (first <= last) {
                const int res = cmp(cache[first], pivot);
                if (res > 0) {
                    break;
                } else if (res == 0) {
                    std::swap(cache[beg_ins++], cache[first]);
                }
                ++first;
            }
            while (first <= last) {
                const int res = cmp(cache[last], pivot);
                if (res < 0) {
                    break;
                } else if (res == 0) {
                    std::swap(cache[end_ins--], cache[last]);
                }
                --last;
            }
            if (first > last)
                break;
            std::swap(cache[first], cache[last]);
            ++first;
            --last;
        }
        // Some calculations to make the code more readable.
        const size_t num_eq_beg = beg_ins;
        const size_t num_eq_end = n-1-end_ins;
                     num_eq     = num_eq_beg+num_eq_end;
                     num_lt     = first-beg_ins;
                     num_gt     = end_ins-last;

        // Swap the equal pointers from the beginning to proper position.
        const size_t size1 = std::min(num_eq_beg, num_lt);
        std::swap_ranges(cache, cache+size1, cache+first-size1);
        // Swap the equal pointers from the end to proper position.
        const size_t size2 = std::min(num_eq_end, num_gt);
        std::swap_ranges(cache+first, cache+first+size2, cache+n-size2);

        // Save offsets for recursive sorting
        this->cache = cache;
        this->n = n;
        this->depth = depth;
        this->idx = 0;
        this->eq_recurse = (pivot.key & 0xFF);
/*
        multikey_cache<false>(cache, num_lt, depth);
        multikey_cache<false>(cache+num_lt+num_eq, num_gt, depth);
        if (pivot.key & 0xFF)
            multikey_cache<true>(cache+num_lt, num_eq, depth + sizeof(key_type));
*/
    }
};

// magic deleter class for tr1::shared_ptr to delete [] array instead of delete
// array.
template<typename T>
struct shared_ptr_array_deleter {
    void operator()(T* p) {
        delete [] p;
    }
};

template <bool CacheDirty>
struct SequentialMKQS : public Job
{
    string*             strings;
    size_t              n, depth;

    BlockQueueType*     block_queue;

    StrCache*           cache;

    // reference counted pointer to base_cache, which might be shared among
    // threads due to work sharing.
    typedef std::tr1::shared_ptr<StrCache> StrCachePtrType;
    StrCachePtrType     cache_base;

    // *** Constructors

    SequentialMKQS(string* _strings, size_t _n, size_t _depth, StrCache* _cache)
        : strings(_strings), n(_n), depth(_depth),
          cache(_cache), cache_base(_cache)
    {
    }

    SequentialMKQS(JobQueue& jobqueue, string* _strings, size_t _n, size_t _depth,
                   BlockQueueType* _block_queue)
        : strings(_strings), n(_n), depth(_depth),
          block_queue(_block_queue), cache(NULL)
    {
        jobqueue.enqueue(this);
    }

    SequentialMKQS(JobQueue& jobqueue, string* _strings, size_t _n, size_t _depth,
                   StrCache* _cache, StrCachePtrType _cache_base)
        : strings(_strings), n(_n), depth(_depth),
          cache(_cache), cache_base(_cache_base)
    {
        jobqueue.enqueue(this);
    }

    // *** Sequential Work

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_seqjobs, "SequentialMKQS for " << n << " strings @ " << this);

        if (!cache)
        {
            if (n == 1) // nothing to sort
            {
                DBG(debug_seqjobs, "copy result to output string ptrs " << srange(strings,n));

                StrCacheBlock* scb = NULL;
                while ( block_queue->try_pop(scb) )
                {
                    assert(scb && (scb->fill == 0 || scb->fill == 1));
                    if (scb->fill == 1) {
                        strings[0] = scb->cache[0].str;
                    }
                    delete scb;
                }

                while ( block_queue->try_pop(scb) )
                {
                    assert(scb && scb->fill == 0);
                    delete scb;
                }

                assert(block_queue->empty());
                delete block_queue;

                return;
            }

            // locally allocate (ptr,cache) array (NUMA-friendly)
            cache = new StrCache[n];
            cache_base = StrCachePtrType(cache, shared_ptr_array_deleter<StrCache>());

            // extract all elements, and maybe update cache

            StrCacheBlock* scb;
            size_t o = 0;

            while ( block_queue->try_pop(scb) )
            {
                for (unsigned int i = 0; i < scb->fill; ++i, ++o)
                {
                    if (CacheDirty) {
                        cache[o].str = scb->cache[i].str;
                        cache[o].key = get_char<uint64_t>(cache[o].str, depth);
                    }
                    else {
                        cache[o] = scb->cache[i];
                    }
                }
                delete scb;
            }

            delete block_queue;
        }

        // sort using cache
        sequential_mkqs(jobqueue);
    }

    inline void
    sequential_mkqs(JobQueue& jobqueue)
    {
        DBG(debug_seqjobs, "SequentialMKQS on area " << srange(strings,n) << " @ job " << this);

        if (n < g_inssort_threshold)
        {
            insertion_sort<true>(cache, n, depth);

            DBG(debug_seqjobs, "copy result to output string ptrs " << srange(strings,n) << " @ job " << this);
            for (size_t i=0; i < n; ++i)
                strings[i] = cache[i].str;

            return;
        }

        // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
        size_t pop_front = 0;
        std::vector<MKQSStep> stack;
        stack.push_back( MKQSStep(cache,n,depth,CacheDirty) );

        // for output, save pointer to finished entry
        StrCache* cache_finished = cache + n;

    jumpout:
        while ( stack.size() > pop_front )
        {
            while ( stack.back().idx < 3 )
            {
                if (use_work_sharing && jobqueue.has_idle())
                {
                    // convert top level of stack into independent jobs

                    MKQSStep& st = stack[pop_front];
                    string* st_strings = strings + (st.cache - cache);

                    DBG(debug_seqjobs, "Queueing front of SequentialMKQS's stack level " << pop_front << ", idx " << st.idx
                        << ", areas lt " << srange(st_strings, st.num_lt)
                        << " eq " << srange(st_strings + st.num_lt, st.num_eq)
                        << " gt " << srange(st_strings + st.num_lt + st.num_eq, st.num_gt)
                        << " @ job " << this);

                    if (st.idx == 0 && st.num_lt != 0)
                    {
                        DBG(debug_seqjobs, "Queueing job for lt-area " << srange(st_strings, st.num_lt) << " @ job " << this);

                        new SequentialMKQS<false>(jobqueue, st_strings, st.num_lt, st.depth,
                                                  st.cache, cache_base);
                    }
                    if (st.idx <= 1 && st.num_eq != 0)
                    {
                        DBG(debug_seqjobs, "Queueing job for eq-area " << srange(st_strings + st.num_lt, st.num_eq) << " @ job " << this);

                        if (st.eq_recurse) {
                            new SequentialMKQS<true>(jobqueue, st_strings + st.num_lt,
                                                     st.num_eq, st.depth + sizeof(key_type),
                                                     st.cache + st.num_lt, cache_base);
                        }
                        else
                        {
                            DBG(debug_seqjobs, "copy result to output string ptrs " << srange(st_strings + st.num_lt, st.num_eq) << " - no recurse equal @ job " << this);

                            // seems overkill to create an extra thread job for this:
                            for (size_t i=st.num_lt; i < st.num_lt+st.num_eq; ++i)
                                st_strings[i] = st.cache[i].str;
                        }
                    }
                    if (st.idx <= 2 && st.num_gt != 0)
                    {
                        DBG(debug_seqjobs, "Queueing job for gt-area " << srange(st_strings + st.num_lt + st.num_eq, st.num_gt) << " @ job " << this);

                        new SequentialMKQS<false>(jobqueue, st_strings + st.num_lt + st.num_eq,
                                                  st.num_gt, st.depth,
                                                  st.cache + st.num_lt + st.num_eq, cache_base);
                    }

                    // recalculate finish pointer for this thread (removes all queued areas)
                    if (st.idx == 0) {
                        cache_finished = st.cache;
                    }
                    else if (st.idx == 1) {
                        cache_finished = st.cache + st.num_lt;
                    }
                    else if (st.idx == 2) {
                        cache_finished = st.cache + st.num_lt + st.num_eq;
                    }

                    // shorten the current stack
                    ++pop_front;
                    goto jumpout;
                }

                MKQSStep& ms = stack.back();
                ++ms.idx; // increment here, because stack may change

                // process the lt-subsequence
                if (ms.idx == 1)
                {
                    if (ms.num_lt == 0)
                        continue;
                    else if (ms.num_lt < g_inssort_threshold)
                        insertion_sort<false>(ms.cache, ms.num_lt, ms.depth);
                    else
                        stack.push_back( MKQSStep(ms.cache, ms.num_lt, ms.depth, false) );
                }
                // process the eq-subsequence
                else if (ms.idx == 2)
                {
                    if (!ms.eq_recurse || ms.num_eq == 0)
                        continue;
                    else if (ms.num_eq < g_inssort_threshold)
                        insertion_sort<true>(ms.cache + ms.num_lt, ms.num_eq,
                                             ms.depth + sizeof(key_type));
                    else
                        stack.push_back( MKQSStep(ms.cache + ms.num_lt, ms.num_eq,
                                                  ms.depth + sizeof(key_type), true) );
                }
                // process the gt-subsequence
                else
                {
                    assert(ms.idx == 3);
                    if (ms.num_gt == 0)
                        continue;
                    else if (ms.num_gt < g_inssort_threshold)
                        insertion_sort<false>(ms.cache + ms.num_lt + ms.num_eq, ms.num_gt,
                                              ms.depth);
                    else
                        stack.push_back( MKQSStep(ms.cache + ms.num_lt + ms.num_eq, ms.num_gt,
                                                  ms.depth, false) );
                }
            }

            stack.pop_back();
        }

        // copy string pointers to output array (for level pop_front)
        {
            size_t n_finished = cache_finished - cache;

            DBG(debug_seqjobs, "copy result to output string ptrs " << srange(strings, n_finished) << " @ job " << this);
            for (size_t i=0; i < n_finished; ++i)
                strings[i] = cache[i].str;
        }
    }
};

static inline void
bingmann_sequential_mkqs_cache8(string* strings, size_t n)
{
    StrCache* cache = new StrCache[n];
    for (size_t i=0; i < n; ++i) {
        cache[i].str = strings[i];
    }
    JobQueue jobqueue;
    SequentialMKQS<true>(strings, n, 0, cache).sequential_mkqs(jobqueue);
}

CONTESTANT_REGISTER(bingmann_sequential_mkqs_cache8,
                    "bingmann/sequential_mkqs_cache8",
                    "multikey_cache with 8byte cache (non-recursive)")

// ****************************************************************************
// *** BlockSource - provide new blocks of unpartitioned input to ParallelMKQS

struct BlockSourceInput
{
    string*             strings;
    size_t              n, depth;

    /// number of block in input
    unsigned int        block_count;

    /// atomic index to currently unprocessed block of input
    std::atomic<unsigned int> block_current;

    BlockSourceInput(string* _strings, size_t _n, size_t _depth)
        : strings(_strings), n(_n), depth(_depth),
          block_count( (n + block_size-1) / block_size ), // round-up
          block_current( 0 )
    {
    }

    key_type get_direct(size_t i) const
    {
        assert(i < n);
        return get_char<uint64_t>(strings[i], depth);
    }

    key_type select_pivot()
    {
        // select pivot from median of 9
        return med3char(
            med3char(get_direct(0),       get_direct(n/8),     get_direct(n/4)),
            med3char(get_direct(n/2-n/8), get_direct(n/2),     get_direct(n/2+n/8)),
            med3char(get_direct(n-1-n/4), get_direct(n-1-n/8), get_direct(n-3))
            );
    }

    StrCacheBlock* get_block(unsigned int &fill)
    {
        unsigned int blk;

        do {
            blk = block_current;
            DBG(debug_blocks, "current input block: " << blk);

            // if block_current == block_count -> no block left.
            if ( blk == block_count ) return NULL;

            // atomic swap in incremented block, which reserves blk
        } while ( !block_current.compare_exchange_weak(blk, blk+1) );

        fill = std::min<size_t>(block_size, n - blk * block_size);

        DBG(debug_blocks, "reserved input block " << blk << " @ " << (strings + blk * block_size) << " fill " << fill);

        StrCacheBlock* scb = new StrCacheBlock;

        // TODO: maybe separate loops?

        // fill cache for processor's part
        string* str = strings + blk * block_size;
        for (StrCache* sc = scb->cache; sc != scb->cache + fill; ++sc, ++str) {
            sc->str = *str;
            sc->key = get_char<uint64_t>(*str, depth);
        }

        return scb;
    }
};

struct BlockSourceQueueUnequal
{
    string*             strings;
    size_t              n, depth;

    BlockQueueType*     block_queue;
    PivotKeyQueueType*  pivot_queue;

    BlockSourceQueueUnequal(string* _strings, size_t _n, size_t _depth,
                            BlockQueueType* _blocks, PivotKeyQueueType* _pivots)
        : strings(_strings), n(_n), depth(_depth),
          block_queue(_blocks), pivot_queue(_pivots)
    {
    }

    ~BlockSourceQueueUnequal()
    {
        delete block_queue;
    }

    key_type select_pivot()
    {
        // extract all cached pivot keys
        std::vector<key_type> pivots;
        pivots.reserve( n / block_size + 1 );

        key_type k;
        while( pivot_queue->try_pop(k) ) {
            pivots.push_back(k);
        }
        delete pivot_queue;

        size_t p = pivots.size();

        return med3char(
            med3char(pivots[0],       pivots[p/8],     pivots[p/4]),
            med3char(pivots[p/2-p/8], pivots[p/2],     pivots[p/2+p/8]),
            med3char(pivots[p-1-p/4], pivots[p-1-p/8], pivots[p-3])
            );
    }

    StrCacheBlock* get_block(unsigned int &fill)
    {
        StrCacheBlock* blk;
        if (!block_queue->try_pop(blk)) return NULL;

        DBG(debug_blocks, "pop()ed input block " << blk << " fill " << blk->fill);
        fill = blk->fill;

        return blk;
    }
};

struct BlockSourceQueueEqual
{
    string*             strings;
    size_t              n, depth;

    BlockQueueType*     block_queue;
    PivotStrQueueType*  pivot_queue;

    BlockSourceQueueEqual(string* _strings, size_t _n, size_t _depth,
                          BlockQueueType* _blocks, PivotStrQueueType* _pivots)
        : strings(_strings), n(_n), depth(_depth),
          block_queue(_blocks), pivot_queue(_pivots)
    {
    }

    ~BlockSourceQueueEqual()
    {
        delete block_queue;
    }

    key_type select_pivot()
    {
        // extract all cached pivot keys
        std::vector<key_type> pivots;
        pivots.reserve( n / block_size + 1 );

        string s;
        while( pivot_queue->try_pop(s) ) {
            pivots.push_back( get_char<key_type>(s, depth) );
        }
        delete pivot_queue;

        size_t p = pivots.size();

        return med3char(
            med3char(pivots[0],       pivots[p/8],     pivots[p/4]),
            med3char(pivots[p/2-p/8], pivots[p/2],     pivots[p/2+p/8]),
            med3char(pivots[p-1-p/4], pivots[p-1-p/8], pivots[p-3])
            );
    }

    StrCacheBlock* get_block(unsigned int &fill)
    {
        StrCacheBlock* blk;
        if (!block_queue->try_pop(blk)) return NULL;

        DBG(debug_blocks, "pop()ed input block " << blk << " fill " << blk->fill);
        fill = blk->fill;

        // refill key cache
        for (unsigned int i = 0; i < fill; ++i)
        {
            blk->cache[i].key = get_char<key_type>(blk->cache[i].str, depth);
        }

        return blk;
    }
};

// ****************************************************************************
// *** ParallelMKQS - split ternary with 8-byte superalphabet

enum { LT, EQ, GT };

template <typename BlockSource>
struct ParallelMKQS
{
    // *** Class Attributes

    BlockSource         blks;
    uint64_t            pivot;

    unsigned int        procs;
    std::atomic<unsigned int> pwork;

    // *** PartitionJob for JobQueue

    struct PartitionJob : public Job
    {
        ParallelMKQS*       step;
        unsigned int        p;

        inline PartitionJob(ParallelMKQS* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual void run(JobQueue& jobqueue)
        {
            step->partition(p, jobqueue);
        }
    };

    // *** Lock-free Queues for Finished Blocks

    // block queues
    BlockQueueType      *oblk_lt, *oblk_eq, *oblk_gt;

    // queue for one potential pivot from each block
    PivotKeyQueueType*  oblk_lt_pivot;
    PivotStrQueueType*  oblk_eq_pivot;
    PivotKeyQueueType*  oblk_gt_pivot;

    // counters for determinting final output position
    std::atomic<size_t> count_lt, count_eq;

    // *** Constructor

    inline ParallelMKQS(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
        : blks(strings, n, depth),
          pivot( blks.select_pivot() ),
          // contruct queues
          oblk_lt( new BlockQueueType() ),
          oblk_eq( new BlockQueueType() ),
          oblk_gt( new BlockQueueType() ),
          oblk_lt_pivot( new PivotKeyQueueType() ),
          oblk_eq_pivot( new PivotStrQueueType() ),
          oblk_gt_pivot( new PivotKeyQueueType() ),
          count_lt(0), count_eq(0)
    {
        enqueue_jobs(jobqueue);
    }

    // for BlockSourceQueueUnequal
    inline ParallelMKQS(JobQueue& jobqueue, string* strings, size_t n, size_t depth,
                        BlockQueueType* iblk, PivotKeyQueueType* iblk_pivot)
        : blks(strings, n, depth, iblk, iblk_pivot),
          pivot( blks.select_pivot() ),
          // contruct queues
          oblk_lt( new BlockQueueType() ),
          oblk_eq( new BlockQueueType() ),
          oblk_gt( new BlockQueueType() ),
          oblk_lt_pivot( new PivotKeyQueueType() ),
          oblk_eq_pivot( new PivotStrQueueType() ),
          oblk_gt_pivot( new PivotKeyQueueType() ),
          count_lt(0), count_eq(0)
    {
        enqueue_jobs(jobqueue);
    }

    //  for BlockSourceQueueEqual
    inline ParallelMKQS(JobQueue& jobqueue, string* strings, size_t n, size_t depth,
                        BlockQueueType* iblk, PivotStrQueueType* iblk_pivot)
        : blks(strings, n, depth, iblk, iblk_pivot),
          pivot( blks.select_pivot() ),
          // contruct queues
          oblk_lt( new BlockQueueType() ),
          oblk_eq( new BlockQueueType() ),
          oblk_gt( new BlockQueueType() ),
          oblk_lt_pivot( new PivotKeyQueueType() ),
          oblk_eq_pivot( new PivotStrQueueType() ),
          oblk_gt_pivot( new PivotKeyQueueType() ),
          count_lt(0), count_eq(0)
    {
        enqueue_jobs(jobqueue);
    }

    void enqueue_jobs(JobQueue& jobqueue)
    {
        procs = blks.n / g_sequential_threshold;
        if (procs == 0) procs = 1;

        DBG(debug_parajobs, "ParallelMKQS on area " << srange(blks.strings, blks.n) << " with " << procs << " threads @ job " << this);

        // create partition jobs
        pwork = procs;
        for (unsigned int p = 0; p < procs; ++p)
            jobqueue.enqueue( new PartitionJob(this, p) );
    }

    // *** Helper to Output to One of the oblk Queues

    inline void oblk_push(const int type, StrCacheBlock*& blk)
    {
        // hopefully this function will be optimized, a templated version is very obfuscated.
        if (type == LT) {
            count_lt += blk->fill;
            oblk_lt_pivot->push( blk->key(blk->fill/2) );
            oblk_lt->push(blk);
        }
        else if (type == EQ) {
            count_eq += blk->fill;
            oblk_eq_pivot->push( blk->str(blk->fill/2) );
            oblk_eq->push(blk);
        }
        else {
            oblk_gt_pivot->push( blk->key(blk->fill/2) );
            oblk_gt->push(blk);
        }
    }

    // *** Main partition() function and collector
    void partition(unsigned int p, JobQueue& jobqueue);
    void partition_finished(JobQueue& jobqueue);
};

// *** Class Representing a Current Block during partition()

template <typename ParallelMKQS>
struct PartitionBlock
{
    unsigned int pos, fill;
    StrCacheBlock* blk;
    bool partial;

    PartitionBlock() : pos(0), fill(0), blk(NULL) { }

    template <int type>
    inline bool has_src_block(ParallelMKQS& mkqs)
    {
        if (pos < fill) return true; // still element in source block

        // try to fetch full block from BlockSource
        unsigned int newfill = 0;
        StrCacheBlock* newblk = mkqs.blks.get_block(newfill);

        if (newblk || fill == block_size)
        {
            // save existing block to correct output queue
            if (blk) {
                assert(pos == fill);
                blk->fill = fill, mkqs.oblk_push(type, blk);
            }
            // switch to new block
            return (pos = 0, fill = newfill, blk = newblk);
        }
        else
        {
            // keep block, as it has free space
            return false;
        }
    }

    template <int type>
    inline void check_partial_block(ParallelMKQS& mkqs)
    {
        if (blk && pos < block_size) return; // blk has free space
        // allocate free partial block
        if (blk) blk->fill = fill, mkqs.oblk_push(type, blk);
        pos = fill = 0;
        blk = new StrCacheBlock;
        partial = true;
    }

    inline StrCache& front_cache()
    {
        assert(blk && pos < fill && pos < block_size);
        return blk->cache[pos];
    }

    inline const key_type& front_key() const
    {
        assert(blk && pos < fill && pos < block_size);
        return blk->key(pos);
    }

    inline StrCache& back_cache()
    {
        assert(blk && fill > 0 && fill-1 < block_size);
        return blk->cache[fill-1];
    }

    template <int type>
    inline void swap_or_move_to(ParallelMKQS& mkqs, PartitionBlock& from)
    {
        if (!partial && pos < block_size) {
            if (pos < fill) {
                DBG(debug_cmp2, "swap with unpartitioned blk[" << pos << "] = " << front_key() << ".");
                std::swap(from.front_cache(), front_cache()), pos++;
            }
            else {
                DBG(debug_cmp2, "move to free-area at blk[" << pos << "].");
                assert(fill < block_size);
                fill++, front_cache() = from.front_cache(), pos++;
                from.front_cache() = from.back_cache(), from.fill--; // move back element in from
            }
        }
        else {
            DBG(debug_cmp2, "move to partial blk[" << pos << "].");
            check_partial_block<type>(mkqs);
            fill++, front_cache() = from.front_cache(), pos++;
            from.front_cache() = from.back_cache(), from.fill--; // move back element in from
        }
    }

    template <int type>
    inline void finish_partial(ParallelMKQS& mkqs, PartitionBlock& lt, PartitionBlock& eq, PartitionBlock& gt)
    {
        if (!blk || partial) return;

        key_type& pivot = mkqs.pivot;

        while ( pos < fill )
        {
            int res = cmp(front_key(), pivot);
            if (res < 0) { // < than pivot
                DBG(debug_cmp2, "blk[" << pos << "] = " << front_key() << " < pivot " << toHex(pivot) << ".");

                if (type == LT) { // good.
                    pos++;
                }
                else { // need to move to LT.
                    lt.swap_or_move_to<LT>(mkqs, *this);
                }
            }
            else if (res == 0) { // = pivot
                DBG(debug_cmp2, "blk[" << pos << "] = " << front_key() << " = pivot " << toHex(pivot) << ".");

                if (type == EQ) { // good.
                    pos++;
                }
                else { // need to move to EQ.
                    eq.swap_or_move_to<EQ>(mkqs, *this);
                }
            }
            else { // > than pivot
                DBG(debug_cmp2, "blk[" << pos << "] = " << front_key() << " > pivot " << toHex(pivot) << ".");

                if (type == GT) { // good.
                    pos++;
                }
                else { // need to move to GT.
                    gt.swap_or_move_to<GT>(mkqs, *this);
                }
            }
        }
    }
};

// *** partition() separating a BlockSource into three Queues lt, eq and gt.

template <typename BlockSource>
void ParallelMKQS<BlockSource>::partition(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_parajobs, "process PartitionJob " << p << " @ " << this << " with pivot " << toHex(pivot));

    // phase 1: partition blocks in-place

    PartitionBlock<ParallelMKQS> lt, eq, gt;

    while (1)
    {
        while ( lt.template has_src_block<LT>(*this) && eq.template has_src_block<EQ>(*this) )
        {
            int res = cmp(lt.front_key(), pivot);
            if (res < 0) { // < than pivot
                DBG(debug_cmp1, "blk_lt[" << lt.pos << "] = " << lt.front_key() << " < pivot " << toHex(pivot) << ", continue.");
                lt.pos++;
            }
            else if (res == 0) { // = pivot
                DBG(debug_cmp1, "blk_lt[" << lt.pos << "] = " << lt.front_key() << " = pivot " << toHex(pivot) << ", swap to blk_eq");
                std::swap(lt.front_cache(), eq.front_cache());
                eq.pos++;
            }
            else { // > than pivot
                DBG(debug_cmp1, "blk_lt[" << lt.pos << "] = " << lt.front_key() << " > pivot " << toHex(pivot) << ", break.");
                goto jump1;
            }
        }
        break;

    jump1:
        while ( gt.template has_src_block<GT>(*this) && eq.template has_src_block<EQ>(*this) )
        {
            int res = cmp(gt.front_key(), pivot);
            if (res < 0) { // < than pivot
                DBG(debug_cmp1, "blk_gt[" << gt.pos << "] = " << gt.front_key() << " < pivot " << toHex(pivot) << ", break.");
                goto jump2;
            }
            else if (res == 0) { // = pivot
                DBG(debug_cmp1, "blk_gt[" << gt.pos << "] = " << gt.front_key() << " = pivot " << toHex(pivot) << ", swap to blk_eq");
                std::swap(gt.front_cache(), eq.front_cache());
                eq.pos++;
            }
            else { // > than pivot
                DBG(debug_cmp1, "blk_gt[" << gt.pos << "] = " << gt.front_key() << " > pivot " << toHex(pivot) << ", continue.");
                gt.pos++;
            }
        }
        break;

    jump2:
        DBG(debug_cmp1, "swap blk_lt[" << lt.pos << "] = " << lt.front_key()
            <<          " and blk_gt[" << gt.pos << "] = " << gt.front_key());

        assert( lt.front_key() > pivot && gt.front_key() < pivot );
        std::swap(lt.front_cache(), gt.front_cache());
        lt.pos++, gt.pos++;
    }

    DBG(debug, "finished full blocks, creating partials @ " << this);
    DBG(debug, "lt " << lt.blk << " eq " << eq.blk << " gt " << gt.blk);

    lt.partial = (lt.blk == NULL), eq.partial = (eq.blk == NULL), gt.partial = (gt.blk == NULL);

    // phase 2: finish partitioning items by allocating extra blocks

    lt.template finish_partial<LT>(*this, lt, eq, gt);
    eq.template finish_partial<EQ>(*this, lt, eq, gt);
    gt.template finish_partial<GT>(*this, lt, eq, gt);

    if (lt.blk) lt.blk->fill = lt.fill, oblk_push(LT, lt.blk);
    if (eq.blk) eq.blk->fill = eq.fill, oblk_push(EQ, eq.blk);
    if (gt.blk) gt.blk->fill = gt.fill, oblk_push(GT, gt.blk);

    if (--pwork == 0)
        partition_finished(jobqueue);
}

template <typename BlockSource>
void ParallelMKQS<BlockSource>::partition_finished(JobQueue& jobqueue)
{
    DBG(debug_parajobs, "finished PartitionJobs @ " << this);
    DBG(debug_parajobs, "finished partitioning - " << count_lt << " lt "
        << count_eq << " eq " << blks.n - count_lt - count_eq << " gt - total " << blks.n);

    if (debug_selfcheck)
    {
        size_t count = 0;

        for(BlockQueueType::iterator it = oblk_lt->unsafe_begin();
            it != oblk_lt->unsafe_end(); ++it)
        {
            std::cout << "fill_lt : " << (*it)->fill << "\n";
            for (unsigned int i = 0; i < (*it)->fill; ++i)
            {
                assert((*it)->key(i) < pivot);
            }
            count += (*it)->fill;
        }
        assert(count == count_lt);

        for(BlockQueueType::iterator it = oblk_eq->unsafe_begin();
            it != oblk_eq->unsafe_end(); ++it)
        {
            std::cout << "fill_eq : " << (*it)->fill << "\n";
            for (unsigned int i = 0; i < (*it)->fill; ++i)
            {
                assert((*it)->key(i) == pivot);
            }
            count += (*it)->fill;
        }
        assert(count == count_lt + count_eq);

        for(BlockQueueType::iterator it = oblk_gt->unsafe_begin();
            it != oblk_gt->unsafe_end(); ++it)
        {
            std::cout << "fill_gt : " << (*it)->fill << "\n";
            for (unsigned int i = 0; i < (*it)->fill; ++i)
            {
                assert((*it)->key(i) > pivot);
            }
            count += (*it)->fill;
        }
        std::cout << "count " << count << " - " << blks.n << "\n";
        assert(count == blks.n);
    }

    // *** Create Recursive Jobs

    // recurse into lt-queue
    if (count_lt == 0) {
    }
    else if (count_lt <= g_sequential_threshold) {
        new SequentialMKQS<false>
            (jobqueue, blks.strings, count_lt, blks.depth,
             oblk_lt);
        delete oblk_lt_pivot;
    }
    else {
        new ParallelMKQS<BlockSourceQueueUnequal>
            (jobqueue, blks.strings, count_lt, blks.depth,
             oblk_lt, oblk_lt_pivot);
    }

    // recurse into eq-queue
    if (count_eq == 0) {
    }
    else if (count_eq <= g_sequential_threshold) {
        new SequentialMKQS<true>
            (jobqueue, blks.strings + count_lt, count_eq, blks.depth + sizeof(key_type),
             oblk_eq);
        delete oblk_eq_pivot;
    }
    else {
        new ParallelMKQS<BlockSourceQueueEqual>
            (jobqueue, blks.strings + count_lt, count_eq, blks.depth + sizeof(key_type),
             oblk_eq, oblk_eq_pivot);
    }

    // recurse into gt-queue
    size_t count_lteq = count_lt + count_eq;
    size_t count_gt = blks.n - count_lteq;
    if (count_gt == 0) {
    }
    else if (count_gt <= g_sequential_threshold) {
        new SequentialMKQS<false>
            (jobqueue, blks.strings + count_lteq, count_gt, blks.depth,
             oblk_gt);
        delete oblk_gt_pivot;
    }
    else {
        new ParallelMKQS<BlockSourceQueueUnequal>
            (jobqueue, blks.strings + count_lteq, count_gt, blks.depth,
             oblk_gt, oblk_gt_pivot);
    }

    delete this;
}

void parallel_mkqs(string* strings, size_t n)
{
    g_strings = strings;
    g_totalsize = n;
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);

    g_statscache >> "block_size" << block_size;

    JobQueue jobqueue;
    new ParallelMKQS<BlockSourceInput>(jobqueue, strings, n, 0);
    jobqueue.loop();
}

CONTESTANT_REGISTER_PARALLEL(parallel_mkqs,
                             "bingmann/parallel_mkqs",
                             "Parallel MKQS with blocks and cache8")

} // namespace bingmann_parallel_mkqs
