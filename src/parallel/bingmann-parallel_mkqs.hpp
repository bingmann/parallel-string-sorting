/*******************************************************************************
 * src/parallel/bingmann-parallel_mkqs.hpp
 *
 * Parallel multikey-quicksort.
 *
 *******************************************************************************
 * Copyright (C) 2013-2015 Timo Bingmann <tb@panthema.net>
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
 ******************************************************************************/

#ifndef PSS_SRC_PARALLEL_BINGMANN_PARALLEL_MKQS_HEADER
#define PSS_SRC_PARALLEL_BINGMANN_PARALLEL_MKQS_HEADER

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "../tools/debug.hpp"
#include "../tools/contest.hpp"
#include "../tools/stringtools.hpp"
#include "../tools/jobqueue.hpp"
#include "../tools/stringset.hpp"

#undef DBGX
#define DBGX DBGX_OMP

#include <tr1/memory>
#include <tbb/concurrent_queue.h>

namespace bingmann_parallel_mkqs {

using stringtools::toHex;
using namespace jobqueue;

static const size_t g_inssort_threshold = 64;

static const unsigned int block_size = 128 * 1024;

template <typename StringSet>
class ParallelMKQS
{
public:
    // *** Debugging Switching ***

    static const bool debug = false;
    static const bool debug_parajobs = false;
    static const bool debug_blocks = false;
    static const bool debug_cmp1 = false;
    static const bool debug_cmp2 = false;
    static const bool debug_seqjobs = false;

    static const bool debug_selfcheck = false;

    static const bool use_work_sharing = true;

    // *** Import StringSet's types ***

    typedef typename StringSet::String String;
    typedef typename StringSet::Iterator Iterator;
    typedef typename StringSet::CharIterator CharIterator;

    //! type used for characters in caches
    typedef uint64_t key_type;

    //! Context class holding global sorting algorithm variables
    class Context
    {
    public:
        //! total size of input
        const StringSet* g_strings;
        //! calculated threshold for sequential sorting
        size_t g_sequential_threshold;
        //! number of threads overall
        size_t g_threadnum;

        //! output string ranges for debugging
        std::string srange(const StringSet& ss)
        {
            std::ostringstream oss;
            oss << "[" << (ss.begin() - g_strings->begin())
                << "," << (ss.end() - g_strings->begin()) << ")=" << ss.size();
            return oss.str();
        }
    };

    /// *** Multi-character cache block for each string ***

    struct StrCache
    {
        uint64_t key;
        String   str;

        friend std::ostream& operator << (std::ostream& os, const StrCache& sc)
        {
            return os << toHex(sc.key);
        }
    };

    struct StrCacheBlock
    {
        unsigned int fill;
        StrCache     cache[block_size];

        String &     str(unsigned int i)
        {
            assert(i < block_size);
            return cache[i].str;
        }

        uint64_t &   key(unsigned int i)
        {
            assert(i < block_size);
            return cache[i].key;
        }
    };

    typedef tbb::concurrent_queue<StrCacheBlock*> BlockQueueType;
    typedef tbb::concurrent_queue<key_type> PivotKeyQueueType;
    typedef tbb::concurrent_queue<String> PivotStrQueueType;

    //! median of three (characters) and return reference to it.
    template <typename CharT>
    static inline const CharT&
    med3char(const CharT& a, const CharT& b, const CharT& c)
    {
        if (a == b) return a;
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

    //! median of three (characters) and return mutable reference to it.
    template <typename CharT>
    static inline CharT&
    med3charRef(CharT& a, CharT& b, CharT& c)
    {
        if (cmp(a, b) == 0) return a;
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

    // *************************************************************************
    // *** Rantala's Multikey-Quicksort with cached 8-byte superalphabet

    //! Insertion sort, ignores any cached characters.
    static inline void
    insertion_sort_nocache(const StringSet& strset, StrCache* cache,
                           size_t n, size_t depth)
    {
        StrCache* pi, * pj;
        String s, t;
        for (pi = cache + 1; --n > 0; ++pi) {
            String tmp = pi->str;
            for (pj = pi; pj > cache; --pj) {
                CharIterator s = strset.get_chars((pj - 1)->str, depth);
                CharIterator t = strset.get_chars(tmp, depth);

                while (*s == *t && *s != 0)
                    ++s, ++t;
                if (*s <= *t)
                    break;
                pj->str = (pj - 1)->str;
            }
            pj->str = tmp;
        }
    }

    //! Insertion sorts the strings only based on the cached characters.
    static inline void
    insertion_sort_cache_block(StrCache* cache, size_t n)
    {
        StrCache* pi, * pj;
        for (pi = cache + 1; --n > 0; ++pi) {
            StrCache tmp = *pi;
            for (pj = pi; pj > cache; --pj) {
                if (cmp((pj - 1)->key, tmp.key) <= 0)
                    break;
                *pj = *(pj - 1);
            }
            *pj = tmp;
        }
    }

    template <bool CacheDirty>
    static inline void
    insertion_sort(const StringSet& strset, StrCache* cache,
                   size_t n, size_t depth)
    {
        if (n == 0) return;
        if (CacheDirty) return insertion_sort_nocache(strset, cache, n, depth);

        insertion_sort_cache_block(cache, n);

        size_t start = 0, cnt = 1;
        for (size_t i = 0; i < n - 1; ++i) {
            if (cmp(cache[i], cache[i + 1]) == 0) {
                ++cnt;
                continue;
            }
            if (cnt > 1 && cache[start].key & 0xFF)
                insertion_sort_nocache(strset, cache + start, cnt,
                                       depth + sizeof(key_type));
            cnt = 1;
            start = i + 1;
        }
        if (cnt > 1 && cache[start].key & 0xFF)
            insertion_sort_nocache(strset, cache + start, cnt,
                                   depth + sizeof(key_type));
    }

    // *************************************************************************
    // *** SequentialMKQS - split ternary with 8-byte superalphabet

    class MKQSStep
    {
    public:
        StrCache* cache;
        size_t num_lt, num_eq, num_gt, n, depth;
        size_t idx;
        bool eq_recurse;

        //! constructor, but also actual algorithm which swaps around cache
        MKQSStep(const StringSet& ss, StrCache* cache,
                 size_t n, size_t depth, bool CacheDirty)
        {
            if (CacheDirty)
            {
                for (size_t i = 0; i < n; ++i)
                    cache[i].key = ss.get_uint64(cache[i].str, depth);
            }
            // Move pivot to first position to avoid wrapping the unsigned values
            // we are using in the main loop from zero to max.
            std::swap(cache[0],
                      med3charRef(
                          med3charRef(cache[0], cache[n / 8], cache[n / 4]),
                          med3charRef(cache[n / 2 - n / 8], cache[n / 2], cache[n / 2 + n / 8]),
                          med3charRef(cache[n - 1 - n / 4], cache[n - 1 - n / 8], cache[n - 3])
                          ));
            StrCache pivot = cache[0];
            size_t first = 1;
            size_t last = n - 1;
            size_t beg_ins = 1;
            size_t end_ins = n - 1;
            while (true) {
                while (first <= last) {
                    const int res = cmp(cache[first], pivot);
                    if (res > 0) {
                        break;
                    }
                    else if (res == 0) {
                        std::swap(cache[beg_ins++], cache[first]);
                    }
                    ++first;
                }
                while (first <= last) {
                    const int res = cmp(cache[last], pivot);
                    if (res < 0) {
                        break;
                    }
                    else if (res == 0) {
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
            const size_t num_eq_end = n - 1 - end_ins;
            num_eq = num_eq_beg + num_eq_end;
            num_lt = first - beg_ins;
            num_gt = end_ins - last;

            // Swap the equal pointers from the beginning to proper position.
            const size_t size1 = std::min(num_eq_beg, num_lt);
            std::swap_ranges(cache, cache + size1, cache + first - size1);
            // Swap the equal pointers from the end to proper position.
            const size_t size2 = std::min(num_eq_end, num_gt);
            std::swap_ranges(cache + first, cache + first + size2, cache + n - size2);

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

    //! magic deleter class for tr1::shared_ptr to delete [] array instead of
    //! delete array.
    template <typename T>
    struct shared_ptr_array_deleter {
        void operator () (T* p)
        { delete[] p; }
    };

    template <bool CacheDirty>
    class SequentialJob : public Job
    {
    public:
        Context& ctx;
        StringSet strset;
        size_t depth;

        BlockQueueType* block_queue;

        StrCache* cache;

        // reference counted pointer to base_cache, which might be shared among
        // threads due to work sharing.
        typedef std::shared_ptr<StrCache> StrCachePtrType;
        StrCachePtrType cache_base;

        // *** Constructors

        //! for testing via algorithm bingmann_sequential_mkqs_cache8
        SequentialJob(Context& _ctx,
                      const StringSet& _strset, size_t _depth,
                      StrCache* _cache)
            : ctx(_ctx), strset(_strset), depth(_depth),
              cache(_cache),
              cache_base(_cache, shared_ptr_array_deleter<StrCache>())
        { }

        SequentialJob(Context& _ctx, JobQueue& jobqueue,
                      const StringSet& _strset, size_t _depth,
                      BlockQueueType* _block_queue)
            : ctx(_ctx), strset(_strset), depth(_depth),
              block_queue(_block_queue), cache(NULL)
        {
            jobqueue.enqueue(this);
        }

        SequentialJob(Context& _ctx, JobQueue& jobqueue,
                      const StringSet& _strset, size_t _depth,
                      StrCache* _cache, StrCachePtrType _cache_base)
            : ctx(_ctx), strset(_strset), depth(_depth),
              cache(_cache), cache_base(_cache_base)
        {
            jobqueue.enqueue(this);
        }

        // *** Sequential Work

        virtual bool run(JobQueue& jobqueue)
        {
            DBG(debug_seqjobs, "SequentialJob "
                "for " << strset.size() << " strings @ " << this);

            if (!cache)
            {
                if (strset.size() == 1) // nothing to sort
                {
                    DBG(debug_seqjobs, "copy result to output string ptrs "
                        << ctx.srange(strset));

                    StrCacheBlock* scb = NULL;
                    while (block_queue->try_pop(scb))
                    {
                        assert(scb && (scb->fill == 0 || scb->fill == 1));
                        if (scb->fill == 1) {
                            strset[strset.begin()] = scb->cache[0].str;
                        }
                        delete scb;
                    }

                    while (block_queue->try_pop(scb))
                    {
                        assert(scb && scb->fill == 0);
                        delete scb;
                    }

                    assert(block_queue->empty());
                    delete block_queue;

                    return true;
                }

                // locally allocate (ptr,cache) array (NUMA-friendly)
                cache = new StrCache[strset.size()];
                cache_base = StrCachePtrType(cache, shared_ptr_array_deleter<StrCache>());

                // extract all elements, and maybe update cache

                StrCacheBlock* scb;
                size_t o = 0;

                while (block_queue->try_pop(scb))
                {
                    for (unsigned int i = 0; i < scb->fill; ++i, ++o)
                    {
                        if (CacheDirty) {
                            cache[o].str = scb->cache[i].str;
                            cache[o].key = strset.get_uint64(cache[o].str, depth);
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

            return true;
        }

        void sequential_mkqs(JobQueue& jobqueue)
        {
            DBG(debug_seqjobs, "SequentialJob on area "
                << ctx.srange(strset) << " @ job " << this);

            if (strset.size() < g_inssort_threshold)
            {
                insertion_sort<true>(strset, cache, strset.size(), depth);

                DBG(debug_seqjobs, "copy result to output string ptrs "
                    << ctx.srange(strset) << " @ job " << this);

                Iterator begin = strset.begin();
                for (size_t i = 0; i < strset.size(); ++i)
                    strset[begin + i] = cache[i].str;

                return;
            }

            // std::deque is much slower than std::vector, so we use an
            // artifical pop_front variable.
            size_t pop_front = 0;
            std::vector<MKQSStep> stack;
            stack.push_back(MKQSStep(strset, cache, strset.size(), depth, CacheDirty));

            // for output, save pointer to finished entry
            StrCache* cache_finished = cache + strset.size();

jumpout:
            while (stack.size() > pop_front)
            {
                while (stack.back().idx < 3)
                {
                    if (use_work_sharing && jobqueue.has_idle())
                    {
                        // convert top level of stack into independent jobs

                        MKQSStep& st = stack[pop_front];
                        Iterator st_strings = strset.begin() + (st.cache - cache);

                        DBG(debug_seqjobs, "Queueing front of SequentialJob's stack level "
                            << pop_front << ", idx " << st.idx
                            << ", areas lt "
                            << ctx.srange(strset.subr(st_strings, st.num_lt))
                            << " eq "
                            << ctx.srange(strset.subr(st_strings + st.num_lt, st.num_eq))
                            << " gt "
                            << ctx.srange(strset.subr(st_strings + st.num_lt + st.num_eq, st.num_gt))
                            << " @ job " << this);

                        if (st.idx == 0 && st.num_lt != 0)
                        {
                            DBG(debug_seqjobs, "Queueing job for lt-area "
                                << ctx.srange(strset.subr(st_strings, st.num_lt))
                                << " @ job " << this);

                            new SequentialJob<false>(
                                ctx, jobqueue,
                                strset.subr(st_strings, st.num_lt),
                                st.depth, st.cache, cache_base);
                        }
                        if (st.idx <= 1 && st.num_eq != 0)
                        {
                            DBG(debug_seqjobs, "Queueing job for eq-area "
                                << ctx.srange(strset.subr(st_strings + st.num_lt, st.num_eq))
                                << " @ job " << this);

                            if (st.eq_recurse) {
                                new SequentialJob<true>(
                                    ctx, jobqueue,
                                    strset.subr(st_strings + st.num_lt, st.num_eq),
                                    st.depth + sizeof(key_type),
                                    st.cache + st.num_lt, cache_base);
                            }
                            else
                            {
                                DBG(debug_seqjobs, "copy result to output string ptrs "
                                    << ctx.srange(strset.subr(st_strings + st.num_lt, st.num_eq))
                                    << " - no recurse equal @ job " << this);

                                // seems overkill to create an extra thread job for this:
                                for (size_t i = st.num_lt; i < st.num_lt + st.num_eq; ++i)
                                    st_strings[i] = st.cache[i].str;
                            }
                        }
                        if (st.idx <= 2 && st.num_gt != 0)
                        {
                            DBG(debug_seqjobs, "Queueing job for gt-area "
                                << ctx.srange(strset.subr(st_strings + st.num_lt + st.num_eq, st.num_gt))
                                << " @ job " << this);

                            new SequentialJob<false>(
                                ctx, jobqueue,
                                strset.subr(st_strings + st.num_lt + st.num_eq, st.num_gt),
                                st.depth,
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
                            insertion_sort<false>(
                                strset, ms.cache, ms.num_lt, ms.depth);
                        else
                            stack.push_back(MKQSStep(
                                                strset, ms.cache, ms.num_lt, ms.depth, false));
                    }
                    // process the eq-subsequence
                    else if (ms.idx == 2)
                    {
                        if (!ms.eq_recurse || ms.num_eq == 0)
                            continue;
                        else if (ms.num_eq < g_inssort_threshold)
                            insertion_sort<true>(
                                strset, ms.cache + ms.num_lt, ms.num_eq,
                                ms.depth + sizeof(key_type));
                        else
                            stack.push_back(MKQSStep(strset, ms.cache + ms.num_lt, ms.num_eq,
                                                     ms.depth + sizeof(key_type), true));
                    }
                    // process the gt-subsequence
                    else
                    {
                        assert(ms.idx == 3);
                        if (ms.num_gt == 0)
                            continue;
                        else if (ms.num_gt < g_inssort_threshold)
                            insertion_sort<false>(
                                strset, ms.cache + ms.num_lt + ms.num_eq, ms.num_gt,
                                ms.depth);
                        else
                            stack.push_back(MKQSStep(strset, ms.cache + ms.num_lt + ms.num_eq, ms.num_gt,
                                                     ms.depth, false));
                    }
                }

                stack.pop_back();
            }

            // copy string pointers to output array (for level pop_front)
            {
                size_t n_finished = cache_finished - cache;

                DBG(debug_seqjobs, "copy result to output string ptrs "
                    << ctx.srange(strset.subi(0, n_finished)) << " @ job " << this);

                Iterator begin = strset.begin();
                for (size_t i = 0; i < n_finished; ++i)
                    strset[begin + i] = cache[i].str;
            }
        }
    };

    // *************************************************************************
    // *** BlockSource - provide new blocks of unpartitioned input to ParallelJob

    class BlockSourceInput
    {
    public:
        StringSet strset;
        size_t depth;

        /// number of block in input
        unsigned int block_count;

        /// atomic index to currently unprocessed block of input
        std::atomic<unsigned int> block_current;

        BlockSourceInput(const StringSet& _strset, size_t _depth)
            : strset(_strset), depth(_depth),
              block_count((strset.size() + block_size - 1) / block_size), // round-up
              block_current(0)
        { }

        key_type get_direct(size_t i) const
        {
            assert(i < strset.size());
            return strset.get_uint64(strset[strset.begin() + i], depth);
        }

        key_type select_pivot()
        {
            size_t n = strset.size();
            // select pivot from median of 9
            return med3char(
                med3char(get_direct(0), get_direct(n / 8), get_direct(n / 4)),
                med3char(get_direct(n / 2 - n / 8), get_direct(n / 2), get_direct(n / 2 + n / 8)),
                med3char(get_direct(n - 1 - n / 4), get_direct(n - 1 - n / 8), get_direct(n - 3))
                );
        }

        StrCacheBlock * get_block(unsigned int& fill)
        {
            unsigned int blk;

            do {
                blk = block_current;
                DBG(debug_blocks, "current input block: " << blk);

                // if block_current == block_count -> no block left.
                if (blk == block_count) return NULL;

                // atomic swap in incremented block, which reserves blk
            } while (!block_current.compare_exchange_weak(blk, blk + 1));

            fill = std::min<size_t>(block_size, strset.size() - blk * block_size);

            DBG(debug_blocks, "reserved input block " << blk <<
                " @ " << (blk * block_size) << " fill " << fill);

            StrCacheBlock* scb = new StrCacheBlock;

            // TODO: maybe separate loops?

            // fill cache for processor's part
            Iterator str = strset.begin() + blk * block_size;
            for (StrCache* sc = scb->cache; sc != scb->cache + fill; ++sc, ++str) {
                sc->str = *str;
                sc->key = strset.get_uint64(*str, depth);
            }

            return scb;
        }
    };

    class BlockSourceQueueUnequal
    {
    public:
        StringSet strset;
        size_t depth;

        BlockQueueType* block_queue;
        PivotKeyQueueType* pivot_queue;

        BlockSourceQueueUnequal(const StringSet& _strset, size_t _depth,
                                BlockQueueType* _blocks, PivotKeyQueueType* _pivots)
            : strset(_strset), depth(_depth),
              block_queue(_blocks), pivot_queue(_pivots)
        { }

        ~BlockSourceQueueUnequal()
        {
            delete block_queue;
        }

        key_type select_pivot()
        {
            // extract all cached pivot keys
            std::vector<key_type> pivots;
            pivots.reserve(strset.size() / block_size + 1);

            key_type k;
            while (pivot_queue->try_pop(k)) {
                pivots.push_back(k);
            }
            delete pivot_queue;

            size_t p = pivots.size();

            return med3char(
                med3char(pivots[0], pivots[p / 8], pivots[p / 4]),
                med3char(pivots[p / 2 - p / 8], pivots[p / 2], pivots[p / 2 + p / 8]),
                med3char(pivots[p - 1 - p / 4], pivots[p - 1 - p / 8], pivots[p - 3])
                );
        }

        StrCacheBlock * get_block(unsigned int& fill)
        {
            StrCacheBlock* blk;
            if (!block_queue->try_pop(blk)) return NULL;

            DBG(debug_blocks, "pop()ed input block " << blk << " fill " << blk->fill);
            fill = blk->fill;

            return blk;
        }
    };

    class BlockSourceQueueEqual
    {
    public:
        StringSet strset;
        size_t depth;

        BlockQueueType* block_queue;
        PivotStrQueueType* pivot_queue;

        BlockSourceQueueEqual(const StringSet& _strset, size_t _depth,
                              BlockQueueType* _blocks, PivotStrQueueType* _pivots)
            : strset(_strset), depth(_depth),
              block_queue(_blocks), pivot_queue(_pivots)
        { }

        ~BlockSourceQueueEqual()
        {
            delete block_queue;
        }

        key_type select_pivot()
        {
            // extract all cached pivot keys
            std::vector<key_type> pivots;
            pivots.reserve(strset.size() / block_size + 1);

            String s;
            while (pivot_queue->try_pop(s)) {
                pivots.push_back(strset.get_uint64(s, depth));
            }
            delete pivot_queue;

            size_t p = pivots.size();

            return med3char(
                med3char(pivots[0], pivots[p / 8], pivots[p / 4]),
                med3char(pivots[p / 2 - p / 8], pivots[p / 2], pivots[p / 2 + p / 8]),
                med3char(pivots[p - 1 - p / 4], pivots[p - 1 - p / 8], pivots[p - 3])
                );
        }

        StrCacheBlock * get_block(unsigned int& fill)
        {
            StrCacheBlock* blk;
            if (!block_queue->try_pop(blk)) return NULL;

            DBG(debug_blocks, "pop()ed input block " << blk << " fill " << blk->fill);
            fill = blk->fill;

            // refill key cache
            for (unsigned int i = 0; i < fill; ++i)
            {
                blk->cache[i].key = strset.get_uint64(blk->cache[i].str, depth);
            }

            return blk;
        }
    };

    // *************************************************************************
    // *** ParallelJob - split ternary with 8-byte superalphabet

    enum { LT, EQ, GT };

    template <typename BlockSource>
    class ParallelJob
    {
    public:
        // *** Class Attributes

        Context& ctx;

        BlockSource blks;
        uint64_t pivot;

        unsigned int procs;
        std::atomic<unsigned int> pwork;

        // *** PartitionJob for JobQueue

        struct PartitionJob : public Job
        {
            ParallelJob  * step;
            unsigned int p;

            PartitionJob(ParallelJob* _step, unsigned int _p)
                : step(_step), p(_p) { }

            virtual bool run(JobQueue& jobqueue)
            {
                step->partition(p, jobqueue);
                return true;
            }
        };

        // *** Lock-free Queues for Finished Blocks

        // block queues
        BlockQueueType* oblk_lt, * oblk_eq, * oblk_gt;

        // queue for one potential pivot from each block
        PivotKeyQueueType* oblk_lt_pivot;
        PivotStrQueueType* oblk_eq_pivot;
        PivotKeyQueueType* oblk_gt_pivot;

        // counters for determinting final output position
        std::atomic<size_t> count_lt, count_eq;

        // *** Constructor

        ParallelJob(Context& _ctx, JobQueue& jobqueue,
                    const StringSet& strset, size_t depth)
            : ctx(_ctx),
              blks(strset, depth),
              pivot(blks.select_pivot()),
              // contruct queues
              oblk_lt(new BlockQueueType()),
              oblk_eq(new BlockQueueType()),
              oblk_gt(new BlockQueueType()),
              oblk_lt_pivot(new PivotKeyQueueType()),
              oblk_eq_pivot(new PivotStrQueueType()),
              oblk_gt_pivot(new PivotKeyQueueType()),
              count_lt(0), count_eq(0)
        {
            enqueue_jobs(jobqueue);
        }

        // for BlockSourceQueueUnequal
        ParallelJob(Context& _ctx, JobQueue& jobqueue,
                    const StringSet& strset, size_t depth,
                    BlockQueueType* iblk, PivotKeyQueueType* iblk_pivot)
            : ctx(_ctx),
              blks(strset, depth, iblk, iblk_pivot),
              pivot(blks.select_pivot()),
              // contruct queues
              oblk_lt(new BlockQueueType()),
              oblk_eq(new BlockQueueType()),
              oblk_gt(new BlockQueueType()),
              oblk_lt_pivot(new PivotKeyQueueType()),
              oblk_eq_pivot(new PivotStrQueueType()),
              oblk_gt_pivot(new PivotKeyQueueType()),
              count_lt(0), count_eq(0)
        {
            enqueue_jobs(jobqueue);
        }

        //  for BlockSourceQueueEqual
        ParallelJob(Context& _ctx, JobQueue& jobqueue,
                    const StringSet& strset, size_t depth,
                    BlockQueueType* iblk, PivotStrQueueType* iblk_pivot,
                    bool /* dummy */)
            : ctx(_ctx),
              blks(strset, depth, iblk, iblk_pivot),
              pivot(blks.select_pivot()),
              // contruct queues
              oblk_lt(new BlockQueueType()),
              oblk_eq(new BlockQueueType()),
              oblk_gt(new BlockQueueType()),
              oblk_lt_pivot(new PivotKeyQueueType()),
              oblk_eq_pivot(new PivotStrQueueType()),
              oblk_gt_pivot(new PivotKeyQueueType()),
              count_lt(0), count_eq(0)
        {
            enqueue_jobs(jobqueue);
        }

        void enqueue_jobs(JobQueue& jobqueue)
        {
            procs = blks.strset.size() / ctx.g_sequential_threshold;
            if (procs == 0) procs = 1;

            DBG(debug_parajobs, "ParallelJob on area "
                << ctx.srange(blks.strset)
                << " with " << procs << " threads @ job " << this);

            // create partition jobs
            pwork = procs;
            for (unsigned int p = 0; p < procs; ++p)
                jobqueue.enqueue(new PartitionJob(this, p));
        }

        // *** Helper to Output to One of the oblk Queues

        void oblk_push(const int type, StrCacheBlock*& blk)
        {
            // hopefully this function will be optimized, a templated version is
            // very obfuscated.
            if (type == LT) {
                count_lt += blk->fill;
                oblk_lt_pivot->push(blk->key(blk->fill / 2));
                oblk_lt->push(blk);
            }
            else if (type == EQ) {
                count_eq += blk->fill;
                oblk_eq_pivot->push(blk->str(blk->fill / 2));
                oblk_eq->push(blk);
            }
            else {
                oblk_gt_pivot->push(blk->key(blk->fill / 2));
                oblk_gt->push(blk);
            }
        }

        // *** Main partition() function and collector

        //! partition() separating a BlockSource into three Queues lt, eq and gt.
        void partition(unsigned int p, JobQueue& jobqueue)
        {
            DBG(debug_parajobs, "process PartitionJob " << p << " @ " << this << " with pivot " << toHex(pivot));

            // phase 1: partition blocks in-place

            PartitionBlock<ParallelJob> lt, eq, gt;

            while (1)
            {
                while (lt.template has_src_block<LT>(*this) && eq.template has_src_block<EQ>(*this))
                {
                    int res = cmp(lt.front_key(), pivot);
                    if (res < 0) {       // < than pivot
                        DBG(debug_cmp1, "blk_lt[" << lt.pos << "] = " <<
                            lt.front_key() << " < pivot " << toHex(pivot) << ", continue.");
                        lt.pos++;
                    }
                    else if (res == 0) { // = pivot
                        DBG(debug_cmp1, "blk_lt[" << lt.pos << "] = " <<
                            lt.front_key() << " = pivot " << toHex(pivot) << ", swap to blk_eq");
                        std::swap(lt.front_cache(), eq.front_cache());
                        eq.pos++;
                    }
                    else {  // > than pivot
                        DBG(debug_cmp1, "blk_lt[" << lt.pos << "] = " <<
                            lt.front_key() << " > pivot " << toHex(pivot) << ", break.");
                        goto jump1;
                    }
                }
                break;

jump1:
                while (gt.template has_src_block<GT>(*this) && eq.template has_src_block<EQ>(*this))
                {
                    int res = cmp(gt.front_key(), pivot);
                    if (res < 0) {       // < than pivot
                        DBG(debug_cmp1, "blk_gt[" << gt.pos << "] = " <<
                            gt.front_key() << " < pivot " << toHex(pivot) << ", break.");
                        goto jump2;
                    }
                    else if (res == 0) { // = pivot
                        DBG(debug_cmp1, "blk_gt[" << gt.pos << "] = " <<
                            gt.front_key() << " = pivot " << toHex(pivot) << ", swap to blk_eq");
                        std::swap(gt.front_cache(), eq.front_cache());
                        eq.pos++;
                    }
                    else {  // > than pivot
                        DBG(debug_cmp1, "blk_gt[" << gt.pos << "] = " <<
                            gt.front_key() << " > pivot " << toHex(pivot) << ", continue.");
                        gt.pos++;
                    }
                }
                break;

jump2:
                DBG(debug_cmp1, "swap blk_lt[" << lt.pos << "] = " <<
                    lt.front_key() << " and blk_gt[" << gt.pos << "] = " << gt.front_key());

                assert(lt.front_key() > pivot && gt.front_key() < pivot);
                std::swap(lt.front_cache(), gt.front_cache());
                lt.pos++, gt.pos++;
            }

            DBG(debug, "finished full blocks, creating partials @ " << this);
            DBG(debug, "lt " << lt.blk << " eq " << eq.blk << " gt " << gt.blk);

            lt.partial = (lt.blk == NULL);
            eq.partial = (eq.blk == NULL);
            gt.partial = (gt.blk == NULL);

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

        void partition_finished(JobQueue& jobqueue)
        {
            DBG(debug_parajobs, "finished PartitionJobs @ " << this);
            DBG(debug_parajobs, "finished partitioning - "
                << count_lt << " lt " << count_eq
                << " eq " << blks.strset.size() - count_lt - count_eq
                << " gt - total " << blks.strset.size());

            if (debug_selfcheck)
            {
                size_t count = 0;

                for (typename BlockQueueType::iterator it = oblk_lt->unsafe_begin();
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

                for (typename BlockQueueType::iterator it = oblk_eq->unsafe_begin();
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

                for (typename BlockQueueType::iterator it = oblk_gt->unsafe_begin();
                     it != oblk_gt->unsafe_end(); ++it)
                {
                    std::cout << "fill_gt : " << (*it)->fill << "\n";
                    for (unsigned int i = 0; i < (*it)->fill; ++i)
                    {
                        assert((*it)->key(i) > pivot);
                    }
                    count += (*it)->fill;
                }
                std::cout << "count " << count << " - " << blks.strset.size() << "\n";
                assert(count == blks.strset.size());
            }

            // *** Create Recursive Jobs

            // recurse into lt-queue
            if (count_lt == 0) { }
            else if (count_lt <= ctx.g_sequential_threshold) {
                new SequentialJob<false>(
                    ctx, jobqueue, blks.strset.subi(0, count_lt), blks.depth,
                    oblk_lt);
                delete oblk_lt_pivot;
            }
            else {
                new ParallelJob<BlockSourceQueueUnequal>(
                    ctx, jobqueue, blks.strset.subi(0, count_lt), blks.depth,
                    oblk_lt, oblk_lt_pivot);
            }

            // recurse into eq-queue
            if (count_eq == 0) { }
            else if (count_eq <= ctx.g_sequential_threshold) {
                new SequentialJob<true>(
                    ctx, jobqueue,
                    blks.strset.subi(count_lt, count_lt + count_eq),
                    blks.depth + sizeof(key_type),
                    oblk_eq);
                delete oblk_eq_pivot;
            }
            else {
                new ParallelJob<BlockSourceQueueEqual>(
                    ctx, jobqueue,
                    blks.strset.subi(count_lt, count_lt + count_eq),
                    blks.depth + sizeof(key_type),
                    oblk_eq, oblk_eq_pivot,
                    true);
            }

            // recurse into gt-queue
            size_t count_lteq = count_lt + count_eq;
            size_t count_gt = blks.strset.size() - count_lteq;
            if (count_gt == 0) { }
            else if (count_gt <= ctx.g_sequential_threshold) {
                new SequentialJob<false>(
                    ctx, jobqueue,
                    blks.strset.subi(count_lteq, count_lteq + count_gt),
                    blks.depth,
                    oblk_gt);
                delete oblk_gt_pivot;
            }
            else {
                new ParallelJob<BlockSourceQueueUnequal>(
                    ctx, jobqueue,
                    blks.strset.subi(count_lteq, count_lteq + count_gt),
                    blks.depth,
                    oblk_gt, oblk_gt_pivot);
            }

            delete this;
        }
    };

    // *** Class Representing a Current Block during partition()

    template <typename ParallelJob>
    class PartitionBlock
    {
    public:
        unsigned int pos, fill;
        StrCacheBlock* blk;
        bool partial;

        PartitionBlock() : pos(0), fill(0), blk(NULL) { }

        template <int type>
        bool has_src_block(ParallelJob& mkqs)
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
        void check_partial_block(ParallelJob& mkqs)
        {
            if (blk && pos < block_size) return; // blk has free space
            // allocate free partial block
            if (blk) blk->fill = fill, mkqs.oblk_push(type, blk);
            pos = fill = 0;
            blk = new StrCacheBlock;
            partial = true;
        }

        StrCache & front_cache()
        {
            assert(blk && pos < fill && pos < block_size);
            return blk->cache[pos];
        }

        const key_type & front_key() const
        {
            assert(blk && pos < fill && pos < block_size);
            return blk->key(pos);
        }

        StrCache & back_cache()
        {
            assert(blk && fill > 0 && fill - 1 < block_size);
            return blk->cache[fill - 1];
        }

        template <int type>
        void swap_or_move_to(ParallelJob& mkqs, PartitionBlock& from)
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
        void finish_partial(ParallelJob& mkqs,
                            PartitionBlock& lt, PartitionBlock& eq, PartitionBlock& gt)
        {
            if (!blk || partial) return;

            key_type& pivot = mkqs.pivot;

            while (pos < fill)
            {
                int res = cmp(front_key(), pivot);
                if (res < 0) {        // < than pivot
                    DBG(debug_cmp2, "blk[" << pos << "] = "
                                           << front_key() << " < pivot " << toHex(pivot) << ".");

                    if (type == LT) { // good.
                        pos++;
                    }
                    else {            // need to move to LT.
                        lt.swap_or_move_to<LT>(mkqs, *this);
                    }
                }
                else if (res == 0) {  // = pivot
                    DBG(debug_cmp2, "blk[" << pos << "] = "
                                           << front_key() << " = pivot " << toHex(pivot) << ".");

                    if (type == EQ) { // good.
                        pos++;
                    }
                    else {            // need to move to EQ.
                        eq.swap_or_move_to<EQ>(mkqs, *this);
                    }
                }
                else {                // > than pivot
                    DBG(debug_cmp2, "blk[" << pos << "] = "
                                           << front_key() << " > pivot " << toHex(pivot) << ".");

                    if (type == GT) { // good.
                        pos++;
                    }
                    else {            // need to move to GT.
                        gt.swap_or_move_to<GT>(mkqs, *this);
                    }
                }
            }
        }
    };
};

template <typename StringSet>
static inline void
bingmann_sequential_mkqs_cache8(const StringSet& ss, size_t depth)
{
    typedef ParallelMKQS<StringSet> MKQS;

    typename MKQS::Context ctx;
    ctx.g_strings = &ss;

    typename MKQS::StrCache * cache = new typename MKQS::StrCache[ss.size()];

    typename StringSet::Iterator begin = ss.begin();
    for (size_t i = 0; i < ss.size(); ++i)
        cache[i].str = ss[begin + i];

    JobQueue jobqueue;
    typename MKQS::template SequentialJob<true>(ctx, ss, depth, cache)
    .sequential_mkqs(jobqueue);
}

template <typename StringSet>
static inline
void bingmann_parallel_mkqs(const StringSet& strset, size_t depth)
{
    typedef ParallelMKQS<StringSet> MKQS;

    typename MKQS::Context ctx;
    ctx.g_strings = &strset;
    ctx.g_threadnum = omp_get_max_threads();
    ctx.g_sequential_threshold
        = std::max(g_inssort_threshold,
                   ctx.g_strings->size() / ctx.g_threadnum);

    g_stats >> "block_size" << block_size;

    JobQueue jobqueue;
    new typename MKQS::template ParallelJob<typename MKQS::BlockSourceInput>(
        ctx, jobqueue, strset, depth);
    jobqueue.loop();
}

} // namespace bingmann_parallel_mkqs

#endif // !PSS_SRC_PARALLEL_BINGMANN_PARALLEL_MKQS_HEADER

/******************************************************************************/
