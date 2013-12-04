/******************************************************************************
 * src/parallel/bingmann-parallel_sample_sort.h
 *
 * Parallel Super Scalar String Sample-Sort, many variant via different
 * Classifier templates.
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

#undef DBGX
#define DBGX DBGX_OMP

#include "../sequential/inssort.h"
#include "../sequential/bingmann-lcp_inssort.h"

namespace bingmann_parallel_sample_sort {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = false;

static const bool debug_splitter = false;
static const bool debug_bucketsize = false;
static const bool debug_recursion = false;
static const bool debug_splitter_tree = false;

static const bool use_work_sharing = true;

static const size_t MAXPROCS = 64+1; // +1 due to round up of processor number

static const size_t l2cache = 256*1024;

static const size_t g_smallsort_threshold = 64*1024;
static const size_t g_inssort_threshold = 64;

typedef uint64_t key_type;
typedef stringtools::StringPtrLcp StringPtr;

// ****************************************************************************
// *** Global Parallel Super Scalar String Sample Sort Context

class Context
{
public:
    //! total size of input
    size_t totalsize;

    //! number of remaining strings to sort
    std::atomic<size_t> restsize;

    //! number of threads overall
    size_t threadnum;

    //! counters
    size_t para_ss_steps, seq_ss_steps, bs_steps;

    //! job queue
    JobQueueT<Context> jobqueue;

    //! return sequential sorting threshold
    size_t sequential_threshold()
    {
        return std::max(g_smallsort_threshold, restsize / threadnum);
    }
};

typedef JobT<Context> Job;

// ****************************************************************************
// *** Classification Variants

unsigned char lcpKeyType(const key_type& a, const key_type& b)
{
    // XOR both values and count the number of zero bytes
    return count_high_zero_bits(a ^ b) / 8;
}

//! Recursive TreeBuilder for full-descent and unrolled variants, constructs a
//! binary tree and the plain splitter array.
template <size_t numsplitters>
struct TreeBuilderFullDescent
{
    key_type*       m_splitter;
    key_type*       m_tree;
    unsigned char*  m_lcp_iter;
    key_type*       m_samples;

    TreeBuilderFullDescent(key_type* splitter, key_type* splitter_tree, unsigned char* splitter_lcp,
                           key_type* samples, size_t samplesize)
        : m_splitter( splitter ),
          m_tree( splitter_tree ),
          m_lcp_iter( splitter_lcp ),
          m_samples( samples )
    {

        key_type sentinel = 0;
        recurse(samples, samples + samplesize, 1, sentinel);

        assert(m_splitter == splitter + numsplitters);
        assert(m_lcp_iter == splitter_lcp + numsplitters);
        splitter_lcp[0] &= 0x80; // overwrite sentinel lcp for first < everything bucket
        splitter_lcp[numsplitters] = 0; // sentinel for > everything bucket
    }

    ptrdiff_t snum(key_type* s) const
    {
        return (ptrdiff_t)(s - m_samples);
    }

    key_type recurse(key_type* lo, key_type* hi, unsigned int treeidx,
                     key_type& rec_prevkey)
    {
        DBG(debug_splitter, "rec_buildtree(" << snum(lo) << "," << snum(hi)
            << ", treeidx=" << treeidx << ")");

        // pick middle element as splitter
        key_type* mid = lo + (ptrdiff_t)(hi - lo) / 2;

        DBG(debug_splitter, "tree[" << treeidx << "] = samples[" << snum(mid) << "] = "
            << toHex(*mid));

        key_type mykey = m_tree[treeidx] = *mid;
#if 1
        key_type* midlo = mid;
        while (lo < midlo && *(midlo-1) == mykey) midlo--;

        key_type* midhi = mid;
        while (midhi+1 < hi && *midhi == mykey) midhi++;

        if (midhi - midlo > 1)
            DBG(0, "key range = [" << snum(midlo) << "," << snum(midhi) << ")");
#else
        key_type *midlo = mid, *midhi = mid+1;
#endif
        if (2 * treeidx < numsplitters)
        {
            key_type prevkey = recurse(lo, midlo, 2 * treeidx + 0, rec_prevkey);

            key_type xorSplit = prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " << toHex(prevkey) << " XOR " << toHex(mykey) << " = "
                << toHex(xorSplit) << " - " << count_high_zero_bits(xorSplit) << " bits = "
                << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            *m_splitter++ = mykey;

            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                | ((mykey & 0xFF) ? 0 : 0x80); // marker for done splitters

            return recurse(midhi, hi, 2 * treeidx + 1, mykey);
        }
        else
        {
            key_type xorSplit = rec_prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " << toHex(rec_prevkey) << " XOR " << toHex(mykey) << " = "
                << toHex(xorSplit) << " - " << count_high_zero_bits(xorSplit) << " bits = "
                << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            *m_splitter++ = mykey;

            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                | ((mykey & 0xFF) ? 0 : 0x80); // marker for done splitters

            return mykey;
        }
    }
};

// ****************************************************************************
// * Classification Subroutines for TreeBuilderFullDescent: rolled, and two
// * unrolled variants

template <size_t treebits>
struct ClassifySimple
{
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter[numsplitters];
    key_type splitter_tree[numsplitters+1];

    /// binary search on splitter array for bucket number
    inline unsigned int
    find_bkt_tree(const key_type& key) const
    {
        // binary tree traversal without left branch

        unsigned int i = 1;

        while ( i <= numsplitters )
        {
#if 0
            // in gcc-4.6.3 this produces a SETA, LEA sequence
            i = 2 * i + (key <= splitter_tree[i] ? 0 : 1);
#else
            // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which is slightly faster
            if (key <= splitter_tree[i])
                i = 2*i + 0;
            else // (key > splitter_tree[i])
                i = 2*i + 1;
#endif
        }

        i -= numsplitters+1;

        size_t b = i * 2;                                   // < bucket
        if (i < numsplitters && splitter[i] == key) b += 1; // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key);
            *bktout++ = b;
        }
    }

    //! return a splitter
    inline key_type get_splitter(unsigned int i) const
    { return splitter[i]; }

    /// build tree and splitter array from sample
    inline void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderFullDescent<numsplitters>(splitter, splitter_tree, splitter_lcp,
                                             samples, samplesize);
    }
};

template <size_t treebits>
struct ClassifyUnrollTree
{
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter[numsplitters];
    key_type splitter_tree[numsplitters+1];

    /// binary search on splitter array for bucket number
    __attribute__((optimize("unroll-all-loops")))
    inline unsigned int
    find_bkt_tree(const key_type& key) const
    {
        // binary tree traversal without left branch

        unsigned int i = 1;

        for (size_t l = 0; l < treebits; ++l)
        {
#if 0
            // in gcc-4.6.3 this produces a SETA, LEA sequence
            i = 2 * i + (key <= splitter_tree[i] ? 0 : 1);
#else
            // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which is slightly faster
            if (key <= splitter_tree[i])
                i = 2*i + 0;
            else // (key > splitter_tree[i])
                i = 2*i + 1;
#endif
        }

        i -= numsplitters+1;

        size_t b = i * 2;                                   // < bucket
        if (i < numsplitters && splitter[i] == key) b += 1; // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key);
            *bktout++ = b;
        }
    }

    //! return a splitter
    inline key_type get_splitter(unsigned int i) const
    { return splitter[i]; }

    /// build tree and splitter array from sample
    inline void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderFullDescent<numsplitters>(splitter, splitter_tree, splitter_lcp,
                                             samples, samplesize);
    }
};

template <size_t treebits>
struct ClassifyUnrollBoth : public ClassifyUnrollTree<treebits>
{
   /// search in splitter tree for bucket number, unrolled for U keys at once.
    template <int U>
    __attribute__((optimize("unroll-all-loops")))
    inline void
    find_bkt_tree_unroll(const key_type key[U], uint16_t obkt[U]) const
    {
        // binary tree traversal without left branch

        static const size_t numsplitters = this->numsplitters;
        const key_type* splitter = this->splitter;
        const key_type* splitter_tree = this->splitter_tree;

        unsigned int i[U];
        std::fill(i+0, i+U, 1);

        for (size_t l = 0; l < treebits; ++l)
        {
            for (int u = 0; u < U; ++u)
            {
#if 0
                // in gcc-4.6.3 this produces a SETA, LEA sequence
                i[u] = 2 * i[u] + (key[u] <= splitter_tree[i[u]] ? 0 : 1);
#else
                // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which is slightly faster
                if (key[u] <= splitter_tree[i[u]])
                    i[u] = 2*i[u] + 0;
                else // (key > splitter_tree[i[u]])
                    i[u] = 2*i[u] + 1;
#endif
            }
        }

        for (int u = 0; u < U; ++u)
            i[u] -= numsplitters+1;

        for (int u = 0; u < U; ++u)
            obkt[u] = i[u] * 2; // < bucket

        for (int u = 0; u < U; ++u)
        {
            if (i[u] < numsplitters && splitter[i[u]] == key[u]) obkt[u] += 1; // equal bucket
        }
    }

    /// classify all strings in area by walking tree and saving bucket id, unrolled loops
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            static const int rollout = 4;
            if (str + rollout < strE)
            {
                key_type key[rollout];
                key[0] = get_char<key_type>(str[0], depth);
                key[1] = get_char<key_type>(str[1], depth);
                key[2] = get_char<key_type>(str[2], depth);
                key[3] = get_char<key_type>(str[3], depth);

                find_bkt_tree_unroll<rollout>(key, bktout);

                str += rollout;
                bktout += rollout;
            }
            else
            {
                // binary search in splitter with equal check
                key_type key = get_char<key_type>(*str++, depth);

                unsigned int b = this->find_bkt_tree(key);
                *bktout++ = b;
            }
        }
    }
};

// ****************************************************************************
// *** Insertion Sort Template-switch

template <typename StringPtrType>
void insertion_sort(const StringPtrType& strptr, size_t n, size_t depth);

template <>
void insertion_sort<stringtools::StringPtr>(const stringtools::StringPtr& strptr, size_t n, size_t depth)
{
    inssort::inssort(strptr.active(), n, depth);
}

template <>
void insertion_sort<stringtools::StringPtrLcp>(const stringtools::StringPtrLcp& strptr, size_t n, size_t depth)
{
    bingmann_lcp_inssort::lcp_insertion_sort(strptr, n, depth);
}

// ****************************************************************************
// *** SampleSort non-recursive in-place sequential sample sort for small sorts

template <template <size_t> class Classify>
void Enqueue(Context& ctx, const StringPtr& strptr, size_t n, size_t depth);

template <template <size_t> class Classify, typename BktSizeType>
struct SmallsortJob : public Job
{
    StringPtr   strptr;
    size_t      n, depth;

    typedef BktSizeType bktsize_type;

    SmallsortJob(Context& ctx, const StringPtr& _strptr, size_t _n, size_t _depth)
        : strptr(_strptr), n(_n), depth(_depth)
    {
        DBG(debug_jobs, "enqueue SmallsortJob n=" << n << " depth=" << depth);

        ctx.jobqueue.enqueue(this);
    }

    struct SeqSampleSortStep
    {
#if 0
        static const size_t numsplitters = 2*16;
#else
        // bounding equations:
        // a) K * key_type splitter_tree size (and maybe K * key_type equal-cmp splitters)
        // b) 2 * K + 1 buckets (bktsize_type) when counting bkt occurances.
        static const size_t numsplitters2 = (l2cache - sizeof(BktSizeType)) / (2 * sizeof(BktSizeType));

        static const size_t treebits = logfloor_<numsplitters2>::value;
        static const size_t numsplitters = (1 << treebits) - 1;
#endif

        static const size_t bktnum = 2 * numsplitters + 1;

        StringPtr strptr;
        size_t idx;
        size_t depth;

        bktsize_type bkt[bktnum+1];

        unsigned char splitter_lcp[numsplitters+1];

        SeqSampleSortStep(Context& ctx, const StringPtr& _strptr, size_t n, size_t _depth, uint16_t* bktcache)
            : strptr(_strptr), idx(0), depth(_depth)
        {
            // step 1: select splitters with oversampling

            const size_t oversample_factor = 2;
            const size_t samplesize = oversample_factor * numsplitters;

            key_type samples[ samplesize ];

            string* strings = strptr.active();
            LCGRandom rng(strings);

            for (unsigned int i = 0; i < samplesize; ++i)
            {
                samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
            }

            std::sort(samples, samples + samplesize);

            Classify<treebits> classifier;

            classifier.build(samples, samplesize, splitter_lcp);

            // step 2: classify all strings

            classifier.classify(strings, strings+n, bktcache, depth);

            // step 2.5: count bucket sizes

            bktsize_type bktsize[bktnum];
            memset(bktsize, 0, bktnum * sizeof(bktsize_type));

            for (size_t si = 0; si < n; ++si)
                ++bktsize[ bktcache[si] ];

            // step 3: prefix sum

            bkt[0] = bktsize[0];
            bktsize_type last_bkt_size = bktsize[0];
            for (unsigned int i=1; i < bktnum; ++i) {
                bkt[i] = bkt[i-1] + bktsize[i];
                if (bktsize[i]) last_bkt_size = bktsize[i];
            }
            assert(bkt[bktnum-1] == n);

            // step 4: premute in-place

            for (size_t i=0, j; i < n - last_bkt_size; )
            {
                string perm = strings[i];
                uint16_t permbkt = bktcache[i];

                while ( (j = --bkt[ permbkt ]) > i )
                {
                    std::swap(perm, strings[j]);
                    std::swap(permbkt, bktcache[j]);
                }

                strings[i] = perm;
                i += bktsize[ permbkt ];
            }

            ++ctx.seq_ss_steps;
        }
    };

    // *** Stack of Recursive Sample Sort Steps

    uint16_t* bktcache;
    size_t bktcache_size;

    size_t ss_pop_front;
    std::vector<SeqSampleSortStep> ss_stack;

    virtual void run(Context& ctx)
    {
        DBG(debug_jobs, "Process SmallsortJob " << this << " of size " << n);

        bktcache = NULL;
        bktcache_size = 0;
        ss_pop_front = 0;
        ms_pop_front = 0;

        if (n >= g_smallsort_threshold && 0) // TODO
        {
            bktcache = new uint16_t[n];
            bktcache_size = n * sizeof(uint16_t);
            sort_sample_sort(ctx);
        }
        else
        {
            sort_mkqs_cache(ctx, strptr, n, depth);
        }
        delete [] bktcache;
    }

    void sort_sample_sort(Context& ctx)
    {
        typedef SeqSampleSortStep Step;

        // sort first level
        ss_pop_front = 0;
        ss_stack.push_back( Step(ctx,strptr,n,depth,bktcache) );

        // step 5: "recursion"

        while ( ss_stack.size() > ss_pop_front )
        {
            while ( ss_stack.back().idx < Step::bktnum )
            {
                Step& s = ss_stack.back();
                size_t i = s.idx++; // process the bucket s.idx

                size_t bktsize = s.bkt[i+1] - s.bkt[i];

                // i is even -> bkt[i] is less-than bucket
                if ((i & 1) == 0)
                {
                    if (bktsize == 0)
                        ;
                    else if (bktsize < g_smallsort_threshold)
                    {
                        assert(i/2 <= Step::numsplitters);
                        if (i == Step::bktnum-1)
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << bktsize << " no lcp");
                        else
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << bktsize << " lcp " << int(s.splitter_lcp[i/2] & 0x7F));

                        sort_mkqs_cache(ctx, s.strptr + s.bkt[i], bktsize, s.depth + (s.splitter_lcp[i/2] & 0x7F));
                    }
                    else
                    {
                        if (i == Step::bktnum-1)
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << bktsize << " no lcp");
                        else
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << bktsize << " lcp " << int(s.splitter_lcp[i/2] & 0x7F));

                        ss_stack.push_back( Step(ctx, s.strptr + s.bkt[i], bktsize, s.depth + (s.splitter_lcp[i/2] & 0x7F), bktcache) );
                    }
                }
                // i is odd -> bkt[i] is equal bucket
                else
                {
                    if (bktsize == 0)
                        ;
                    else if ( s.splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key, done.
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " is done!");
                        (s.strptr + s.bkt[i]).copy_back(bktsize);
                        ctx.restsize -= bktsize;
                    }
                    else if (bktsize < g_smallsort_threshold)
                    {
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " lcp keydepth!");

                        sort_mkqs_cache(ctx, s.strptr + s.bkt[i], bktsize, s.depth + sizeof(key_type));
                    }
                    else
                    {
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " lcp keydepth!");

                        ss_stack.push_back( Step(ctx, s.strptr + s.bkt[i], bktsize, s.depth + sizeof(key_type), bktcache) );
                    }
                }

                if (use_work_sharing && ctx.jobqueue.has_idle())
                    sample_sort_free_work(ctx);
            }

            ss_stack.pop_back();
        }
    }

    void sample_sort_free_work(Context& ctx)
    {
        assert(ss_stack.size() >= ss_pop_front);

        if (ss_stack.size() == ss_pop_front) {
            // ss_stack is empty, check other stack
            //return radix_sort_free_work(ctx);
            return mkqs_free_work(ctx);
        }

        // convert top level of stack into independent jobs
        DBG(debug_jobs, "Freeing top level of SmallsortJob's sample_sort stack");

        typedef SeqSampleSortStep Step;
        Step& s = ss_stack[ss_pop_front];

        while ( s.idx < Step::bktnum )
        {
            size_t i = s.idx++; // process the bucket s.idx

            size_t bktsize = s.bkt[i+1] - s.bkt[i];

            // i is even -> bkt[i] is less-than bucket
            if ((i & 1) == 0)
            {
                if (bktsize == 0)
                    ;
                else
                {
                    if (i == Step::bktnum-1)
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << bktsize << " no lcp");
                    else
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << bktsize << " lcp " << int(s.splitter_lcp[i/2] & 0x7F));

                    Enqueue<Classify>( ctx, s.strptr + s.bkt[i], bktsize, s.depth + (s.splitter_lcp[i/2] & 0x7F) );
                }
            }
            // i is odd -> bkt[i] is equal bucket
            else
            {
                if (bktsize == 0)
                    ;
                else if ( s.splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key, done.
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " is done!");
                    (s.strptr + s.bkt[i]).copy_back(bktsize);
                    ctx.restsize -= bktsize;
                }
                else
                {
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " lcp keydepth!");

                    Enqueue<Classify>( ctx, s.strptr + s.bkt[i], bktsize, s.depth + sizeof(key_type) );
                }
            }
        }

        // shorten the current stack
        ++ss_pop_front;
    }

    // *********************************************************************************
    // *** Stack of Recursive MKQS Steps

    static inline int cmp(const key_type& a, const key_type& b)
    {
        return (a > b) ? 1 : (a < b) ? -1 : 0;
    }

    template <typename Type>
    static inline size_t
    med3(Type* A, size_t i, size_t j, size_t k)
    {
        if (A[i] == A[j])                 return i;
        if (A[k] == A[i] || A[k] == A[j]) return k;
        if (A[i] < A[j]) {
            if (A[j] < A[k]) return j;
            if (A[i] < A[k]) return k;
            return i;
        }
        else {
            if (A[j] > A[k]) return j;
            if (A[i] < A[k]) return i;
            return k;
        }
    }

    // Insertion sort the strings only based on the cached characters.
    static inline void
    insertion_sort_cache_block(const StringPtr& strptr, key_type* cache, int n)
    {
        string* strings = strptr.active();
        unsigned int pi, pj;
        for (pi = 1; --n > 0; ++pi) {
            string tmps = strings[pi];
            key_type tmpc = cache[pi];
            for (pj = pi; pj > 0; --pj) {
                if (cmp(cache[pj-1], tmpc) <= 0)
                    break;
                strings[pj] = strings[pj-1];
                cache[pj] = cache[pj-1];
            }
            strings[pj] = tmps;
            cache[pj] = tmpc;
        }
    }

    // Insertion sort, but use cached characters if possible.
    template <bool CacheDirty>
    static inline void
    insertion_sort_cache(const StringPtr& strptr, key_type* cache, int n, size_t depth)
    {
        if (n <= 1) return;
        if (CacheDirty) return insertion_sort(strptr, n, depth);

        insertion_sort_cache_block(strptr, cache, n);

        size_t start = 0, bktsize = 1;
        for (int i = 0; i < n - 1; ++i) {
            if (cache[i] == cache[i+1]) {
                ++bktsize;
                continue;
            }
            if (bktsize > 1 && cache[start] & 0xFF)
                insertion_sort(strptr + start, bktsize, depth + sizeof(key_type));
            if (start != 0)
                strptr.lcp(start) = depth + lcpKeyType(cache[start-1], cache[start]);
            bktsize = 1;
            start = i+1;
        }
        if (bktsize > 1 && cache[start] & 0xFF)
            insertion_sort(strptr + start, bktsize, depth + sizeof(key_type));
        if (start != 0)
            strptr.lcp(start) = depth + lcpKeyType(cache[start-1], cache[start]);
    }

    struct MKQSStep
    {
        StringPtr strptr;
        key_type* cache;
        size_t num_lt, num_eq, num_gt, n, depth;
        size_t idx;
        bool eq_recurse;

        MKQSStep(Context& ctx, const StringPtr& _strptr, key_type* _cache, size_t _n, size_t _depth, bool CacheDirty)
            : strptr(_strptr), cache(_cache), n(_n), depth(_depth), idx(0)
        {
            string* strings = strptr.active();

            if (CacheDirty) {
                for (size_t i = 0; i < n; ++i) {
                    cache[i] = get_char<key_type>(strings[i], depth);
                }
            }
            // select median of 9
            size_t p = med3(cache,
                            med3(cache, 0,       n/8,     n/4),
                            med3(cache, n/2-n/8, n/2,     n/2+n/8),
                            med3(cache, n-1-n/4, n-1-n/8, n-3)
                );
            // swap pivot to first position
            std::swap(strings[0], strings[p]);
            std::swap(cache[0], cache[p]);
            // save the pivot value
            key_type pivot = cache[0];
            // for immediate LCP calculation
            key_type max_lt = 0, min_gt = std::numeric_limits<key_type>::max();
            // indexes into array: 0 [pivot] 1 [===] leq [<<<] llt [???] rgt [>>>] req [===] n-1
            size_t leq = 1, llt = 1, rgt = n-1, req = n-1;
            while (true)
            {
                while (llt <= rgt)
                {
                    int r = cmp(cache[llt], pivot);
                    if (r > 0) {
                        min_gt = std::min(min_gt, cache[llt]);
                        break;
                    }
                    else if (r == 0) {
                        std::swap(strings[leq], strings[llt]);
                        std::swap(cache[leq], cache[llt]);
                        leq++;
                    }
                    else {
                        max_lt = std::max(max_lt, cache[llt]);
                    }
                    ++llt;
                }
                while (llt <= rgt)
                {
                    int r = cmp(cache[rgt], pivot);
                    if (r < 0) {
                        max_lt = std::max(max_lt, cache[rgt]);
                        break;
                    }
                    else if (r == 0) {
                        std::swap(strings[req], strings[rgt]);
                        std::swap(cache[req], cache[rgt]);
                        req--;
                    }
                    else {
                        min_gt = std::min(min_gt, cache[rgt]);
                    }
                    --rgt;
                }
                if (llt > rgt)
                    break;
                std::swap(strings[llt], strings[rgt]);
                std::swap(cache[llt], cache[rgt]);
                ++llt;
                --rgt;
            }
            // calculate size of areas = < and >, save into struct
            size_t num_leq = leq, num_req = n-1-req;
            num_eq = num_leq + num_req;
            num_lt = llt - leq;
            num_gt = req - rgt;
            assert(num_eq > 0);
            assert(num_lt + num_eq + num_gt == n);

            // swap equal values from left to center
            const size_t size1 = std::min(num_leq, num_lt);
            std::swap_ranges(strings, strings+size1, strings+llt-size1);
            std::swap_ranges(cache, cache+size1, cache+llt-size1);

            // swap equal values from right to center
            const size_t size2 = std::min(num_req, num_gt);
            std::swap_ranges(strings+llt, strings+llt+size2, strings+n-size2);
            std::swap_ranges(cache+llt, cache+llt+size2, cache+n-size2);

            // No recursive sorting if pivot has a zero byte
            this->eq_recurse = (pivot & 0xFF);

            // for immediate LCP calculation, max and min could also be
            // calculated in the loop above.
            if (num_lt > 0)
            {
                //key_type _max_lt = *std::max_element(cache + 0, cache + num_lt);
                //assert(max_lt == _max_lt);

                unsigned int lcp = depth + lcpKeyType(max_lt, pivot);
                //std::cout << "lcp with pivot: " << lcp << "\n";

                strptr.lcp(num_lt) = lcp;
            }
            if (num_gt > 0)
            {
                //key_type _min_gt = *std::min_element(cache + num_lt + num_eq, cache + n);
                //assert(min_gt == _min_gt);

                unsigned int lcp = depth + lcpKeyType(pivot, min_gt);
                //std::cout << "lcp with pivot: " << lcp << "\n";

                strptr.lcp(num_lt + num_eq) = lcp;
            }

            ++ctx.bs_steps;
        }
    };

    size_t ms_pop_front;
    size_t ms_depth; // mkqs sorting base depth
    std::vector<MKQSStep> ms_stack;

    void sort_mkqs_cache(Context& ctx, const StringPtr& strptr, size_t n, size_t depth)
    {
        if (n < g_inssort_threshold) {
            insertion_sort(strptr, n, depth);
            ctx.restsize -= n;
            strptr.copy_back(n);
            return;
        }

        if (bktcache_size < n * sizeof(key_type)) {
            delete [] bktcache;
            bktcache = (uint16_t*)new key_type[n];
            bktcache_size = n * sizeof(key_type);
        }

        key_type* cache = (key_type*)bktcache;

        // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
        ms_pop_front = 0;
        ms_depth = depth;
        ms_stack.clear();
        ms_stack.push_back( MKQSStep(ctx, strptr, cache, n, ms_depth, true) );

    jumpout:
        while ( ms_stack.size() > ms_pop_front )
        {
            while ( ms_stack.back().idx < 3 )
            {
                if (use_work_sharing && ctx.jobqueue.has_idle() && 0) { // TODO
                    sample_sort_free_work(ctx);
                    goto jumpout;
                }

                MKQSStep& ms = ms_stack.back();
                ++ms.idx; // increment here, because stack may change

                // process the lt-subsequence
                if (ms.idx == 1)
                {
                    if (ms.num_lt == 0)
                        ;
                    else if (ms.num_lt < g_inssort_threshold) {
                        insertion_sort_cache<false>(ms.strptr, ms.cache,
                                                    ms.num_lt, ms.depth);
                        ctx.restsize -= ms.num_lt;
                        ms.strptr.copy_back(ms.num_lt);
                    }
                    else {
                        ms_stack.push_back(
                            MKQSStep(ctx,
                                     ms.strptr, ms.cache,
                                     ms.num_lt, ms.depth, false) );
                    }
                }
                // process the eq-subsequence
                else if (ms.idx == 2)
                {
                    StringPtr sp = ms.strptr + ms.num_lt;

                    if (!ms.eq_recurse) {
                        sp.copy_back(ms.num_eq);
                        ctx.restsize -= ms.num_eq;
                    }
                    else if (ms.num_eq < g_inssort_threshold) {
                        insertion_sort_cache<true>(sp, ms.cache + ms.num_lt,
                                                   ms.num_eq, ms.depth + sizeof(key_type));
                        ctx.restsize -= ms.num_eq;
                        sp.copy_back(ms.num_eq);
                    }
                    else {
                        ms_stack.push_back(
                            MKQSStep(ctx, sp,
                                     ms.cache + ms.num_lt,
                                     ms.num_eq, ms.depth + sizeof(key_type), true) );
                    }
                }
                // process the gt-subsequence
                else
                {
                    StringPtr sp = ms.strptr + ms.num_lt + ms.num_eq;

                    assert(ms.idx == 3);
                    if (ms.num_gt == 0)
                        ;
                    else if (ms.num_gt < g_inssort_threshold) {
                        insertion_sort_cache<false>(sp, ms.cache + ms.num_lt + ms.num_eq,
                                                    ms.num_gt, ms.depth);
                        ctx.restsize -= ms.num_gt;
                        sp.copy_back(ms.num_gt);
                    }
                    else {
                        ms_stack.push_back(
                            MKQSStep(ctx, sp,
                                     ms.cache + ms.num_lt + ms.num_eq,
                                     ms.num_gt, ms.depth, false) );
                    }
                }
            }

            ms_stack.pop_back();
        }

        ms_pop_front = 0; // clear stack
        ms_stack.clear();
    }

    void mkqs_free_work(Context& ctx)
    {
        assert(ms_stack.size() >= ms_pop_front);

        if (ms_stack.size() == ms_pop_front) {
            return;
        }

        DBG(debug_jobs, "Freeing top level of SmallsortJob's mkqs stack");

        // convert top level of stack into independent jobs

        MKQSStep& st = ms_stack[ms_pop_front];

        if (st.idx == 0 && st.num_lt != 0)
        {
            Enqueue<Classify>(ctx, st.strptr, st.num_lt, st.depth);
        }
        if (st.idx <= 1) // st.num_eq > 0 always
        {
            if (st.eq_recurse) {
                Enqueue<Classify>(ctx, st.strptr + st.num_lt,
                                  st.num_eq, st.depth + sizeof(key_type));
            }
            else {
                (st.strptr + st.num_lt).copy_back(st.num_eq);
                ctx.restsize -= st.num_eq;
            }
        }
        if (st.idx <= 2 && st.num_gt != 0)
        {
            Enqueue<Classify>(ctx, st.strptr + st.num_lt + st.num_eq,
                              st.num_gt, st.depth);
        }

        // shorten the current stack
        ++ms_pop_front;
    }
};

// ****************************************************************************
// *** SampleSortStep out-of-place parallel sample sort with separate Jobs

template <template <size_t> class Classify>
struct SampleSortStep
{
#if 0
    static const size_t numsplitters = 2*16;
#else
    // bounding equation: 2*K+1 buckets when counting bkt occurances.
    static const size_t numsplitters2 = (l2cache - sizeof(size_t)) / (2 * sizeof(size_t));

    static const size_t treebits = logfloor_<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;
#endif

    static const size_t bktnum = 2 * numsplitters + 1;

    //! string pointers, size, and current sorting depth
    StringPtr           strptr;
    size_t              n, depth;

    //! number of parts into which the strings were split
    unsigned int        parts;
    //! size of all parts except the last
    size_t              psize;
    //! number of threads still working
    std::atomic<unsigned int> pwork;

    //! classifier instance and variables (contains splitter tree
    Classify<treebits>  classifier;
    //! LCPs of splitters, needed for recursive calls
    unsigned char       splitter_lcp[numsplitters+1];

    //! individual bucket array of threads, keep bkt[0] for DistributeJob
    size_t*             bkt[MAXPROCS];
    //! bucket ids cache, created by classifier and later counted
    uint16_t*           bktcache[MAXPROCS];

    // *** Classes for JobQueue

    struct SampleJob : public Job
    {
        SampleSortStep*     step;

        SampleJob(SampleSortStep* _step)
            : step(_step) { }

        virtual void run(Context& ctx)
        {
            step->sample(ctx);
        }
    };

    struct CountJob : public Job
    {
        SampleSortStep*     step;
        unsigned int        p;

        CountJob(SampleSortStep* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual void run(Context& ctx)
        {
            step->count(p, ctx);
        }
    };

    struct DistributeJob : public Job
    {
        SampleSortStep*     step;
        unsigned int        p;

        DistributeJob(SampleSortStep* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual void run(Context& ctx)
        {
            step->distribute(p, ctx);
        }
    };

    // *** Constructor

    SampleSortStep(Context& ctx, const StringPtr& _strptr, size_t _n, size_t _depth)
        : strptr(_strptr), n(_n), depth(_depth)
    {
        parts = (n + ctx.sequential_threshold()-1) / ctx.sequential_threshold();
        if (parts == 0) parts = 1;

        psize = (n + parts-1) / parts;
        DBG(debug_jobs, "enqueue SampleSortStep n=" << n << " parts=" << parts << " psize=" << psize);

        ctx.jobqueue.enqueue( new SampleJob(this) );
        ++ctx.para_ss_steps;
    }

    // *** Sample Step

    void sample(Context& ctx)
    {
        //DBG(debug_jobs, "Process SampleJob @ " << this);

        const size_t oversample_factor = 2;
        size_t samplesize = oversample_factor * numsplitters;

        string* strings = strptr.active();
        key_type samples[ samplesize ];

        LCGRandom rng(&samples);

        for (unsigned int i = 0; i < samplesize; ++i)
        {
            samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
        }

        std::sort(samples, samples + samplesize);

        classifier.build(samples, samplesize, splitter_lcp);

        // create new jobs
        pwork = parts;
        for (unsigned int p = 0; p < parts; ++p)
            ctx.jobqueue.enqueue( new CountJob(this, p) );
    }

    // *** Counting Step

    void count(unsigned int p, Context& ctx)
    {
        //DBG(debug_jobs, "Process CountJob " << p << " @ " << this);

        string* strB = strptr.active() + p * psize;
        string* strE = strptr.active() + std::min((p+1) * psize, n);
        if (strE < strB) strE = strB;

        uint16_t* mybktcache = bktcache[p] = new uint16_t[strE-strB];
        classifier.classify(strB, strE, mybktcache, depth);

        size_t* mybkt = bkt[p] = new size_t[bktnum + (p == 0 ? 1 : 0)];
        memset(mybkt, 0, bktnum * sizeof(size_t));

        for (uint16_t* bc = mybktcache; bc != mybktcache + (strE-strB); ++bc)
            ++mybkt[ *bc ];

        if (--pwork == 0)
            count_finished(ctx);
    }

    void count_finished(Context& ctx)
    {
        //DBG(debug_jobs, "Finishing CountJob " << this << " with prefixsum");

        // inclusive prefix sum over bkt
        size_t sum = 0;
        for (unsigned int i = 0; i < bktnum; ++i)
        {
            for (unsigned int p = 0; p < parts; ++p)
            {
                bkt[p][i] = (sum += bkt[p][i]);
            }
        }
        assert(sum == n);

        // create new jobs
        pwork = parts;
        for (unsigned int p = 0; p < parts; ++p)
            ctx.jobqueue.enqueue( new DistributeJob(this, p) );
    }

    // *** Distribute Step

    void distribute(unsigned int p, Context& ctx)
    {
        //DBG(debug_jobs, "Process DistributeJob " << p << " @ " << this);

        string* strB = strptr.active() + p * psize;
        string* strE = strptr.active() + std::min((p+1) * psize, n);
        if (strE < strB) strE = strB;

        string* sorted = strptr.shadow(); // get alternative shadow pointer array

        uint16_t* mybktcache = bktcache[p];
        size_t* mybkt = bkt[p];

        for (string* str = strB; str != strE; ++str, ++mybktcache)
            sorted[ --mybkt[ *mybktcache ] ] = *str;

        if (p != 0) // p = 0 is needed for recursion into bkts
            delete [] bkt[p];

        delete [] bktcache[p];

        if (--pwork == 0)
            distribute_finished(ctx);
    }

    void distribute_finished(Context& ctx)
    {
        //DBG(debug_jobs, "Finishing DistributeJob " << this << " with enqueuing subjobs");

        size_t* bkt = this->bkt[0];
        assert(bkt);

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
            else if (bktsize == 1) { // just one string pointer, copyback
                strptr.flip(bkt[i]).copy_back(1);
                ctx.restsize -= 1;
            }
            else
            {
                DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bkt[i] << " size " << bktsize << " lcp " << int(splitter_lcp[i/2] & 0x7F));
                Enqueue<Classify>(ctx, strptr.flip(bkt[i]), bktsize, depth + (splitter_lcp[i/2] & 0x7F));
            }
            ++i;
            // i is odd -> bkt[i] is equal bucket
            bktsize = bkt[i+1] - bkt[i];
            if (bktsize == 0)
                ;
            else if (bktsize == 1) { // just one string pointer, copyback
                strptr.flip(bkt[i]).copy_back(1);
                ctx.restsize -= 1;
            }
            else
            {
                if ( splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key, done.
                    DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " is done!");
                    strptr.flip(bkt[i]).copy_back(bktsize);
                    ctx.restsize -= bktsize;
                }
                else {
                    DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " lcp keydepth!");
                    Enqueue<Classify>(ctx, strptr.flip(bkt[i]), bktsize, depth + sizeof(key_type));
                }
            }
            ++i;
        }

        size_t bktsize = bkt[i+1] - bkt[i];

        if (bktsize == 0)
            ;
        else if (bktsize == 1) { // just one string pointer, copyback
            strptr.flip(bkt[i]).copy_back(1);
            ctx.restsize -= 1;
        }
        else
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bkt[i] << " size " << bktsize << " no lcp");
            Enqueue<Classify>(ctx, strptr.flip(bkt[i]), bktsize, depth);
        }

        delete [] bkt;
        delete this;
    }

    static inline void put_stats()
    {
        g_statscache >> "l2cache" << size_t(l2cache)
                     >> "splitter_treebits" << size_t(treebits)
                     >> "numsplitters" << size_t(numsplitters);
    }
};

template <template <size_t> class Classify>
void Enqueue(Context& ctx, const StringPtr& strptr, size_t n, size_t depth)
{
    if (n > ctx.sequential_threshold() && 0) { // TODO
        new SampleSortStep<Classify>(ctx, strptr, n, depth);
    }
    else {
        if (n < ((uint64_t)1 << 32))
            new SmallsortJob<Classify,uint32_t>(ctx, strptr, n, depth);
        else
            new SmallsortJob<Classify,uint64_t>(ctx, strptr, n, depth);
    }
}

template <template <size_t> class Classify>
void parallel_sample_sort_base(string* strings, size_t n, size_t depth)
{
    Context ctx;
    ctx.totalsize = ctx.restsize = n;
    ctx.threadnum = omp_get_max_threads();
    ctx.para_ss_steps = ctx.seq_ss_steps = ctx.bs_steps = 0;

    SampleSortStep<Classify>::put_stats();

    string* shadow = new string[n]; // allocate shadow pointer array
    size_t* lcp = new size_t[n];
    lcp[0] = -1;

    Enqueue<Classify>(ctx, StringPtr(strings, shadow, false, lcp), n, depth);
    ctx.jobqueue.loop(ctx);

    delete [] shadow;

    for (size_t i = 1; i < n; ++i)
    {
        string s1 = strings[i-1], s2 = strings[i];
        size_t h = 0;
        while (*s1 != 0 && *s1 == *s2)
            ++h, ++s1, ++s2;

        if (h != lcp[i]) {
            std::cout << "lcp[" << i << "] mismatch " << h << " != " << lcp[i] << std::endl;
        }
    }

    delete [] lcp;

    //assert(ctx.restsize == 0);

    g_statscache >> "steps_para_sample_sort" << ctx.para_ss_steps
                 >> "steps_seq_sample_sort" << ctx.seq_ss_steps
                 >> "steps_base_sort" << ctx.bs_steps;
}

void parallel_sample_sortBTC(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifySimple>(strings, n, 0);
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTC,
                             "bingmann/parallel_sample_sortBTC",
                             "bingmann/parallel_sample_sortBTC: binary tree, bktcache")

void parallel_sample_sortBTCU1(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyUnrollTree>(strings, n, 0);
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTCU1,
                             "bingmann/parallel_sample_sortBTCU1",
                             "bingmann/parallel_sample_sortBTCU1: binary tree, bktcache, unroll tree")

void parallel_sample_sortBTCU2(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyUnrollBoth>(strings, n, 0);
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTCU2,
                             "bingmann/parallel_sample_sortBTCU2",
                             "bingmann/parallel_sample_sortBTCU2: binary tree, bktcache, unroll tree and strings")

////////////////////////////////////////////////////////////////////////////////

// ****************************************************************************
// *** Classification with Binary Search in Splitter Array (old BSC variant)

template <size_t treebits>
struct ClassifyBinarySearch
{
    // NOTE: for binary search numsplitters need not be 2^k-1, any size will
    // do, but the tree implementations are always faster, so we keep this only
    // for historical reasons.
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter[numsplitters];

    /// binary search on splitter array for bucket number
    inline unsigned int
    find_bkt_tree(const key_type& key) const
    {
        unsigned int lo = 0, hi = numsplitters;

        while ( lo < hi )
        {
            unsigned int mid = (lo + hi) >> 1;
            assert(mid < numsplitters);

            if (key <= splitter[mid])
                hi = mid;
            else // (key > splitter[mid])
                lo = mid + 1;
        }

        size_t b = lo * 2;                                    // < bucket
        if (lo < numsplitters && splitter[lo] == key) b += 1; // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key);
            *bktout++ = b;
        }
    }

    //! return a splitter
    inline key_type get_splitter(unsigned int i) const
    { return splitter[i]; }

    /// build tree and splitter array from sample
    inline void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        const size_t oversample_factor = samplesize / numsplitters;
        DBG(debug_splitter, "oversample_factor: " << oversample_factor);

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
    }
};

void parallel_sample_sortBSC(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyBinarySearch>(strings, n, 0);
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBSC,
                             "bingmann/parallel_sample_sortBSC",
                             "bingmann/parallel_sample_sortBSC: binary search, bktcache")

////////////////////////////////////////////////////////////////////////////////

// ****************************************************************************
// *** Classification with Equality Checking at Each Node

template <size_t numsplitters>
struct TreeBuilderEqual
{
    key_type*       m_tree;
    unsigned char*  m_lcp_iter;
    key_type*       m_samples;

    TreeBuilderEqual(key_type* splitter_tree, unsigned char* splitter_lcp,
                     key_type* samples, size_t samplesize)
        : m_tree( splitter_tree ),
          m_lcp_iter( splitter_lcp ),
          m_samples( samples )
    {

        key_type sentinel = 0;
        recurse(samples, samples + samplesize, 1, sentinel);

        assert(m_lcp_iter == splitter_lcp + numsplitters);
        splitter_lcp[0] &= 0x80; // overwrite sentinel lcp for first < everything bucket
        splitter_lcp[numsplitters] = 0; // sentinel for > everything bucket
    }

    ptrdiff_t snum(key_type* s) const
    {
        return (ptrdiff_t)(s - m_samples);
    }

    key_type recurse(key_type* lo, key_type* hi, unsigned int treeidx,
                     key_type& rec_prevkey)
    {
        DBG(debug_splitter, "rec_buildtree(" << snum(lo) << "," << snum(hi)
            << ", treeidx=" << treeidx << ")");

        // pick middle element as splitter
        key_type* mid = lo + (ptrdiff_t)(hi - lo) / 2;

        DBG(debug_splitter, "tree[" << treeidx << "] = samples[" << snum(mid) << "] = "
            << toHex(*mid));

        key_type mykey = m_tree[treeidx] = *mid;
#if 1
        key_type* midlo = mid;
        while (lo < midlo && *(midlo-1) == mykey) midlo--;

        key_type* midhi = mid;
        while (midhi+1 < hi && *midhi == mykey) midhi++;

        if (midhi - midlo > 1)
            DBG(0, "key range = [" << snum(midlo) << "," << snum(midhi) << ")");
#else
        key_type *midlo = mid, *midhi = mid+1;
#endif
        if (2 * treeidx < numsplitters)
        {
            key_type prevkey = recurse(lo, midlo, 2 * treeidx + 0, rec_prevkey);

            key_type xorSplit = prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " << toHex(prevkey) << " XOR " << toHex(mykey) << " = "
                << toHex(xorSplit) << " - " << count_high_zero_bits(xorSplit) << " bits = "
                << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                | ((mykey & 0xFF) ? 0 : 0x80); // marker for done splitters

            return recurse(midhi, hi, 2 * treeidx + 1, mykey);
        }
        else
        {
            key_type xorSplit = rec_prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " << toHex(rec_prevkey) << " XOR " << toHex(mykey) << " = "
                << toHex(xorSplit) << " - " << count_high_zero_bits(xorSplit) << " bits = "
                << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                | ((mykey & 0xFF) ? 0 : 0x80); // marker for done splitters

            return mykey;
        }
    }
};

// ****************************************************************************
// *** Classification Subroutines: rolled, and two unrolled variants

template <size_t treebits>
struct ClassifyEqual
{
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters+1];

    /// binary search on splitter array for bucket number
    inline unsigned int
    find_bkt_tree(const key_type& key) const
    {
        unsigned int i;

#if 1
        // hand-coded assembler binary tree traversal with equality, using CMOV
        asm("mov    $1, %%rax \n"             // rax = i
            // body of while loop
            "1: \n"
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            "cmp    %[numsplitters1], %%rax \n"  // i < numsplitters+1
            "jb     1b \n"
            "sub    %[numsplitters1], %%rax \n"  // i -= numsplitters+1;
            "lea    (%%rax,%%rax), %%rax \n"     // i = i*2
            "jmp    3f \n"
            "2: \n"
            "bsr    %%rax, %%rdx \n"             // dx = bit number of highest one
            "mov    %[treebits], %%rcx \n"
            "sub    %%rdx, %%rcx \n"             // cx = treebits - highest
            "shl    %%cl, %%rax \n"              // shift ax to left
            "and    %[numsplitters], %%rax \n"   // mask off other bits
            "lea    -1(%%rcx), %%rcx \n"
            "mov    $1, %%rdx \n"                // dx = (1 << (hi-1))
            "shl    %%cl, %%rdx \n"              //
            "or     %%rdx, %%rax \n"             // ax = OR of both
            "lea    -1(%%rax,%%rax), %%rax \n"    // i = i * 2 - 1
            "3: \n"
            : "=&a" (i)
            : [key] "r" (key), [splitter_tree] "r" (splitter_tree),
              [numsplitters1] "g" (numsplitters+1),
              [treebits] "g" (treebits),
              [numsplitters] "g" (numsplitters)
            : "rcx", "rdx");
#else
        // hand-coded assembler binary tree traversal with equality, using SETA
        // this version is slightly slower than the CMOV above.
        asm("mov    $1, %%rax \n"             // rax = i
            // body of while loop
            "1: \n"
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "seta   %%cl \n"                     // cl = 1 if key > splitter_tree
            "movzb  %%cl, %%rcx \n"              // pad cl with zeros
            "je     2f \n"
            "lea    (%%rcx,%%rax,2), %%rax \n"   // rax += 2 * i + 0/1
            "cmp    %[numsplitters1], %%rax \n"  // i < numsplitters+1
            "jb     1b \n"
            "sub    %[numsplitters1], %%rax \n"  // i -= numsplitters+1;
            "lea    (%%rax,%%rax), %%rax \n"     // i = i*2
            "jmp    3f \n"
            "2: \n"
            "bsr    %%rax, %%rdx \n"             // dx = bit number of highest one
            "mov    %[treebits], %%rcx \n"
            "sub    %%rdx, %%rcx \n"             // cx = treebits - highest
            "shl    %%cl, %%rax \n"              // shift ax to left
            "and    %[numsplitters], %%rax \n"   // mask off other bits
            "lea    -1(%%rcx), %%rcx \n"
            "mov    $1, %%rdx \n"                // dx = (1 << (hi-1))
            "shl    %%cl, %%rdx \n"              //
            "or     %%rdx, %%rax \n"             // ax = OR of both
            "lea    -1(%%rax,%%rax), %%rax \n"    // i = i * 2 - 1
            "3: \n"
            : "=&a" (i)
            : [key] "r" (key), [splitter_tree] "r" (splitter_tree),
              [numsplitters1] "g" (numsplitters+1),
              [treebits] "g" (treebits),
              [numsplitters] "g" (numsplitters)
            : "rcx", "rdx");
#endif
        return i;
    }

    /// classify all strings in area by walking tree and saving bucket id
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key);
            *bktout++ = b;
        }
    }

    //! return a splitter
    inline key_type get_splitter(unsigned int i) const
    {
        return splitter_tree[ TreeCalculations<treebits>::in_to_levelorder(i) ];
    }

    /// build tree and splitter array from sample
    inline void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderEqual<numsplitters>(splitter_tree, splitter_lcp,
                                       samples, samplesize);
    }
};

template <size_t TreeBits>
struct ClassifyEqualUnrollTree
{
    static const size_t treebits = TreeBits;
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters+1];

    /// specialized implementation of this find_bkt_tree are below
    inline unsigned int
    find_bkt_tree(const key_type& key) const;

    /// classify all strings in area by walking tree and saving bucket id
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key);
            *bktout++ = b;
        }
    }

    //! return a splitter
    inline key_type get_splitter(unsigned int i) const
    {
        return splitter_tree[ TreeCalculations<treebits>::in_to_levelorder(i) ];
    }

    /// build tree and splitter array from sample
    inline void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderEqual<numsplitters>(splitter_tree, splitter_lcp,
                                       samples, samplesize);
    }

#define SPLITTER_TREE_STEP                                              \
    /* inside the tree */                                               \
    "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"                      \
    "je     2f \n"                                                      \
    "lea    (%%rax,%%rax), %%rax \n"                                    \
    "lea    1(%%rax), %%rcx \n"                                         \
    "cmova  %%rcx, %%rax \n"             /* CMOV rax = 2 * i + 1 */

#define SPLITTER_TREE_END                                               \
    /* at leaf level */                                                 \
    "sub    %[numsplitters1], %%rax \n"  /* i -= numsplitters+1; */     \
    "lea    (%%rax,%%rax), %%rax \n"     /* i = i*2 */                  \
    "jmp    3f \n"                                                      \
    "2: \n"                                                             \
    "bsr    %%rax, %%rdx \n"             /* dx = bit number of highest one */ \
    "mov    %[treebits], %%rcx \n"                                      \
    "sub    %%rdx, %%rcx \n"             /* cx = treebits - highest */  \
    "shl    %%cl, %%rax \n"              /* shift ax to left */         \
    "and    %[numsplitters], %%rax \n"   /* mask off other bits */      \
    "lea    -1(%%rcx), %%rcx \n"                                        \
    "mov    $1, %%rdx \n"                /* dx = (1 << (hi-1)) */       \
    "shl    %%cl, %%rdx \n"                                             \
    "or     %%rdx, %%rax \n"             /* ax = OR of both */          \
    "lea    -1(%%rax,%%rax), %%rax \n"   /* i = i * 2 - 1 */            \
            "3: \n"
};

/// binary search on splitter array for bucket number
template <> inline unsigned int
ClassifyEqualUnrollTree<11>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm("mov    $1, %%rax \n"             // rax = i
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_STEP

        SPLITTER_TREE_END

        : "=&a" (i)
        : [key] "r" (key), [splitter_tree] "r" (splitter_tree),
          [numsplitters1] "g" (numsplitters+1),
          [treebits] "g" (treebits),
          [numsplitters] "g" (numsplitters)
        : "rcx", "rdx");

    return i;
}

template <> inline unsigned int
ClassifyEqualUnrollTree<12>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm("mov    $1, %%rax \n"             // rax = i
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_END

        : "=&a" (i)
        : [key] "r" (key), [splitter_tree] "r" (splitter_tree),
          [numsplitters1] "g" (numsplitters+1),
          [treebits] "g" (treebits),
          [numsplitters] "g" (numsplitters)
        : "rcx", "rdx");

    return i;
}

template <> inline unsigned int
ClassifyEqualUnrollTree<13>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm("mov    $1, %%rax \n"             // rax = i
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_END

        : "=&a" (i)
        : [key] "r" (key), [splitter_tree] "r" (splitter_tree),
          [numsplitters1] "g" (numsplitters+1),
          [treebits] "g" (treebits),
          [numsplitters] "g" (numsplitters)
        : "rcx", "rdx");

    return i;
}

template <> inline unsigned int
ClassifyEqualUnrollTree<14>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm("mov    $1, %%rax \n"             // rax = i
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP
        SPLITTER_TREE_STEP

        SPLITTER_TREE_END

        : "=&a" (i)
        : [key] "r" (key), [splitter_tree] "r" (splitter_tree),
          [numsplitters1] "g" (numsplitters+1),
          [treebits] "g" (treebits),
          [numsplitters] "g" (numsplitters)
        : "rcx", "rdx");

    return i;
}

#undef SPLITTER_TREE_STEP
#undef SPLITTER_TREE_END

void parallel_sample_sortBTCE(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyEqual>(strings, n, 0);

}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTCE,
		"bingmann/parallel_sample_sortBTCE",
		"bingmann/parallel_sample_sortBTCE: binary tree, equality, bktcache")

void parallel_sample_sortBTCEU1(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyEqualUnrollTree>(strings, n, 0);
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTCEU1,
        "bingmann/parallel_sample_sortBTCEU1",
        "bingmann/parallel_sample_sortBTCEU1: binary tree, equality, bktcache, unroll tree")

//! Call for NUMA aware parallel sorting
void parallel_sample_sort_numa(string * strings, size_t n,
                               int numaNode, int numberOfThreads)
{
    Context ctx;
    ctx.totalsize = ctx.restsize = n;
    ctx.threadnum = numberOfThreads;
    ctx.para_ss_steps = ctx.seq_ss_steps = ctx.bs_steps = 0;

    SampleSortStep<ClassifyUnrollBoth>::put_stats();

    string* shadow = new string[n]; // allocate shadow pointer array

    Enqueue<ClassifyUnrollBoth>(ctx, StringPtr(strings, shadow), n, 0);
    ctx.jobqueue.numaLoop(numaNode, numberOfThreads, ctx);

    delete [] shadow;

    assert(ctx.restsize == 0);

    g_statscache >> "steps_para_sample_sort" << ctx.para_ss_steps
                 >> "steps_seq_sample_sort" << ctx.seq_ss_steps
                 >> "steps_base_sort" << ctx.bs_steps;
}

} // namespace bingmann_parallel_sample_sort
