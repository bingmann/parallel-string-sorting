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
#undef DBGX
#define DBGX DBGX_OMP

#include "../tools/lcgrandom.h"
#include "../tools/contest.h"
#include "../tools/stringtools.h"
#include "../tools/jobqueue.h"
#include "../tools/logfloor.h"

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
static const bool debug_lcp = false;

static const bool use_work_sharing = false;

static const size_t MAXPROCS = 64+1; // +1 due to round up of processor number

static const size_t l2cache = 256*1024;

static const size_t g_smallsort_threshold = 2*1024; // TODO
static const size_t g_inssort_threshold = 64;

typedef uint64_t key_type;

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
// *** SortStep to Keep Track of Substeps

class SortStep
{
private:

    //! Number of substeps still running
    std::atomic<unsigned int> m_substep_working;

    //! Pure virtual function called by substep when all substeps are done.
    virtual void substep_all_done() = 0;

protected:
    SortStep()
        : m_substep_working(0)
    { }

    virtual ~SortStep()
    {
        assert(m_substep_working == 0);
    }

    //! Register new substep
    void substep_add()
    {
        ++m_substep_working;
    }

public:
    //! Notify superstep that the currently substep is done.
    void substep_notify_done()
    {
        assert(m_substep_working > 0);
        if (--m_substep_working == 0)
            substep_all_done();
    }
};

// ****************************************************************************
// *** Classification Variants

unsigned char lcpKeyType(const key_type& a, const key_type& b)
{
    // XOR both values and count the number of zero bytes
    return count_high_zero_bits(a ^ b) / 8;
}

unsigned char lcpKeyDepth(const key_type& a)
{
    // count number of non-zero bytes
    return sizeof(key_type) - (count_low_zero_bits(a) / 8);
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

void insertion_sort(const StringPtr& strptr, size_t n, size_t depth)
{
    assert(!strptr.flipped());
    //inssort::inssort(strptr.active(), n, depth);
    bingmann_lcp_inssort::lcp_insertion_sort(strptr, n, depth);
}

// ****************************************************************************
// *** LCP Calculation for finished Sample Sort Steps

template <size_t bktnum, typename Classify, typename BktSizeType>
void sample_sort_lcp(const Classify& classifier, const StringPtr& strptr, size_t depth, const BktSizeType* bkt)
{
    size_t b = 0; // current bucket number
    key_type prevkey; // previous key

    // the following while loops only check b < bktnum when b is odd,
    // because bktnum is always odd. We need a goto to jump into the loop,
    // as b == 0 start even.
    goto even_first;

    // find first non-empty bucket
    while (b < bktnum)
    {
        // odd bucket: = bkt
        if (bkt[b] != bkt[b+1])
        {
            prevkey = classifier.get_splitter(b/2);
            assert( prevkey == get_char<key_type>(strptr.str(bkt[b+1]-1), depth) );
            break;
        }
        ++b;
    even_first:
        // even bucket: <, << or > bkt
        if (bkt[b] != bkt[b+1])
        {
            prevkey = get_char<key_type>( strptr.str(bkt[b+1]-1), depth );
            break;
        }
        ++b;
    }
    ++b;

    // goto depends on whether the first non-empty bucket was odd or
    // even. the while loop below encodes this in the program counter.
    if (b < bktnum && b % 2 == 0) goto even_bucket;

    // find next non-empty bucket
    while (b < bktnum)
    {
        // odd bucket: = bkt
        if (bkt[b] != bkt[b+1])
        {
            key_type thiskey = classifier.get_splitter(b/2);
            if (thiskey != get_char<key_type>(strptr.str(bkt[b]), depth))
                std::cout << "mismath: " << thiskey << " != " << get_char<key_type>(strptr.str(bkt[b]), depth) << "\n";

            if (thiskey != get_char<key_type>(strptr.str(bkt[b]), depth))
                DBG(1, toHex(thiskey) << " != " << toHex(get_char<key_type>(strptr.str(bkt[b]), depth)));

            assert( thiskey == get_char<key_type>(strptr.str(bkt[b]), depth) );

            strptr.set_lcp(bkt[b], depth + lcpKeyType(prevkey, thiskey));
            DBG(debug_lcp, "LCP at odd-bucket " << b << " [" << bkt[b] << "," << bkt[b+1]
                << ") is " << strptr.lcp(bkt[b]));

            prevkey = thiskey;

            if (prevkey != get_char<key_type>(strptr.str(bkt[b+1]-1), depth))
                DBG(1, toHex(prevkey) << " != " << toHex(get_char<key_type>(strptr.str(bkt[b+1]-1), depth)));
            assert( prevkey == get_char<key_type>(strptr.str(bkt[b+1]-1), depth) );
        }
        ++b;
    even_bucket:
        // even bucket: <, << or > bkt
        if (bkt[b] != bkt[b+1])
        {
            key_type thiskey = get_char<key_type>( strptr.str(bkt[b]), depth );

            strptr.set_lcp(bkt[b], depth + lcpKeyType(prevkey, thiskey));
            DBG(debug_lcp, "LCP at even-bucket " << b << " [" << bkt[b] << "," << bkt[b+1]
                << ") is " << strptr.lcp(bkt[b]));

            prevkey = get_char<key_type>( strptr.str(bkt[b+1]-1), depth );
        }
        ++b;
    }

}

// ****************************************************************************
// *** SampleSort non-recursive in-place sequential sample sort for small sorts

template <template <size_t> class Classify>
void Enqueue(Context& ctx, SortStep* sstep, const StringPtr& strptr, size_t n, size_t depth);

template <template <size_t> class Classify, typename BktSizeType>
struct SmallsortJob : public Job, public SortStep
{
    //! parent sort step
    SortStep*   pstep;

    StringPtr   in_strptr;
    size_t      in_n, in_depth;

    typedef BktSizeType bktsize_type;

    SmallsortJob(Context& ctx, SortStep* _pstep,
                 const StringPtr& strptr, size_t n, size_t depth)
        : pstep(_pstep), in_strptr(strptr), in_n(n), in_depth(depth)
    {
        DBG(debug_jobs, "enqueue SmallsortJob n=" << in_n << " depth=" << in_depth);

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

        Classify<treebits> classifier;
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

            classifier.build(samples, samplesize, splitter_lcp);

            // step 2: classify all strings

            classifier.classify(strings, strings + n, bktcache, depth);

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
            bkt[bktnum] = n;

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

            // fix prefix sum??? TODO

            bkt[0] = 0;
            for (unsigned int i=1; i <= bktnum; ++i) {
                bkt[i] = bkt[i-1] + bktsize[i-1];
            }
            assert(bkt[bktnum] == n);

            ++ctx.seq_ss_steps;
        }

        void calculate_lcp()
        {
            DBG(debug_lcp, "Calculate LCP after sample sort step");
            sample_sort_lcp<bktnum>(classifier, strptr.front(), depth, bkt);
        }
    };

    // *** Stack of Recursive Sample Sort Steps

    uint16_t* bktcache;
    size_t bktcache_size;

    size_t ss_pop_front;
    std::vector<SeqSampleSortStep> ss_stack;

    virtual bool run(Context& ctx)
    {
        DBG(debug_jobs, "Process SmallsortJob " << this << " of size " << in_n);

        // create anonymous wrapper job
        substep_add();

        bktcache = NULL;
        bktcache_size = 0;
        ss_pop_front = 0;
        ms_pop_front = 0;

        if (in_n >= g_smallsort_threshold)
        {
            bktcache = new uint16_t[in_n];
            bktcache_size = in_n * sizeof(uint16_t);
            sort_sample_sort(ctx, in_strptr, in_n, in_depth);
        }
        else
        {
            sort_mkqs_cache(ctx, in_strptr, in_n, in_depth);
        }

        delete [] bktcache;

        // finish wrapper job, handler delete's this
        substep_notify_done();

        return false;
    }

    void sort_sample_sort(Context& ctx, const StringPtr& strptr, size_t n, size_t depth)
    {
        typedef SeqSampleSortStep Step;

        assert(ss_pop_front == 0);
        assert(ss_stack.size() == 0);

        // sort first level
        ss_stack.push_back( Step(ctx,strptr,n,depth,bktcache) );

        // step 5: "recursion"

    jumpout:
        while ( ss_stack.size() > ss_pop_front )
        {
            while ( ss_stack.back().idx < Step::bktnum )
            {
                Step& s = ss_stack.back();
                size_t i = s.idx++; // process the bucket s.idx

                size_t bktsize = s.bkt[i+1] - s.bkt[i];

                // i is even -> bkt[i] is less-than bucket
                if (i % 2 == 0)
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
                        StringPtr sp = (s.strptr + s.bkt[i]).copy_back(bktsize);
                        sp.fill_lcp(bktsize, s.depth + lcpKeyDepth(s.classifier.get_splitter(i/2)));
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

                if (use_work_sharing && ctx.jobqueue.has_idle()) { // TODO
                    sample_sort_free_work(ctx);
                    goto jumpout;
                }
            }
            // after full sort: calculate LCPs at this level
            ss_stack.back().calculate_lcp();

            ss_stack.pop_back();
        }
    }

    void sample_sort_free_work(Context& ctx)
    {
        assert(ss_stack.size() >= ss_pop_front);

        if (ss_stack.size() == ss_pop_front) {
            // ss_stack is empty, check other stack
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
            if (i % 2 == 0)
            {
                if (bktsize == 0)
                    ;
                else
                {
                    if (i == Step::bktnum-1)
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << bktsize << " no lcp");
                    else
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << bktsize << " lcp " << int(s.splitter_lcp[i/2] & 0x7F));

                    substep_add();
                    Enqueue<Classify>( ctx, this, s.strptr + s.bkt[i], bktsize, s.depth + (s.splitter_lcp[i/2] & 0x7F) );
                }
            }
            // i is odd -> bkt[i] is equal bucket
            else
            {
                if (bktsize == 0)
                    ;
                else if ( s.splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key, done.
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " is done!");
                    StringPtr sp = (s.strptr + s.bkt[i]).copy_back(bktsize);
                    sp.fill_lcp(bktsize, s.depth + lcpKeyDepth(s.classifier.get_splitter(i/2)));
                    ctx.restsize -= bktsize;
                }
                else
                {
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " lcp keydepth!");

                    substep_add();
                    Enqueue<Classify>( ctx, this, s.strptr + s.bkt[i], bktsize, s.depth + sizeof(key_type) );
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
                if (cache[pj-1] <= tmpc)
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
    insertion_sort_cache(const StringPtr& _strptr, key_type* cache, int n, size_t depth)
    {
        StringPtr strptr = _strptr.copy_back(n);

        if (n <= 1) return;
        if (CacheDirty) return insertion_sort(strptr, n, depth);

        insertion_sort_cache_block(strptr, cache, n);

        size_t start = 0, bktsize = 1;
        for (int i = 0; i < n - 1; ++i) {
            // group areas with equal cache values
            if (cache[i] == cache[i+1]) {
                ++bktsize;
                continue;
            }
            // calculate LCP between group areas
            if (start != 0) {
                strptr.set_lcp(start, depth + lcpKeyType(cache[start-1], cache[start]));
            }
            // sort group areas deeper if needed
            if (bktsize > 1) {
                if (cache[start] & 0xFF) // need deeper sort
                    insertion_sort(strptr + start, bktsize, depth + sizeof(key_type));
                else // cache contains NULL-termination
                    (strptr + start).fill_lcp(bktsize, depth + lcpKeyDepth(cache[start]));
            }
            bktsize = 1;
            start = i+1;
        }
        // tail of loop for last item
        if (start != 0) {
            strptr.set_lcp(start, depth + lcpKeyType(cache[start-1], cache[start]));
        }
        if (bktsize > 1) {
            if (cache[start] & 0xFF) // need deeper sort
                insertion_sort(strptr + start, bktsize, depth + sizeof(key_type));
            else // cache contains NULL-termination
                (strptr + start).fill_lcp(bktsize, depth + lcpKeyDepth(cache[start]));
        }
    }

    struct MKQSStep
    {
        StringPtr strptr;
        key_type* cache;
        size_t num_lt, num_eq, num_gt, depth;
        size_t idx;
        bool eq_recurse;
        unsigned char lcp_lt, lcp_eq, lcp_gt;

        MKQSStep(Context& ctx, const StringPtr& _strptr, key_type* _cache, size_t n, size_t _depth, bool CacheDirty)
            : strptr(_strptr), cache(_cache), depth(_depth), idx(0)
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

            // save LCP values for writing into LCP array after sorting further
            if (num_lt > 0)
            {
                assert(max_lt == *std::max_element(cache + 0, cache + num_lt));

                lcp_lt = lcpKeyType(max_lt, pivot);
                DBG(debug_lcp, "LCP lt with pivot: " << depth + lcp_lt);
            }

            // calculate equal area lcp: +1 for the equal zero termination byte
            lcp_eq = lcpKeyDepth(pivot);

            if (num_gt > 0)
            {
                assert(min_gt == *std::min_element(cache + num_lt + num_eq, cache + n));

                lcp_gt = lcpKeyType(pivot, min_gt);
                DBG(debug_lcp, "LCP pivot with gt: " << depth + lcp_gt);
             }

            ++ctx.bs_steps;
        }

        void calculate_lcp()
        {
#if 1
            if (num_lt > 0)
                strptr.front().set_lcp(num_lt, depth + lcp_lt);

            if (num_gt > 0)
                strptr.front().set_lcp(num_lt + num_eq, depth + lcp_gt);
#else
            key_type pivot = get_char<key_type>(strptr.str(num_lt), depth);

            if (num_lt > 0)
            {
                key_type max_lt = get_char<key_type>(strptr.str(num_lt-1), depth);

                unsigned int lcp = depth + lcpKeyType(max_lt, pivot);
                DBG(debug_lcp, "LCP lt with pivot: " << lcp);

                assert(lcp == depth + lcp_lt);

                strptr.set_lcp(num_lt, lcp);
            }
            if (num_gt > 0)
            {
                key_type min_gt = get_char<key_type>(strptr.str(num_lt + num_eq), depth);

                unsigned int lcp = depth + lcpKeyType(pivot, min_gt);
                DBG(debug_lcp, "LCP pivot with gt: " << lcp);

                assert(lcp == depth + lcp_gt);

                strptr.set_lcp(num_lt + num_eq, lcp);
            }
#endif
        }
    };

    size_t ms_pop_front;
    std::vector<MKQSStep> ms_stack;

    void sort_mkqs_cache(Context& ctx, const StringPtr& strptr, size_t n, size_t depth)
    {
        if (n < g_inssort_threshold) { // TODO
            insertion_sort(strptr.copy_back(n), n, depth);
            ctx.restsize -= n;
            return;
        }

        if (bktcache_size < n * sizeof(key_type)) {
            delete [] bktcache;
            bktcache = (uint16_t*)new key_type[n];
            bktcache_size = n * sizeof(key_type);
        }

        key_type* cache = (key_type*)bktcache;

        assert(ms_pop_front == 0);
        assert(ms_stack.size() == 0);

        // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
        ms_stack.push_back( MKQSStep(ctx, strptr, cache, n, depth, true) );

    jumpout:
        while ( ms_stack.size() > ms_pop_front )
        {
            while ( ms_stack.back().idx < 3 )
            {
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

                    assert(ms.num_eq > 0);

                    if (!ms.eq_recurse) {
                        sp = sp.copy_back(ms.num_eq);
                        sp.fill_lcp(ms.num_eq, ms.depth + ms.lcp_eq);
                        ctx.restsize -= ms.num_eq;
                    }
                    else if (ms.num_eq < g_inssort_threshold) {
                        insertion_sort_cache<true>(sp, ms.cache + ms.num_lt,
                                                   ms.num_eq, ms.depth + sizeof(key_type));
                        ctx.restsize -= ms.num_eq;
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
                    assert(ms.idx == 3);

                    StringPtr sp = ms.strptr + ms.num_lt + ms.num_eq;

                    if (ms.num_gt == 0)
                        ;
                    else if (ms.num_gt < g_inssort_threshold) {
                        insertion_sort_cache<false>(sp, ms.cache + ms.num_lt + ms.num_eq,
                                                    ms.num_gt, ms.depth);
                        ctx.restsize -= ms.num_gt;
                    }
                    else {
                        ms_stack.push_back(
                            MKQSStep(ctx, sp,
                                     ms.cache + ms.num_lt + ms.num_eq,
                                     ms.num_gt, ms.depth, false) );
                    }
                }

                if (use_work_sharing && ctx.jobqueue.has_idle()) { // TODO
                    sample_sort_free_work(ctx);
                    goto jumpout;
                }
            }
            // calculate LCP after the three parts are sorted
            ms_stack.back().calculate_lcp();

            ms_stack.pop_back();
        }
    }

    void mkqs_free_work(Context& ctx)
    {
        assert(ms_stack.size() >= ms_pop_front);

        if (ms_stack.size() == ms_pop_front) {
            return;
        }

        DBG(debug_jobs, "Freeing top level of SmallsortJob's mkqs stack");

        // convert top level of stack into independent jobs

        MKQSStep& ms = ms_stack[ms_pop_front];

        if (ms.idx == 0 && ms.num_lt != 0)
        {
            substep_add();
            Enqueue<Classify>(ctx, this, ms.strptr, ms.num_lt, ms.depth);
        }
        if (ms.idx <= 1) // st.num_eq > 0 always
        {
            assert(ms.num_eq > 0);

            StringPtr sp = ms.strptr + ms.num_lt;

            if (ms.eq_recurse) {
                substep_add();
                Enqueue<Classify>(ctx, this, sp,
                                  ms.num_eq, ms.depth + sizeof(key_type));
            }
            else {
                sp = sp.copy_back(ms.num_eq);
                sp.fill_lcp(ms.num_eq, ms.depth + ms.lcp_eq);
                ctx.restsize -= ms.num_eq;
            }
        }
        if (ms.idx <= 2 && ms.num_gt != 0)
        {
            substep_add();
            Enqueue<Classify>(ctx, this, ms.strptr + ms.num_lt + ms.num_eq,
                              ms.num_gt, ms.depth);
        }

        // shorten the current stack
        ++ms_pop_front;
    }

    virtual void substep_all_done()
    {
        DBG(1 || debug_recursion, "SmallSort[" << in_depth << "] all substeps done -> LCP calculation");

        while (ms_pop_front > 0) {
            DBG(1 || debug_lcp, "SmallSort[" << in_depth << "] ms_pop_front: " << ms_pop_front);
            ms_stack[--ms_pop_front].calculate_lcp();
        }

        while (ss_pop_front > 0) {
            DBG(1 || debug_lcp, "SmallSort[" << in_depth << "] ss_pop_front: " << ss_pop_front);
            ss_stack[--ss_pop_front].calculate_lcp();
        }

        if (pstep) pstep->substep_notify_done();

        delete this;
    }
};

// ****************************************************************************
// *** SampleSortStep out-of-place parallel sample sort with separate Jobs

template <template <size_t> class Classify>
struct SampleSortStep : public SortStep
{
#if 0
    static const size_t numsplitters2 = 16;
#else
    // bounding equation: 2*K+1 buckets when counting bkt occurances.
    static const size_t numsplitters2 = (l2cache - sizeof(size_t)) / (2 * sizeof(size_t));
#endif
    static const size_t treebits = logfloor_<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;

    static const size_t bktnum = 2 * numsplitters + 1;

    //! parent sort step notification
    SortStep*           pstep;

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

        virtual bool run(Context& ctx)
        {
            step->sample(ctx);
            return true;
        }
    };

    struct CountJob : public Job
    {
        SampleSortStep*     step;
        unsigned int        p;

        CountJob(SampleSortStep* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual bool run(Context& ctx)
        {
            step->count(p, ctx);
            return true;
        }
    };

    struct DistributeJob : public Job
    {
        SampleSortStep*     step;
        unsigned int        p;

        DistributeJob(SampleSortStep* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual bool run(Context& ctx)
        {
            step->distribute(p, ctx);
            return true;
        }
    };

    // *** Constructor

    SampleSortStep(Context& ctx, SortStep* _pstep,
                   const StringPtr& _strptr, size_t _n, size_t _depth)
        : pstep(_pstep), strptr(_strptr), n(_n), depth(_depth)
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
        DBG(debug_jobs, "Process SampleJob @ " << this);

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
        DBG(debug_jobs, "Process CountJob " << p << " @ " << this);

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
        DBG(debug_jobs, "Finishing CountJob " << this << " with prefixsum");

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
        DBG(debug_jobs, "Process DistributeJob " << p << " @ " << this);

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
        DBG(debug_jobs, "Finishing DistributeJob " << this << " with enqueuing subjobs");

        size_t* bkt = this->bkt[0];
        assert(bkt);

        // first processor's bkt pointers are boundaries between bkts, just add sentinel:
        assert(bkt[0] == 0);
        bkt[bktnum] = n;

        // keep anonymous subjob handle while creating subjobs
        substep_add();

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
                substep_add();
                Enqueue<Classify>(ctx, this, strptr.flip(bkt[i]), bktsize, depth + (splitter_lcp[i/2] & 0x7F));
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
                    StringPtr sp = strptr.flip(bkt[i]).copy_back(bktsize);
                    sp.fill_lcp(bktsize, depth + lcpKeyDepth(classifier.get_splitter(i/2)));
                    ctx.restsize -= bktsize;
                }
                else {
                    DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " lcp keydepth!");
                    substep_add();
                    Enqueue<Classify>(ctx, this, strptr.flip(bkt[i]), bktsize, depth + sizeof(key_type));
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
            substep_add();
            Enqueue<Classify>(ctx, this, strptr.flip(bkt[i]), bktsize, depth);
        }

        substep_notify_done(); // release anonymous subjob handle
    }

    // *** After Recursive Sorting

    virtual void substep_all_done()
    {
        DBG(1 || debug_recursion, "pSampleSortStep[" << depth << "]: all substeps done.");

        sample_sort_lcp<bktnum>(classifier, strptr.front(), depth, bkt[0]);
        delete [] bkt[0];

        if (pstep) pstep->substep_notify_done();
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
void Enqueue(Context& ctx, SortStep* pstep,
             const StringPtr& strptr, size_t n, size_t depth)
{
    if (n > ctx.sequential_threshold()) { // TODO
        new SampleSortStep<Classify>(ctx, pstep, strptr, n, depth);
    }
    else {
        if (n < ((uint64_t)1 << 32))
            new SmallsortJob<Classify,uint32_t>(ctx, pstep, strptr, n, depth);
        else
            new SmallsortJob<Classify,uint64_t>(ctx, pstep, strptr, n, depth);
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

    StringPtr strptr(strings, shadow);

    Enqueue<Classify>(ctx, NULL, strptr, n, depth);
    ctx.jobqueue.loop(ctx);

    if (1)
    {
        unsigned int err = 0;
        for (size_t i = 1; i < n; ++i)
        {
            string s1 = strptr.str(i-1), s2 = strptr.str(i);
            size_t h = 0;
            while (*s1 != 0 && *s1 == *s2)
                ++h, ++s1, ++s2;

            if (h != strptr.lcp(i)) {
                std::cout << "lcp[" << i << "] mismatch " << h << " != " << strptr.lcp(i) << std::endl;
                ++err;
            }
        }
        std::cout << err << " lcp errors" << std::endl;
    }

    assert(ctx.restsize == 0);

    delete [] shadow;

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

template <> inline unsigned int
ClassifyEqualUnrollTree<4>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm("mov    $1, %%rax \n"             // rax = i
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
ClassifyEqualUnrollTree<5>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm("mov    $1, %%rax \n"             // rax = i
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

    StringPtr strptr(strings, shadow);

    Enqueue<ClassifyUnrollBoth>(ctx, NULL, strptr, n, 0);
    ctx.jobqueue.numaLoop(numaNode, numberOfThreads, ctx);

    delete [] shadow;

    assert(ctx.restsize == 0);

    g_statscache >> "steps_para_sample_sort" << ctx.para_ss_steps
                 >> "steps_seq_sample_sort" << ctx.seq_ss_steps
                 >> "steps_base_sort" << ctx.bs_steps;
}

} // namespace bingmann_parallel_sample_sort
