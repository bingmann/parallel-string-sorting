/*******************************************************************************
 * src/parallel/bingmann-parallel_sample_sort.hpp
 *
 * Parallel Super Scalar String Sample-Sort, many variant via different
 * Classifier templates.
 *
 *******************************************************************************
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
 ******************************************************************************/

#ifndef PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_HEADER
#define PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_HEADER

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>
#include <vector>
#include <algorithm>

#include "../tools/debug.hpp"
#undef DBGX
#define DBGX DBGX_OMP

#include "../tools/lcgrandom.hpp"
#include "../tools/contest.hpp"
#include "../tools/stringtools.hpp"
#include "../tools/jobqueue.hpp"
#include "../tools/lockfree.hpp"

#include "../sequential/inssort.hpp"
#include "../sequential/bingmann-lcp_inssort.hpp"

#include <tlx/string/hexdump.hpp>
#include <tlx/meta/log2.hpp>

#ifndef CALC_LCP
#define CALC_LCP 0
#endif

// first MKQS lcp variant: keep min/max during ternary split
#if CALC_LCP
#define CALC_LCP_MKQS 1
#endif

#ifdef PSS_CONTESTANT_PARALLEL_LCP
#undef PSS_CONTESTANT_PARALLEL_LCP
#endif

#if !CALC_LCP
namespace bingmann_parallel_sample_sort {

#define PSS_CONTESTANT_PARALLEL_LCP(func, name, desc) \
    PSS_CONTESTANT_PARALLEL(func, name, desc)

//#define PSS_CONTESTANT_PARALLEL_LCP CONTESTANT_REGISTER_PARALLEL

#else
namespace bingmann_parallel_sample_sort_lcp {

#define PSS_CONTESTANT_PARALLEL_LCP(func, name, desc) \
    PSS_CONTESTANT_PARALLEL(func, name "_lcp", desc "_lcp")

#endif

using namespace stringtools;
using namespace parallel_string_sorting;
using namespace jobqueue;

static const bool debug_steps = false;
static const bool debug_jobs = false;

static const bool debug_splitter = false;
static const bool debug_bucketsize = false;
static const bool debug_recursion = false;
static const bool debug_splitter_tree = false;
static const bool debug_lcp = false;

//! enable work freeing
static const bool use_work_sharing = true;

//! whether the base sequential_threshold() on the remaining unsorted string
//! set or on the whole string set.
#define PS5_ENABLE_RESTSIZE false

//! use LCP insertion sort for non-LCP pS5 ?
static const bool use_lcp_inssort = false;

//! terminate sort after first parallel sample sort step
#ifndef PS5_SINGLE_STEP
#define PS5_SINGLE_STEP false
#endif
static const bool use_only_first_sortstep = PS5_SINGLE_STEP;

//! maximum number of threads, used in a few static arrays
static const size_t MAXPROCS = 2 * 64 + 1; // +1 due to round up of processor number

//! L2 cache size, used to calculate classifier tree sizes
#ifndef PS5_L2CACHE
#define PS5_L2CACHE     256 * 1024
#endif
static const size_t l2cache = PS5_L2CACHE;

static const size_t g_smallsort_threshold = 1024 * 1024;
static const size_t g_inssort_threshold = 64;

typedef uint64_t key_type;

//! step timer ids for different sorting steps
enum { TM_WAITING, TM_PARA_SS, TM_SEQ_SS, TM_MKQS, TM_INSSORT };

//! replace TimerArrayMT with a no-op implementation
#ifndef TIMERARRAY_REAL
typedef ::TimerArrayDummy TimerArrayMT;
typedef ::ScopedTimerKeeperDummy ScopedTimerKeeperMT;
#endif

// ****************************************************************************
// *** Global Parallel Super Scalar String Sample Sort Context

template <template <typename> class JobQueueGroupType = DefaultJobQueueGroup>
class Context
{
public:
    //! total size of input
    size_t totalsize;

#if PS5_ENABLE_RESTSIZE
    //! number of remaining strings to sort
    lockfree::lazy_counter<MAXPROCS> restsize;
#endif

    //! number of threads overall
    size_t threadnum;

    //! counters
    size_t para_ss_steps, seq_ss_steps, bs_steps;

    //! timers for individual sorting steps
    TimerArrayMT timers;

    //! type of job queue group (usually a No-Op Class)
    typedef JobQueueGroupType<Context> jobqueuegroup_type;

    //! type of job queue
    typedef typename jobqueuegroup_type::jobqueue_type jobqueue_type;

    //! typedef of compatible job type
    typedef typename jobqueue_type::job_type job_type;

    //! job queue
    jobqueue_type jobqueue;

    //! context constructor
    Context(jobqueuegroup_type* jqg = NULL)
        : para_ss_steps(0), seq_ss_steps(0), bs_steps(0),
          timers(16),
          jobqueue(*this, jqg)
    { }

    //! return sequential sorting threshold
    size_t sequential_threshold()
    {
#if PS5_ENABLE_RESTSIZE
        size_t wholesize = restsize.update().get();
#else
        size_t wholesize = totalsize;
#endif
        return std::max(g_smallsort_threshold, wholesize / threadnum);
    }

    //! decrement number of unordered strings
    void donesize(size_t n, size_t tid)
    {
#if PS5_ENABLE_RESTSIZE
        restsize.add(-n, tid);
#else
        (void)n;
        (void)tid;
#endif
    }
};

// ****************************************************************************
// *** SortStep to Keep Track of Substeps

#if CALC_LCP
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
#else
class SortStep
{
protected:
    SortStep()
    { }

    void substep_add()
    { }

public:
    void substep_notify_done()
    { }
};
#endif

// ****************************************************************************
// *** Classification Variants

static inline unsigned char
lcpKeyType(const key_type& a, const key_type& b)
{
    // XOR both values and count the number of zero bytes
    return count_high_zero_bits(a ^ b) / 8;
}

static inline unsigned char
lcpKeyDepth(const key_type& a)
{
    // count number of non-zero bytes
    return sizeof(key_type) - (count_low_zero_bits(a) / 8);
}

//! return the d-th character in the (swapped) key
static inline unsigned char
getCharAtDepth(const key_type& a, unsigned char d)
{
    return static_cast<unsigned char>(a >> (8 * (sizeof(key_type) - 1 - d)));
}

//! Recursive TreeBuilder for full-descent and unrolled variants, constructs a
//! binary tree and the plain splitter array.
template <size_t numsplitters>
class TreeBuilderFullDescent
{
public:
    key_type* m_splitter;
    key_type* m_tree;
    unsigned char* m_lcp_iter;
    key_type* m_samples;

    TreeBuilderFullDescent(
        key_type* splitter, key_type* splitter_tree,
        unsigned char* splitter_lcp,
        key_type* samples, size_t samplesize)
        : m_splitter(splitter),
          m_tree(splitter_tree),
          m_lcp_iter(splitter_lcp),
          m_samples(samples)
    {
        key_type sentinel = 0;
        recurse(samples, samples + samplesize, 1, sentinel);

        assert(m_splitter == splitter + numsplitters);
        assert(m_lcp_iter == splitter_lcp + numsplitters);
        // overwrite sentinel lcp for first < everything bucket
        splitter_lcp[0] &= 0x80;
        // sentinel for > everything bucket
        splitter_lcp[numsplitters] = 0;
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

        DBG(debug_splitter, "tree[" << treeidx << "] = samples[" << snum(mid) <<
            "] = " << tlx::hexdump_type(*mid));

        key_type mykey = m_tree[treeidx] = *mid;
#if 1
        key_type* midlo = mid;
        while (lo < midlo && *(midlo - 1) == mykey) midlo--;

        key_type* midhi = mid;
        while (midhi + 1 < hi && *midhi == mykey) midhi++;

        if (midhi - midlo > 1)
            DBG(0, "key range = [" << snum(midlo) << "," << snum(midhi) << ")");
#else
        key_type* midlo = mid, * midhi = mid + 1;
#endif
        if (2 * treeidx < numsplitters)
        {
            key_type prevkey = recurse(lo, midlo, 2 * treeidx + 0, rec_prevkey);

            key_type xorSplit = prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " <<
                tlx::hexdump_type(prevkey) << " XOR " << tlx::hexdump_type(mykey) << " = " <<
                tlx::hexdump_type(xorSplit) << " - " << count_high_zero_bits(xorSplit) <<
                " bits = " << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            *m_splitter++ = mykey;

            // marker for done splitters
            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                            | ((mykey & 0xFF) ? 0 : 0x80);

            return recurse(midhi, hi, 2 * treeidx + 1, mykey);
        }
        else
        {
            key_type xorSplit = rec_prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " <<
                tlx::hexdump_type(rec_prevkey) << " XOR " << tlx::hexdump_type(mykey) << " = " <<
                tlx::hexdump_type(xorSplit) << " - " << count_high_zero_bits(xorSplit) <<
                " bits = " << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            *m_splitter++ = mykey;

            // marker for done splitters
            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                            | ((mykey & 0xFF) ? 0 : 0x80);

            return mykey;
        }
    }
};

// ****************************************************************************
// * Classification Subroutines for TreeBuilderFullDescent: rolled, and two
// * unrolled variants

template <size_t treebits>
class ClassifySimple
{
public:
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter[numsplitters];
    key_type splitter_tree[numsplitters + 1];

    /// binary search on splitter array for bucket number
    unsigned int find_bkt_tree(const key_type& key) const
    {
        // binary tree traversal without left branch

        unsigned int i = 1;

        while (i <= numsplitters)
        {
#if 0
            // in gcc-4.6.3 this produces a SETA, LEA sequence
            i = 2 * i + (key <= splitter_tree[i] ? 0 : 1);
#else
            // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which is
            // slightly faster
            if (key <= splitter_tree[i])
                i = 2 * i + 0;
            else    // (key > splitter_tree[i])
                i = 2 * i + 1;
#endif
        }

        i -= numsplitters + 1;

        size_t b = i * 2;                                    // < bucket
        if (i < numsplitters && splitter[i] == key) b += 1;  // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    template <typename StringSet>
    void classify(
        const StringSet& strset,
        typename StringSet::Iterator begin, typename StringSet::Iterator end,
        uint16_t* bktout, size_t depth) const
    {
        while (begin != end)
        {
            key_type key = strset.get_uint64(*begin++, depth);
            *bktout++ = find_bkt_tree(key);
        }
    }

    //! return a splitter
    key_type get_splitter(unsigned int i) const
    { return splitter[i]; }

    /// build tree and splitter array from sample
    void build(key_type* samples, size_t samplesize,
               unsigned char* splitter_lcp)
    {
        TreeBuilderFullDescent<numsplitters>(
            splitter, splitter_tree, splitter_lcp,
            samples, samplesize);
    }
};

template <size_t treebits>
class ClassifyUnrollTree
{
public:
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter[numsplitters];
    key_type splitter_tree[numsplitters + 1];

    /// binary search on splitter array for bucket number
    __attribute__ ((optimize("unroll-all-loops")))
    unsigned int
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
            // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which is
            // slightly faster
            if (key <= splitter_tree[i])
                i = 2 * i + 0;
            else    // (key > splitter_tree[i])
                i = 2 * i + 1;
#endif
        }

        i -= numsplitters + 1;

        size_t b = i * 2;                                    // < bucket
        if (i < numsplitters && splitter[i] == key) b += 1;  // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    template <typename StringSet>
    void classify(
        const StringSet& strset,
        typename StringSet::Iterator begin, typename StringSet::Iterator end,
        uint16_t* bktout, size_t depth) const
    {
        while (begin != end)
        {
            key_type key = strset.get_uint64(*begin++, depth);
            *bktout++ = find_bkt_tree(key);
        }
    }

    //! return a splitter
    key_type get_splitter(unsigned int i) const
    { return splitter[i]; }

    /// build tree and splitter array from sample
    void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderFullDescent<numsplitters>(
            splitter, splitter_tree, splitter_lcp,
            samples, samplesize);
    }
};

template <size_t treebits>
class ClassifyUnrollBoth : public ClassifyUnrollTree<treebits>
{
public:
    /// search in splitter tree for bucket number, unrolled for U keys at once.
    template <int U>
    __attribute__ ((optimize("unroll-all-loops")))
    void
    find_bkt_tree_unroll(const key_type key[U], uint16_t obkt[U]) const
    {
        // binary tree traversal without left branch

        static const size_t numsplitters = this->numsplitters;
        const key_type* splitter = this->splitter;
        const key_type* splitter_tree = this->splitter_tree;

        unsigned int i[U];
        std::fill(i + 0, i + U, 1);

        for (size_t l = 0; l < treebits; ++l)
        {
            for (int u = 0; u < U; ++u)
            {
#if 0
                // in gcc-4.6.3 this produces a SETA, LEA sequence
                i[u] = 2 * i[u] + (key[u] <= splitter_tree[i[u]] ? 0 : 1);
#else
                // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which
                // is slightly faster
                if (key[u] <= splitter_tree[i[u]])
                    i[u] = 2 * i[u] + 0;
                else    // (key > splitter_tree[i[u]])
                    i[u] = 2 * i[u] + 1;
#endif
            }
        }

        for (int u = 0; u < U; ++u)
            i[u] -= numsplitters + 1;

        for (int u = 0; u < U; ++u) {
            // < bucket
            obkt[u] = i[u] * 2;
        }

        for (int u = 0; u < U; ++u)
        {
            // equal bucket
            if (i[u] < numsplitters && splitter[i[u]] == key[u]) obkt[u] += 1;
        }
    }

    /// classify all strings in area by walking tree and saving bucket id,
    /// unrolled loops
    template <typename StringSet>
    void classify(
        const StringSet& strset,
        typename StringSet::Iterator begin, typename StringSet::Iterator end,
        uint16_t* bktout, size_t depth) const
    {
        while (begin != end)
        {
            static const int rollout = 4;
            if (begin + rollout < end)
            {
                key_type key[rollout];
                key[0] = strset.get_uint64(begin[0], depth);
                key[1] = strset.get_uint64(begin[1], depth);
                key[2] = strset.get_uint64(begin[2], depth);
                key[3] = strset.get_uint64(begin[3], depth);

                find_bkt_tree_unroll<rollout>(key, bktout);

                begin += rollout;
                bktout += rollout;
            }
            else
            {
                // binary search in splitter with equal check
                key_type key = strset.get_uint64(*begin++, depth);
                *bktout++ = this->find_bkt_tree(key);
            }
        }
    }
};

// ****************************************************************************
// *** Insertion Sort Type-Switch

template <typename StringSet>
static inline void
insertion_sort(const stringtools::StringShadowPtr<StringSet>& strptr,
               size_t depth)
{
    assert(!strptr.flipped());

    if (!use_lcp_inssort)
        inssort::inssort_generic(strptr.output(), depth);
    else
        bingmann_lcp_inssort::lcp_insertion_sort_nolcp(strptr.output(), depth);
}

template <typename StringSet>
static inline void
insertion_sort(const stringtools::StringShadowLcpPtr<StringSet>& strptr,
               size_t depth)
{
    assert(!strptr.flipped());

    bingmann_lcp_inssort::lcp_insertion_sort</* SaveCache */ false, StringSet>(
        strptr.output(), strptr.lcparray(), nullptr, depth);
}

template <typename StringSet>
static inline void
insertion_sort(const stringtools::StringShadowOutPtr<StringSet>& strptr,
               size_t depth)
{
    assert(!strptr.flipped());

    if (!use_lcp_inssort)
        inssort::inssort_generic(strptr.output(), depth);
    else
        bingmann_lcp_inssort::lcp_insertion_sort_nolcp(strptr.output(), depth);
}

template <typename StringSet>
static inline void
insertion_sort(const stringtools::StringShadowLcpOutPtr<StringSet>& strptr,
               size_t depth)
{
    assert(!strptr.flipped());

    bingmann_lcp_inssort::lcp_insertion_sort(
        strptr.output(), strptr.lcparray(), depth);
}

template <typename StringSet>
static inline void
insertion_sort(const stringtools::StringShadowLcpCacheOutPtr<StringSet>& strptr,
               size_t depth)
{
    assert(!strptr.flipped());

    bingmann_lcp_inssort::lcp_insertion_sort</* SaveCache */ true>(
        strptr.output(), strptr.lcparray(), strptr.cache(), depth);
}

// ****************************************************************************
// *** LCP Calculation for finished Sample Sort Steps

template <size_t bktnum, typename Classify, typename StringPtr, typename BktSizeType>
void sample_sort_lcp(const Classify& classifier,
                     const StringPtr& strptr, size_t depth,
                     const BktSizeType* bkt)
{
    assert(!strptr.flipped());
    assert(strptr.check());

    const typename StringPtr::StringSet& strset = strptr.output();

    size_t b = 0;         // current bucket number
    key_type prevkey = 0; // previous key

    // the following while loops only check b < bktnum when b is odd,
    // because bktnum is always odd. We need a goto to jump into the loop,
    // as b == 0 start even.
    goto even_first;

    // find first non-empty bucket
    while (b < bktnum)
    {
        // odd bucket: = bkt
        if (bkt[b] != bkt[b + 1])
        {
            prevkey = classifier.get_splitter(b / 2);
            assert(prevkey == strset.get_uint64(strset.at(bkt[b + 1] - 1), depth));
            break;
        }
        ++b;
even_first:
        // even bucket: <, << or > bkt
        if (bkt[b] != bkt[b + 1])
        {
            prevkey = strset.get_uint64(strset.at(bkt[b + 1] - 1), depth);
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
        if (bkt[b] != bkt[b + 1])
        {
            key_type thiskey = classifier.get_splitter(b / 2);
            assert(thiskey == strset.get_uint64(strset.at(bkt[b]), depth));

            int rlcp = lcpKeyType(prevkey, thiskey);
            strptr.set_lcp(bkt[b], depth + rlcp);
            strptr.set_cache(bkt[b], getCharAtDepth(thiskey, rlcp));

            DBG(debug_lcp, "LCP at odd-bucket " << b <<
                " [" << bkt[b] << "," << bkt[b + 1] << ")" <<
                " is " << strptr.lcp(bkt[b]));

            prevkey = thiskey;
            assert(prevkey == strset.get_uint64(strset.at(bkt[b + 1] - 1), depth));
        }
        ++b;
even_bucket:
        // even bucket: <, << or > bkt
        if (bkt[b] != bkt[b + 1])
        {
            key_type thiskey = strset.get_uint64(strset.at(bkt[b]), depth);

            int rlcp = lcpKeyType(prevkey, thiskey);
            strptr.set_lcp(bkt[b], depth + rlcp);
            strptr.set_cache(bkt[b], getCharAtDepth(thiskey, rlcp));

            DBG(debug_lcp, "LCP at even-bucket " << b <<
                " [" << bkt[b] << "," << bkt[b + 1] << ")" <<
                " is " << strptr.lcp(bkt[b]));

            prevkey = strset.get_uint64(strset.at(bkt[b + 1] - 1), depth);
        }
        ++b;
    }
}

// ****************************************************************************
// *** SampleSort non-recursive in-place sequential sample sort for small sorts

template <template <size_t> class Classify, typename Context, typename StringPtr>
void Enqueue(Context& ctx, SortStep* sstep, const StringPtr& strptr, size_t depth);

template <typename Context, template <size_t> class Classify,
          typename StringPtr, typename BktSizeType>
class SmallsortJob : public Context::job_type, public SortStep
{
public:
    //! parent sort step
    SortStep* pstep;

    size_t thrid;

    StringPtr in_strptr;
    size_t in_depth;

    typedef typename StringPtr::StringSet StringSet;

    typedef BktSizeType bktsize_type;

    SmallsortJob(Context& ctx, SortStep* _pstep,
                 const StringPtr& strptr, size_t depth)
        : pstep(_pstep), in_strptr(strptr), in_depth(depth)
    {
        DBG(debug_steps, "enqueue depth=" << in_depth <<
            " size=" << in_strptr.size() << " flip=" << in_strptr.flipped());

        ctx.jobqueue.enqueue(this);
    }

    class SeqSampleSortStep
    {
    public:
#if 0
        static const size_t numsplitters2 = 2 * 16;
#else
        // bounding equations:
        // a) K * key_type splitter_tree size (and maybe K * key_type equal-cmp splitters)
        // b) 2 * K + 1 buckets (bktsize_type) when counting bkt occurances.
        static const size_t numsplitters2
            = (l2cache - sizeof(size_t)) / (2 * sizeof(size_t));
#endif

        static const size_t treebits = tlx::Log2Floor<numsplitters2>::value;
        static const size_t numsplitters = (1 << treebits) - 1;

        static const size_t bktnum = 2 * numsplitters + 1;

        StringPtr strptr;
        size_t idx;
        size_t depth;

        bktsize_type bkt[bktnum + 1];

#if CALC_LCP
        Classify<treebits> classifier;
#endif
        unsigned char splitter_lcp[numsplitters + 1];

        SeqSampleSortStep(Context& ctx, const StringPtr& _strptr, size_t _depth,
                          uint16_t* bktcache)
            : strptr(_strptr), idx(0), depth(_depth)
        {
            size_t n = strptr.size();

            // step 1: select splitters with oversampling

            const size_t oversample_factor = 2;
            const size_t samplesize = oversample_factor * numsplitters;

            key_type samples[samplesize];

            const StringSet& strset = strptr.active();
            typename StringSet::Iterator begin = strset.begin();

            LCGRandom rng(&samples);

            for (unsigned int i = 0; i < samplesize; ++i)
            {
                samples[i] = strset.get_uint64(strset[begin + rng() % n], depth);
            }

            std::sort(samples, samples + samplesize);

#if !CALC_LCP
            Classify<treebits> classifier;
#endif
            classifier.build(samples, samplesize, splitter_lcp);

            // step 2: classify all strings

            classifier.classify(
                strset, strset.begin(), strset.end(), bktcache, depth);

            // step 2.5: count bucket sizes

            bktsize_type bktsize[bktnum];
            memset(bktsize, 0, bktnum * sizeof(bktsize_type));

            for (size_t si = 0; si < n; ++si)
                ++bktsize[bktcache[si]];

            // step 3: inclusive prefix sum

            bkt[0] = bktsize[0];
            for (unsigned int i = 1; i < bktnum; ++i) {
                bkt[i] = bkt[i - 1] + bktsize[i];
            }
            assert(bkt[bktnum - 1] == n);
            bkt[bktnum] = n;

            // step 4: premute out-of-place

            const StringSet& strB = strptr.active();
            // get alternative shadow pointer array
            const StringSet& sorted = strptr.shadow();
            typename StringSet::Iterator sbegin = sorted.begin();

            for (typename StringSet::Iterator str = strB.begin();
                 str != strB.end(); ++str, ++bktcache)
                *(sbegin + --bkt[*bktcache]) = std::move(*str);

            // bkt is afterwards the exclusive prefix sum of bktsize

            // statistics

            ++ctx.seq_ss_steps;
        }

        void calculate_lcp()
        {
#if CALC_LCP
            //DBG(debug_lcp, "Calculate LCP after sample sort step " << strptr);
            sample_sort_lcp<bktnum>(classifier, strptr.original(), depth, bkt);
#endif
        }
    };

    // *** Stack of Recursive Sample Sort Steps

    uint16_t* bktcache;
    size_t bktcache_size;

    size_t ss_pop_front;
    std::vector<SeqSampleSortStep> ss_stack;

    virtual bool run(Context& ctx)
    {
        ScopedTimerKeeperMT tm_seq_ss(ctx.timers, TM_SEQ_SS);

        size_t n = in_strptr.size();

        DBG(debug_jobs, "Process SmallsortJob " << this << " of size " << n);

        thrid = PS5_ENABLE_RESTSIZE ? omp_get_thread_num() : 0;

        // create anonymous wrapper job
        substep_add();

        bktcache = NULL;
        bktcache_size = 0;
        ss_pop_front = 0;
        ms_pop_front = 0;

        if (n >= g_smallsort_threshold)
        {
            bktcache = new uint16_t[n];
            bktcache_size = n * sizeof(uint16_t);
            sort_sample_sort(ctx, in_strptr, in_depth);
        }
        else
        {
            sort_mkqs_cache(ctx, in_strptr, in_depth);
        }

        delete[] bktcache;

        // finish wrapper job, handler delete's this
        substep_notify_done();

        return false;
    }

    void sort_sample_sort(Context& ctx, const StringPtr& strptr, size_t depth)
    {
        typedef SeqSampleSortStep Step;

        assert(ss_pop_front == 0);
        assert(ss_stack.size() == 0);

        // sort first level
        ss_stack.push_back(Step(ctx, strptr, depth, bktcache));

        // step 5: "recursion"

        while (ss_stack.size() > ss_pop_front)
        {
            Step& s = ss_stack.back();
            size_t i = s.idx++; // process the bucket s.idx

            if (i < Step::bktnum)
            {
                size_t bktsize = s.bkt[i + 1] - s.bkt[i];

                StringPtr sp = s.strptr.flip(s.bkt[i], bktsize);

                // i is even -> bkt[i] is less-than bucket
                if (i % 2 == 0)
                {
                    if (bktsize == 0)
                        ;
                    else if (bktsize < g_smallsort_threshold)
                    {
                        assert(i / 2 <= Step::numsplitters);
                        if (i == Step::bktnum - 1)
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << bktsize << " no lcp");
                        else
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << bktsize << " lcp " << int(s.splitter_lcp[i / 2] & 0x7F));

                        sort_mkqs_cache(ctx, sp, s.depth + (s.splitter_lcp[i / 2] & 0x7F));
                    }
                    else
                    {
                        if (i == Step::bktnum - 1)
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << bktsize << " no lcp");
                        else
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << bktsize << " lcp " << int(s.splitter_lcp[i / 2] & 0x7F));

                        ss_stack.push_back(Step(ctx, sp, s.depth + (s.splitter_lcp[i / 2] & 0x7F), bktcache));
                    }
                }
                // i is odd -> bkt[i] is equal bucket
                else
                {
                    if (bktsize == 0)
                        ;
                    else if (s.splitter_lcp[i / 2] & 0x80) { // equal-bucket has NULL-terminated key, done.
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " is done!");
                        StringPtr spb = sp.copy_back();
#if CALC_LCP
                        spb.fill_lcp(s.depth + lcpKeyDepth(s.classifier.get_splitter(i / 2)));
#else
                        UNUSED(spb);
#endif
                        ctx.donesize(bktsize, thrid);
                    }
                    else if (bktsize < g_smallsort_threshold)
                    {
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " lcp keydepth!");

                        sort_mkqs_cache(ctx, sp, s.depth + sizeof(key_type));
                    }
                    else
                    {
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " lcp keydepth!");

                        ss_stack.push_back(Step(ctx, sp, s.depth + sizeof(key_type), bktcache));
                    }
                }
            }
            else
            {
                // finished sort
                assert(ss_stack.size() > ss_pop_front);

                // after full sort: calculate LCPs at this level
                ss_stack.back().calculate_lcp();

                ss_stack.pop_back();
            }

            if (use_work_sharing && ctx.jobqueue.has_idle()) {
                sample_sort_free_work(ctx);
            }
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

        while (s.idx < Step::bktnum)
        {
            size_t i = s.idx++; // process the bucket s.idx

            size_t bktsize = s.bkt[i + 1] - s.bkt[i];

            StringPtr sp = s.strptr.flip(s.bkt[i], bktsize);

            // i is even -> bkt[i] is less-than bucket
            if (i % 2 == 0)
            {
                if (bktsize == 0)
                    ;
                else
                {
                    if (i == Step::bktnum - 1)
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << bktsize << " no lcp");
                    else
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << bktsize << " lcp " << int(s.splitter_lcp[i / 2] & 0x7F));

                    substep_add();
                    Enqueue<Classify>(ctx, this, sp, s.depth + (s.splitter_lcp[i / 2] & 0x7F));
                }
            }
            // i is odd -> bkt[i] is equal bucket
            else
            {
                if (bktsize == 0)
                    ;
                else if (s.splitter_lcp[i / 2] & 0x80) { // equal-bucket has NULL-terminated key, done.
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " is done!");
                    StringPtr spb = sp.copy_back();
#if CALC_LCP
                    spb.fill_lcp(s.depth + lcpKeyDepth(s.classifier.get_splitter(i / 2)));
#else
                    UNUSED(spb);
#endif
                    ctx.donesize(bktsize, thrid);
                }
                else
                {
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << bktsize << " lcp keydepth!");

                    substep_add();
                    Enqueue<Classify>(ctx, this, sp, s.depth + sizeof(key_type));
                }
            }
        }

        // shorten the current stack
        ++ss_pop_front;
    }

    // *************************************************************************
    // *** Stack of Recursive MKQS Steps

    static inline int cmp(const key_type& a, const key_type& b)
    {
        return (a > b) ? 1 : (a < b) ? -1 : 0;
    }

    template <typename Type>
    static inline size_t
    med3(Type* A, size_t i, size_t j, size_t k)
    {
        if (A[i] == A[j]) return i;
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
    insertion_sort_cache_block(const StringPtr& strptr, key_type* cache)
    {
        const StringSet& strings = strptr.output();
        unsigned int n = strptr.size();
        unsigned int pi, pj;
        for (pi = 1; --n > 0; ++pi) {
            typename StringSet::String tmps = std::move(strings.at(pi));
            key_type tmpc = cache[pi];
            for (pj = pi; pj > 0; --pj) {
                if (cache[pj - 1] <= tmpc)
                    break;
                strings.at(pj) = std::move(strings.at(pj - 1));
                cache[pj] = cache[pj - 1];
            }
            strings.at(pj) = std::move(tmps);
            cache[pj] = tmpc;
        }
    }

    // Insertion sort, but use cached characters if possible.
    template <bool CacheDirty>
    static inline void
    insertion_sort_cache(const StringPtr& _strptr, key_type* cache, size_t depth)
    {
        StringPtr strptr = _strptr.copy_back();

        if (strptr.size() <= 1) return;
        if (CacheDirty) return insertion_sort(strptr, depth);

        insertion_sort_cache_block(strptr, cache);

        size_t start = 0, bktsize = 1;
        for (size_t i = 0; i < strptr.size() - 1; ++i) {
            // group areas with equal cache values
            if (cache[i] == cache[i + 1]) {
                ++bktsize;
                continue;
            }
            // calculate LCP between group areas
            if (start != 0) {
                int rlcp = lcpKeyType(cache[start - 1], cache[start]);
                strptr.set_lcp(start, depth + rlcp);
                strptr.set_cache(start, getCharAtDepth(cache[start], rlcp));
            }
            // sort group areas deeper if needed
            if (bktsize > 1) {
                if (cache[start] & 0xFF) // need deeper sort
                    insertion_sort(strptr.sub(start, bktsize), depth + sizeof(key_type));
                else {                   // cache contains NULL-termination
                    strptr.sub(start, bktsize).fill_lcp(depth + lcpKeyDepth(cache[start]));
                }
            }
            bktsize = 1;
            start = i + 1;
        }
        // tail of loop for last item
        if (start != 0) {
            int rlcp = lcpKeyType(cache[start - 1], cache[start]);
            strptr.set_lcp(start, depth + rlcp);
            strptr.set_cache(start, getCharAtDepth(cache[start], rlcp));
        }
        if (bktsize > 1) {
            if (cache[start] & 0xFF) // need deeper sort
                insertion_sort(strptr.sub(start, bktsize), depth + sizeof(key_type));
            else                     // cache contains NULL-termination
                strptr.sub(start, bktsize).fill_lcp(depth + lcpKeyDepth(cache[start]));
        }
    }

    class MKQSStep
    {
    public:
        StringPtr strptr;
        key_type* cache;
        size_t num_lt, num_eq, num_gt, depth;
        size_t idx;
        unsigned char eq_recurse;
#if CALC_LCP_MKQS == 1
        char_type dchar_eq, dchar_gt;
        unsigned char lcp_lt, lcp_eq, lcp_gt;
#elif CALC_LCP_MKQS == 2
        key_type pivot;
#endif

        MKQSStep(Context& ctx, const StringPtr& _strptr,
                 key_type* _cache, size_t _depth, bool CacheDirty)
            : strptr(_strptr), cache(_cache), depth(_depth), idx(0)
        {
            size_t n = strptr.size();

            const StringSet& strset = strptr.active();

            if (CacheDirty) {
                typename StringSet::Iterator it = strset.begin();
                for (size_t i = 0; i < n; ++i, ++it) {
                    cache[i] = strset.get_uint64(*it, depth);
                }
            }
            // select median of 9
            size_t p = med3(cache,
                            med3(cache, 0, n / 8, n / 4),
                            med3(cache, n / 2 - n / 8, n / 2, n / 2 + n / 8),
                            med3(cache, n - 1 - n / 4, n - 1 - n / 8, n - 3)
                            );
            // swap pivot to first position
            std::swap(strset.at(0), strset.at(p));
            std::swap(cache[0], cache[p]);
            // save the pivot value
            key_type pivot = cache[0];
#if CALC_LCP_MKQS == 1
            // for immediate LCP calculation
            key_type max_lt = 0, min_gt = std::numeric_limits<key_type>::max();
#elif CALC_LCP_MKQS == 2
            this->pivot = pivot;
#endif
            // indexes into array: 0 [pivot] 1 [===] leq [<<<] llt [???] rgt [>>>] req [===] n-1
            size_t leq = 1, llt = 1, rgt = n - 1, req = n - 1;
            while (true)
            {
                while (llt <= rgt)
                {
                    int r = cmp(cache[llt], pivot);
                    if (r > 0) {
#if CALC_LCP_MKQS == 1
                        min_gt = std::min(min_gt, cache[llt]);
#endif
                        break;
                    }
                    else if (r == 0) {
                        std::swap(strset.at(leq), strset.at(llt));
                        std::swap(cache[leq], cache[llt]);
                        leq++;
                    }
                    else {
#if CALC_LCP_MKQS == 1
                        max_lt = std::max(max_lt, cache[llt]);
#endif
                    }
                    ++llt;
                }
                while (llt <= rgt)
                {
                    int r = cmp(cache[rgt], pivot);
                    if (r < 0) {
#if CALC_LCP_MKQS == 1
                        max_lt = std::max(max_lt, cache[rgt]);
#endif
                        break;
                    }
                    else if (r == 0) {
                        std::swap(strset.at(req), strset.at(rgt));
                        std::swap(cache[req], cache[rgt]);
                        req--;
                    }
                    else {
#if CALC_LCP_MKQS == 1
                        min_gt = std::min(min_gt, cache[rgt]);
#endif
                    }
                    --rgt;
                }
                if (llt > rgt)
                    break;
                std::swap(strset.at(llt), strset.at(rgt));
                std::swap(cache[llt], cache[rgt]);
                ++llt;
                --rgt;
            }
            // calculate size of areas = < and >, save into struct
            size_t num_leq = leq, num_req = n - 1 - req;
            num_eq = num_leq + num_req;
            num_lt = llt - leq;
            num_gt = req - rgt;
            assert(num_eq > 0);
            assert(num_lt + num_eq + num_gt == n);

            // swap equal values from left to center
            const size_t size1 = std::min(num_leq, num_lt);
            std::swap_ranges(strset.begin(), strset.begin() + size1,
                             strset.begin() + llt - size1);
            std::swap_ranges(cache, cache + size1, cache + llt - size1);

            // swap equal values from right to center
            const size_t size2 = std::min(num_req, num_gt);
            std::swap_ranges(strset.begin() + llt, strset.begin() + llt + size2,
                             strset.begin() + n - size2);
            std::swap_ranges(cache + llt, cache + llt + size2,
                             cache + n - size2);

            // No recursive sorting if pivot has a zero byte
            this->eq_recurse = (pivot & 0xFF);

#if CALC_LCP_MKQS == 1
            // save LCP values for writing into LCP array after sorting further
            if (num_lt > 0)
            {
                assert(max_lt == *std::max_element(cache + 0, cache + num_lt));

                lcp_lt = lcpKeyType(max_lt, pivot);
                dchar_eq = getCharAtDepth(pivot, lcp_lt);
                DBG(debug_lcp, "LCP lt with pivot: " << depth + lcp_lt);
            }

            // calculate equal area lcp: +1 for the equal zero termination byte
            lcp_eq = lcpKeyDepth(pivot);

            if (num_gt > 0)
            {
                assert(min_gt == *std::min_element(cache + num_lt + num_eq, cache + n));

                lcp_gt = lcpKeyType(pivot, min_gt);
                dchar_gt = getCharAtDepth(min_gt, lcp_gt);
                DBG(debug_lcp, "LCP pivot with gt: " << depth + lcp_gt);
            }
#endif
            ++ctx.bs_steps;
        }

        void calculate_lcp()
        {
#if CALC_LCP_MKQS == 1
            if (num_lt > 0)
            {
                strptr.original().set_lcp(num_lt, depth + lcp_lt);
                strptr.original().set_cache(num_lt, dchar_eq);
            }

            if (num_gt > 0)
            {
                strptr.original().set_lcp(num_lt + num_eq, depth + lcp_gt);
                strptr.original().set_cache(num_lt + num_eq, dchar_gt);
            }
#elif CALC_LCP_MKQS == 2
            if (num_lt > 0)
            {
                key_type max_lt = strptr.original().output().get_uint64(
                    strptr.original().out(num_lt - 1), depth);

                unsigned int rlcp = lcpKeyType(max_lt, pivot);
                DBG(debug_lcp, "LCP lt with pivot: " << depth + rlcp);

                strptr.original().set_lcp(num_lt, depth + rlcp);
                strptr.original().set_cache(num_lt, getCharAtDepth(pivot, rlcp));
            }
            if (num_gt > 0)
            {
                key_type min_gt = strptr.original().output().get_uint64(
                    strptr.original().out(num_lt + num_eq), depth);

                unsigned int rlcp = lcpKeyType(pivot, min_gt);
                DBG(debug_lcp, "LCP pivot with gt: " << depth + rlcp);

                strptr.original().set_lcp(num_lt + num_eq, depth + rlcp);
                strptr.original().set_cache(num_lt + num_eq, getCharAtDepth(min_gt, rlcp));
            }
#endif
        }
    };

    size_t ms_pop_front;
    std::vector<MKQSStep> ms_stack;

    void sort_mkqs_cache(Context& ctx, const StringPtr& strptr, size_t depth)
    {
        ScopedTimerKeeperMT tm_mkqs(ctx.timers, TM_MKQS);

        if (strptr.size() < g_inssort_threshold) {
            ScopedTimerKeeperMT tm_inssort(ctx.timers, TM_INSSORT);
            insertion_sort(strptr.copy_back(), depth);
            ctx.donesize(strptr.size(), thrid);
            return;
        }

        if (bktcache_size < strptr.size() * sizeof(key_type)) {
            delete[] bktcache;
            bktcache = (uint16_t*)new key_type[strptr.size()];
            bktcache_size = strptr.size() * sizeof(key_type);
        }

        key_type* cache = (key_type*)bktcache; // reuse bktcache as keycache

        assert(ms_pop_front == 0);
        assert(ms_stack.size() == 0);

        // std::deque is much slower than std::vector, so we use an artificial
        // pop_front variable.
        ms_stack.push_back(MKQSStep(ctx, strptr, cache, depth, true));

        while (ms_stack.size() > ms_pop_front)
        {
            MKQSStep& ms = ms_stack.back();
            ++ms.idx; // increment here, because stack may change

            // process the lt-subsequence
            if (ms.idx == 1)
            {
                if (ms.num_lt == 0)
                    ;
                else if (ms.num_lt < g_inssort_threshold) {
                    ScopedTimerKeeperMT tm_inssort(ctx.timers, TM_INSSORT);
                    insertion_sort_cache<false>(ms.strptr.sub(0, ms.num_lt),
                                                ms.cache, ms.depth);
                    ctx.donesize(ms.num_lt, thrid);
                }
                else {
                    ms_stack.push_back(
                        MKQSStep(ctx,
                                 ms.strptr.sub(0, ms.num_lt),
                                 ms.cache, ms.depth, false));
                }
            }
            // process the eq-subsequence
            else if (ms.idx == 2)
            {
                StringPtr sp = ms.strptr.sub(ms.num_lt, ms.num_eq);

                assert(ms.num_eq > 0);

                if (!ms.eq_recurse) {
                    StringPtr spb = sp.copy_back();
#if CALC_LCP_MKQS == 1
                    spb.fill_lcp(ms.depth + ms.lcp_eq);
#elif CALC_LCP_MKQS == 2
                    spb.fill_lcp(ms.depth + lcpKeyDepth(ms.pivot));
#endif
                    ctx.donesize(spb.size(), thrid);
                }
                else if (ms.num_eq < g_inssort_threshold) {
                    ScopedTimerKeeperMT tm_inssort(ctx.timers, TM_INSSORT);
                    insertion_sort_cache<true>(sp, ms.cache + ms.num_lt,
                                               ms.depth + sizeof(key_type));
                    ctx.donesize(ms.num_eq, thrid);
                }
                else {
                    ms_stack.push_back(
                        MKQSStep(ctx, sp,
                                 ms.cache + ms.num_lt,
                                 ms.depth + sizeof(key_type), true));
                }
            }
            // process the gt-subsequence
            else if (ms.idx == 3)
            {
                StringPtr sp = ms.strptr.sub(ms.num_lt + ms.num_eq, ms.num_gt);

                if (ms.num_gt == 0)
                    ;
                else if (ms.num_gt < g_inssort_threshold) {
                    ScopedTimerKeeperMT tm_inssort(ctx.timers, TM_INSSORT);
                    insertion_sort_cache<false>(sp, ms.cache + ms.num_lt + ms.num_eq,
                                                ms.depth);
                    ctx.donesize(ms.num_gt, thrid);
                }
                else {
                    ms_stack.push_back(
                        MKQSStep(ctx, sp,
                                 ms.cache + ms.num_lt + ms.num_eq,
                                 ms.depth, false));
                }
            }
            // calculate lcps
            else
            {
                // finished sort
                assert(ms_stack.size() > ms_pop_front);

                // calculate LCP after the three parts are sorted
                ms_stack.back().calculate_lcp();

                ms_stack.pop_back();
            }

            if (use_work_sharing && ctx.jobqueue.has_idle()) {
                sample_sort_free_work(ctx);
            }
        }
    }

    void mkqs_free_work(Context& ctx)
    {
        assert(ms_stack.size() >= ms_pop_front);

        for (unsigned int fl = 0; fl < 8; ++fl)
        {
            if (ms_stack.size() == ms_pop_front) {
                return;
            }

            DBG(debug_jobs, "Freeing top level of SmallsortJob's mkqs stack - size " << ms_stack.size());

            // convert top level of stack into independent jobs

            MKQSStep& ms = ms_stack[ms_pop_front];

            if (ms.idx == 0 && ms.num_lt != 0)
            {
                substep_add();
                Enqueue<Classify>(ctx, this, ms.strptr.sub(0, ms.num_lt), ms.depth);
            }
            if (ms.idx <= 1) // st.num_eq > 0 always
            {
                assert(ms.num_eq > 0);

                StringPtr sp = ms.strptr.sub(ms.num_lt, ms.num_eq);

                if (ms.eq_recurse) {
                    substep_add();
                    Enqueue<Classify>(ctx, this, sp,
                                      ms.depth + sizeof(key_type));
                }
                else {
                    StringPtr spb = sp.copy_back();
#if CALC_LCP_MKQS == 1
                    spb.fill_lcp(ms.depth + ms.lcp_eq);
#elif CALC_LCP_MKQS == 2
                    spb.fill_lcp(ms.depth + lcpKeyDepth(ms.pivot));
#else
                    UNUSED(spb);
#endif
                    ctx.donesize(ms.num_eq, thrid);
                }
            }
            if (ms.idx <= 2 && ms.num_gt != 0)
            {
                substep_add();
                Enqueue<Classify>(ctx, this,
                                  ms.strptr.sub(ms.num_lt + ms.num_eq, ms.num_gt), ms.depth);
            }

            // shorten the current stack
            ++ms_pop_front;
        }
    }

    virtual void substep_all_done()
    {
#if CALC_LCP
        DBG(debug_recursion, "SmallSort[" << in_depth << "] all substeps done -> LCP calculation");

        while (ms_pop_front > 0) {
            DBG(debug_lcp, "SmallSort[" << in_depth << "] ms_pop_front: " << ms_pop_front);
            ms_stack[--ms_pop_front].calculate_lcp();
        }

        while (ss_pop_front > 0) {
            DBG(debug_lcp, "SmallSort[" << in_depth << "] ss_pop_front: " << ss_pop_front);
            ss_stack[--ss_pop_front].calculate_lcp();
        }

        if (pstep) pstep->substep_notify_done();

        delete this;
#endif
    }
};

// ****************************************************************************
// *** SampleSortStep out-of-place parallel sample sort with separate Jobs

template <typename Context, template <size_t> class Classify, typename StringPtr>
class SampleSortStep : public SortStep
{
public:
#if 0
    static const size_t numsplitters2 = 16;
#else
    // bounding equation: 2*K+1 buckets when counting bkt occurances.
    static const size_t numsplitters2 = (l2cache - sizeof(size_t)) / (2 * sizeof(size_t));
#endif
    static const size_t treebits = tlx::Log2Floor<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;

    static const size_t bktnum = 2 * numsplitters + 1;

    typedef typename StringPtr::StringSet StringSet;
    typedef typename StringSet::Iterator StrIterator;

    //! type of Job
    typedef typename Context::job_type job_type;

    //! parent sort step notification
    SortStep* pstep;

    //! string pointers, size, and current sorting depth
    StringPtr strptr;
    size_t depth;

    //! number of parts into which the strings were split
    unsigned int parts;
    //! size of all parts except the last
    size_t psize;
    //! number of threads still working
    std::atomic<unsigned int> pwork;

    //! classifier instance and variables (contains splitter tree
    Classify<treebits> classifier;
    //! LCPs of splitters, needed for recursive calls
    unsigned char splitter_lcp[numsplitters + 1];

    //! individual bucket array of threads, keep bkt[0] for DistributeJob
    size_t* bkt[MAXPROCS];
    //! bucket ids cache, created by classifier and later counted
    uint16_t* bktcache[MAXPROCS];

    // *** Classes for JobQueue

    struct SampleJob : public job_type
    {
        SampleSortStep* step;

        SampleJob(SampleSortStep* _step)
            : step(_step) { }

        virtual bool run(Context& ctx)
        {
            ScopedTimerKeeperMT tm_seq_ss(ctx.timers, TM_PARA_SS);

            step->sample(ctx);
            return true;
        }
    };

    struct CountJob : public job_type
    {
        SampleSortStep* step;
        unsigned int  p;

        CountJob(SampleSortStep* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual bool run(Context& ctx)
        {
            ScopedTimerKeeperMT tm_seq_ss(ctx.timers, TM_PARA_SS);

            step->count(p, ctx);
            return true;
        }
    };

    struct DistributeJob : public job_type
    {
        SampleSortStep* step;
        unsigned int  p;

        DistributeJob(SampleSortStep* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual bool run(Context& ctx)
        {
            ScopedTimerKeeperMT tm_seq_ss(ctx.timers, TM_PARA_SS);

            step->distribute(p, ctx);
            return true;
        }
    };

    // *** Constructor

    SampleSortStep(Context& ctx, SortStep* _pstep,
                   const StringPtr& _strptr, size_t _depth)
        : pstep(_pstep), strptr(_strptr), depth(_depth)
    {
        parts = strptr.size() / ctx.sequential_threshold() * 2;
        if (parts == 0) parts = 1;
        if (parts > MAXPROCS) parts = MAXPROCS;

        psize = (strptr.size() + parts - 1) / parts;

        DBG(debug_steps, "enqueue depth=" << depth << " size=" << strptr.size() << " parts=" << parts << " psize=" << psize << " flip=" << strptr.flipped());

        ctx.jobqueue.enqueue(new SampleJob(this));
        ++ctx.para_ss_steps;
    }

    // *** Sample Step

    void sample(Context& ctx)
    {
        DBG(debug_jobs, "Process SampleJob @ " << this);

        const size_t oversample_factor = 2;
        size_t samplesize = oversample_factor * numsplitters;

        const StringSet& strset = strptr.active();
        StrIterator begin = strset.begin();
        size_t n = strset.size();

        key_type samples[samplesize];

        LCGRandom rng(&samples);

        for (unsigned int i = 0; i < samplesize; ++i)
        {
            samples[i] = strset.get_uint64(strset[begin + rng() % n], depth);
        }

        std::sort(samples, samples + samplesize);

        classifier.build(samples, samplesize, splitter_lcp);

        // create new jobs
        pwork = parts;
        for (unsigned int p = 0; p < parts; ++p)
            ctx.jobqueue.enqueue(new CountJob(this, p));
    }

    // *** Counting Step

    void count(unsigned int p, Context& ctx)
    {
        DBG(debug_jobs, "Process CountJob " << p << " @ " << this);

        const StringSet& strset = strptr.active();

        StrIterator strB = strset.begin() + p * psize;
        StrIterator strE = strset.begin() + std::min((p + 1) * psize, strptr.size());
        if (strE < strB) strE = strB;

        uint16_t* mybktcache = bktcache[p] = new uint16_t[strE - strB];
        classifier.classify(strset, strB, strE, mybktcache, depth);

        size_t* mybkt = bkt[p] = new size_t[bktnum + (p == 0 ? 1 : 0)];
        memset(mybkt, 0, bktnum * sizeof(size_t));

        for (uint16_t* bc = mybktcache; bc != mybktcache + (strE - strB); ++bc)
            ++mybkt[*bc];

        if (--pwork == 0)
            count_finished(ctx);
    }

    void count_finished(Context& ctx)
    {
        DBG(debug_jobs, "Finishing CountJob " << this << " with prefixsum");

        // abort sorting if we're measuring only the top level
        if (use_only_first_sortstep)
            return;

        // inclusive prefix sum over bkt
        size_t sum = 0;
        for (unsigned int i = 0; i < bktnum; ++i)
        {
            for (unsigned int p = 0; p < parts; ++p)
            {
                bkt[p][i] = (sum += bkt[p][i]);
            }
        }
        assert(sum == strptr.size());

        // create new jobs
        pwork = parts;
        for (unsigned int p = 0; p < parts; ++p)
            ctx.jobqueue.enqueue(new DistributeJob(this, p));
    }

    // *** Distribute Step

    void distribute(unsigned int p, Context& ctx)
    {
        DBG(debug_jobs, "Process DistributeJob " << p << " @ " << this);

        const StringSet& strset = strptr.active();

        StrIterator strB = strset.begin() + p * psize;
        StrIterator strE = strset.begin() + std::min((p + 1) * psize, strptr.size());
        if (strE < strB) strE = strB;

        const StringSet& sorted = strptr.shadow(); // get alternative shadow pointer array
        typename StringSet::Iterator sbegin = sorted.begin();

        uint16_t* mybktcache = bktcache[p];
        size_t* mybkt = bkt[p];

        for (StrIterator str = strB; str != strE; ++str, ++mybktcache)
            *(sbegin + --mybkt[*mybktcache]) = std::move(*str);

        if (p != 0) // p = 0 is needed for recursion into bkts
            delete[] bkt[p];

        delete[] bktcache[p];

        if (--pwork == 0)
            distribute_finished(ctx);
    }

    void distribute_finished(Context& ctx)
    {
        DBG(debug_jobs, "Finishing DistributeJob " << this << " with enqueuing subjobs");

        size_t thrid = PS5_ENABLE_RESTSIZE ? omp_get_thread_num() : 0;

        size_t* bkt = this->bkt[0];
        assert(bkt);

        // first processor's bkt pointers are boundaries between bkts, just add sentinel:
        assert(bkt[0] == 0);
        bkt[bktnum] = strptr.size();

        // keep anonymous subjob handle while creating subjobs
        substep_add();

        size_t i = 0;
        while (i < bktnum - 1)
        {
            // i is even -> bkt[i] is less-than bucket
            size_t bktsize = bkt[i + 1] - bkt[i];
            if (bktsize == 0)
                ;
            else if (bktsize == 1) { // just one string pointer, copyback
                strptr.flip(bkt[i], 1).copy_back();
                ctx.donesize(1, thrid);
            }
            else
            {
                DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bkt[i] << " size " << bktsize << " lcp " << int(splitter_lcp[i / 2] & 0x7F));
                substep_add();
                Enqueue<Classify>(ctx, this, strptr.flip(bkt[i], bktsize), depth + (splitter_lcp[i / 2] & 0x7F));
            }
            ++i;
            // i is odd -> bkt[i] is equal bucket
            bktsize = bkt[i + 1] - bkt[i];
            if (bktsize == 0)
                ;
            else if (bktsize == 1) { // just one string pointer, copyback
                strptr.flip(bkt[i], 1).copy_back();
                ctx.donesize(1, thrid);
            }
            else
            {
                if (splitter_lcp[i / 2] & 0x80) { // equal-bucket has NULL-terminated key, done.
                    DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " is done!");
                    StringPtr sp = strptr.flip(bkt[i], bktsize).copy_back();
                    sp.fill_lcp(depth + lcpKeyDepth(classifier.get_splitter(i / 2)));
                    ctx.donesize(bktsize, thrid);
                }
                else {
                    DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bkt[i] << " size " << bktsize << " lcp keydepth!");
                    substep_add();
                    Enqueue<Classify>(ctx, this, strptr.flip(bkt[i], bktsize), depth + sizeof(key_type));
                }
            }
            ++i;
        }

        size_t bktsize = bkt[i + 1] - bkt[i];

        if (bktsize == 0)
            ;
        else if (bktsize == 1) { // just one string pointer, copyback
            strptr.flip(bkt[i], 1).copy_back();
            ctx.donesize(1, thrid);
        }
        else
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bkt[i] << " size " << bktsize << " no lcp");
            substep_add();
            Enqueue<Classify>(ctx, this, strptr.flip(bkt[i], bktsize), depth);
        }

        substep_notify_done(); // release anonymous subjob handle

#if !CALC_LCP
        delete[] bkt;
        delete this;
#endif
    }

    // *** After Recursive Sorting

#if CALC_LCP
    virtual void substep_all_done()
    {
        DBG(debug_steps, "pSampleSortStep[" << depth << "]: all substeps done.");

        sample_sort_lcp<bktnum>(classifier, strptr.original(), depth, bkt[0]);
        delete[] bkt[0];

        if (pstep) pstep->substep_notify_done();
        delete this;
    }
#endif

    static inline void put_stats()
    {
        g_stats >> "l2cache" << size_t(l2cache)
            >> "splitter_treebits" << size_t(treebits)
            >> "numsplitters" << size_t(numsplitters)
            >> "use_work_sharing" << use_work_sharing
            >> "use_restsize" << PS5_ENABLE_RESTSIZE
            >> "use_lcp_inssort" << use_lcp_inssort;
    }
};

template <template <size_t> class Classify, typename Context, typename StringPtr>
void Enqueue(Context& ctx, SortStep* pstep,
             const StringPtr& strptr, size_t depth)
{
    if (strptr.size() > ctx.sequential_threshold() || use_only_first_sortstep) {
        new SampleSortStep<Context, Classify, StringPtr>(ctx, pstep, strptr, depth);
    }
    else {
        if (strptr.size() < ((uint64_t)1 << 32))
            new SmallsortJob<Context, Classify, StringPtr, uint32_t>
                (ctx, pstep, strptr, depth);
        else
            new SmallsortJob<Context, Classify, StringPtr, uint64_t>
                (ctx, pstep, strptr, depth);
    }
}

template <template <size_t> class Classify = ClassifyUnrollBoth,
          typename StringPtr>
void parallel_sample_sort(const StringPtr& strptr, size_t depth)
{
    Context<> ctx;
    ctx.totalsize = strptr.size();
#if PS5_ENABLE_RESTSIZE
    ctx.restsize = strptr.size();
#endif
    ctx.threadnum = omp_get_max_threads();

    SampleSortStep<Context<>, Classify, StringPtr>::put_stats();

    ctx.timers.start(ctx.threadnum);

    Enqueue<Classify>(ctx, NULL, strptr, depth);
    ctx.jobqueue.loop();

    ctx.timers.stop();

#if PS5_ENABLE_RESTSIZE
    assert(!PS5_ENABLE_RESTSIZE || ctx.restsize.update().get() == 0);
#endif

    g_stats >> "steps_para_sample_sort" << ctx.para_ss_steps
        >> "steps_seq_sample_sort" << ctx.seq_ss_steps
        >> "steps_base_sort" << ctx.bs_steps;

    if (ctx.timers.is_real)
    {
        g_stats >> "tm_waiting" << ctx.timers.get(TM_WAITING)
            >> "tm_jq_work" << ctx.jobqueue.m_timers.get(ctx.jobqueue.TM_WORK)
            >> "tm_jq_idle" << ctx.jobqueue.m_timers.get(ctx.jobqueue.TM_IDLE)
            >> "tm_para_ss" << ctx.timers.get(TM_PARA_SS)
            >> "tm_seq_ss" << ctx.timers.get(TM_SEQ_SS)
            >> "tm_mkqs" << ctx.timers.get(TM_MKQS)
            >> "tm_inssort" << ctx.timers.get(TM_INSSORT)
            >> "tm_sum" << ctx.timers.get_sum();
    }
}

template <template <size_t> class Classify = ClassifyUnrollBoth,
          typename StringSet>
void parallel_sample_sort_base(const StringSet& strset, size_t depth)
{
    typedef stringtools::StringShadowPtr<StringSet> StringPtr;
    typedef typename StringSet::Container Container;

    // allocate shadow pointer array
    Container shadow = strset.allocate(strset.size());

    StringPtr strptr(strset, StringSet(shadow));

    parallel_sample_sort<Classify, StringPtr>(strptr, depth);

    StringSet::deallocate(shadow);
}

template <template <size_t> class Classify>
void parallel_sample_sort_base(string* strings, size_t n, size_t depth)
{
    return parallel_sample_sort_base<Classify>(
        UCharStringSet(strings, strings + n), depth);
}

template <template <size_t> class Classify = ClassifyUnrollBoth,
          typename StringSet>
void parallel_sample_sort_out_base(
    const StringSet& strset, const StringSet& output, size_t depth)
{
    typedef stringtools::StringShadowOutPtr<StringSet> StringOutPtr;
    typedef typename StringSet::Container Container;

    // allocate shadow pointer array
    Container shadow = strset.allocate(strset.size());

    StringOutPtr strptr(strset, StringSet(shadow), output);

    parallel_sample_sort<Classify, StringOutPtr>(strptr, depth);

    StringSet::deallocate(shadow);
}

template <template <size_t> class Classify>
void parallel_sample_sort_out_base(
    string* strings, string* output, size_t n, size_t depth)
{
    return parallel_sample_sort_out_base<Classify>(
        UCharStringSet(strings, strings + n),
        UCharStringSet(output, output + n),
        depth);
}

template <template <size_t> class Classify = ClassifyUnrollBoth,
          typename StringSet>
void parallel_sample_sort_out_test(const StringSet& strset, size_t depth)
{
    typename StringSet::Container out = strset.allocate(strset.size());
    StringSet output(out);
    parallel_sample_sort_out_base(strset, output, depth);

    // move strings back to strset
    std::move(output.begin(), output.end(), strset.begin());

    StringSet::deallocate(out);
}

/******************************************************************************/

#if CALC_LCP

template <template <size_t> class Classify = ClassifyUnrollBoth,
          typename StringSet>
void parallel_sample_sort_lcp_base(
    const StringSet& strset, uintptr_t* lcp, size_t depth)
{
    typedef stringtools::StringShadowLcpPtr<StringSet> StringPtr;
    typedef typename StringSet::Container Container;

    // allocate shadow pointer array
    Container shadow = strset.allocate(strset.size());

    StringPtr strptr(strset, StringSet(shadow), lcp);

    parallel_sample_sort<Classify, StringPtr>(strptr, depth);

    StringSet::deallocate(shadow);
}

template <template <size_t> class Classify = ClassifyUnrollBoth,
          typename StringSet>
void parallel_sample_sort_lcp_verify(const StringSet& strset, size_t depth)
{
    std::vector<uintptr_t> tmp_lcp(strset.size());
    tmp_lcp[0] = 42;                 // must keep lcp[0] unchanged
    std::fill(tmp_lcp.begin() + 1, tmp_lcp.end(), -1);
    parallel_sample_sort_lcp_base<Classify>(strset, tmp_lcp.data(), depth);
    assert(stringtools::verify_lcp(strset, tmp_lcp.data(), 42));
}

//! Call for NUMA aware parallel sorting
void parallel_sample_sort_numa(string* strings, size_t n,
                               int numaNode, int numberOfThreads,
                               const LcpCacheStringPtr& output)
{
    // tie thread to a NUMA node
    numa_run_on_node(numaNode);
    numa_set_preferred(numaNode);

    Context<> ctx;
    ctx.totalsize = n;
#if PS5_ENABLE_RESTSIZE
    ctx.restsize = n;
#endif
    ctx.threadnum = numberOfThreads;

    typedef UCharStringSet StringSet;
    StringSet strset(strings, strings + n);
    StringSet shadow((string*)output.lcps, (string*)output.lcps + n);
    StringSet outputss((string*)output.strings, (string*)output.strings + n);

    StringShadowLcpCacheOutPtr<StringSet>
    strptr(strset, shadow, outputss, output.lcps, output.cachedChars);

    Enqueue<ClassifyUnrollBoth>(ctx, NULL, strptr, 0);
    ctx.jobqueue.numaLoop(numaNode, numberOfThreads);

    // fixup first entry of LCP and charcache
    output.firstLcp() = 0;
    output.firstCached() = output.firstString()[0];

#if PS5_ENABLE_RESTSIZE
    assert(ctx.restsize.update().get() == 0);
#endif

    g_stats >> "steps_para_sample_sort" << ctx.para_ss_steps
        >> "steps_seq_sample_sort" << ctx.seq_ss_steps
        >> "steps_base_sort" << ctx.bs_steps;
}

//! Call for NUMA aware parallel sorting
void parallel_sample_sort_numa2(const UCharStringShadowLcpCacheOutPtr* strptr,
                                unsigned numInputs)
{
    typedef Context<NumaJobQueueGroup> context_type;

    context_type::jobqueuegroup_type group;

    // construct one Context per input
    context_type* ctx[numInputs];

    for (unsigned i = 0; i < numInputs; ++i)
    {
        ctx[i] = new context_type(&group);

        ctx[i]->totalsize = strptr[i].size();
#if PS5_ENABLE_RESTSIZE
        ctx[i]->restsize = strptr[i].size();
#endif
        ctx[i]->threadnum = group.calcThreadNum(i, numInputs);
        if (ctx[i]->threadnum == 0)
            ctx[i]->threadnum = 1;

        Enqueue<ClassifyUnrollBoth>(*ctx[i], NULL, strptr[i], 0);

        group.add_jobqueue(&ctx[i]->jobqueue);
    }

    group.numaLaunch();

    for (unsigned i = 0; i < numInputs; ++i)
    {
        // fixup first entry of LCP and charcache
        strptr[i].lcparray()[0] = 0;
        strptr[i].set_cache(0, strptr[i].out(0)[0]);

#if PS5_ENABLE_RESTSIZE
        assert(ctx[i].restsize.update().get() == 0);
#endif

        g_stats >> "steps_para_sample_sort" << ctx[i]->para_ss_steps
            >> "steps_seq_sample_sort" << ctx[i]->seq_ss_steps
            >> "steps_base_sort" << ctx[i]->bs_steps;

        delete ctx[i];
    }
}

#endif // CALC_LCP

/******************************************************************************/

static inline void
parallel_sample_sortBTC(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifySimple>(strings, n, 0);
}

PSS_CONTESTANT_PARALLEL_LCP(
    parallel_sample_sortBTC,
    "bingmann/parallel_sample_sortBTC",
    "pS5: binary tree, bktcache")

static inline void
parallel_sample_sortBTCU1(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyUnrollTree>(strings, n, 0);
}

PSS_CONTESTANT_PARALLEL_LCP(
    parallel_sample_sortBTCU1,
    "bingmann/parallel_sample_sortBTCU1",
    "pS5: binary tree, bktcache, unroll tree")

static inline void
parallel_sample_sortBTCU2(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyUnrollBoth>(strings, n, 0);
}

PSS_CONTESTANT_PARALLEL_LCP(
    parallel_sample_sortBTCU2,
    "bingmann/parallel_sample_sortBTCU2",
    "pS5: binary tree, bktcache, unroll tree and strings")

static inline void
parallel_sample_sortBTCU2_out(string* strings, size_t n)
{
    string* output = new string[n];

    parallel_sample_sort_out_base<ClassifyUnrollBoth>(strings, output, n, 0);

    // copy back for verification
    memcpy(strings, output, n * sizeof(string));
    delete[] output;
}

PSS_CONTESTANT_PARALLEL_LCP(
    parallel_sample_sortBTCU2_out,
    "bingmann/parallel_sample_sortBTCU2_out",
    "pS5: binary tree, bktcache, unroll tree and strings, separate output")

////////////////////////////////////////////////////////////////////////////////

// ****************************************************************************
// *** Classification with Binary Search in Splitter Array (old BSC variant)

template <size_t treebits>
class ClassifyBinarySearch
{
public:
    // NOTE: for binary search numsplitters need not be 2^k-1, any size will
    // do, but the tree implementations are always faster, so we keep this only
    // for historical reasons.
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter[numsplitters];

    /// binary search on splitter array for bucket number
    unsigned int find_bkt_tree(const key_type& key) const
    {
        unsigned int lo = 0, hi = numsplitters;

        while (lo < hi)
        {
            unsigned int mid = (lo + hi) >> 1;
            assert(mid < numsplitters);

            if (key <= splitter[mid])
                hi = mid;
            else                                               // (key > splitter[mid])
                lo = mid + 1;
        }

        size_t b = lo * 2;                                     // < bucket
        if (lo < numsplitters && splitter[lo] == key) b += 1;  // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    template <typename StringSet>
    void classify(
        const StringSet& strset,
        typename StringSet::Iterator begin, typename StringSet::Iterator end,
        uint16_t* bktout, size_t depth) const
    {
        while (begin != end)
        {
            key_type key = strset.get_uint64(*begin++, depth);
            *bktout++ = find_bkt_tree(key);
        }
    }

    //! return a splitter
    key_type get_splitter(unsigned int i) const
    { return splitter[i]; }

    /// build tree and splitter array from sample
    void build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        const size_t oversample_factor = samplesize / numsplitters;
        DBG(debug_splitter, "oversample_factor: " << oversample_factor);

        DBG(debug_splitter, "splitter:");
        splitter_lcp[0] = 0; // sentinel for first < everything bucket
        for (size_t i = 0, j = oversample_factor / 2; i < numsplitters; ++i)
        {
            splitter[i] = samples[j];
            DBG(debug_splitter, "key " << tlx::hexdump_type(splitter[i]));

            if (i != 0) {
                key_type xorSplit = splitter[i - 1] ^ splitter[i];

                DBG1(debug_splitter, "    XOR -> " << tlx::hexdump_type(xorSplit) << " - ");

                DBG3(debug_splitter, count_high_zero_bits(xorSplit) <<
                     " bits = " << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

                splitter_lcp[i] = count_high_zero_bits(xorSplit) / 8;
            }

            j += oversample_factor;
        }
    }
};

static inline void
parallel_sample_sortBSC(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyBinarySearch>(strings, n, 0);
}

PSS_CONTESTANT_PARALLEL_LCP(
    parallel_sample_sortBSC,
    "bingmann/parallel_sample_sortBSC",
    "pS5: binary search, bktcache")

////////////////////////////////////////////////////////////////////////////////

// Different Classifier variants are moved to other files
#include "bingmann-parallel_sample_sort_equal.hpp"
#include "bingmann-parallel_sample_sort_treecalc.hpp"

} // namespace bingmann_parallel_sample_sort(_lcp)

#endif // !PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_HEADER

/******************************************************************************/
