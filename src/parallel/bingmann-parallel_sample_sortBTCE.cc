/******************************************************************************
 * src/parallel/bingmann-parallel_sample_sortBTCE.h
 *
 * Parallel Super Scalar String Sample-Sort with work-balancing, variant BTCE:
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

#undef DBGX
#define DBGX DBGX_OMP

#include "../sequential/inssort.h"

extern size_t g_smallsort;

namespace bingmann_parallel_sample_sortBTCE {

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

size_t g_totalsize;             // total size of input
size_t g_sequential_threshold;  // calculated threshold for sequential sorting
size_t g_threadnum;             // number of threads overall

size_t g_para_ss_steps, g_seq_ss_steps, g_bs_steps;  // counters

typedef uint64_t key_type;

// ****************************************************************************
// *** Recursive TreeBuilder

template <size_t numsplitters>
struct TreeBuilder
{
    key_type*       m_tree;
    unsigned char*  m_lcp_iter;
    key_type*       m_samples;

    TreeBuilder(key_type* splitter_tree, unsigned char* splitter_lcp, key_type* samples, size_t samplesize)
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
struct ClassifySimple
{
    static const size_t numsplitters = (1 << treebits) - 1;

    /// binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt_tree(const key_type& key, const key_type* splitter_tree)
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
    static inline void
    classify(string* strB, string* strE, uint16_t* bktout,
             const key_type* splitter_tree, size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key, splitter_tree);
            *bktout++ = b;
        }
    }
};

template <size_t TreeBits>
struct ClassifyUnrollTree
{
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

template <>
struct ClassifyUnrollTree<11>
{
    static const size_t treebits = 11;
    static const size_t numsplitters = (1 << treebits) - 1;

    /// binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt_tree(const key_type& key, const key_type* splitter_tree)
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

    /// classify all strings in area by walking tree and saving bucket id
    static inline void
    classify(string* strB, string* strE, uint16_t* bktout,
             const key_type* splitter_tree, size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key, splitter_tree);
            *bktout++ = b;
        }
    }
};

template <>
struct ClassifyUnrollTree<12>
{
    static const size_t treebits = 12;
    static const size_t numsplitters = (1 << treebits) - 1;

    /// binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt_tree(const key_type& key, const key_type* splitter_tree)
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

    /// classify all strings in area by walking tree and saving bucket id
    static inline void
    classify(string* strB, string* strE, uint16_t* bktout,
             const key_type* splitter_tree, size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key, splitter_tree);
            *bktout++ = b;
        }
    }
};

template <>
struct ClassifyUnrollTree<13>
{
    static const size_t treebits = 13;
    static const size_t numsplitters = (1 << treebits) - 1;

    /// binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt_tree(const key_type& key, const key_type* splitter_tree)
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

    /// classify all strings in area by walking tree and saving bucket id
    static inline void
    classify(string* strB, string* strE, uint16_t* bktout,
             const key_type* splitter_tree, size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key, splitter_tree);
            *bktout++ = b;
        }
    }
};

#undef SPLITTER_TREE_STEP
#undef SPLITTER_TREE_END

// ****************************************************************************
// *** SampleSort non-recursive in-place sequential sample sort for small sorts

template <template <size_t> class Classify>
struct SmallsortJobBTCE : public Job
{
    StringPtr   strptr;
    size_t      n, depth;

    typedef uint32_t bktsize_type;

    SmallsortJobBTCE(JobQueue& jobqueue, const StringPtr& _strptr, size_t _n, size_t _depth)
        : strptr(_strptr), n(_n), depth(_depth)
    {
        jobqueue.enqueue(this);
    }

    struct SeqSampleSortStep
    {
#if 0
        static const size_t numsplitters = 2*16;
#else
        // bounding equations:
        // a) K * key_type splitter_tree size (and maybe K * key_type equal-cmp splitters)
        // b) 2 * K + 1 buckets (bktsize_type) when counting bkt occurances.
        static const size_t numsplitters2 = (l2cache - sizeof(size_t)) / (2 * sizeof(size_t));

        static const size_t treebits = logfloor_<numsplitters2>::value;
        static const size_t numsplitters = (1 << treebits) - 1;
#endif

        static const size_t bktnum = 2 * numsplitters + 1;

        string* str;
        size_t idx;
        size_t depth;
        bktsize_type bktsize[bktnum];

        unsigned char splitter_lcp[numsplitters+1];

        SeqSampleSortStep(string* strings, size_t n, size_t depth, uint16_t* bktcache)
        {
            // step 1: select splitters with oversampling

            const size_t oversample_factor = 1;
            const size_t samplesize = oversample_factor * numsplitters;

            key_type samples[ samplesize ];

            LCGRandom rng(&strings);

            for (unsigned int i = 0; i < samplesize; ++i)
            {
                samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
            }

            std::sort(samples, samples + samplesize);

            key_type splitter_tree[numsplitters+1];

            TreeBuilder<numsplitters>(splitter_tree, splitter_lcp, samples, samplesize);

            // step 2: classify all strings

            Classify<treebits>::classify(strings, strings+n, bktcache,
                                         splitter_tree, depth);

            // step 2.5: count bucket sizes

            memset(bktsize, 0, bktnum * sizeof(bktsize_type));

            for (size_t si = 0; si < n; ++si)
                ++bktsize[ bktcache[si] ];

            // step 3: prefix sum

            bktsize_type bktindex[bktnum];
            bktindex[0] = bktsize[0];
            bktsize_type last_bkt_size = bktsize[0];
            for (unsigned int i=1; i < bktnum; ++i) {
                bktindex[i] = bktindex[i-1] + bktsize[i];
                if (bktsize[i]) last_bkt_size = bktsize[i];
            }
            assert(bktindex[bktnum-1] == n);

            // step 4: premute in-place

            for (size_t i=0, j; i < n - last_bkt_size; )
            {
                string perm = strings[i];
                uint16_t permbkt = bktcache[i];

                while ( (j = --bktindex[ permbkt ]) > i )
                {
                    std::swap(perm, strings[j]);
                    std::swap(permbkt, bktcache[j]);
                }

                strings[i] = perm;
                i += bktsize[ permbkt ];
            }

            str = strings;
            idx = 0;
            this->depth = depth;

            ++g_seq_ss_steps;
        }
    };

    // *** Stack of Recursive Sample Sort Steps

    uint16_t* bktcache;
    size_t bktcache_size;

    size_t ss_pop_front;
    std::vector<SeqSampleSortStep> ss_stack;

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process SmallsortJobBTCE " << this << " of size " << n);

        ASSERT(n < (((uint64_t)1) << 32)); // must fit in uint32_t

        bktcache = NULL;
        bktcache_size = 0;
        ss_pop_front = 0;
        ms_pop_front = 0;

        if (n >= g_smallsort_threshold)
        {
            bktcache = new uint16_t[n];
            bktcache_size = n * sizeof(uint16_t);
            sort_sample_sort(jobqueue);
        }
        else
        {
            sort_mkqs_cache(jobqueue, strptr.to_original(n), n, depth);
        }
        delete [] bktcache;
    }

    void sort_sample_sort(JobQueue& jobqueue)
    {
        string* strings = strptr.to_original(n);

        typedef SeqSampleSortStep Step;

        // sort first level
        ss_pop_front = 0;
        ss_stack.push_back( Step(strings,n,depth,bktcache) );

        // step 5: "recursion"

        while ( ss_stack.size() > ss_pop_front )
        {
            while ( ss_stack.back().idx < Step::bktnum )
            {
                Step& s = ss_stack.back();
                size_t i = s.idx++; // process the bucket s.idx

                // i is even -> bkt[i] is less-than bucket
                if ((i & 1) == 0)
                {
                    if (s.bktsize[i] == 0)
                        ;
                    else if (s.bktsize[i] < g_smallsort_threshold)
                    {
                        assert(i/2 <= Step::numsplitters);
                        if (i == Step::bktnum-1)
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << s.bktsize[i] << " no lcp");
                        else
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << s.bktsize[i] << " lcp " << int(s.splitter_lcp[i/2] & 0x7F));

                        s.str += s.bktsize[i];
                        sort_mkqs_cache(jobqueue, s.str - s.bktsize[i], s.bktsize[i], s.depth + (s.splitter_lcp[i/2] & 0x7F));
                    }
                    else
                    {
                        if (i == Step::bktnum-1)
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << s.bktsize[i] << " no lcp");
                        else
                            DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << s.bktsize[i] << " lcp " << int(s.splitter_lcp[i/2] & 0x7F));

                        s.str += s.bktsize[i];
                        ss_stack.push_back( Step(s.str - s.bktsize[i], s.bktsize[i], s.depth + (s.splitter_lcp[i/2] & 0x7F), bktcache) );
                        // cannot add here, because s may have invalidated
                    }
                }
                // i is odd -> bkt[i] is equal bucket
                else
                {
                    if (s.bktsize[i] == 0)
                        ;
                    else if ( s.splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key, done.
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << s.bktsize[i] << " is done!");
                        s.str += s.bktsize[i];
                    }
                    else if (s.bktsize[i] < g_smallsort_threshold)
                    {
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << s.bktsize[i] << " lcp keydepth!");

                        s.str += s.bktsize[i];
                        sort_mkqs_cache(jobqueue, s.str - s.bktsize[i], s.bktsize[i], s.depth + sizeof(key_type));
                    }
                    else
                    {
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << s.bktsize[i] << " lcp keydepth!");

                        s.str += s.bktsize[i];
                        ss_stack.push_back( Step(s.str - s.bktsize[i], s.bktsize[i], s.depth + sizeof(key_type), bktcache) );
                        // cannot add here, because s may have invalidated
                    }
                }

                if (use_work_sharing && jobqueue.has_idle())
                    sample_sort_free_work(jobqueue);
            }

            ss_stack.pop_back();
        }
    }

    void sample_sort_free_work(JobQueue& jobqueue)
    {
        assert(ss_stack.size() >= ss_pop_front);

        if (ss_stack.size() == ss_pop_front) {
            // ss_stack is empty, check other stack
            //return radix_sort_free_work(jobqueue);
            return mkqs_free_work(jobqueue);
        }

        // convert top level of stack into independent jobs
        DBG(debug_jobs, "Freeing top level of SmallsortJobBTCE's sample_sort stack");

        typedef SeqSampleSortStep Step;
        Step& s = ss_stack[ss_pop_front];

        while ( s.idx < Step::bktnum )
        {
            size_t i = s.idx++; // process the bucket s.idx

            // i is even -> bkt[i] is less-than bucket
            if ((i & 1) == 0)
            {
                if (s.bktsize[i] == 0)
                    ;
                else
                {
                    if (i == Step::bktnum-1)
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << s.bktsize[i] << " no lcp");
                    else
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << s.bktsize[i] << " lcp " << int(s.splitter_lcp[i/2] & 0x7F));

                    new SmallsortJobBTCE<Classify>( jobqueue, StringPtr(s.str), s.bktsize[i], s.depth + (s.splitter_lcp[i/2] & 0x7F) );
                }
                s.str += s.bktsize[i];
            }
            // i is odd -> bkt[i] is equal bucket
            else
            {
                if (s.bktsize[i] == 0)
                    ;
                else if ( s.splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key, done.
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << s.bktsize[i] << " is done!");
                }
                else
                {
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << s.bktsize[i] << " lcp keydepth!");

                    new SmallsortJobBTCE<Classify>( jobqueue, StringPtr(s.str), s.bktsize[i], s.depth + sizeof(key_type) );
                }
                s.str += s.bktsize[i];
            }
        }

        // shorten the current stack
        ++ss_pop_front;
    }

#if 0
    // *********************************************************************************
    // *** Stack of Recursive Radix Sort Steps

    struct RadixStep8_CI
    {
        string* str;
        size_t idx;
        bktsize_type bktsize[256];

        RadixStep8_CI(string* strings, size_t n, size_t depth, uint8_t* charcache)
        {
            // count character occurances
#if 1
            // DONE: first variant to first fill charcache and then count. This is
            // 2x as fast as the second variant.
            for (size_t i=0; i < n; ++i)
                charcache[i] = strings[i][depth];

            memset(bktsize, 0, sizeof(bktsize));
            for (size_t i=0; i < n; ++i)
                ++bktsize[ charcache[i] ];
#else
            memset(bktsize, 0, sizeof(bktsize));
            for (size_t i=0; i < n; ++i) {
                ++bktsize[ (charcache[i] = strings[i][depth]) ];
            }
#endif
            // inclusive prefix sum
            size_t bkt[256];
            bkt[0] = bktsize[0];
            size_t last_bkt_size = bktsize[0];
            for (unsigned int i=1; i < 256; ++i) {
                bkt[i] = bkt[i-1] + bktsize[i];
                if (bktsize[i]) last_bkt_size = bktsize[i];
            }

            // premute in-place
            for (size_t i=0, j; i < n - last_bkt_size; )
            {
                string perm = strings[i];
                uint8_t permch = charcache[i];
                while ( (j = --bkt[ permch ]) > i )
                {
                    std::swap(perm, strings[j]);
                    std::swap(permch, charcache[j]);
                }
                strings[i] = perm;
                i += bktsize[ permch ];
            }

            str = strings + bktsize[0];
            idx = 0; // will increment to 1 on first process, bkt 0 is not sorted further

            ++g_bs_steps;
        }
    };

    size_t rs_pop_front;
    size_t rs_depth; // radix sorting base depth
    std::vector<RadixStep8_CI> radixstack;

    void sort_radix_sort(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
    {
        uint8_t* charcache = (uint8_t*)bktcache;

        // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
        rs_pop_front = 0;
        rs_depth = depth;
        radixstack.clear();
        radixstack.push_back( RadixStep8_CI(strings,n,rs_depth,charcache) );

        while ( radixstack.size() > rs_pop_front )
        {
            while ( radixstack.back().idx < 255 )
            {
                RadixStep8_CI& rs = radixstack.back();
                size_t i = ++rs.idx; // process the bucket rs.idx

                if (rs.bktsize[i] == 0)
                    continue;
                else if (rs.bktsize[i] < g_inssort_threshold)
                {
                    inssort::inssort(rs.str, rs.bktsize[i], depth + radixstack.size());
                    rs.str += rs.bktsize[i];
                }
                else
                {
                    rs.str += rs.bktsize[i];
                    radixstack.push_back( RadixStep8_CI(rs.str - rs.bktsize[i],
                                                        rs.bktsize[i],
                                                        depth + radixstack.size(),
                                                        charcache) );
                    // cannot add to rs.str here, because rs may have invalidated
                }

                if (use_work_sharing && jobqueue.has_idle())
                    sample_sort_free_work(jobqueue);
            }
            radixstack.pop_back();
        }

        rs_pop_front = 0; // clear stack
        radixstack.clear();
    }

    void radix_sort_free_work(JobQueue& jobqueue)
    {
        assert(radixstack.size() >= rs_pop_front);

        if (radixstack.size() == rs_pop_front) {
            return;
        }

        // convert top level of stack into independent jobs
        DBG(debug_jobs, "Freeing top level of SmallsortJobBTCE's radixsort stack");

        RadixStep8_CI& rt = radixstack[rs_pop_front];

        while ( rt.idx < 255 )
        {
            ++rt.idx; // enqueue the bucket rt.idx

            if (rt.bktsize[rt.idx] == 0) continue;
            new SmallsortJobBTCE<Classify>(jobqueue, StringPtr(rt.str), rt.bktsize[rt.idx], rs_depth + rs_pop_front);
            rt.str += rt.bktsize[rt.idx];
        }

        // shorten the current stack
        ++rs_pop_front;
    }
#endif

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
        if (cmp(A[i], A[j]) == 0)                         return i;
        if (cmp(A[k], A[i]) == 0 || cmp(A[k], A[j]) == 0) return k;
        if (cmp(A[i], A[j]) < 0) {
            if (cmp(A[j], A[k]) < 0) return j;
            if (cmp(A[i], A[k]) < 0) return k;
            return i;
        }
        if (cmp(A[j], A[k]) > 0) return j;
        if (cmp(A[i], A[k]) < 0) return i;
        return k;
    }

    // Insertion sort the strings only based on the cached characters.
    static inline void
    insertion_sort_cache_block(string* strings, key_type* cache, int n)
    {
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

    template <bool CacheDirty>
    static inline void
    insertion_sort(string* strings, key_type* cache, int n, size_t depth)
    {
        if (n == 0) return;
        if (CacheDirty) return inssort::inssort(strings, n, depth);

        insertion_sort_cache_block(strings, cache, n);

        size_t start=0, cnt=1;
        for (int i=0; i < n-1; ++i) {
            if (cmp(cache[i], cache[i+1]) == 0) {
                ++cnt;
                continue;
            }
            if (cnt > 1 && cache[start] & 0xFF)
                inssort::inssort(strings + start, cnt, depth + sizeof(key_type));
            cnt = 1;
            start = i+1;
        }
        if (cnt > 1 && cache[start] & 0xFF)
            inssort::inssort(strings+start, cnt, depth + sizeof(key_type));
    }

    struct MKQSStep
    {
        string* strings;
        key_type* cache;
        size_t num_lt, num_eq, num_gt, n, depth;
        size_t idx;
        bool eq_recurse;

        MKQSStep(string* strings, key_type* cache, size_t n, size_t depth, bool CacheDirty)
        {
            if (CacheDirty) {
                for (size_t i=0; i < n; ++i) {
                    cache[i] = get_char<key_type>(strings[i], depth);
                }
            }
            // Move pivot to first position to avoid wrapping the unsigned values
            // we are using in the main loop from zero to max.
            size_t p = med3(cache,
                            med3(cache, 0,       n/8,     n/4),
                            med3(cache, n/2-n/8, n/2,     n/2+n/8),
                            med3(cache, n-1-n/4, n-1-n/8, n-3)
                );
            std::swap(strings[0], strings[p]);
            std::swap(cache[0], cache[p]);
            key_type pivot = cache[0];
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
                        std::swap(strings[beg_ins], strings[first]);
                        std::swap(cache[beg_ins], cache[first]);
                        beg_ins++;
                    }
                    ++first;
                }
                while (first <= last) {
                    const int res = cmp(cache[last], pivot);
                    if (res < 0) {
                        break;
                    } else if (res == 0) {
                        std::swap(strings[end_ins], strings[last]);
                        std::swap(cache[end_ins], cache[last]);
                        end_ins--;
                    }
                    --last;
                }
                if (first > last)
                    break;
                std::swap(strings[first], strings[last]);
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
            assert(num_lt + num_eq + num_gt == n);

            // Swap the equal pointers from the beginning to proper position.
            const size_t size1 = std::min(num_eq_beg, num_lt);
            std::swap_ranges(strings, strings+size1, strings+first-size1);
            std::swap_ranges(cache, cache+size1, cache+first-size1);
            // Swap the equal pointers from the end to proper position.
            const size_t size2 = std::min(num_eq_end, num_gt);
            std::swap_ranges(strings+first, strings+first+size2, strings+n-size2);
            std::swap_ranges(cache+first, cache+first+size2, cache+n-size2);

            // Save offsets for recursive sorting
            this->strings = strings;
            this->cache = cache;
            this->n = n;
            this->depth = depth;
            this->idx = 0;
            this->eq_recurse = (pivot & 0xFF);

            ++g_bs_steps;
        }
    };

    size_t ms_pop_front;
    size_t ms_depth; // mkqs sorting base depth
    std::vector<MKQSStep> ms_stack;

    void sort_mkqs_cache(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
    {
        if (n < g_inssort_threshold)
            return inssort::inssort(strings, n, depth);

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
        ms_stack.push_back( MKQSStep(strings,cache,n,ms_depth,true) );

    jumpout:
        while ( ms_stack.size() > ms_pop_front )
        {
            while ( ms_stack.back().idx < 3 )
            {
                if (use_work_sharing && jobqueue.has_idle()) {
                    sample_sort_free_work(jobqueue);
                    goto jumpout;
                }

                MKQSStep& ms = ms_stack.back();
                ++ms.idx; // increment here, because stack may change

                // process the lt-subsequence
                if (ms.idx == 1)
                {
                    if (ms.num_lt == 0)
                        ;
                    else if (ms.num_lt < g_inssort_threshold)
                        insertion_sort<false>(ms.strings, ms.cache, ms.num_lt, ms.depth);
                    else
                        ms_stack.push_back( MKQSStep(ms.strings, ms.cache, ms.num_lt, ms.depth, false) );
                }
                // process the eq-subsequence
                else if (ms.idx == 2)
                {
                    if (!ms.eq_recurse || ms.num_eq == 0)
                        ;
                    else if (ms.num_eq < g_inssort_threshold)
                        insertion_sort<true>(ms.strings + ms.num_lt,
                                             ms.cache + ms.num_lt,
                                             ms.num_eq, ms.depth + sizeof(key_type));
                    else
                        ms_stack.push_back( MKQSStep(ms.strings + ms.num_lt,
                                                     ms.cache + ms.num_lt,
                                                     ms.num_eq, ms.depth + sizeof(key_type), true) );
                }
                // process the gt-subsequence
                else
                {
                    assert(ms.idx == 3);
                    if (ms.num_gt == 0)
                        ;
                    else if (ms.num_gt < g_inssort_threshold)
                        insertion_sort<false>(ms.strings + ms.num_lt + ms.num_eq,
                                              ms.cache + ms.num_lt + ms.num_eq,
                                              ms.num_gt,
                                              ms.depth);
                    else
                        ms_stack.push_back( MKQSStep(ms.strings + ms.num_lt + ms.num_eq,
                                                     ms.cache + ms.num_lt + ms.num_eq,
                                                     ms.num_gt, ms.depth, false) );
                }
            }

            ms_stack.pop_back();
        }

        ms_pop_front = 0; // clear stack
        ms_stack.clear();
    }

    void mkqs_free_work(JobQueue& jobqueue)
    {
        assert(ms_stack.size() >= ms_pop_front);

        if (ms_stack.size() == ms_pop_front) {
            return;
        }

        DBG(debug_jobs, "Freeing top level of SmallsortJobBTCE's mkqs stack");

        // convert top level of stack into independent jobs

        MKQSStep& st = ms_stack[ms_pop_front];

        if (st.idx == 0 && st.num_lt != 0)
        {
            new SmallsortJobBTCE(jobqueue, StringPtr(st.strings), st.num_lt, st.depth);
        }
        if (st.idx <= 1 && st.num_eq != 0 && st.eq_recurse)
        {
            new SmallsortJobBTCE(jobqueue, StringPtr(st.strings + st.num_lt),
                                st.num_eq, st.depth + sizeof(key_type));
        }
        if (st.idx <= 2 && st.num_gt != 0)
        {
            new SmallsortJobBTCE(jobqueue, StringPtr(st.strings + st.num_lt + st.num_eq),
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
    static const size_t numsplitters = 2*16;       // TODO: calculate to match L2 cache
#else
    // bounding equation: 2*K+1 buckets when counting bkt occurances.
    static const size_t numsplitters2 = (l2cache - sizeof(size_t)) / (2 * sizeof(size_t));

    static const size_t treebits = logfloor_<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;
#endif

    static const size_t bktnum = 2 * numsplitters + 1;

    StringPtr           strptr;
    size_t              n, depth;

    unsigned int        parts;
    size_t              psize;
    std::atomic<unsigned int> pwork;
    key_type            splitter_tree[numsplitters+1];
    unsigned char       splitter_lcp[numsplitters+1];
    size_t*             bkt;

    uint16_t*           bktcache[MAXPROCS];

    // *** Classes for JobQueue

    struct SampleJob : public Job
    {
        SampleSortStep*     step;

        SampleJob(SampleSortStep* _step)
            : step(_step) { }

        virtual void run(JobQueue& jobqueue)
        {
            step->sample(jobqueue);
        }
    };

    struct CountJob : public Job
    {
        SampleSortStep*     step;
        unsigned int        p;

        CountJob(SampleSortStep* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual void run(JobQueue& jobqueue)
        {
            step->count(p, jobqueue);
        }
    };

    struct DistributeJob : public Job
    {
        SampleSortStep*     step;
        unsigned int        p;

        DistributeJob(SampleSortStep* _step, unsigned int _p)
            : step(_step), p(_p) { }

        virtual void run(JobQueue& jobqueue)
        {
            step->distribute(p, jobqueue);
        }
    };

    // *** Constructor

    SampleSortStep(JobQueue& jobqueue, const StringPtr& _strptr, size_t _n, size_t _depth)
        : strptr(_strptr), n(_n), depth(_depth)
    {
        parts = (n + g_sequential_threshold-1) / g_sequential_threshold;
        if (parts == 0) parts = 1;

        psize = (n + parts-1) / parts;

        jobqueue.enqueue( new SampleJob(this) );
        ++g_para_ss_steps;
    }

    // *** Sample Step

    void sample(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process SampleJob @ " << this);

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

        TreeBuilder<numsplitters>(splitter_tree, splitter_lcp, samples, samplesize);

        bkt = new size_t[bktnum * parts + 1];

        // create new jobs
        pwork = parts;
        for (unsigned int p = 0; p < parts; ++p)
            jobqueue.enqueue( new CountJob(this, p) );
    }

    // *** Counting Step

    void count(unsigned int p, JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process CountJob " << p << " @ " << this);

        string* strB = strptr.active() + p * psize;
        string* strE = strptr.active() + std::min((p+1) * psize, n);
        if (strE < strB) strE = strB;

        uint16_t* mybktcache = bktcache[p] = new uint16_t[strE-strB];
        uint16_t* bktout = mybktcache;

        Classify<treebits>::classify(strB, strE, bktout,
                                     splitter_tree, depth);

        size_t* mybkt = bkt + p * bktnum;
        memset(mybkt, 0, bktnum * sizeof(size_t));

        for (uint16_t* bc = mybktcache; bc != mybktcache + (strE-strB); ++bc)
            ++mybkt[ *bc ];

        if (--pwork == 0)
            count_finished(jobqueue);
    }

    void count_finished(JobQueue& jobqueue)
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

    // *** Distrbute Step

    void distribute(unsigned int p, JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process DistributeJob " << p << " @ " << this);

        string* strB = strptr.active() + p * psize;
        string* strE = strptr.active() + std::min((p+1) * psize, n);
        if (strE < strB) strE = strB;

        string* sorted = strptr.shadow(); // get alternative shadow pointer array

        uint16_t* mybktcache = bktcache[p];
        size_t mybkt[bktnum]; // copy bkt array to local stack
        memcpy(mybkt, bkt + p * bktnum, sizeof(mybkt));

        for (string* str = strB; str != strE; ++str, ++mybktcache)
            sorted[ --mybkt[ *mybktcache ] ] = *str;

        if (p == 0) // these are needed for recursion into bkts
            memcpy(bkt, mybkt, sizeof(mybkt));

        delete [] bktcache[p];

        if (--pwork == 0)
            distribute_finished(jobqueue);
    }

    void distribute_finished(JobQueue& jobqueue)
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
                DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bkt[i] << " size " << bktsize << " lcp " << int(splitter_lcp[i/2] & 0x7F));
                Enqueue(jobqueue, strptr.flip_ptr(bkt[i]), bktsize, depth + (splitter_lcp[i/2] & 0x7F));
            }
            ++i;
            // i is odd -> bkt[i] is equal bucket
            bktsize = bkt[i+1] - bkt[i];
            if (bktsize == 0)
                ;
            else if (bktsize == 1) // just one string pointer, copyback
                strptr.flip_ptr(bkt[i]).to_original(1);
            else
            {
                if ( splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key
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

    // *** Enqueue Deeper Sort Jobs

    static void Enqueue(JobQueue& jobqueue, const StringPtr& strptr, size_t n, size_t depth)
    {
        if (n > g_sequential_threshold) {
            new SampleSortStep(jobqueue, strptr, n, depth);
        }
        else {
            new SmallsortJobBTCE<Classify>(jobqueue, strptr, n, depth);
        }
   }

    static inline void put_stats()
    {
        g_statscache >> "l2cache" << size_t(l2cache)
                     >> "splitter_treebits" << size_t(treebits)
                     >> "numsplitters" << size_t(numsplitters);
    }
};

void parallel_sample_sortBTCE(string* strings, size_t n)
{
    g_totalsize = n;
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);
    g_para_ss_steps = g_seq_ss_steps = g_bs_steps = 0;

    SampleSortStep<ClassifySimple>::put_stats();

    string* shadow = new string[n]; // allocate shadow pointer array

    JobQueue jobqueue;
    SampleSortStep<ClassifySimple>::Enqueue(jobqueue, StringPtr(strings, shadow), n, 0);
    jobqueue.loop();

    delete [] shadow;

    g_statscache >> "steps_para_sample_sort" << g_para_ss_steps
                 >> "steps_seq_sample_sort" << g_seq_ss_steps
                 >> "steps_base_sort" << g_bs_steps;
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTCE,
                             "bingmann/parallel_sample_sortBTCE",
                             "bingmann/parallel_sample_sortBTCE: binary tree, equality, bktcache")

void parallel_sample_sortBTCEU1(string* strings, size_t n)
{
    g_totalsize = n;
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);
    g_para_ss_steps = g_seq_ss_steps = g_bs_steps = 0;

    SampleSortStep<ClassifyUnrollTree>::put_stats();

    string* shadow = new string[n]; // allocate shadow pointer array

    JobQueue jobqueue;
    SampleSortStep<ClassifyUnrollTree>::Enqueue(jobqueue, StringPtr(strings, shadow), n, 0);
    jobqueue.loop();

    delete [] shadow;

    g_statscache >> "steps_para_sample_sort" << g_para_ss_steps
                 >> "steps_seq_sample_sort" << g_seq_ss_steps
                 >> "steps_base_sort" << g_bs_steps;
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTCEU1,
                             "bingmann/parallel_sample_sortBTCEU1",
                             "bingmann/parallel_sample_sortBTCEU1: binary tree, equality, bktcache, unroll tree")

} // namespace bingmann_parallel_sample_sortBTCE
