/******************************************************************************
 * src/parallel/bingmann-parallel_sample_sortBTCE.h
 *
 * Parallel Super Scalar String Sample-Sort with work-balancing, variant BTCE:
 * with binary splitter tree, equality checking and cache.
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

namespace bingmann_parallel_sample_sortBTCE {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = false;

static const bool debug_splitter = false;
static const bool debug_bucketsize = false;
static const bool debug_recursion = false;
static const bool debug_splitter_tree = false;

static const size_t MAXPROCS = 64+1; // +1 due to round up of processor number

static const size_t g_inssort_threshold = 64;

size_t g_totalsize;             // total size of input
size_t g_sequential_threshold;  // calculated threshold for sequential sorting
size_t g_threadnum;             // number of threads overall

// ****************************************************************************
// *** SampleSortStep out-of-place parallel sample sort with separate Jobs

typedef uint64_t key_type;

// *** rolled, and two unrolled classification subroutines

struct ClassifySimple
{
    /// binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt_tree(const key_type& key, const key_type* splitter_tree, size_t treebits, size_t numsplitters)
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
             const key_type* splitter_tree, size_t numsplitters,
             size_t treebits, size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key, splitter_tree, treebits, numsplitters);
            *bktout++ = b;
        }
    }
};

struct ClassifyUnrollTree
{
    /// binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt_tree(const key_type& key, const key_type* splitter_tree, const size_t treebits, const size_t numsplitters)
    {
        unsigned int i;

        assert(treebits == 11);

        // hand-coded assembler binary tree traversal with equality, using CMOV
        asm("mov    $1, %%rax \n"             // rax = i
            // level 1
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 2
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 3
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 4
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 5
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 6
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 7
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 8
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 9
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 10
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // level 11
            "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
            "je     2f \n"
            "lea    (%%rax,%%rax), %%rax \n"
            "lea    1(%%rax), %%rcx \n"
            "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
            // at leaf level
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

        return i;
    }

    /// classify all strings in area by walking tree and saving bucket id
    static inline void
    classify(string* strB, string* strE, uint16_t* bktout,
             const key_type* splitter_tree, const size_t numsplitters,
             const size_t treebits, size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            unsigned int b = find_bkt_tree(key, splitter_tree, treebits, numsplitters);
            *bktout++ = b;
        }
    }
};

template <typename Classify>
struct SampleSortStep
{
#if 0
    static const size_t numsplitters = 2*16;       // TODO: calculate to match L2 cache
#else
    static const size_t l2cache = 64*1024;

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
    key_type            splitter_tree[numsplitters+1];
    unsigned char       splitter_lcp[numsplitters];
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
    }

    // *** Sample Step

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

        TreeBuilder(splitter_tree, splitter_lcp, samples, samplesize);

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

        Classify::classify(strB, strE, bktout,
                           splitter_tree, numsplitters,
                           treebits, depth);

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
        if (n > g_sequential_threshold)
            new SampleSortStep(jobqueue, strptr, n, depth);
        else
            bingmann_parallel_radix_sort::EnqueueSmall(jobqueue, strptr, n, depth);
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

    SampleSortStep<ClassifySimple>::put_stats();

    string* shadow = new string[n]; // allocate shadow pointer array

    JobQueue jobqueue;
    SampleSortStep<ClassifySimple>::Enqueue(jobqueue, StringPtr(strings, shadow), n, 0);
    jobqueue.loop();

    delete [] shadow;
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTCE,
                             "bingmann/parallel_sample_sortBTCE",
                             "bingmann/parallel_sample_sortBTCE: binary tree, bktcache")

void parallel_sample_sortBTCEU1(string* strings, size_t n)
{
    g_totalsize = n;
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);

    SampleSortStep<ClassifyUnrollTree>::put_stats();

    string* shadow = new string[n]; // allocate shadow pointer array

    JobQueue jobqueue;
    SampleSortStep<ClassifyUnrollTree>::Enqueue(jobqueue, StringPtr(strings, shadow), n, 0);
    jobqueue.loop();

    delete [] shadow;
}

CONTESTANT_REGISTER_PARALLEL(parallel_sample_sortBTCEU1,
                             "bingmann/parallel_sample_sortBTCEU1",
                             "bingmann/parallel_sample_sortBTCEU1: binary tree, bktcache, unroll tree")

} // namespace bingmann_parallel_sample_sortBTC
