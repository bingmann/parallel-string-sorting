/*******************************************************************************
 * src/parallel/bingmann-parallel_radix_sort.cc
 *
 * Parallel radix sort with work-balancing.
 *
 * The set of strings is sorted using a 8- or 16-bit radix sort
 * algorithm. Recursive sorts are processed in parallel using a lock-free job
 * queue and OpenMP threads. Two radix sort implementations are used:
 * sequential in-place and parallelized out-of-place.
 *
 * The sequential radix sort is implemented using an explicit recursion stack,
 * which enables threads to "free up" work from the top of the stack when other
 * threads become idle. This variant uses in-place permuting and a character
 * cache (oracle).
 *
 * To parallelize sorting of large buckets an out-of-place variant is
 * implemented using a sequence of three Jobs: count, distribute and
 * copyback. All threads work on a job-specific context called RadixStepCE,
 * which encapsules all variables of an 8-bit or 16-bit radix sort (templatized
 * with key_type).
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

#include <cstdlib>
#include <cstring>

#include <iostream>
#include <vector>

#include "../tools/debug.hpp"
#include "../tools/contest.hpp"
#include "../tools/stringtools.hpp"
#include "../tools/jobqueue.hpp"

#undef DBGX
#define DBGX DBGX_OMP

#include "../sequential/inssort.hpp"

extern size_t g_smallsort;

namespace bingmann_parallel_radix_sort {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = false;

static const bool use_work_sharing = true;

typedef unsigned char* string;
typedef StringShadowPtr StringPtr;

static const size_t g_inssort_threshold = 64;

size_t g_totalsize;             // total size of input
size_t g_sequential_threshold;  // calculated threshold for sequential sorting
size_t g_threadnum;             // number of threads overall

/// Prototype called to schedule deeper sorts
template <typename bigsort_key_type>
void Enqueue(JobQueue& jobqueue, const StringPtr& strings, size_t depth);

// ****************************************************************************
// *** SmallsortJob8 - sort 8-bit radix in-place with explicit stack-based recursion

void EnqueueSmallsortJob8(JobQueue& jobqueue, const StringPtr& strptr, size_t depth);

template <typename BktSizeType>
struct SmallsortJob8 : public Job
{
    StringPtr strptr;
    size_t    depth;

    typedef BktSizeType bktsize_type;

    SmallsortJob8(JobQueue& jobqueue, const StringPtr& _strptr, size_t _depth)
        : strptr(_strptr), depth(_depth)
    {
        jobqueue.enqueue(this);
    }

    struct RadixStep8_CI
    {
        StringPtr    strptr;
        size_t       idx;
        bktsize_type bkt[256 + 1];

        RadixStep8_CI(const StringPtr& _strptr, size_t depth, uint8_t* charcache)
            : strptr(_strptr)
        {
            string* strings = strptr.active();
            size_t n = strptr.size();

            // count character occurances
            bktsize_type bktsize[256];
#if 1
            // DONE: first variant to first fill charcache and then count. This is
            // 2x as fast as the second variant.
            for (size_t i = 0; i < n; ++i)
                charcache[i] = strings[i][depth];

            memset(bktsize, 0, sizeof(bktsize));
            for (size_t i = 0; i < n; ++i)
                ++bktsize[charcache[i]];
#else
            memset(bktsize, 0, sizeof(bktsize));
            for (size_t i = 0; i < n; ++i) {
                ++bktsize[(charcache[i] = strings[i][depth])];
            }
#endif
            // inclusive prefix sum
            bkt[0] = bktsize[0];
            bktsize_type last_bkt_size = bktsize[0];
            for (unsigned int i = 1; i < 256; ++i) {
                bkt[i] = bkt[i - 1] + bktsize[i];
                if (bktsize[i]) last_bkt_size = bktsize[i];
            }

            // premute in-place
            for (size_t i = 0, j; i < n - last_bkt_size; )
            {
                string perm = strings[i];
                uint8_t permch = charcache[i];
                while ((j = --bkt[permch]) > i)
                {
                    std::swap(perm, strings[j]);
                    std::swap(permch, charcache[j]);
                }
                strings[i] = perm;
                i += bktsize[permch];
            }

            // fix prefix sum
            bkt[0] = 0;
            for (unsigned int i = 1; i <= 256; ++i) {
                bkt[i] = bkt[i - 1] + bktsize[i - 1];
            }
            assert(bkt[256] == n);

            idx = 0; // will increment to 1 on first process, bkt 0 is not sorted further
        }
    };

    virtual bool run(JobQueue& jobqueue)
    {
        size_t n = strptr.size();

        DBG(debug_jobs, "Process SmallsortJob8 " << this << " of size " << n);

        strptr = strptr.copy_back();

        if (n < g_inssort_threshold) {
            inssort::inssort(strptr.active(), n, depth);
            return true;
        }

        uint8_t* charcache = new uint8_t[n];

        // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
        size_t pop_front = 0;
        std::vector<RadixStep8_CI> radixstack;
        radixstack.push_back(RadixStep8_CI(strptr, depth, charcache));

        while (radixstack.size() > pop_front)
        {
            while (radixstack.back().idx < 255)
            {
                RadixStep8_CI& rs = radixstack.back();
                size_t b = ++rs.idx; // process the bucket rs.idx

                size_t bktsize = rs.bkt[b + 1] - rs.bkt[b];

                if (bktsize == 0)
                    continue;
                else if (bktsize < g_inssort_threshold)
                {
                    inssort::inssort(rs.strptr.active() + rs.bkt[b], bktsize,
                                     depth + radixstack.size());
                }
                else
                {
                    radixstack.push_back(RadixStep8_CI(rs.strptr.sub(rs.bkt[b], bktsize),
                                                       depth + radixstack.size(),
                                                       charcache));
                }

                if (use_work_sharing && jobqueue.has_idle())
                {
                    // convert top level of stack into independent jobs
                    DBG(debug_jobs, "Freeing top level of SmallsortJob8's radixsort stack");

                    RadixStep8_CI& rt = radixstack[pop_front];

                    while (rt.idx < 255)
                    {
                        b = ++rt.idx; // enqueue the bucket rt.idx

                        size_t bktsize = rt.bkt[b + 1] - rt.bkt[b];

                        if (bktsize == 0) continue;
                        EnqueueSmallsortJob8(jobqueue, rt.strptr.sub(rt.bkt[b], bktsize),
                                             depth + pop_front);
                    }

                    // shorten the current stack
                    ++pop_front;
                }
            }
            radixstack.pop_back();
        }

        delete[] charcache;

        return true;
    }
};

void EnqueueSmallsortJob8(JobQueue& jobqueue, const StringPtr& strptr, size_t depth)
{
    if (strptr.size() < ((uint64_t)1 << 32))
        new SmallsortJob8<uint32_t>(jobqueue, strptr, depth);
    else
        new SmallsortJob8<uint64_t>(jobqueue, strptr, depth);
}

// ****************************************************************************
// *** SmallsortJob16 - sort 16-bit radix in-place with explicit stack-based recursion

void EnqueueSmallsortJob16(JobQueue& jobqueue, const StringPtr& strptr, size_t depth);

template <typename BktSizeType>
struct SmallsortJob16 : public Job
{
    StringPtr strptr;
    size_t    depth;

    typedef BktSizeType bktsize_type;

    SmallsortJob16(JobQueue& jobqueue, const StringPtr& _strptr, size_t _depth)
        : strptr(_strptr), depth(_depth)
    {
        jobqueue.enqueue(this);
    }

    struct RadixStep16_CI
    {
        typedef uint16_t key_type;
        static const size_t numbkts = key_traits<key_type>::radix;

        StringPtr           strptr;
        size_t              idx;
        bktsize_type        bkt[numbkts + 1];

        RadixStep16_CI(const StringPtr& _strptr, size_t depth, key_type* charcache)
            : strptr(_strptr)
        {
            string* strings = strptr.active();
            size_t n = strptr.size();

            // fill character cache
            for (size_t i = 0; i < n; ++i)
                charcache[i] = get_char<key_type>(strings[i], depth);

            // count character occurances
            bktsize_type bktsize[numbkts];
            memset(bktsize, 0, sizeof(bktsize));
            for (size_t i = 0; i < n; ++i)
                ++bktsize[charcache[i]];

            // inclusive prefix sum
            bkt[0] = bktsize[0];
            bktsize_type last_bkt_size = bktsize[0];
            for (unsigned int i = 1; i < numbkts; ++i) {
                bkt[i] = bkt[i - 1] + bktsize[i];
                if (bktsize[i]) last_bkt_size = bktsize[i];
            }

            // premute in-place
            for (size_t i = 0, j; i < n - last_bkt_size; )
            {
                string perm = strings[i];
                key_type permch = charcache[i];
                while ((j = --bkt[permch]) > i)
                {
                    std::swap(perm, strings[j]);
                    std::swap(permch, charcache[j]);
                }
                strings[i] = perm;
                i += bktsize[permch];
            }

            // fix prefix sum
            bkt[0] = 0;
            for (unsigned int i = 1; i <= numbkts; ++i) {
                bkt[i] = bkt[i - 1] + bktsize[i - 1];
            }
            assert(bkt[numbkts] == n);

            idx = 0x0100; // will increment to 0x0101 on first process, bkt 0 is not sorted further
        }
    };

    virtual bool run(JobQueue& jobqueue)
    {
        size_t n = strptr.size();

        DBG(debug_jobs, "Process SmallsortJob16 " << this << " of size " << n);

        strptr = strptr.copy_back();

        if (n < g_inssort_threshold) {
            inssort::inssort(strptr.active(), n, depth);
            return true;
        }

        typedef uint16_t key_type;

        static const size_t numbkts = key_traits<key_type>::radix;

        key_type* charcache = new key_type[n];

        // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
        size_t pop_front = 0;
        std::vector<RadixStep16_CI> radixstack;
        radixstack.push_back(RadixStep16_CI(strptr, depth, charcache));

        while (radixstack.size() > pop_front)
        {
            while (radixstack.back().idx < numbkts - 1)
            {
                RadixStep16_CI& rs = radixstack.back();
                size_t b = ++rs.idx;   // process the bucket rs.idx

                size_t bktsize = rs.bkt[b + 1] - rs.bkt[b];

                if ((b & 0xFF) == 0) { // skip over finished 0x??00 buckets
                    continue;
                }

                if (bktsize <= 1)
                    continue;
                else if (bktsize < g_inssort_threshold)
                {
                    inssort::inssort(rs.strptr.active() + rs.bkt[b], bktsize,
                                     depth + 2* radixstack.size());
                }
                else
                {
                    radixstack.push_back(RadixStep16_CI(rs.strptr.sub(rs.bkt[b], bktsize),
                                                        depth + 2 * radixstack.size(),
                                                        charcache));
                }

                if (use_work_sharing && jobqueue.has_idle())
                {
                    // convert top level of stack into independent jobs
                    DBG(debug_jobs, "Freeing top level of SmallsortJob16's radixsort stack");

                    RadixStep16_CI& rt = radixstack[pop_front];

                    while (rt.idx < numbkts - 1)
                    {
                        b = ++rt.idx; // enqueue the bucket rt.idx

                        size_t bktsize = rt.bkt[b + 1] - rt.bkt[b];

                        if (bktsize == 0) continue;
                        EnqueueSmallsortJob16(jobqueue, rt.strptr.sub(rt.bkt[b], bktsize),
                                              depth + 2 * pop_front);
                    }

                    // shorten the current stack
                    ++pop_front;
                }
            }
            radixstack.pop_back();
        }

        delete[] charcache;

        return true;
    }
};

void EnqueueSmallsortJob16(JobQueue& jobqueue, const StringPtr& strptr, size_t depth)
{
    if (strptr.size() < ((uint64_t)1 << 32))
        new SmallsortJob16<uint32_t>(jobqueue, strptr, depth);
    else
        new SmallsortJob16<uint64_t>(jobqueue, strptr, depth);
}

// ****************************************************************************
// *** RadixStepCE out-of-place 8- or 16-bit parallel radix sort with Jobs

template <typename KeyType>
struct RadixStepCE
{
    typedef KeyType key_type;

    static const size_t       numbkts = key_traits<key_type>::radix; // 256 or 65536

    StringPtr                 strptr;
    size_t                    depth;

    unsigned int              parts;
    size_t                    psize;
    std::atomic<unsigned int> pwork;
    size_t                    * bkt;

    key_type                  * charcache;

    RadixStepCE(JobQueue& jobqueue, const StringPtr& _strptr, size_t _depth);

    void                      count(unsigned int p, JobQueue& jobqueue);
    void                      count_finished(JobQueue& jobqueue);

    void                      distribute(unsigned int p, JobQueue& jobqueue);
    void                      distribute_finished(JobQueue& jobqueue);
};

template <typename key_type>
struct CountJob : public Job
{
    RadixStepCE<key_type>* step;
    unsigned int         p;

    CountJob(RadixStepCE<key_type>* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual bool run(JobQueue& jobqueue)
    {
        step->count(p, jobqueue);
        return true;
    }
};

template <typename key_type>
struct DistributeJob : public Job
{
    RadixStepCE<key_type>* step;
    unsigned int         p;

    DistributeJob(RadixStepCE<key_type>* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual bool run(JobQueue& jobqueue)
    {
        step->distribute(p, jobqueue);
        return true;
    }
};

template <typename key_type>
RadixStepCE<key_type>::RadixStepCE(JobQueue& jobqueue, const StringPtr& _strptr, size_t _depth)
    : strptr(_strptr), depth(_depth)
{
    size_t n = strptr.size();

    parts = (n + g_sequential_threshold - 1) / g_sequential_threshold;
    if (parts == 0) parts = 1;

    psize = (n + parts - 1) / parts;

    DBG(debug_jobs, "Area split into " << parts << " parts of size " << psize);

    bkt = new size_t[numbkts * parts + 1];
    charcache = new key_type[n];

    // create worker jobs
    pwork = parts;
    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue(new CountJob<key_type>(this, p));
}

template <typename key_type>
void RadixStepCE<key_type>::count(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process CountJob " << p << " @ " << this);

    string* strB = strptr.active() + p * psize;
    string* strE = strptr.active() + std::min((p + 1) * psize, strptr.size());
    if (strE < strB) strE = strB;

    key_type* mycache = charcache + p * psize;
    key_type* mycacheE = mycache + (strE - strB);

    // DONE: check if processor-local stack + copy is faster. On 48-core AMD
    // Opteron it is not faster to copy first. On 32-core Intel Xeon the second
    // loop is slightly faster.
#if 0
    for (string* str = strB; str != strE; ++str, ++mycache)
        *mycache = get_char<key_type>(*str, depth);

    size_t* mybkt = bkt + p * numbkts;
    memset(mybkt, 0, numbkts * sizeof(size_t));
    for (mycache = charcache + p * psize; mycache != mycacheE; ++mycache)
        ++mybkt[*mycache];
#else
    for (string* str = strB; str != strE; ++str, ++mycache)
        *mycache = get_char<key_type>(*str, depth);

    size_t mybkt[numbkts] = { 0 };
    for (mycache = charcache + p * psize; mycache != mycacheE; ++mycache)
        ++mybkt[*mycache];
    memcpy(bkt + p * numbkts, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        count_finished(jobqueue);
}

template <typename key_type>
void RadixStepCE<key_type>::count_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing CountJob " << this << " with prefixsum");

    // inclusive prefix sum over bkt
    size_t sum = 0;
    for (unsigned int i = 0; i < numbkts; ++i)
    {
        for (unsigned int p = 0; p < parts; ++p)
        {
            bkt[p * numbkts + i] = (sum += bkt[p * numbkts + i]);
        }
    }
    assert(sum == strptr.size());

    // create new jobs
    pwork = parts;
    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue(new DistributeJob<key_type>(this, p));
}

template <typename key_type>
void RadixStepCE<key_type>::distribute(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process DistributeJob " << p << " @ " << this);

    string* strB = strptr.active() + p * psize;
    string* strE = strptr.active() + std::min((p + 1) * psize, strptr.size());
    if (strE < strB) strE = strB;

    string* sorted = strptr.shadow(); // get alternative shadow pointer array
    key_type* mycache = charcache + p * psize;

    // DONE: check if processor-local stack + copy is faster. On 48-core AMD
    // Opteron it is not faster to copy first. On 32-core Intel Xeon the second
    // loop is slightly faster.
#if 0
    size_t* mybkt = bkt + p * numbkts;

    for (string* str = strB; str != strE; ++str, ++mycache)
        sorted[--mybkt[*mycache]] = *str;
#else
    size_t mybkt[numbkts];
    memcpy(mybkt, bkt + p * numbkts, sizeof(mybkt));

    for (string* str = strB; str != strE; ++str, ++mycache)
        sorted[--mybkt[*mycache]] = *str;

    if (p == 0) // these are needed for recursion into bkts
        memcpy(bkt, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        distribute_finished(jobqueue);
}

template <>
void RadixStepCE<uint8_t>::distribute_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing DistributeJob " << this << " with enqueuing subjobs");

    delete[] charcache;

    // first p's bkt pointers are boundaries between bkts, just add sentinel:
    assert(bkt[0] == 0);
    bkt[numbkts] = strptr.size();

    for (unsigned int i = 1; i < numbkts; ++i)
    {
        if (bkt[i] == bkt[i + 1])
            continue;
        else if (bkt[i] + 1 == bkt[i + 1]) // just one string pointer, copyback
            strptr.flip(bkt[i], 1).copy_back();
        else
            Enqueue<key_type>(jobqueue, strptr.flip(bkt[i], bkt[i + 1] - bkt[i]),
                              depth + key_traits<key_type>::add_depth);
    }

    // maybe copy back finished pointers from shadow
    strptr.flip(0, bkt[1] - bkt[0]).copy_back();

    delete[] bkt;
    delete this;
}

template <>
void RadixStepCE<uint16_t>::distribute_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing DistributeJob " << this << " with enqueuing subjobs");

    delete[] charcache;

    // first p's bkt pointers are boundaries between bkts, just add sentinel:
    assert(bkt[0] == 0);
    bkt[numbkts] = strptr.size();

    for (unsigned int i = 0x0101; i < numbkts; ++i)
    {
        if ((i & 0xFF) == 0) ++i;          // skip over finished buckets 0x??00

        if (bkt[i] == bkt[i + 1])
            continue;
        else if (bkt[i] + 1 == bkt[i + 1]) // just one string pointer, copyback
            strptr.flip(bkt[i], 1).copy_back();
        else
            Enqueue<key_type>(jobqueue, strptr.flip(bkt[i], bkt[i + 1] - bkt[i]),
                              depth + key_traits<key_type>::add_depth);
    }

    // maybe copy back finished pointers from shadow
    if (!strptr.flipped()) // if sorted[] = shadow
    {
        strptr.flip(0, bkt[0x0100] - bkt[0]).copy_back();

        for (unsigned int i = 0x0100; i < numbkts; i += 0x0100)
        {
            if (bkt[i] == bkt[i + 1]) continue;
            strptr.flip(bkt[i], bkt[i + 1] - bkt[i]).copy_back();
        }
    }

    delete[] bkt;
    delete this;
}

void EnqueueSmall(JobQueue& jobqueue, const StringPtr& strptr, size_t depth)
{
    EnqueueSmallsortJob8(jobqueue, strptr, depth);
}

template <typename bigsort_key_type>
void Enqueue(JobQueue& jobqueue, const StringPtr& strptr, size_t depth)
{
    if (strptr.size() > g_sequential_threshold)
        new RadixStepCE<bigsort_key_type>(jobqueue, strptr, depth);
    else
        EnqueueSmallsortJob8(jobqueue, strptr, depth);
}

static void parallel_radix_sort_8bit(string* strings, size_t n)
{
    g_totalsize = n;
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);

    string* shadow = new string[n]; // allocate shadow pointer array

    JobQueue jobqueue;
    Enqueue<uint8_t>(jobqueue, StringPtr(strings, shadow, n), 0);
    jobqueue.loop();

    delete[] shadow;
}

PSS_CONTESTANT_PARALLEL(parallel_radix_sort_8bit,
                        "bingmann/parallel_radix_sort_8bit",
                        "Parallel MSD Radix sort with load balancing, 8-bit BigSorts")

static void parallel_radix_sort_16bit(string* strings, size_t n)
{
    g_totalsize = n;
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);

    string* shadow = new string[n]; // allocate shadow pointer array

    JobQueue jobqueue;
    Enqueue<uint16_t>(jobqueue, StringPtr(strings, shadow, n), 0);
    jobqueue.loop();

    delete[] shadow;
}

PSS_CONTESTANT_PARALLEL(parallel_radix_sort_16bit,
                        "bingmann/parallel_radix_sort_16bit",
                        "Parallel MSD Radix sort with load balancing, 16-bit BigSorts")

} // namespace bingmann_parallel_radix_sort

/******************************************************************************/
