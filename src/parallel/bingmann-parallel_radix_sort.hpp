/*******************************************************************************
 * src/parallel/bingmann-parallel_radix_sort.hpp
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

#ifndef PSS_SRC_PARALLEL_BINGMANN_PARALLEL_RADIX_SORT_HEADER
#define PSS_SRC_PARALLEL_BINGMANN_PARALLEL_RADIX_SORT_HEADER

#include <cstdlib>
#include <cstring>

#include <iostream>
#include <vector>

#include "../tools/contest.hpp"
#include "../tools/globals.hpp"
#include "../tools/stringtools.hpp"
#include "../tools/stringset.hpp"
#include "../tools/jobqueue.hpp"

#include "../sequential/inssort.hpp"

#include <tlx/die.hpp>

namespace bingmann_parallel_radix_sort {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = false;

static const bool use_work_sharing = true;

static const size_t g_inssort_threshold = 32;

size_t g_totalsize;             // total size of input
size_t g_sequential_threshold;  // calculated threshold for sequential sorting
size_t g_threadnum;             // number of threads overall

/// Prototype called to schedule deeper sorts
template <typename bigsort_key_type, typename StringPtr>
void Enqueue(JobQueue& job_queue, const StringPtr& strings, size_t depth);

// ****************************************************************************
// *** SmallsortJob8 - sort 8-bit radix in-place with explicit stack-based recursion

template <typename StringPtr>
void EnqueueSmallsortJob8(JobQueue& job_queue, const StringPtr& strptr, size_t depth);

template <typename BktSizeType, typename StringPtr>
struct SmallsortJob8 final : public Job
{
    StringPtr strptr;
    size_t    depth;

    typedef BktSizeType bktsize_type;

    typedef typename StringPtr::StringSet StringSet;
    typedef typename StringSet::String String;
    typedef typename StringSet::Char Char;
    typedef typename StringSet::Iterator Iterator;

    SmallsortJob8(JobQueue& job_queue, const StringPtr& strptr, size_t _depth)
        : strptr(strptr), depth(_depth)
    {
        job_queue.enqueue(this);
    }

    struct RadixStep8_CI
    {
        StringPtr    strptr;
        size_t       idx;
        bktsize_type bkt[256 + 1];

        RadixStep8_CI(const StringPtr& strptr, size_t depth, Char* charcache)
            : strptr(strptr)
        {
            StringSet ss = strptr.active();
            size_t n = ss.size();

            // count character occurances
            bktsize_type bktsize[256];
#if 1
            // DONE: first variant to first fill charcache and then count. This is
            // 2x as fast as the second variant.
            Char* cc = charcache;
            for (Iterator i = ss.begin(); i != ss.end(); ++i, ++cc)
                *cc = ss.get_uint8(ss[i], depth);

            memset(bktsize, 0, sizeof(bktsize));
            for (cc = charcache; cc != charcache + n; ++cc)
                ++bktsize[static_cast<size_t>(*cc)];
#else
            memset(bktsize, 0, sizeof(bktsize));
            Iterator it = ss.begin();
            for (size_t i = 0; i < n; ++i, ++it) {
                ++bktsize[(charcache[i] = ss.get_uint8(ss[it], depth)];
                           }
#endif
            // inclusive prefix sum
            bkt[0] = bktsize[0];
            bktsize_type last_bkt_size = bktsize[0];
            for (size_t i = 1; i < 256; ++i) {
                bkt[i] = bkt[i - 1] + bktsize[i];
                if (bktsize[i]) last_bkt_size = bktsize[i];
            }

            // premute in-place
            for (size_t i = 0, j; i < n - last_bkt_size; )
            {
                String perm = std::move(ss[ss.begin() + i]);
                Char permch = charcache[i];
                while ((j = --bkt[static_cast<size_t>(permch)]) > i)
                {
                    std::swap(perm, ss[ss.begin() + j]);
                    std::swap(permch, charcache[j]);
                }
                ss[ss.begin() + i] = std::move(perm);
                i += bktsize[static_cast<size_t>(permch)];
            }

            // fix prefix sum
            bkt[0] = 0;
            for (size_t i = 1; i <= 256; ++i) {
                bkt[i] = bkt[i - 1] + bktsize[i - 1];
            }
            assert(bkt[256] == n);

            idx = 0; // will increment to 1 on first process, bkt 0 is not sorted further
        }
    };

    virtual bool run(JobQueue& job_queue)
    {
        size_t n = strptr.size();

        LOGC(debug_jobs)
            << "Process SmallsortJob8 " << this << " of size " << n;

        strptr = strptr.copy_back();

        if (n < g_inssort_threshold) {
            inssort::inssort_generic(strptr.copy_back().active(), depth);
            return true;
        }

        Char* charcache = new Char[n];

        // std::deque is much slower than std::vector, so we use an artificial
        // pop_front variable.
        size_t pop_front = 0;
        std::vector<RadixStep8_CI> radixstack;
        radixstack.emplace_back(strptr, depth, charcache);

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
                    inssort::inssort_generic(
                        rs.strptr.sub(rs.bkt[b], bktsize).copy_back().active(),
                        depth + radixstack.size());
                }
                else
                {
                    radixstack.emplace_back(rs.strptr.sub(rs.bkt[b], bktsize),
                                            depth + radixstack.size(),
                                            charcache);
                }

                if (use_work_sharing && job_queue.has_idle())
                {
                    // convert top level of stack into independent jobs
                    LOGC(debug_jobs)
                        << "Freeing top level of SmallsortJob8's radixsort stack";

                    RadixStep8_CI& rt = radixstack[pop_front];

                    while (rt.idx < 255)
                    {
                        b = ++rt.idx; // enqueue the bucket rt.idx

                        size_t bktsize = rt.bkt[b + 1] - rt.bkt[b];

                        if (bktsize == 0) continue;
                        EnqueueSmallsortJob8(
                            job_queue, rt.strptr.sub(rt.bkt[b], bktsize),
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

template <typename StringPtr>
void EnqueueSmallsortJob8(JobQueue& job_queue, const StringPtr& strptr, size_t depth)
{
    if (strptr.size() < ((uint64_t)1 << 32))
        new SmallsortJob8<uint32_t, StringPtr>(job_queue, strptr, depth);
    else
        new SmallsortJob8<uint64_t, StringPtr>(job_queue, strptr, depth);
}

// ****************************************************************************
// *** SmallsortJob16 - sort 16-bit radix in-place with explicit stack-based recursion

template <typename StringPtr>
void EnqueueSmallsortJob16(JobQueue& job_queue, const StringPtr& strptr, size_t depth);

template <typename BktSizeType, typename StringPtr>
struct SmallsortJob16 final : public Job
{
    StringPtr strptr;
    size_t    depth;

    typedef BktSizeType bktsize_type;

    SmallsortJob16(JobQueue& job_queue, const StringPtr& strptr, size_t depth)
        : strptr(strptr), depth(depth)
    {
        job_queue.enqueue(this);
    }

    struct RadixStep16_CI
    {
        typedef uint16_t key_type;
        static const size_t numbkts = key_traits<key_type>::radix;

        StringPtr           strptr;
        size_t              idx;
        bktsize_type        bkt[numbkts + 1];

        RadixStep16_CI(const StringPtr& strptr, size_t depth, key_type* charcache)
            : strptr(strptr)
        {
            string* strings = strptr.active().begin();
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
            for (size_t i = 1; i < numbkts; ++i) {
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
            for (size_t i = 1; i <= numbkts; ++i) {
                bkt[i] = bkt[i - 1] + bktsize[i - 1];
            }
            assert(bkt[numbkts] == n);

            idx = 0x0100; // will increment to 0x0101 on first process, bkt 0 is not sorted further
        }
    };

    virtual bool run(JobQueue& job_queue)
    {
        size_t n = strptr.size();

        LOGC(debug_jobs)
            << "Process SmallsortJob16 " << this << " of size " << n;

        strptr = strptr.copy_back();

        if (n < g_inssort_threshold) {
            inssort::inssort_generic(strptr.copy_back().active(), depth);
            return true;
        }

        typedef uint16_t key_type;

        static const size_t numbkts = key_traits<key_type>::radix;

        key_type* charcache = new key_type[n];

        // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
        size_t pop_front = 0;
        std::vector<RadixStep16_CI> radixstack;
        radixstack.emplace_back(strptr, depth, charcache);

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
                    inssort::inssort_generic(
                        rs.strptr.sub(rs.bkt[b], bktsize).copy_back().active(),
                        depth + 2 * radixstack.size());
                }
                else
                {
                    radixstack.emplace_back(
                        rs.strptr.sub(rs.bkt[b], bktsize),
                        depth + 2 * radixstack.size(),
                        charcache);
                }

                if (use_work_sharing && job_queue.has_idle())
                {
                    // convert top level of stack into independent jobs
                    LOGC(debug_jobs)
                        << "Freeing top level of SmallsortJob16's radixsort stack";

                    RadixStep16_CI& rt = radixstack[pop_front];

                    while (rt.idx < numbkts - 1)
                    {
                        b = ++rt.idx; // enqueue the bucket rt.idx

                        size_t bktsize = rt.bkt[b + 1] - rt.bkt[b];

                        if (bktsize == 0) continue;
                        EnqueueSmallsortJob16(
                            job_queue, rt.strptr.sub(rt.bkt[b], bktsize),
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

template <typename StringPtr>
void EnqueueSmallsortJob16(JobQueue& job_queue, const StringPtr& strptr, size_t depth)
{
    if (strptr.size() < ((uint64_t)1 << 32))
        new SmallsortJob16<uint32_t, StringPtr>(job_queue, strptr, depth);
    else
        new SmallsortJob16<uint64_t, StringPtr>(job_queue, strptr, depth);
}

// ****************************************************************************
// *** RadixStepCE out-of-place 8- or 16-bit parallel radix sort with Jobs

template <typename KeyType, typename StringPtr>
struct RadixStepCE
{
    typedef KeyType key_type;

    typedef typename StringPtr::StringSet StringSet;
    typedef typename StringSet::String String;
    typedef typename StringSet::Char Char;
    typedef typename StringSet::Iterator Iterator;

    static const size_t numbkts = key_traits<key_type>::radix;       // 256 or 65536

    StringPtr           strptr;
    size_t              depth;

    size_t              parts;
    size_t              psize;
    std::atomic<size_t> pwork;
    size_t              * bkt;

    key_type            * charcache;

    RadixStepCE(JobQueue& job_queue, const StringPtr& strptr, size_t _depth);

    void                count(size_t p, JobQueue& job_queue);
    void                count_finished(JobQueue& job_queue);

    void                distribute(size_t p, JobQueue& job_queue);
    void                distribute_finished(JobQueue& job_queue);
};

template <typename key_type, typename StringPtr>
struct CountJob final : public Job
{
    RadixStepCE<key_type, StringPtr>* step;
    size_t                          p;

    CountJob(RadixStepCE<key_type, StringPtr>* step, size_t p)
        : step(step), p(p) { }

    virtual bool run(JobQueue& job_queue)
    {
        step->count(p, job_queue);
        return true;
    }
};

template <typename key_type, typename StringPtr>
struct DistributeJob final : public Job
{
    RadixStepCE<key_type, StringPtr>* step;
    size_t                          p;

    DistributeJob(RadixStepCE<key_type, StringPtr>* step, size_t p)
        : step(step), p(p) { }

    virtual bool run(JobQueue& job_queue)
    {
        step->distribute(p, job_queue);
        return true;
    }
};

template <typename key_type, typename StringPtr>
RadixStepCE<key_type, StringPtr>::RadixStepCE(
    JobQueue& job_queue, const StringPtr& strptr, size_t _depth)
    : strptr(strptr), depth(_depth)
{
    size_t n = strptr.size();

    parts = (n + g_sequential_threshold - 1) / g_sequential_threshold;
    if (parts == 0) parts = 1;

    psize = (n + parts - 1) / parts;

    LOGC(debug_jobs)
        << "Area split into " << parts << " parts of size " << psize;

    bkt = new size_t[numbkts * parts + 1];
    charcache = new key_type[n];

    // create worker jobs
    pwork = parts;
    for (size_t p = 0; p < parts; ++p)
        job_queue.enqueue(new CountJob<key_type, StringPtr>(this, p));
}

template <typename key_type, typename StringPtr>
void RadixStepCE<key_type, StringPtr>::count(size_t p, JobQueue& job_queue)
{
    LOGC(debug_jobs)
        << "Process CountJob " << p << " @ " << this;

    const StringSet& ss = strptr.active();
    Iterator strB = ss.begin() + p * psize;
    Iterator strE = ss.begin() + std::min((p + 1) * psize, strptr.size());
    if (strE < strB) strE = strB;

    key_type* mycache = charcache + p * psize;
    key_type* mycacheE = mycache + (strE - strB);

    // DONE: check if processor-local stack + copy is faster. On 48-core AMD
    // Opteron it is not faster to copy first. On 32-core Intel Xeon the second
    // loop is slightly faster.
#if 0
    for (Iterator str = strB; str != strE; ++str, ++mycache)
        *mycache = ss.get_key<key_type>(*str, depth);

    size_t* mybkt = bkt + p * numbkts;
    memset(mybkt, 0, numbkts * sizeof(size_t));
    for (mycache = charcache + p * psize; mycache != mycacheE; ++mycache)
        ++mybkt[*mycache];
#else
    for (Iterator str = strB; str != strE; ++str, ++mycache)
        *mycache = parallel_string_sorting::get_key<key_type>(ss, *str, depth);

    size_t mybkt[numbkts] = { 0 };
    for (mycache = charcache + p * psize; mycache != mycacheE; ++mycache)
        ++mybkt[*mycache];
    memcpy(bkt + p * numbkts, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        count_finished(job_queue);
}

template <typename key_type, typename StringPtr>
void RadixStepCE<key_type, StringPtr>::count_finished(JobQueue& job_queue)
{
    LOGC(debug_jobs)
        << "Finishing CountJob " << this << " with prefixsum";

    // inclusive prefix sum over bkt
    size_t sum = 0;
    for (size_t i = 0; i < numbkts; ++i)
    {
        for (size_t p = 0; p < parts; ++p)
        {
            bkt[p * numbkts + i] = (sum += bkt[p * numbkts + i]);
        }
    }
    assert(sum == strptr.size());

    // create new jobs
    pwork = parts;
    for (size_t p = 0; p < parts; ++p)
        job_queue.enqueue(new DistributeJob<key_type, StringPtr>(this, p));
}

template <typename key_type, typename StringPtr>
void RadixStepCE<key_type, StringPtr>::distribute(size_t p, JobQueue& job_queue)
{
    LOGC(debug_jobs)
        << "Process DistributeJob " << p << " @ " << this;

    const StringSet& ss = strptr.active();
    Iterator strB = ss.begin() + p * psize;
    Iterator strE = ss.begin() + std::min((p + 1) * psize, strptr.size());
    if (strE < strB) strE = strB;

    Iterator sorted = strptr.shadow().begin(); // get alternative shadow pointer array
    key_type* mycache = charcache + p * psize;

    // DONE: check if processor-local stack + copy is faster. On 48-core AMD
    // Opteron it is not faster to copy first. On 32-core Intel Xeon the second
    // loop is slightly faster.
#if 0
    size_t* mybkt = bkt + p * numbkts;

    for (Iterator str = strB; str != strE; ++str, ++mycache)
        sorted[--mybkt[*mycache]] = *str;
#else
    size_t mybkt[numbkts];
    memcpy(mybkt, bkt + p * numbkts, sizeof(mybkt));

    for (Iterator str = strB; str != strE; ++str, ++mycache)
        sorted[--mybkt[*mycache]] = std::move(*str);

    if (p == 0) // these are needed for recursion into bkts
        memcpy(bkt, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        distribute_finished(job_queue);
}

template <typename key_type, typename StringPtr>
void RadixStepCE<key_type, StringPtr>::distribute_finished(JobQueue& job_queue)
{
    LOGC(debug_jobs)
        << "Finishing DistributeJob " << this << " with enqueuing subjobs";

    delete[] charcache;

    if (sizeof(key_type) == 1)
    {
        // first p's bkt pointers are boundaries between bkts, just add sentinel:
        assert(bkt[0] == 0);
        bkt[numbkts] = strptr.size();

        for (size_t i = 1; i < numbkts; ++i)
        {
            if (bkt[i] == bkt[i + 1])
                continue;
            else if (bkt[i] + 1 == bkt[i + 1]) // just one string pointer, copyback
                strptr.flip(bkt[i], 1).copy_back();
            else
                Enqueue<key_type>(job_queue, strptr.flip(bkt[i], bkt[i + 1] - bkt[i]),
                                  depth + key_traits<key_type>::add_depth);
        }

        // maybe copy back finished pointers from shadow
        strptr.flip(0, bkt[1] - bkt[0]).copy_back();
    }
    else if (sizeof(key_type) == 2)
    {
        // first p's bkt pointers are boundaries between bkts, just add sentinel:
        assert(bkt[0] == 0);
        bkt[numbkts] = strptr.size();

        for (size_t i = 0x0101; i < numbkts; ++i)
        {
            // skip over finished buckets 0x??00
            if ((i & 0xFF) == 0) ++i;

            if (bkt[i] == bkt[i + 1])
                continue;
            else if (bkt[i] + 1 == bkt[i + 1]) // just one string pointer, copyback
                strptr.flip(bkt[i], 1).copy_back();
            else
                Enqueue<key_type>(job_queue, strptr.flip(bkt[i], bkt[i + 1] - bkt[i]),
                                  depth + key_traits<key_type>::add_depth);
        }

        // maybe copy back finished pointers from shadow
        if (!strptr.flipped()) // if sorted[] = shadow
        {
            strptr.flip(0, bkt[0x0100] - bkt[0]).copy_back();

            for (size_t i = 0x0100; i < numbkts; i += 0x0100)
            {
                if (bkt[i] == bkt[i + 1]) continue;
                strptr.flip(bkt[i], bkt[i + 1] - bkt[i]).copy_back();
            }
        }
    }
    else
        die("impossible");

    delete[] bkt;
    delete this;
}

template <typename StringPtr>
void EnqueueSmall(JobQueue& job_queue, const StringPtr& strptr, size_t depth)
{
    EnqueueSmallsortJob8(job_queue, strptr, depth);
}

template <typename bigsort_key_type, typename StringPtr>
void Enqueue(JobQueue& job_queue, const StringPtr& strptr, size_t depth)
{
    if (strptr.size() > g_sequential_threshold && 0)
        new RadixStepCE<bigsort_key_type, StringPtr>(job_queue, strptr, depth);
    else
        EnqueueSmallsortJob8(job_queue, strptr, depth);
}

/******************************************************************************/
// Frontends

template <typename StringSet>
void parallel_radix_sort_8bit_generic(const StringSet& ss, size_t depth)
{
    g_totalsize = ss.size();
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);

    // allocate shadow pointer array
    typename StringSet::Container shadow = ss.allocate(ss.size());

    JobQueue job_queue;
    Enqueue<uint8_t>(
        job_queue, StringShadowPtr<StringSet>(ss, StringSet(shadow)), depth);
    job_queue.loop();

    StringSet::deallocate(shadow);
}

template <typename StringSet>
void parallel_radix_sort_16bit_generic(const StringSet& ss, size_t depth)
{
    g_totalsize = ss.size();
    g_threadnum = omp_get_max_threads();
    g_sequential_threshold = std::max(g_inssort_threshold, g_totalsize / g_threadnum);

    // allocate shadow pointer array
    typename StringSet::Container shadow = ss.allocate(ss.size());

    JobQueue job_queue;
    Enqueue<uint16_t>(
        job_queue, StringShadowPtr<StringSet>(ss, StringSet(shadow)), depth);
    job_queue.loop();

    StringSet::deallocate(shadow);
}

} // namespace bingmann_parallel_radix_sort

#endif // !PSS_SRC_PARALLEL_BINGMANN_PARALLEL_RADIX_SORT_HEADER

/******************************************************************************/
