/******************************************************************************
 * src/parallel/bingmann-parallel_radix_sort3.h
 *
 * Parallel radix sort with work-balancing, variant 3.
 *
 * The set of strings is sorted using a 8- or 16-bit radix sort
 * algorithm. Recursive sorts are processed in parallel using a lock-free job
 * queue and OpenMP threads. Two radix sort implementations are used:
 * sequential in-place and parallelized out-of-place.
 *
 * The sequential radix sort is implemented using an explicit recursion stack,
 * which enables threads to "free up" work from the top of the stack when other
 * threads become idle. This variant uses in-place permuting and does _not_ use
 * a character cache (oracle).
 *
 * To parallelize sorting of large buckets an out-of-place variant is
 * implemented using a sequence of three Jobs: count, distribute and
 * copyback. All threads work on a job-specific context called RadixStepCE,
 * which encapsules all variables of a 16-bit radix sort (templatized with
 * key_type).
 *
 * Difference to variant 3: all radix sorts now use a character cache (also
 * called oracle in rantala sources).
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

namespace bingmann_parallel_radix_sort3 {

using namespace stringtools;
using namespace jobqueue;

static const bool debug_jobs = false;

static const bool use_work_stealing = true;

typedef unsigned char* string;

static const size_t g_inssort_threshold = 64;

/// Prototype called to schedule deeper sorts
void Enqueue(JobQueue& jobqueue, string* strings, size_t n, size_t depth);

// ****************************************************************************
// *** SmallsortJob - sort 8-bit radix in-place with explicit stack-based recursion

struct SmallsortJob : public Job
{
    string*             strings;
    size_t              n;
    size_t              depth;

    SmallsortJob(string* _strings, size_t _n, size_t _depth)
        : strings(_strings), n(_n), depth(_depth) { }

    virtual void run(JobQueue& jobqueue);
};

struct RadixStepCI
{
    string* str;
    size_t idx;
    size_t bktsize[256];

    RadixStepCI(string* strings, size_t n, size_t depth, uint8_t* charcache)
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
    }
};

void SmallsortJob::run(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process SmallsortJob " << this << " of size " << n);

    if (n < g_inssort_threshold)
        return inssort::inssort(strings,n,depth);

    uint8_t* charcache = new uint8_t[n];

    // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
    size_t pop_front = 0;
    std::vector<RadixStepCI> radixstack;
    radixstack.push_back( RadixStepCI(strings,n,depth,charcache) );

    while ( radixstack.size() > pop_front )
    {
        while ( radixstack.back().idx < 255 )
        {
            RadixStepCI& rs = radixstack.back();
            ++rs.idx; // process the bucket rs.idx

            if (rs.bktsize[rs.idx] == 0)
                continue;
            else if (rs.bktsize[rs.idx] < g_inssort_threshold)
            {
                inssort::inssort(rs.str, rs.bktsize[rs.idx], depth + radixstack.size());
                rs.str += rs.bktsize[ rs.idx ];
            }
            else
            {
                rs.str += rs.bktsize[rs.idx];
                radixstack.push_back( RadixStepCI(rs.str - rs.bktsize[rs.idx],
                                                  rs.bktsize[rs.idx],
                                                  depth + radixstack.size(),
                                                  charcache) );
                // cannot add to rs.str here, because rs may have invalidated
            }

            if (use_work_stealing && jobqueue.has_idle())
            {
                // convert top level of stack into independent jobs
                DBG(debug_jobs, "Freeing top level of SmallsortJob's radixsort stack");

                RadixStepCI& rt = radixstack[pop_front];

                while ( rt.idx < 255 )
                {
                    ++rt.idx; // enqueue the bucket rt.idx

                    if (rt.bktsize[rt.idx] == 0) continue;
                    jobqueue.enqueue( new SmallsortJob(rt.str, rt.bktsize[rt.idx], depth + pop_front) );
                    rt.str += rt.bktsize[rt.idx];
                }

                // shorten the current stack
                ++pop_front;
            }
        }
        radixstack.pop_back();
    }

    delete [] charcache;
}

// ****************************************************************************
// *** SmallsortJob16 - sort 16-bit radix in-place with explicit stack-based recursion

struct SmallsortJob16 : public Job
{
    string*             strings;
    size_t              n;
    size_t              depth;

    SmallsortJob16(string* _strings, size_t _n, size_t _depth)
        : strings(_strings), n(_n), depth(_depth) { }

    virtual void run(JobQueue& jobqueue);
};

struct RadixStep16
{
    typedef uint16_t key_type;
    static const size_t numbkts = key_traits<key_type>::radix;

    string* str;
    size_t idx;
    size_t bktsize[numbkts];

    RadixStep16(string* strings, size_t n, size_t depth, key_type* charcache)
    {
        // fill character cache
        for (size_t i=0; i < n; ++i)
            charcache[i] = get_char<key_type>(strings[i], depth);

        // count character occurances
        memset(bktsize, 0, sizeof(bktsize));
        for (size_t i=0; i < n; ++i)
            ++bktsize[ charcache[i] ];

        // inclusive prefix sum
        size_t bkt[numbkts];
        bkt[0] = bktsize[0];
        size_t last_bkt_size = bktsize[0];
        for (unsigned int i=1; i < numbkts; ++i) {
            bkt[i] = bkt[i-1] + bktsize[i];
            if (bktsize[i]) last_bkt_size = bktsize[i];
        }

        // premute in-place
        for (size_t i=0, j; i < n - last_bkt_size; )
        {
            string perm = strings[i];
            key_type permch = charcache[i];
            while ( (j = --bkt[ permch ]) > i )
            {
                std::swap(perm, strings[j]);
                std::swap(permch, charcache[j]);
            }
            strings[i] = perm;
            i += bktsize[ permch ];
        }

        idx = 0x0100; // will increment to 0x0101 on first process, bkt 0 is not sorted further
        str = strings + bkt[idx+1];
    }
};

void SmallsortJob16::run(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process SmallsortJob16 " << this << " of size " << n);

    if (n < g_inssort_threshold)
        return inssort::inssort(strings,n,depth);

    typedef uint16_t key_type;

    static const size_t numbkts = key_traits<key_type>::radix;

    key_type* charcache = new key_type[n];

    // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
    size_t pop_front = 0;
    std::vector<RadixStep16> radixstack;
    radixstack.push_back( RadixStep16(strings,n,depth,charcache) );

    while ( radixstack.size() > pop_front )
    {
        while ( radixstack.back().idx < numbkts-1 )
        {
            RadixStep16& rs = radixstack.back();
            ++rs.idx; // process the bucket rs.idx

            if ( (rs.idx & 0xFF) == 0 ) { // skip over finished 0x??00 buckets
                rs.str += rs.bktsize[rs.idx];
                ++rs.idx;
            }

            if (rs.bktsize[rs.idx] == 0)
                continue;
            else if (rs.bktsize[rs.idx] < g_inssort_threshold)
            {
                inssort::inssort(rs.str, rs.bktsize[rs.idx], depth + 2*radixstack.size());
                rs.str += rs.bktsize[rs.idx];
            }
            else
            {
                rs.str += rs.bktsize[rs.idx];
                radixstack.push_back( RadixStep16(rs.str - rs.bktsize[rs.idx],
                                                  rs.bktsize[rs.idx],
                                                  depth + 2*radixstack.size(),
                                                  charcache) );
                // cannot add to rs.str here, because rs may have invalidated
            }

            if (use_work_stealing && jobqueue.has_idle())
            {
                // convert top level of stack into independent jobs
                DBG(debug_jobs, "Freeing top level of SmallsortJob16's radixsort stack");

                RadixStep16& rt = radixstack[pop_front];

                while ( rt.idx < numbkts-1 )
                {
                    ++rt.idx; // enqueue the bucket rt.idx

                    if (rt.bktsize[rt.idx] == 0) continue;
                    jobqueue.enqueue( new SmallsortJob16(rt.str, rt.bktsize[rt.idx], depth + 2*pop_front) );
                    rt.str += rt.bktsize[rt.idx];
                }

                // shorten the current stack
                ++pop_front;
            }
        }
        radixstack.pop_back();
    }

    delete [] charcache;
}

void EnqueueSmall(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
{
    jobqueue.enqueue( new SmallsortJob(strings,n,depth) );
}

// ****************************************************************************
// *** RadixStepCE out-of-place 8- or 16-bit parallel radix sort with Jobs

template <typename KeyType>
struct RadixStepCE
{
    typedef KeyType key_type;

    static const size_t numbkts = key_traits<key_type>::radix;

    string*             strings;
    size_t              n;
    size_t              depth;

    unsigned int        parts;
    size_t              psize;
    std::atomic<unsigned int> pwork;
    size_t*             bkt;

    string*             sorted;

    key_type*           charcache;

    RadixStepCE() {}

    void count(unsigned int p, JobQueue& jobqueue);
    void count_finished(JobQueue& jobqueue);

    void distribute(unsigned int p, JobQueue& jobqueue);
    void distribute_finished(JobQueue& jobqueue);

    void copyback(unsigned int p, JobQueue& jobqueue);
    void copyback_finished(JobQueue& jobqueue);
};

template <typename key_type>
struct CountJob : public Job
{
    RadixStepCE<key_type>*      step;
    unsigned int                p;

    CountJob(RadixStepCE<key_type>* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        step->count(p, jobqueue);
    }
};

template <typename key_type>
struct DistributeJob : public Job
{
    RadixStepCE<key_type>*      step;
    unsigned int                p;

    DistributeJob(RadixStepCE<key_type>* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        step->distribute(p, jobqueue);
    }
};

template <typename key_type>
struct CopybackJob : public Job
{
    RadixStepCE<key_type>*      step;
    unsigned int                p;

    CopybackJob(RadixStepCE<key_type>* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        step->copyback(p, jobqueue);
    }
};

template <typename key_type>
void RadixStepCE<key_type>::count(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process CountJob " << p << " @ " << this);

    string* strB = strings + p * psize;
    string* strE = strings + std::min((p+1) * psize, n);
    if (strE < strB) strE = strB;

    key_type* mycache = charcache + p * psize;
    key_type* mycacheE = mycache + (strE - strB);

    // DONE: check if processor-local stack + copy is faster. On 48-core AMD
    // Opteron it is not faster to copy first.
#if 1
    for (string* str = strB; str != strE; ++str, ++mycache)
        *mycache = get_char<key_type>(*str, depth);

    size_t* mybkt = bkt + p * numbkts;
    memset(mybkt, 0, numbkts * sizeof(size_t));
    for (mycache = charcache + p * psize; mycache != mycacheE; ++mycache)
        ++mybkt[ *mycache ];
#else
    for (string* str = strB; str != strE; ++str, ++mycache)
        *mycache = get_char<key_type>(*str, depth);

    size_t mybkt[numbkts] = { 0 };
    for (mycache = charcache + p * psize; mycache != mycacheE; ++mycache)
        ++mybkt[ *mycache ];
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
    assert(sum == n);

    // allocate out-of-place memory
    sorted = (string*)malloc(n * sizeof(string));

    // create new jobs
    pwork = parts;

    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new DistributeJob<key_type>(this, p) );
}

template <typename key_type>
void RadixStepCE<key_type>::distribute(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process DistributeJob " << p << " @ " << this);

    string* strB = strings + p * psize;
    string* strE = strings + std::min((p+1) * psize, n);
    if (strE < strB) strE = strB;

    key_type* mycache = charcache + p * psize;

    // DONE: check if processor-local stack + copy is faster. On 48-core AMD
    // Opteron it is not faster to copy first.
#if 1
    size_t* mybkt = bkt + p * numbkts;

    for (string* str = strB; str != strE; ++str, ++mycache)
        sorted[ --mybkt[ *mycache ] ] = *str;
#else
    size_t mybkt[numbkts];
    memcpy(mybkt, bkt + p * numbkts, sizeof(mybkt));

    for (string* str = strB; str != strE; ++str, ++mycache)
        sorted[ --mybkt[ *mycache ] ] = *str;

    if (p == 0) // these are needed for recursion into bkts
        memcpy(bkt, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        distribute_finished(jobqueue);
}

template <typename key_type>
void RadixStepCE<key_type>::distribute_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing DistributeJob " << this << " with copy to original");

    // create new jobs
    pwork = parts;

    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new CopybackJob<key_type>(this, p) );
}

template <typename key_type>
void RadixStepCE<key_type>::copyback(unsigned int p, JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process CopybackJob " << p << " @ " << this);

    size_t strB = p * psize;
    size_t strE = std::min((p+1) * psize, n);
    if (strE < strB) strE = strB;

    memcpy(strings + strB, sorted + strB, (strE - strB) * sizeof(string));

    if (--pwork == 0)
        copyback_finished(jobqueue);
}

template <>
void RadixStepCE<uint8_t>::copyback_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing CopybackJob " << this << " with copy to original");

    free(sorted);
    delete [] charcache;

    // first p's bkt pointers are boundaries between bkts, just add sentinel:
    bkt[numbkts] = n;

    for (unsigned int i = 1; i < numbkts; ++i)
    {
        if (bkt[i] + 2 > bkt[i+1]) continue;
        Enqueue(jobqueue, strings + bkt[i], bkt[i+1] - bkt[i], depth + key_traits<key_type>::add_depth);
    }

    delete [] bkt;
    delete this;
}

template <>
void RadixStepCE<uint16_t>::copyback_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing CopybackJob " << this << " with copy to original");

    free(sorted);
    delete [] charcache;

    // first p's bkt pointers are boundaries between bkts, just add sentinel:
    bkt[numbkts] = n;

    for (unsigned int i = 0x0101; i < numbkts; ++i)
    {
        if ((i & 0xFF) == 0) ++i; // skip over finished buckets 0x??00

        if (bkt[i] + 2 > bkt[i+1]) continue;
        Enqueue(jobqueue, strings + bkt[i], bkt[i+1] - bkt[i], depth + key_traits<key_type>::add_depth);
    }

    delete [] bkt;
    delete this;
}

size_t totalsize;

void EnqueueBig(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
{
    size_t parts = omp_get_max_threads() * n / totalsize;
    if (parts == 0) parts = 1;
    size_t psize = (n + parts-1) / parts;

    typedef uint16_t key_type;
    typedef RadixStepCE<key_type> RadixStep_type;

    RadixStep_type* step = new RadixStep_type;
    step->strings = strings;
    step->n = n;
    step->depth = depth;

    step->parts = parts;
    step->psize = psize;
    step->pwork = parts;
    step->bkt = new size_t[RadixStep_type::numbkts * parts + 1];

    step->charcache = new key_type[n];

    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new CountJob<key_type>(step, p) );
}

void Enqueue(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
{
    // TODO: tune parameter
    if (n > 128*1024)
        return EnqueueBig(jobqueue, strings, n, depth);
    else
        return EnqueueSmall(jobqueue, strings, n, depth);
}

void parallel_radix_sort3(string* strings, size_t n)
{
    totalsize = n;

    JobQueue jobqueue;
    Enqueue(jobqueue, strings, n, 0);
    jobqueue.loop();
}

CONTESTANT_REGISTER_PARALLEL(parallel_radix_sort3,
                             "bingmann/parallel_radix_sort3",
                             "bingmann/parallel_radix_sort3")

} // namespace bingmann_parallel_radix_sort3
