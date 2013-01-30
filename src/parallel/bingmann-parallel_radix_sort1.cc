/******************************************************************************
 * src/parallel/bingmann-parallel_radix_sort1.h
 *
 * Parallel radix sort with work-balancing, variant 1.
 *
 * The set of strings is sorted using a 8-bit radix sort algorithm. Recursive
 * sorts are processed in parallel using a lock-free job queue and OpenMP
 * threads. Two radix sort implementations are used: sequential in-place and
 * parallelized out-of-place.
 *
 * The sequential radix sort is implemented using an explicit recursion stack,
 * which enables threads to "free up" work from the top of the stack when other
 * threads become idle. This variant uses in-place permuting and does _not_ use
 * a character cache (oracle).
 *
 * To parallelize sorting of large buckets an out-of-place variant is
 * implemented using a sequence of three Jobs: count, distribute and
 * copyback. All threads work on a job-specific context called RadixStepCE.
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

#include <omp.h>

#include <iostream>
#include <vector>
#include <atomic>

#include "../tools/debug.h"
#include "../tools/contest.h"
#include "../tools/jobqueue.h"

namespace bingmann_parallel_radix_sort1 {

using namespace jobqueue;

static const bool debug_jobs = false;
static const bool use_work_stealing = true;

typedef unsigned char* string;

/// Prototype called to schedule deeper sorts
void Enqueue(JobQueue& jobqueue, string* strings, size_t n, size_t depth);

// ****************************************************************************
// *** SmallsortJob - sort in-place with explicit stack-based recursion

struct SmallsortJob : public Job
{
    string*             strings;
    size_t              n;
    size_t              depth;

    SmallsortJob(string* _strings, size_t _n, size_t _depth)
        : strings(_strings), n(_n), depth(_depth) { }

    virtual void run(JobQueue& jobqueue);
};

static inline void
insertion_sort(string* strings, int n, size_t depth)
{
    for (string* i = strings + 1; --n > 0; ++i) {
        string* j = i;
        string tmp = *i;
        while (j > strings) {
            string s = *(j-1)+depth;
            string t = tmp+depth;
            while (*s == *t && *s != 0) ++s, ++t;
            if (*s <= *t) break;
            *j = *(j-1);
            --j;
        }
        *j = tmp;
    }
}

void SmallsortJob::run(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Process SmallsortJob " << this << " of size " << n);

    if (n < 32)
        return insertion_sort(strings,n,depth);

    struct RadixStep
    {
        string* str;
        size_t idx;
        size_t bktsize[256];

        RadixStep(string* strings, size_t n, size_t depth)
        {
            // count character occurances
            memset(bktsize, 0, sizeof(bktsize));
            for (size_t i=0; i < n; ++i)
                ++bktsize[ strings[i][depth] ];

            // inclusive prefix sum
            size_t bkt[256];
            bkt[0] = bktsize[0];
            size_t last_bkt_size = bktsize[0];
            for (unsigned int i=1; i < 256; ++i) {
                bkt[i] = bkt[i-1] + bktsize[i];
                if (bktsize[i]) last_bkt_size = bktsize[i];
            }

            // premute in-place
            for (size_t i=0, j; i < n-last_bkt_size; )
            {
                string perm = strings[i];
                while ( (j = --bkt[ perm[depth] ]) > i )
                {
                    std::swap(perm, strings[j]);
                }
                strings[i] = perm;
                i += bktsize[ perm[depth] ];
            }

            str = strings + bktsize[0];
            idx = 0; // will increment to 1 on first process, bkt 0 is not sorted further
        }
    };

    // std::deque is much slower than std::vector, so we use an artifical pop_front variable.
    size_t pop_front = 0;
    std::vector<RadixStep> radixstack;
    radixstack.push_back( RadixStep(strings,n,depth) );

    while ( radixstack.size() > pop_front )
    {
        while ( radixstack.back().idx < 255 )
        {
            RadixStep& rs = radixstack.back();
            ++rs.idx; // process the bucket rs.idx

            if (rs.bktsize[rs.idx] == 0)
                continue;
            else if (rs.bktsize[rs.idx] < 32)
            {
                insertion_sort(rs.str, rs.bktsize[rs.idx], depth + radixstack.size());
                rs.str += rs.bktsize[ rs.idx ];
            }
            else
            {
                rs.str += rs.bktsize[rs.idx];
                radixstack.push_back( RadixStep(rs.str - rs.bktsize[rs.idx],
                                                rs.bktsize[rs.idx],
                                                depth + radixstack.size()) );
                // cannot add to rs.str here, because rs may have invalidated
            }

            if (use_work_stealing && jobqueue.has_idle())
            {
                // convert top level of stack into independent jobs
                DBG(debug_jobs, "Freeing top level of SmallsortJob's radixsort stack");

                RadixStep& rt = radixstack[pop_front];

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
}

void EnqueueSmall(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
{
    jobqueue.enqueue( new SmallsortJob(strings,n,depth) );
}

// ****************************************************************************
// *** RadixStepCE out-of-place parallel radix sort with separate Jobs

struct RadixStepCE
{
    string*             strings;
    size_t              n;
    size_t              depth;
    
    unsigned int        parts;
    size_t              psize;
    std::atomic<unsigned int> pwork;
    size_t*             bkt;

    string*             sorted;

    void count(unsigned int p, JobQueue& jobqueue);
    void count_finished(JobQueue& jobqueue);

    void distribute(unsigned int p, JobQueue& jobqueue);
    void distribute_finished(JobQueue& jobqueue);

    void copyback(unsigned int p, JobQueue& jobqueue);
    void copyback_finished(JobQueue& jobqueue);
};

struct CountJob : public Job
{
    RadixStepCE*        step;
    unsigned int        p;

    CountJob(RadixStepCE* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process CountJob " << p << " @ " << step);
        step->count(p, jobqueue);
    }
};

struct DistributeJob : public Job
{
    RadixStepCE*        step;
    unsigned int        p;

    DistributeJob(RadixStepCE* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process DistributeJob " << p << " @ " << step);
        step->distribute(p, jobqueue);
    }
};

struct CopybackJob : public Job
{
    RadixStepCE*        step;
    unsigned int        p;

    CopybackJob(RadixStepCE* _step, unsigned int _p)
        : step(_step), p(_p) { }

    virtual void run(JobQueue& jobqueue)
    {
        DBG(debug_jobs, "Process CopybackJob " << p << " @ " << step);
        step->copyback(p, jobqueue);
    }
};

void RadixStepCE::count(unsigned int p, JobQueue& jobqueue)
{
    string* strB = strings + p * psize;
    string* strE = strings + std::min((p+1) * psize, n);

    // TODO: check if processor-local stack + copy is faster
#if 1
    size_t* mybkt = bkt + p * 256;

    memset(mybkt, 0, 256 * sizeof(size_t));
    for (string* str = strB; str != strE; ++str)
        ++mybkt[ (*str)[depth] ];
#else
    size_t mybkt[256] = { 0 };

    for (string* str = strB; str != strE; ++str)
        ++mybkt[ (*str)[depth] ];

    memcpy(bkt + p * 256, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        count_finished(jobqueue);
}

void RadixStepCE::count_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing CountJob " << this << " with prefixsum");

    // inclusive prefix sum over bkt
    size_t sum = 0;
    for (unsigned int i = 0; i < 256; ++i)
    {
        for (unsigned int p = 0; p < parts; ++p)
        {
            bkt[p * 256 + i] = (sum += bkt[p * 256 + i]);
        }
    }
    assert(sum == n);

    // allocate out-of-place memory
    sorted = (string*)malloc(n * sizeof(string));

    // create new jobs
    pwork = parts;

    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new DistributeJob(this, p) );
}

void RadixStepCE::distribute(unsigned int p, JobQueue& jobqueue)
{
    string* strB = strings + p * psize;
    string* strE = strings + std::min((p+1) * psize, n);

    // TODO: check if processor-local stack + copy is faster
#if 1
    size_t* mybkt = bkt + p * 256;

    for (string* str = strB; str != strE; ++str)
        sorted[ --mybkt[(*str)[depth]] ] = *str;
#else
    size_t mybkt[256];
    memcpy(mybkt, bkt + p * 256, sizeof(mybkt));

    for (string* str = strB; str != strE; ++str)
        sorted[ --mybkt[(*str)[depth]] ] = *str;

    if (p == 0) // these are needed for recursion into bkts
        memcpy(bkt, mybkt, sizeof(mybkt));
#endif

    if (--pwork == 0)
        distribute_finished(jobqueue);
}

void RadixStepCE::distribute_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing DistributeJob " << this << " with copy to original");

    // create new jobs
    pwork = parts;

    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new CopybackJob(this, p) );
}

void RadixStepCE::copyback(unsigned int p, JobQueue& jobqueue)
{
    size_t strB = p * psize;
    size_t strE = std::min((p+1) * psize, n);

    memcpy(strings + strB, sorted + strB, (strE - strB) * sizeof(string));

    if (--pwork == 0)
        copyback_finished(jobqueue);
}

void RadixStepCE::copyback_finished(JobQueue& jobqueue)
{
    DBG(debug_jobs, "Finishing CopybackJob " << this << " with copy to original");

    free(sorted);

    // first p's bkt pointers are boundaries between bkts, just add sentinel:
    bkt[256] = n;

    for (unsigned int i = 1; i < 256; ++i)
    {
        if (bkt[i] == bkt[i+1]) continue;
        Enqueue(jobqueue, strings + bkt[i], bkt[i+1] - bkt[i], depth+1);
    }

    delete [] bkt;
    delete this;
}

void EnqueueBig(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
{
    unsigned int parts = omp_get_max_threads();
    size_t psize = (n + parts-1) / parts;

    RadixStepCE* step = new RadixStepCE;
    step->strings = strings;
    step->n = n;
    step->depth = depth;

    step->parts = parts;
    step->psize = psize;
    step->pwork = parts;
    step->bkt = new size_t[256 * parts + 1];

    for (unsigned int p = 0; p < parts; ++p)
        jobqueue.enqueue( new CountJob(step, p) );
}

void Enqueue(JobQueue& jobqueue, string* strings, size_t n, size_t depth)
{
    // TODO: tune parameter
    if (n > 1024*1024)
        return EnqueueBig(jobqueue, strings, n, depth);
    else
        return EnqueueSmall(jobqueue, strings, n, depth);
}

void bingmann_parallel_radix_sort1(string* strings, size_t n)
{
    JobQueue jobqueue;
    Enqueue(jobqueue, strings,n,0);
    jobqueue.loop();
}

//CONTESTANT_REGISTER_UCARRAY(bingmann_parallel_radix_sort1, "bingmann/parallel_radix_sort1")

CONTESTANT_REGISTER_UCARRAY_PARALLEL(bingmann_parallel_radix_sort1, "bingmann/parallel_radix_sort1")

} // namespace bingmann_parallel_radix_sort1
