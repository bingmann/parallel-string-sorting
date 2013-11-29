/******************************************************************************
 * src/tools/jobqueue.h
 *
 * Job queue class for work-balancing parallel string sorting algorithms.
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

#ifndef JOBQUEUE_H_
#define JOBQUEUE_H_

#include <iostream>
#include <assert.h>
#include <omp.h>
#include <numa.h>

#include "src/config.h"

#if defined(HAVE_ATOMIC_H)
#include <atomic>
#elif defined(HAVE_CSTDATOMIC_H)
#include <cstdatomic>
#endif

#include <tbb/concurrent_queue.h>

#include "../tools/debug.h"

namespace jobqueue {

static const bool debug_queue = false;

// ****************************************************************************
// *** Job and JobQueue system with lock-free queue and OpenMP threads

template <typename CookieType>
class JobT
{
public:
    virtual ~JobT()
    { }

    /// local typedef of cookie
    typedef CookieType cookie_type;

    /// virtual function that is called by the JobQueue
    virtual void run(cookie_type& cookie) = 0;
};

template <typename CookieType>
class JobQueueT
{
public:
    /// local typedef of cookie
    typedef CookieType cookie_type;

    /// typedef of compatible Job
    typedef JobT<CookieType> job_type;

private:

    /// lock-free data structure containing pointers to Job objects.
    tbb::concurrent_queue<job_type*> m_queue;

    /// number of threads idle
    std::atomic<int> m_idle_count;

public:

    JobQueueT()
        : m_queue(), m_idle_count(0)
    {
    }

    bool has_idle() const
    {
        return (m_idle_count != 0);
    }

    void enqueue(job_type* job) {
        m_queue.push(job);
    }

    inline void executeThreadWork(cookie_type& cookie)
    {
        job_type* job = NULL;

        while (1)
        {
            if (m_queue.try_pop(job))
            {
            RUNJOB:
                job->run(cookie);
                delete job;
            }
            else
            {
                DBG(debug_queue, "Queue is empty");
                ++m_idle_count;

                while (m_idle_count != omp_get_num_threads())
                {
                    DBG(debug_queue, "Idle thread - m_idle_count: " << m_idle_count);
                    if (m_queue.try_pop(job))
                    {
                        --m_idle_count;
                        goto RUNJOB;
                    }
                }

                assert(m_idle_count == omp_get_num_threads());
                break;
            }
        }
    }

    void loop(cookie_type& cookie)
    {
#pragma omp parallel
        {
            executeThreadWork(cookie);
        } // end omp parallel
    }

    void numaLoop(int numaNode, int numberOfThreads, cookie_type& cookie)
    {
#pragma omp parallel num_threads(numberOfThreads)
        {
            numa_run_on_node(numaNode);
            numa_set_preferred(numaNode);

            executeThreadWork(cookie);
        } // end omp parallel
    }
};

class JobQueue : public JobQueueT<JobQueue>
{
public:
    typedef JobQueueT<JobQueue> super_type;

    void loop()
    {
        return super_type::loop(*this);
    }

    void numaLoop(int numaNode, int numberOfThreads)
    {
        return super_type::numaLoop(numaNode, numberOfThreads, *this);
    }
};

typedef JobT<JobQueue> Job;

} // namespace jobqueue

#endif // JOBQUEUE_H_
