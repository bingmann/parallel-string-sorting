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

class Job
{
public:
    virtual ~Job()
    { }

    /// virtual function that is called by the JobQueue
    virtual void run(class JobQueue& jobqueue) = 0;
};

class JobQueue
{
private:

    /// lock-free data structure containing pointers to Job objects.
    tbb::concurrent_queue<Job*> m_queue;

    /// number of threads idle
    std::atomic<int> m_idle_count;

public:

    JobQueue()
        : m_queue(), m_idle_count(0)
    {
    }

    bool has_idle() const
    {
        return (m_idle_count != 0);
    }

    void enqueue(class Job* job) {
        m_queue.push(job);
    }

    inline void executeThreadWork()
    {
        Job* job = NULL;

        while (1)
        {
            if (m_queue.try_pop(job))
            {
            RUNJOB:
                job->run(*this);
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

    void loop()
    {
#pragma omp parallel
        {
            executeThreadWork();
        } // end omp parallel
    }

    void numaLoop(int numaNode, int numberOfThreads)
    {
#pragma omp parallel num_threads(numberOfThreads)
        {
            numa_run_on_node(numaNode);
            numa_set_preferred(numaNode);

            executeThreadWork();
        } // end omp parallel
    }
};

} // namespace jobqueue

#endif // JOBQUEUE_H_
