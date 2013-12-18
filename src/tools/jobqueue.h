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
#include "../tools/agglogger.h"

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

    /// virtual function that is called by the JobQueue, delete object if run()
    /// returns true.
    virtual bool run(cookie_type& cookie) = 0;
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
    std::atomic<unsigned int> m_idle_count;

    //! SizeLogger or a dummy class
    typedef AggregateLogger<unsigned int> IntLogger;
    //typedef IntLogger::LockingAverageLogger logger_type;
    //typedef IntLogger::LockingMinLogger logger_type;
    typedef IntLogger::DummyLogger logger_type;

    logger_type m_logger, m_work_logger;

public:
    //typedef TimerArrayDummy TimerArrayMT;

    enum { TM_WORK, TM_IDLE };

    //! TimerArray for measing working and idle time (or a dummy class)
    TimerArrayMT m_timers;

public:

    JobQueueT()
        : m_queue(), m_idle_count(0),
          m_logger("jobqueue.txt", 0.01, 100),
          m_work_logger("worker_count.txt", 0.01, 4),
          m_timers(2)
    {
    }

    bool has_idle() const
    {
        return (m_idle_count != 0);
    }

    void enqueue(job_type* job)
    {
        m_queue.push(job);
        m_logger << m_queue.unsafe_size();
    }

    inline void executeThreadWork(cookie_type& cookie)
    {
        job_type* job = NULL;
        unsigned int numthrs = omp_get_num_threads();

        m_timers.change(TM_WORK);

        while (m_idle_count != numthrs)
        {
            if (m_queue.try_pop(job))
            {
                m_logger << m_queue.unsafe_size();

                if (job->run(cookie))
                    delete job;
            }
            else
            {
                DBG(debug_queue, "Queue is empty");

                ++m_idle_count;
                m_timers.change(TM_IDLE);
                m_work_logger << (numthrs - m_idle_count);

                while (m_idle_count != numthrs)
                {
                    DBG(debug_queue, "Idle thread - m_idle_count: " << m_idle_count);
                    if (m_queue.try_pop(job))
                    {
                        // got a new job -> not idle anymore
                        --m_idle_count;
                        m_timers.change(TM_WORK);

                        m_logger << m_queue.unsafe_size();
                        m_work_logger << (numthrs - m_idle_count);

                        if (job->run(cookie))
                            delete job;

                        break;
                    }
                }
            }

            m_work_logger << (numthrs - m_idle_count);
        }
    }

    void loop(cookie_type& cookie)
    {
        m_timers.start(omp_get_max_threads());

#pragma omp parallel
        {
            executeThreadWork(cookie);
        } // end omp parallel

        m_timers.stop();

        assert(m_queue.unsafe_size() == 0);
    }

    void numaLoop(int numaNode, int numberOfThreads, cookie_type& cookie)
    {
        m_timers.start(omp_get_max_threads());

#pragma omp parallel num_threads(numberOfThreads)
        {
            // tie thread to a NUMA node
            numa_run_on_node(numaNode);
            numa_set_preferred(numaNode);

            executeThreadWork(cookie);
        } // end omp parallel

        m_timers.stop();

        assert(m_queue.unsafe_size() == 0);
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
