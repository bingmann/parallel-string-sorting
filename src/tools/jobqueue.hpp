/*******************************************************************************
 * src/tools/jobqueue.hpp
 *
 * Job queue class for work-balancing parallel string sorting algorithms.
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

#ifndef PSS_SRC_TOOLS_JOBQUEUE_HEADER
#define PSS_SRC_TOOLS_JOBQUEUE_HEADER

#include <iostream>
#include <cassert>

#include <omp.h>
#include <numa.h>

#include "src/config.h"

#if defined(HAVE_ATOMIC_H)
#include <atomic>
#elif defined(HAVE_CSTDATOMIC_H)
#include <cstdatomic>
#endif

#include <tbb/concurrent_queue.h>

#include "../tools/agglogger.hpp"
#include "../tools/timer.hpp"
#include "../tools/timer_array.hpp"
#include "../tools/globals.hpp"

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
class DefaultJobQueueGroup;

//! Define JobQueue with templatized cookie parameter, which is passed to run()
//! functions.
template <typename CookieType,
          template <typename> class JobQueueGroupType = DefaultJobQueueGroup>
class JobQueueT
{
public:
    /// local typedef of cookie
    typedef CookieType cookie_type;

    /// typedef of compatible Job
    typedef JobT<CookieType> job_type;

    /// typedef of JobQueueGroup
    typedef JobQueueGroupType<CookieType> jobqueuegroup_type;

private:
    /// lock-free data structure containing pointers to Job objects.
    tbb::concurrent_queue<job_type*> m_queue;

    //! number of threads working on queue
    unsigned m_numthrs;

    /// number of threads idle
    std::atomic<unsigned int> m_idle_count;

    //! reference to the cookie of this JobQueue
    cookie_type& m_cookie;

    //! reference to associated JobQueueGroup
    jobqueuegroup_type* m_group;

    //! number of this JobQueue in JobQueueGroup
    unsigned m_id;

    //! SizeLogger or a dummy class
    typedef AggregateLogger<unsigned int> IntLogger;
    //typedef IntLogger::LockingAverageLogger logger_type;
    typedef IntLogger::DummyLogger logger_type;

    logger_type m_logger, m_work_logger;

public:
    typedef TimerArrayDummy TimerArrayMT;

    enum { TM_WORK, TM_IDLE };

    //! TimerArray for measing working and idle time (or a dummy class)
    TimerArrayMT m_timers;

public:
    JobQueueT(cookie_type& cookie,
              jobqueuegroup_type* group)
        : m_queue(),
          m_numthrs(0),
          m_idle_count(0),
          m_cookie(cookie),
          m_group(group),
          m_logger("jobqueue.txt", 0.005, 10000),
          m_work_logger("worker_count.txt", 0.005, 10000),
          m_timers(2)
    { }

    bool has_idle() const
    {
        return (m_idle_count.load(std::memory_order_relaxed) != 0);
    }

    void enqueue(job_type* job)
    {
        m_queue.push(job);
        m_logger << m_queue.unsafe_size();
    }

    void set_id(unsigned id)
    {
        m_id = id;
    }

    //! try to run one jobs from the queue, returns false if queue is finished,
    //! true if ran jobs or queue not finished.
    bool try_run()
    {
        job_type* job = NULL;

        if (!m_queue.try_pop(job))
            return (m_idle_count != m_numthrs);

        m_logger << m_queue.unsafe_size();

        if (job->run(m_cookie))
            delete job;

        return true;
    }

    inline void executeThreadWork()
    {
        job_type* job = NULL;
        m_numthrs = omp_get_num_threads();

        m_timers.change(TM_WORK);
        m_logger.start();
        m_work_logger.start();

        while (true)
        {
            while (m_queue.try_pop(job))
            {
                m_logger << m_queue.unsafe_size();

                if (job->run(m_cookie))
                    delete job;
            }

            LOGC(debug_queue) << "Queue" << m_id << " is empty";

            // no more jobs -> switch to idle
            m_timers.change(TM_IDLE);
            ++m_idle_count;

            m_logger << m_queue.unsafe_size();
            m_work_logger << (m_numthrs - m_idle_count);

            while (!m_queue.try_pop(job))
            {
                LOGC(debug_queue)
                    << "Idle thread - m_idle_count: " << m_idle_count;

                if ( //!m_group->assist(m_id) &&
                    m_idle_count == m_numthrs)
                {
                    // assist other JobQueues before terminating.
                    while (m_group->assist(m_id)) { }
                    return;
                }
            }

            // got a new job -> not idle anymore
            m_timers.change(TM_WORK);
            --m_idle_count;

            m_logger << m_queue.unsafe_size();
            m_work_logger << (m_numthrs - m_idle_count);

            if (job->run(m_cookie))
                delete job;
        }
    }

    void loop()
    {
        m_timers.start(omp_get_max_threads());
        m_idle_count = 0;

#pragma omp parallel
        {
            if (gopt_memory_type == "mmap_node0")
            {
                // tie thread to first NUMA node
                numa_run_on_node(0);
                numa_set_preferred(0);
            }

            executeThreadWork();
        }   // end omp parallel

        m_timers.stop();

        assert(m_queue.unsafe_size() == 0);
    }

    void numaLoop(int numaNode, int numberOfThreads)
    {
        m_timers.start(omp_get_max_threads());

#pragma omp parallel num_threads(numberOfThreads)
        {
            // tie thread to a NUMA node
            numa_run_on_node(numaNode);
            numa_set_preferred(numaNode);

            executeThreadWork();
        }   // end omp parallel

        m_timers.stop();

        assert(m_queue.unsafe_size() == 0);
    }
};

//! Define no-operation JobQueueGroup for "standard" JobQueue
template <typename CookieType>
class DefaultJobQueueGroup
{
public:
    /// typedef of compatible JobQueue
    typedef JobQueueT<CookieType, DefaultJobQueueGroup> jobqueue_type;

    /// typedef of compatible Job
    typedef JobT<CookieType> job_type;

public:
    static inline bool assist(unsigned)
    {
        return false;
    }
};

//! Define NumaJobQueueGroup to group JobQueue which assist each other when idle.
template <typename CookieType>
class NumaJobQueueGroup
{
public:
    /// typedef of compatible JobQueue
    typedef JobQueueT<CookieType, NumaJobQueueGroup> jobqueue_type;

    /// typedef of compatible Job
    typedef JobT<CookieType> job_type;

protected:
    //! List of managed JobQueues.
    std::vector<jobqueue_type*> m_queues;

public:
    //! Register a JobQueue in the group, this function is NOT THREAD-SAFE.
    void add_jobqueue(jobqueue_type* jq)
    {
        jq->set_id(m_queues.size());
        m_queues.push_back(jq);
    }

    //! Calculate number of the threads the k-th JobQueue of numJobQueues gets.
    static unsigned calcThreadNum(int k, int numJobQueues)
    {
        int realNumaNodes = numa_num_configured_nodes();
        if (realNumaNodes < 1) realNumaNodes = 1;

        int numThreadsPerNode = omp_get_max_threads() / numJobQueues;
        int remainThreads = omp_get_max_threads() % numJobQueues;

        int nodeThreads = numThreadsPerNode;
        if (k < remainThreads) nodeThreads++; // distribute extra threads

        LOG1 << "JobQueue[" << k << "] prospective " << nodeThreads << " threads";

        return nodeThreads;
    }

    //! Launch all threads divided evenly among JobQueues.
    void numaLaunch()
    {
        // check how many real NUMA nodes there are
        int realNumaNodes = numa_num_configured_nodes();
        if (realNumaNodes < 1) realNumaNodes = 1;

        if (realNumaNodes == 1) {
            LOG1 << "No or just one NUMA nodes detected on the system.";
            LOG1 << "Continuing anyway, at your own peril!";
        }

        g_stats >> "num_real_numa_nodes" << realNumaNodes;

        if ((int)m_queues.size() != realNumaNodes || g_numa_nodes == 0)
        {
            LOG1 << "!!! WARNING !!! emulating NUMA nodes! "
                 << "Remove --numa-nodes for REAL EXPERIMENTS.";
        }

        g_stats >> "num_jobqueues" << m_queues.size();

        // distribute threads among NUMA job queues
        int numJobQueues = m_queues.size();
        int numThreadsPerNode = omp_get_max_threads() / numJobQueues;
        int remainThreads = omp_get_max_threads() % numJobQueues;

        if (numThreadsPerNode == 0)
        {
            LOG1 << "Fewer threads than NUMA nodes detected.";
            LOG1 << "Strange things may happen, continuing anyway, at your own peril!";

            // We will start fewer threads than JobQueues, and wait for the
            // first to finish, which will then assist the JobQueues without
            // threads.

            //abort();
        }

        int runThreads = std::min(omp_get_max_threads(), numJobQueues);

        omp_set_nested(true); // enable nested parallel regions

        ClockTimer timer;

#pragma omp parallel for num_threads(runThreads) schedule(dynamic)
        for (int k = 0; k < numJobQueues; k++)
        {
            int nodeThreads = numThreadsPerNode;
            int numaNode = k % realNumaNodes;
            if (k < remainThreads) nodeThreads++; // distribute extra threads

            LOG1 << "JobQueue[" << k << "] gets " << nodeThreads << " threads";

            if (nodeThreads == 0) nodeThreads = 1;

            ClockTimer timer;

            m_queues[k]->numaLoop(numaNode, nodeThreads);

            LOG1 << "JobQueue[" << k << "] took : " << timer.elapsed() << " s";
        }
    }

    //! called by JobQueue's when they want to assist other queues.
    bool assist(unsigned qid)
    {
        unsigned id = qid;

        for (unsigned i = 1; i < m_queues.size(); ++i)
        {
            // go through queues round-robin starting at own
            if (++id >= m_queues.size()) id = 0;

            if (m_queues[id]->try_run())
            {
                LOGC(debug_queue) << "JobQueue[" << qid << "] assisted " << id;
                return true;
            }
        }

        return false;
    }
};

//! Define "standard" JobQueue, which passes a reference to itself as cookie
//! parameter to each run() call.
class JobQueue : public JobQueueT<JobQueue>
{
public:
    typedef JobQueueT<JobQueue> super_type;

    //! Constructor, set JobQueue cookie to ourselves.
    JobQueue()
        : super_type(*this, NULL)
    { }
    void loop()
    {
        return super_type::loop();
    }

    void numaLoop(int numaNode, int numberOfThreads)
    {
        return super_type::numaLoop(numaNode, numberOfThreads);
    }
};

//! Define "standard" Job for "standard" JobQueue.
typedef JobT<JobQueue> Job;

} // namespace jobqueue

#endif // !PSS_SRC_TOOLS_JOBQUEUE_HEADER

/******************************************************************************/
