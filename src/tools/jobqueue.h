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

#include <atomic>
#include <iostream>
#include <omp.h>

extern "C" {
#include <liblfds611.h>
}

#include "../tools/debug.h"

namespace jobqueue {

static const bool debug_queue = false;

// ****************************************************************************
// *** Job and JobQueue system with lock-free queue and OpenMP threads

class Job
{
public:
    virtual ~Job() {};

    /// virtual function that is called by the JobQueue
    virtual void run(class JobQueue& jobqueue) = 0;
};

class JobQueue
{
private:

    /// lock-free data structure containing void*s to Job objects.
    struct lfds611_queue_state*  m_queue;

    /// number of threads idle
    std::atomic<int>            m_idle_count;
    
public:

    JobQueue()
    {
        lfds611_queue_new(&m_queue, 2 * omp_get_num_threads());
        m_idle_count = 0;
    }

    ~JobQueue()
    {
        lfds611_queue_delete(m_queue, NULL, NULL);
    }

    bool has_idle() const
    {
        return (m_idle_count != 0);
    }

    void enqueue(class Job* job)
    {
        if (!lfds611_queue_enqueue(m_queue, job))
        {
            if (!lfds611_queue_guaranteed_enqueue(m_queue, job))
            {
                DBG(1, "Error with queue_guaranteed_enqueue!");
                abort();
            }
        }
    }
    
    void loop()
    {
#pragma omp parallel
        {
            void* job_data;
            lfds611_queue_use(m_queue);

            while(1)
            {
                if (lfds611_queue_dequeue(m_queue, &job_data))
                {
                RUNJOB:
                    Job* j = reinterpret_cast<Job*>(job_data);
                    j->run(*this);
                    delete j;
                }
                else
                {
                    DBG(debug_queue, "Queue is empty");
                    ++m_idle_count;

                    while (m_idle_count != omp_get_num_threads())
                    {
                        DBG(debug_queue, "Idle thread - m_idle_count: " << m_idle_count);
                        if (lfds611_queue_dequeue(m_queue, &job_data))
                        {
                            --m_idle_count;
                            goto RUNJOB;
                        }
                    }

                    assert(m_idle_count == omp_get_num_threads());
                    break;
                }
            }
        } // end omp parallel
    }
};

} // namespace jobqueue
