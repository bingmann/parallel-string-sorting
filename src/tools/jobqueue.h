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

#include <iostream>
#include <assert.h>
#include <omp.h>

#include "src/config.h"

#if defined(HAVE_ATOMIC_H)
#include <atomic>
#elif defined(HAVE_CSTDATOMIC_H)
#include <cstdatomic>
#endif

extern "C" {
#include <liblfds611.h>
}

#include <tbb/concurrent_queue.h>

#include "../tools/debug.h"

namespace jobqueue {

static const bool debug_queue = false;

// ****************************************************************************
// *** C++ frontend to lfds611_queue

template <typename Type>
class LockfreeQueue
{
private:
    /// lock-free data structure containing void*s.
    struct lfds611_queue_state* m_queue;

public:
    /// This function instantiates a queue. After instantiation any thread
    /// (except the instantiating thread) must before using the queue first
    /// call lfds611_queue_use.
    inline LockfreeQueue(lfds611_atom_t number_elements=0)
    {
        lfds611_queue_new(&m_queue, number_elements);
    }

    inline ~LockfreeQueue()
    {
        lfds611_queue_delete(m_queue, NULL, NULL);
/*
  lfds611_queue_delete( struct lfds611_queue_state *qs,
  void (*user_data_delete_function)(void *user_data, void *user_state),
  void *user_state );
*/
    }

    /// After a queue has been instantiated by calling lfds611_queue_new, any
    /// thread (except the instantiating thread) wishing to use that queue must
    /// first call this function.
    inline void use()
    {
        lfds611_queue_use(m_queue);
    }

    /// Enqueuing only fails if the queue has exhausted its supply of freelist
    /// elements. In this case, the user can call
    /// lfds611_queue_guaranteed_enqueue, which will allocate a new element and
    /// enqueue using that. Remember however that the queue can never shrink,
    /// so any such call will permanently increase the size of the queue by one
    /// element.
    inline bool enqueue(const Type& data)
    {
        return lfds611_queue_enqueue(m_queue, data);
    }

    /// The function lfds611_queue_enqueue fails only when the queue's freelist
    /// is empty. In this event, lfds611_queue_guaranteed_push can be called,
    /// which allocates a new element and enqueues using that new element, thus
    /// guaranteeing an enqueue, barring the event of malloc failure.
    inline bool guaranteed_enqueue(const Type& data)
    {
        return lfds611_queue_guaranteed_enqueue(m_queue, data);
    }

    /// Combination of enqueue() and guaranteed_enqueue().
    inline void push(const Type& data)
    {
        if (!enqueue(data) && !guaranteed_enqueue(data))
        {
            DBG(1, "Error with queue_guaranteed_enqueue!");
            abort();
        }
    }

    /// The queue being empty is the only situation in which dequeuing does not
    /// occur. Note that on an unsuccessful dequeue, *user_data is untouched;
    /// it is not set to NULL, since NULL is a valid user data value. Only the
    /// return value indicates whether or not the dequeue was successful.
    inline bool try_pop(Type& data)
    {
        return lfds611_queue_dequeue(m_queue, reinterpret_cast<void**>(&data));
    }

    /// Return approximate size of queue.
    inline size_t size()
    {
        lfds611_atom_t size;
        lfds611_queue_query(m_queue, LFDS611_QUEUE_QUERY_ELEMENT_COUNT, NULL, &size);
        return size;
    }

    /// swap pointers with other object
    inline void swap(LockfreeQueue& other)
    {
        std::swap(m_queue, other.m_queue);
    }
};


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

    /// lock-free data structure containing pointers to Job objects.
    //LockfreeQueue<Job*> m_queue;
    tbb::concurrent_queue<Job*>   m_queue;

    /// number of threads idle
    std::atomic<int>            m_idle_count;

public:

    JobQueue()
        : m_queue(),
          m_idle_count( 0 )
    {
    }

    bool has_idle() const
    {
        return (m_idle_count != 0);
    }

    void enqueue(class Job* job)
    {
        m_queue.push(job);
    }

    void loop()
    {
#pragma omp parallel
        {
            Job* job;

            while(1)
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
        } // end omp parallel
    }
};

} // namespace jobqueue
