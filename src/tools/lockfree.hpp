/*******************************************************************************
 * src/tools/lockfree.hpp
 *
 * Lock-free auxiliary data structures
 *
 *******************************************************************************
 * Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
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

#ifndef PSS_SRC_TOOLS_LOCKFREE_HEADER
#define PSS_SRC_TOOLS_LOCKFREE_HEADER

namespace lockfree {

template <unsigned MaxThreads>
class lazy_counter
{
protected:
    //! per thread information, spaced out onto cache-lines
    struct PerThread
    {
        std::atomic<ssize_t> delta;
        char                 filled[64 - sizeof(ssize_t)];
    };

    //! per thread aggregated delta.
    struct PerThread m_thr[MaxThreads];

    //! lazy counter value, updated only on request
    std::atomic<ssize_t> m_counter;

    //! boolean whether any thread changed a value
    std::atomic<bool> m_changed;

    //! number of valid threads
    size_t m_nthr;

public:
    lazy_counter(ssize_t initial = 0, size_t nthr = omp_get_num_procs())
        : m_counter(initial),
          m_changed(false),
          m_nthr(nthr)
    {
        memset(m_thr, 0, sizeof(m_thr));
    }

    //! direct assignment (unsafe)
    lazy_counter& operator = (ssize_t initial)
    {
        m_counter = initial;
        m_changed = false;
        memset(m_thr, 0, sizeof(m_thr));
        return *this;
    }

    //! lazy adding of a delta
    lazy_counter & add(ssize_t delta, size_t thr)
    {
        m_thr[thr].delta.fetch_add(delta, std::memory_order_relaxed);
        m_changed.store(true, std::memory_order_relaxed);
        return *this;
    }

    //! lazy adding of a delta
    lazy_counter & add(ssize_t delta)
    {
        return add(delta, omp_get_thread_num());
    }

    //! update the lazy couter
    lazy_counter & update()
    {
        if (!m_changed) return *this;

        ssize_t agg_delta = 0;

        for (size_t t = 0; t < m_nthr; ++t)
        {
            // get delta and decrement from thread's counter

            ssize_t delta = m_thr[t].delta;
            if (delta == 0) continue;

            m_thr[t].delta.fetch_sub(delta, std::memory_order_relaxed);

            agg_delta += delta;
        }

        // add aggregated delta
        m_counter.fetch_add(agg_delta, std::memory_order_relaxed);

        m_changed.store(false, std::memory_order_relaxed);
        return *this;
    }

    //! get lazy counter value
    ssize_t get() const
    {
        return m_counter;
    }
};

} // namespace lockfree

#endif // !PSS_SRC_TOOLS_LOCKFREE_HEADER

/******************************************************************************/
