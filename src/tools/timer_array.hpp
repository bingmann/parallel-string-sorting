/*******************************************************************************
 * src/tools/timer_array.hpp
 *
 * Class to output statistics in a flexible text file as key=value pairs.
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

#ifndef PSS_SRC_TOOLS_TIMER_ARRAY_HEADER
#define PSS_SRC_TOOLS_TIMER_ARRAY_HEADER

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <unistd.h>

#include <omp.h>

/// Class to measure different parts of a funciton by switching between
/// different aggregating timers. Immediately start with timer 0. Use enums in
/// your code to give the timer numbers names.
class TimerArray
{
private:
    //! clock time of last call
    struct timespec m_tplast;

    //! currently running timer
    unsigned int m_tmcurr;

    //! array of timers (usually preallocated)
    std::vector<struct timespec> m_tpvector;

public:
    //! real or dummy implementation
    static const bool is_real = true;

    TimerArray(unsigned int timers)
    {
        m_tplast.tv_sec = m_tplast.tv_nsec = 0;
        m_tpvector.resize(timers, m_tplast);

        clear();
    }

    /// clear all timers and start counting in timer 0.
    void clear()
    {
        m_tplast.tv_sec = m_tplast.tv_nsec = 0;
        std::fill(m_tpvector.begin(), m_tpvector.end(), m_tplast);
        m_tmcurr = 0;

        if (clock_gettime(CLOCK_MONOTONIC, &m_tplast)) {
            perror("Could not clock_gettime(CLOCK_MONOTONIC)");
        }
    }

    /// switch to other timer
    inline void change(unsigned int tm)
    {
        assert(tm < m_tpvector.size());
        struct timespec tpnow;

        // get current time
        if (clock_gettime(CLOCK_MONOTONIC, &tpnow)) {
            perror("Could not clock_gettime(CLOCK_MONOTONIC)");
        }

        // add difference to current timer
        m_tpvector[m_tmcurr].tv_sec += tpnow.tv_sec - m_tplast.tv_sec;
        m_tpvector[m_tmcurr].tv_nsec += tpnow.tv_nsec - m_tplast.tv_nsec;

        m_tplast = tpnow;
        m_tmcurr = tm;
    }

    /// return amount of time spent in a timer
    inline double get(unsigned int tm)
    {
        assert(tm < m_tpvector.size());
        return (m_tpvector[tm].tv_sec + m_tpvector[tm].tv_nsec / 1e9);
    }
};

/// Dummy class to replace TimerArray calls with no-ops.
class TimerArrayDummy
{
public:
    //! real or dummy implementation
    static const bool is_real = false;

    inline TimerArrayDummy(unsigned int /* timers */)
    { }

    //! clear all timers and start counting in timer 0.
    inline void clear()
    { }

    //! clear timers and start nthreads timers
    inline void start(unsigned /* nthreads */)
    { }

    //! stop all timers, add remaining time
    inline void stop()
    { }

    //! switch to other timer
    inline void change(unsigned int /* tm */)
    { }

    //! return amount of time spent in a timer
    inline double get(unsigned int /* tm */)
    {
        return 0;
    }

    //! return sum over all timers
    inline double get_sum() const
    {
        return 0;
    }
};

/// Class to measure different parts of a funciton by switching between
/// different aggregating timers. Immediately start with timer 0. Use enums in
/// your code to give the timer numbers names. Multi-threading aware version.
class TimerArrayMT
{
private:
    //! struct of information per thread
    struct ThreadInfo
    {
        //! clock time of last call per thread
        struct timespec m_tplast;

        //! currently running timer per thread
        unsigned int    m_tmcurr;

        //! filler to put info into different cache lines
        unsigned char   m_filler[64 - sizeof(struct timespec) - sizeof(unsigned int)];
    } __attribute__ ((packed));

    //! array of timers (preallocated)
    std::vector<struct timespec> m_timers;

    //! array of thread info
    std::vector<ThreadInfo> m_thread;

public:
    //! real or dummy implementation?
    static const bool is_real = true;

    //! clear all timers at construction
    TimerArrayMT(unsigned ntimers)
        : m_timers(ntimers)
    { }

    //! clear all timers and start counting in timer 0.
    void start(unsigned nthreads)
    {
        struct timespec tp;
        memset(&tp, 0, sizeof(tp));

        std::fill(m_timers.begin(), m_timers.end(), tp);

        if (clock_gettime(CLOCK_MONOTONIC, &tp)) {
            perror("Could not clock_gettime(CLOCK_MONOTONIC)");
        }

        m_thread.resize(nthreads);

        for (size_t i = 0; i < m_thread.size(); ++i) {
            m_thread[i].m_tplast = tp;
            m_thread[i].m_tmcurr = 0;
        }
    }

    //! stop all timers, add remaining time
    void stop()
    {
        for (size_t t = 0; t < m_thread.size(); ++t) {
            change(0, t);
        }
    }

    //! switch to other timer tm for the thread tid
    inline void change(unsigned int tm, unsigned int tid)
    {
        assert(tm < m_timers.size());
        assert(tid < m_thread.size());

        // get current time
        struct timespec tpnow;
        if (clock_gettime(CLOCK_MONOTONIC, &tpnow)) {
            perror("Could not clock_gettime(CLOCK_MONOTONIC)");
        }

        // thread info
        ThreadInfo& ti = m_thread[tid];

        // add difference to current timer
#pragma omp atomic
        m_timers[ti.m_tmcurr].tv_sec += (tpnow.tv_sec - ti.m_tplast.tv_sec);
#pragma omp atomic
        m_timers[ti.m_tmcurr].tv_nsec += tpnow.tv_nsec - ti.m_tplast.tv_nsec;

        ti.m_tplast = tpnow;
        ti.m_tmcurr = tm;
    }

    //! switch to other timer tm for the current thread
    inline void change(unsigned int tm)
    {
        return change(tm, omp_get_thread_num());
    }

    //! return current timer of a thread
    inline unsigned int get_tmcurr(unsigned int tid)
    {
        assert(tid < m_thread.size());
        return m_thread[tid].m_tmcurr;
    }

    //! return current timer of the current thread
    inline unsigned int get_tmcurr()
    {
        return get_tmcurr(omp_get_thread_num());
    }

    //! return amount of time spent in a timer
    inline double get(unsigned int tm) const
    {
        assert(tm < m_timers.size());
        return (m_timers[tm].tv_sec + m_timers[tm].tv_nsec / 1e9);
    }

    //! return sum over all timers
    inline double get_sum() const
    {
        double sum = 0;
        for (size_t i = 0; i < m_timers.size(); ++i)
            sum += get(i);
        return sum;
    }
};

//! Scoped timer keeper: changes timer array to new timer on construction and
//! back to old one on destruction.
class ScopedTimerKeeperMT
{
protected:
    //! reference to timer array
    TimerArrayMT& m_ta;

    //! previous timer identifier
    unsigned int m_tmprev;

public:
    //! construct and change timer to tm
    ScopedTimerKeeperMT(TimerArrayMT& ta, unsigned int tm)
        : m_ta(ta),
          m_tmprev(ta.get_tmcurr())
    {
        m_ta.change(tm);
    }

    //! change back timer to previous timer.
    ~ScopedTimerKeeperMT()
    {
        m_ta.change(m_tmprev);
    }
};

//! Dummy implementation of scoped TimerKeeper
class ScopedTimerKeeperDummy
{
public:
    template <typename Whatever>
    ScopedTimerKeeperDummy(Whatever /* ta */, unsigned int /* tm */)
    { }
};

#endif // !PSS_SRC_TOOLS_TIMER_ARRAY_HEADER

/******************************************************************************/
