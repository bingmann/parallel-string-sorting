/******************************************************************************
 * src/tools/timer.h
 *
 * Class to output statistics in a flexible text file as key=value pairs.
 *
 ******************************************************************************
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
 *****************************************************************************/

#ifndef TOOLS_TIMER_H
#define TOOLS_TIMER_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <unistd.h>

#include <assert.h>
#include <omp.h>

//! very simple class to measure runtime of function using clock_gettime.
template <clockid_t clk_id>
class ClockTimerBase
{
protected:
    //! time of start
    struct timespec     m_tstart;

    //! call gettime
    static inline void get_time(struct timespec& ts)
    {
        if (clock_gettime(clk_id, &ts)) {
            perror("Could not clock_gettime()");
        }
    }

public:
    //! Initialize and (usually) start clock
    ClockTimerBase(bool do_start = true)
    {
        if (do_start) start();
    }

    //! return the resolution of the clock used
    static inline double resolution()
    {
        struct timespec tp_res;

        if (clock_getres(clk_id, &tp_res)) {
            perror("Could not clock_getres()");
            return -1;
        }

        return tp_res.tv_sec + tp_res.tv_nsec / 1e9;
    }

    //! Start timing
    inline void start()
    {
        get_time(m_tstart);
    }

    //! Return time elapsed in seconds between start() and now.
    inline double elapsed() const
    {
        struct timespec now;
        get_time(now);

        return (now.tv_sec - m_tstart.tv_sec)
            +  (now.tv_nsec - m_tstart.tv_nsec) / 1e9;
    }
};

//! Most simple ClockTimer instance for easy measurements
typedef ClockTimerBase<CLOCK_MONOTONIC> ClockTimer;

//! Extended class to measure runtime of function using clock_gettime, with
//! start() - stop() semantics.
template <clockid_t clk_id>
class ClockIntervalBase : public ClockTimerBase<clk_id>
{
protected:
    //! type of super class
    typedef ClockTimerBase<clk_id> super_type;

    //! time of stop
    struct timespec     m_tstop;

public:
    //! Initialize, but do not start the clock
    ClockIntervalBase()
        : super_type(false)
    {
    }

    //! Stop timing
    inline void stop()
    {
        super_type::get_time(m_tstop);
    }

    //! Return delta in seconds between start() and stop().
    inline double delta() const
    {
        return (m_tstop.tv_sec - super_type::m_tstart.tv_sec)
            +  (m_tstop.tv_nsec - super_type::m_tstart.tv_nsec) / 1e9;
    }

    //! Retrun delta in seconds between start() and now.
    inline double delta_now() const
    {
        return super_type::delta();
    }
};

#endif // TOOLS_TIMER_H
