/******************************************************************************
 * src/tools/statsfile.h
 *
 * Class to output statistics in a flexible text file as key=value pairs.
 *
 ******************************************************************************
 * Copyright (C) 2012 Timo Bingmann <tb@panthema.net>
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

#ifndef TOOLS_STATSFILE_H
#define TOOLS_STATSFILE_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <assert.h>
#include <omp.h>

/// Cache of key=value stats during run of algorithm
class StatsCache
{
public:

    typedef std::pair<std::string, std::string> strpair_type;
    typedef std::vector<strpair_type>  statsvec_type;

private:

    statsvec_type               m_statsvec;

    std::string                 m_thiskey;

    std::ostringstream          m_curritem;

public:

    /// Clear all data
    void clear()
    {
        m_statsvec.clear();
        m_thiskey.clear();
        m_curritem.str("");
    }

    /// Append a substring to the current or new key.
    template <typename Type>
    StatsCache& operator>> (const Type& t)
    {
        if (m_thiskey.size())
        {
            m_statsvec.push_back( strpair_type(m_thiskey, m_curritem.str()) );
            m_thiskey.clear();
            m_curritem.str("");
        }

        m_curritem << t;
        return *this;
    }

    /// Append the value to a previously >>-ed key.
    template <typename Type>
    StatsCache& operator<< (const Type& t)
    {
        if (m_thiskey.size() == 0)
        {
            m_thiskey = m_curritem.str();
            assert(m_thiskey.size() && "Key is empty!");
            m_curritem.str("");
        }

        m_curritem << t;
        return *this;
    }

    /// Return vector for inclusion in a StatsWriter.
    const statsvec_type& get_statsvec()
    {
        if (m_thiskey.size())
        {
            m_statsvec.push_back( strpair_type(m_thiskey, m_curritem.str()) );
            m_thiskey.clear();
            m_curritem.str("");
        }

        return m_statsvec;
    }
};

/// Simple writer of statistic files containing key=value pairs per line.
class StatsWriter
{
private:

    std::ofstream       m_out;

    unsigned int        m_firstfield;

    std::ostringstream  m_line;

public:

    StatsWriter(const char* filename)
    {
        m_out.open(filename, std::ios::app);

        m_line << "RESULT\t";

        // output date, time and hostname to m_line

        char datetime[64];
        time_t tnow = time(NULL);

        strftime(datetime,sizeof(datetime),"%Y-%m-%d %H:%M:%S", localtime(&tnow));
        m_line << "datetime=" << datetime;

        char hostname[128];
        gethostname(hostname, sizeof(hostname));

        m_line << "\thost=" << hostname;
    }

    ~StatsWriter()
    {
        m_out << m_line.str() << "\n";
        std::cout << m_line.str() << "\n";
    }

    // Append a key
    template <typename Type>
    StatsWriter& operator>> (const Type& t)
    {
        m_firstfield = 1;
        m_line << "\t" << t;

        return *this;
    }

    // Append a value
    template <typename Type>
    StatsWriter& operator<< (const Type& t)
    {
        if (m_firstfield) {
            m_line << "=";
            m_firstfield = 0;
        }

        m_line << t;

        return *this;
    }

    // Append a stats map
    void append_stats(StatsCache& sc)
    {
        const StatsCache::statsvec_type& sm = sc.get_statsvec();

        for (StatsCache::statsvec_type::const_iterator si = sm.begin();
             si != sm.end(); ++si)
        {
            m_line << "\t" << si->first << "=" << si->second;
        }
    }
};

class SizeLogger
{
private:

    /// log output file
    std::ofstream       m_logfile;

    /// begin timestamp of current averaging
    double              m_begintime;

    /// end timestamp of current averaging
    double              m_endtime;

    /// count of current averaging
    double              m_avgcount;

    /// sum of current averaging
    double              m_avgsum;

    /// timestamp function
    static inline double timestamp() {
        return omp_get_wtime();
    }

public:

    SizeLogger(const char* logname)
        : m_logfile(logname, std::ios::app),
          m_begintime(0)
    {
    }

    SizeLogger& operator<< (unsigned long value)
    {
        double thistime = timestamp();

        if (m_begintime == 0)
        {
            m_begintime = m_endtime = thistime;
            m_avgcount = 1;
            m_avgsum = value;
        }
        else if (m_begintime - thistime > 0.01 || m_avgcount >= 1000)
        {
            m_logfile << std::setprecision(16) << ((m_begintime + m_endtime) / 2.0) << " "
                      << std::setprecision(16) << (m_avgsum / m_avgcount) << " " << m_avgcount << "\n";

            m_begintime = m_endtime = thistime;
            m_avgcount = 1;
            m_avgsum = value;
        }
        else
        {
            m_endtime = thistime;
            m_avgcount++;
            m_avgsum += value;
        }

        return *this;
    }

    ~SizeLogger()
    {
        if (m_begintime != 0)
        {
            m_logfile << std::setprecision(16) << ((m_begintime + m_endtime) / 2.0) << " "
                      << std::setprecision(16) << (m_avgsum / m_avgcount) << " " << m_avgcount << "\n";
        }
    }
};

/// Very simple class to measure runtime of function using clock_gettime.
template <clockid_t clk_id>
class MeasureTime
{
private:

    struct timespec     m_tp1, m_tp2;

public:
    /// return the resolution of the clock used
    inline double resolution() const
    {
        struct timespec tp_res;

        if (clock_getres(clk_id, &tp_res)) {
            perror("Could not clock_getres()");
            return -1;
        }

        return tp_res.tv_sec + tp_res.tv_nsec / 1e9;
    }

    /// Start timing
    inline void start()
    {
        if (clock_gettime(clk_id, &m_tp1)) {
            perror("Could not clock_gettime()");
        }
    }

    /// End timing
    inline void stop()
    {
        if (clock_gettime(clk_id, &m_tp2)) {
            perror("Could not clock_gettime()");
        }
    }

    /// Return delta in seconds between start() and stop().
    inline double delta()
    {
        return (m_tp2.tv_sec - m_tp1.tv_sec)
            +  (m_tp2.tv_nsec - m_tp1.tv_nsec) / 1e9;
    }
};

/// Class to measure different parts of a funciton by switching between
/// different aggregating timers. Immediately start with timer 0. Use enums in
/// your code to give the timer numbers names.
class TimerArray
{
private:

    // clock time of last call
    struct timespec     m_tplast;

    // currently running timer
    unsigned int        m_tmcurr;

    // array of timers (usually preallocated)
    std::vector<struct timespec> m_tpvector;

public:

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
        m_tpvector[ m_tmcurr ].tv_sec += tpnow.tv_sec - m_tplast.tv_sec;
        m_tpvector[ m_tmcurr ].tv_nsec += tpnow.tv_nsec - m_tplast.tv_nsec;

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

/// Dummy class to replace TimerArray with no-ops.
class TimerArrayDummy
{
public:

    TimerArrayDummy(unsigned int /* timers */)
    {
    }

    /// clear all timers and start counting in timer 0.
    void clear()
    {
    }

    /// switch to other timer
    inline void change(unsigned int /* tm */)
    {
    }

    /// return amount of time spent in a timer
    inline double get(unsigned int /* tm */)
    {
        return 0;
    }
};

/// Class to read /proc/<pid>/smaps for memory usage
struct SMapsInfo
{
    size_t      size, rss, pss;
    size_t      referenced, anonymous, locked;

    void read()
    {
        size = rss = pss = 0;
        referenced = anonymous = locked = 0;

        std::ifstream in("/proc/self/smaps");
        std::string line;

        while ( std::getline(in, line) )
        {
            unsigned long mem_from, mem_to, mem_size;
            char mem_info[65];

            if (sscanf(line.c_str(), "%lx-%lx", &mem_from, &mem_to) == 2) {
                //std::cout << "new area " << mem_from << " - " << mem_to << std::endl;
            }
            else if (sscanf(line.c_str(), "%64[^:]: %lu kB", mem_info, &mem_size) == 2)
            {
                std::string info = mem_info;
                if (info == "Size") size += mem_size;
                else if (info == "Rss") rss += mem_size;
                else if (info == "Pss") pss += mem_size;
                else if (info == "Referenced") referenced += mem_size;
                else if (info == "Anonymous") anonymous += mem_size;
                else if (info == "Locked") locked += mem_size;
                else {
                    //std::cout << "type " << mem_info << " - " << mem_size << std::endl;
                }
            }
        }
    }
};

static inline size_t smaps_delta(const size_t& start, const size_t& end)
{
    if (end < start) return 0;
    else return end - start;
}

static inline void smaps_delta_stats(StatsCache& stats, const SMapsInfo& start, const SMapsInfo& end)
{
    stats >> "mem_size" << smaps_delta(start.size, end.size)
          >> "mem_rss" << smaps_delta(start.rss, end.rss)
          >> "mem_pss" << smaps_delta(start.pss, end.pss)
          >> "mem_referenced" << smaps_delta(start.referenced, end.referenced)
          >> "mem_anonymous" << smaps_delta(start.anonymous, end.anonymous)
          >> "mem_locked" << smaps_delta(start.locked, end.locked);
}

#endif // TOOLS_STATSFILE_H
