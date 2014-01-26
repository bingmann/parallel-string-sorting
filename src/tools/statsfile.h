/******************************************************************************
 * src/tools/statsfile.h
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

#ifndef TOOLS_STATSFILE_H
#define TOOLS_STATSFILE_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <unistd.h>

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

protected:

    class Entry
    {
    public:
        class StatsCache& m_sc;
        std::string m_key;

        Entry(StatsCache& sc, const std::string& key)
            : m_sc(sc), m_key(key)
        {
        }

        StatsCache& operator << (const std::string& v)
        {
            return m_sc.put(m_key, v);
        }

        template <typename ValueType>
        StatsCache& operator << (const ValueType& v)
        {
            std::ostringstream vstr;
            vstr << v;
            return operator << (vstr.str());
        }
    };

public:

    /// Clear all data
    void clear()
    {
        m_statsvec.clear();
    }

    // Append a (key,value) pair
    StatsCache& put(const std::string& k, const std::string& v)
    {
#pragma omp critical
        m_statsvec.push_back( strpair_type(k, v) );
        return *this;
    }

    // Append a (key,value) pair
    template <typename KeyType, typename ValueType>
    StatsCache& put(const KeyType& k, const ValueType& v)
    {
        std::ostringstream kstr, vstr;
        kstr << k; vstr << v;
        return put(kstr.str(), vstr.str());
    }

    // Append a (key,value) pair as ">> key << value"
    Entry operator >> (const std::string& k)
    {
        return Entry(*this, k);
    }

    // Append a (key,value) pair as ">> key << value"
    template <typename KeyType>
    Entry operator >> (const KeyType& k)
    {
        std::ostringstream kstr;
        kstr << k;
        return operator >> (kstr.str());
    }

    /// Return vector for inclusion in a StatsWriter.
    const statsvec_type& get_statsvec() const
    {
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
        m_line << '\t' << t;

        return *this;
    }

    // Append a value
    template <typename Type>
    StatsWriter& operator<< (const Type& t)
    {
        if (m_firstfield) {
            m_line << '=';
            m_firstfield = 0;
        }

        m_line << t;

        return *this;
    }

    // Append a (key,value) pair
    template <typename KeyType, typename ValueType>
    StatsWriter& put(const KeyType& k, const ValueType& v)
    {
        assert(m_firstfield == 0);
        m_line << '\t' << k << '=' << v;
        return *this;
    }

    // Append a stats map
    void append_stats(StatsCache& sc)
    {
        const StatsCache::statsvec_type& sm = sc.get_statsvec();

        for (StatsCache::statsvec_type::const_iterator si = sm.begin();
             si != sm.end(); ++si)
        {
            m_line << '\t' << si->first << '=' << si->second;
        }
    }
};

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

/// Class to measure different parts of a funciton by switching between
/// different aggregating timers. Immediately start with timer 0. Use enums in
/// your code to give the timer numbers names.
class TimerArray
{
private:

    //! clock time of last call
    struct timespec     m_tplast;

    //! currently running timer
    unsigned int        m_tmcurr;

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
        struct timespec     m_tplast;

        //! currently running timer per thread
        unsigned int        m_tmcurr;

        //! filler to put info into different cache lines
        unsigned char       m_filler[64 - sizeof(struct timespec) - sizeof(unsigned int)];
    } __attribute__((packed));

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
        m_timers[ ti.m_tmcurr ].tv_sec += (tpnow.tv_sec - ti.m_tplast.tv_sec);
#pragma omp atomic
        m_timers[ ti.m_tmcurr ].tv_nsec += tpnow.tv_nsec - ti.m_tplast.tv_nsec;

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
    TimerArrayMT&       m_ta;

    //! previous timer identifier
    unsigned int        m_tmprev;

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
    {
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
