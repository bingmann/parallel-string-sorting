/******************************************************************************
 * src/tools/agglogger.h
 *
 * Class to write aggregated measurement values into gnuplottable text format.
 *
 * When logging measured amounts in experimetal programs, one often was to keep
 * track of a continously changing value. But that the same time one does not
 * want to write out the value each time it changes, because this would result
 * in a huge log file and also severely impact the experiemtal program's
 * performance.
 *
 * This is where this aggregating logger comes into use: it will aggregate
 * multiple samples and output only one value. The value is written either when
 * the absolute value changed sufficiently or if a configured amount of time
 * has elapsed since the last write. The outputted sample value is either the
 * average, the minimum or maximum of all measured values since the last
 * output.
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

#ifndef TOOLS_AGGLOGGER_H
#define TOOLS_AGGLOGGER_H

#include <fstream>
#include <iomanip>
#include <limits>

#include <omp.h>

template <typename ValueType>
class AggregateLogger
{
public:
    //! type of measured values
    typedef ValueType value_type;

protected:
    //! timestamp function
    static inline double timestamp()
    {
        return omp_get_wtime();
    }

protected:
    /**
     * Output aggregate over all sampled values since the last written value.
     */
    template <typename AggFunctor>
    class TemplateLogger
    {
    protected:

        //! log output file
        std::ofstream       m_logfile;

        //! begin timestamp of current time range
        double              m_begintime;

        //! end timestamp of current time range
        double              m_endtime;

        //! timestamp defined as zero time, output is relative to this
        double              m_zerotime;

        //! count of samples in current time range
        size_t              m_count;

        //! current aggregate value
        value_type          m_aggregate;

        //! maximum duration between two outputted values
        double              m_max_interval;

        //! maximum number of values over which an output is aggregated
        size_t              m_max_count;

        //! output aggregated value
        inline void output()
        {
            m_logfile << std::setprecision(16)
                      << ((m_begintime + m_endtime) / 2.0) - m_zerotime << " "
                      << AggFunctor::output(m_aggregate, m_count) << " "
                      << m_count << std::endl;
        }

    public:

        //! initialize logger with logname, and maximum interval/count between
        //! writes
        TemplateLogger(const char* logname, double max_interval = 0.01,
                       size_t max_count = 1000, bool append = false)
            : m_logfile(logname, append ? std::ios::app : std::ios::out),
              m_begintime(0),
              m_zerotime(timestamp()),
              m_max_interval(max_interval),
              m_max_count(max_count)
        {
        }

        //! Define current timestamp as zero.
        TemplateLogger& start()
        {
            m_zerotime = timestamp();
            return *this;
        }

        //! Put a value into the logger
        TemplateLogger& operator << (const value_type& value)
        {
            double thistime = timestamp();

            if (m_begintime == 0) // first value
            {
                m_begintime = m_endtime = thistime;
                m_count = 1;
                m_aggregate = AggFunctor::aggregate(AggFunctor::initial(), value);
            }
            else if (thistime - m_begintime > m_max_interval ||
                     m_count >= m_max_count)
            {
                // output an average value
                output();

                m_begintime = m_endtime = thistime;
                m_count = 1;
                m_aggregate = AggFunctor::aggregate(AggFunctor::initial(), value);
            }
            else
            {
                // add to running average
                m_endtime = thistime;
                m_count++;
                m_aggregate = AggFunctor::aggregate(m_aggregate, value);
            }

            return *this;
        }

        ~TemplateLogger()
        {
            if (m_begintime != 0)
                output();
        }
    };

    struct MaxLoggerFunctor
    {
        static inline value_type initial()
        {
            return std::numeric_limits<value_type>::min();
        }

        static inline value_type aggregate(const value_type& prev, const value_type& sample)
        {
            return std::max(prev, sample);
        }

        static inline const value_type& output(const value_type& aggregate, size_t /* count */)
        {
            return aggregate;
        }
    };

    struct MinLoggerFunctor
    {
        static inline value_type initial()
        {
            return std::numeric_limits<value_type>::max();
        }

        static inline value_type aggregate(const value_type& prev, const value_type& sample)
        {
            return std::min(prev, sample);
        }

        static inline const value_type& output(const value_type& aggregate, size_t /* count */)
        {
            return aggregate;
        }
    };

    struct AvgLoggerFunctor
    {
        static inline value_type initial()
        {
            return 0;
        }

        static inline value_type aggregate(const value_type& prev, const value_type& sample)
        {
            return prev + sample;
        }

        static inline double output(const value_type& aggregate, size_t count)
        {
            return aggregate / (double)count;
        }
    };

public:

    typedef TemplateLogger<MaxLoggerFunctor> MaximumLogger;
    typedef TemplateLogger<MaxLoggerFunctor> MaxLogger;

    typedef TemplateLogger<MinLoggerFunctor> MinimumLogger;
    typedef TemplateLogger<MinLoggerFunctor> MinLogger;

    typedef TemplateLogger<AvgLoggerFunctor> AverageLogger;
    typedef TemplateLogger<AvgLoggerFunctor> AvgLogger;

protected:
    //! Thread-safe facade template class for Loggers
    template <typename BaseLogger>
    class LockingLogger : protected BaseLogger
    {
    public:

        LockingLogger(const char* logname, double max_interval = 0.01,
                      size_t max_count = 1000, bool append = false)
            : BaseLogger(logname, max_interval, max_count, append)
        {
        }

        LockingLogger& start()
        {
            BaseLogger::start();
            return *this;
        }

        LockingLogger& operator << (const value_type& value)
        {
#pragma omp critical
            BaseLogger::operator << (value);
            return *this;
        }
    };

public:
    typedef LockingLogger<MinLogger> LockingMinLogger;
    typedef LockingLogger<MinLogger> LockingMinimumLogger;

    typedef LockingLogger<MaxLogger> LockingMaxLogger;
    typedef LockingLogger<MaxLogger> LockingMaximumLogger;

    typedef LockingLogger<AvgLogger> LockingAvgLogger;
    typedef LockingLogger<AvgLogger> LockingAverageLogger;

protected:
    //! Class to replace SizeLogger with no-ops
    class DummyLogger
    {
    public:

        DummyLogger(const char* /* logname */, double /* max_interval */ = 0,
                    double /* max_count */ = 0, bool /* append */ = false)
        { }

        DummyLogger& start()
        {
            return *this;
        }

        DummyLogger& operator<< (const value_type& /* value */)
        {
            return *this;
        }
    };

}; // class AggregateLogger

#endif // TOOLS_AGGLOGGER_H
