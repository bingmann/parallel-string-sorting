/******************************************************************************
 * examples/stats_writer.h
 *
 * Class to collect and output statistics as key=value pairs.
 *
 * The usual method to use this stats collector is to make a global object
 * stats_writer g_stats, which is filled by algorithms using sequences as.
 *
 * g_stats >> "key" << "value " << 42;
 *
 * After the program was run, the stats are formatted as a RESULT line using
 * get(), which can be outputted to a file or stdout.
 *
 ******************************************************************************
 * Copyright (C) 2012-2014 Timo Bingmann <tb@panthema.net>
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

#ifndef SQLPLOTS_STATS_WRITER_H
#define SQLPLOTS_STATS_WRITER_H

#include <string>
#include <sstream>

#include <unistd.h>
#include <time.h>

/*!
 * Collect key=value pairs, which are given by operator >> and operator <<
 * sequences. After all stats are set, the final output line can be fetched.
 */
class stats_writer
{
protected:

    //! All collected key=value values.
    std::ostringstream  m_line;

    //! An internal class to collect key=value pairs as a sequence of >> and <<
    //! operator calls.
    class entry
    {
    protected:
        //! Reference to parent stats writer object
        class stats_writer& m_sw;

        //! Collected key and value string items
        std::string m_key, m_value;

    public:

        //! Start entry collection for the given key
        entry(stats_writer& sw, const std::string& key)
            : m_sw(sw), m_key(key)
        { }

        //! Collect more information about the value
        entry& operator << (const std::string& v)
        {
            m_value += v;
            return *this;
        }

        //! Collect more information about the value
        template <typename ValueType>
        entry& operator << (const ValueType& v)
        {
            std::ostringstream vstr;
            vstr << v;
            return operator << (vstr.str());
        }

        //! Start another key= entry via parent
        template <typename ValueType>
        entry operator >> (const ValueType& v)
        {
            // put key=value into writer before returning next entry
            m_sw.put(m_key, m_value);
            m_key.clear(); m_value.clear();
            return m_sw.operator >> (v);
        }

        //! Output key=value to stats writer for the last entry
        ~entry()
        {
            if (m_key.size() || m_value.size())
                m_sw.put(m_key, m_value);
        }
    };

public:

    //! Clear all data in the stats writer.
    void clear()
    {
        m_line.str("");
    }

    //! Append a (key,value) pair as ">> key << value << more"
    entry operator >> (const std::string& k)
    {
        return entry(*this, k);
    }

    //! Append a (key,value) pair as ">> key << value << more"
    template <typename KeyType>
    entry operator >> (const KeyType& k)
    {
        std::ostringstream kstr;
        kstr << k;
        return operator >> (kstr.str());
    }

    //! Append a (key,value) pair as strings
    stats_writer& put(const std::string& k, const std::string& v)
    {
#if _OPENMP
#pragma omp critical
#endif
        m_line << '\t' << k << '=' << v;
        return *this;
    }

    //! Append a (key,value) pair with automatic conversion to strings
    template <typename KeyType, typename ValueType>
    stats_writer& put(const KeyType& k, const ValueType& v)
    {
        std::ostringstream kstr, vstr;
        kstr << k; vstr << v;
        return put(kstr.str(), vstr.str());
    }

    //! Return RESULT string for outputting.
    std::string get() const
    {
        std::ostringstream out;
        out << "RESULT";

        // output date, time and hostname

        char datetime[64];
        time_t tnow = time(NULL);

        strftime(datetime,sizeof(datetime),"%Y-%m-%d %H:%M:%S", localtime(&tnow));
        out << "\tdatetime=" << datetime;

        char hostname[128];
        gethostname(hostname, sizeof(hostname));

        out << "\thost=" << hostname;

        // output collected key=values

        out << m_line.str();

        return out.str();
    }

    //! Return RESULT string for outputting.
    friend std::ostream& operator << (std::ostream& os, const stats_writer& sw)
    {
        return os << sw.get();
    }
};

#endif // SQLPLOTS_STATS_WRITER_H
