/******************************************************************************
 * src/tools/memprofile.h
 *
 * Class to write the datafile for a memory profile plot.
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

#ifndef _MEM_PROFILE_H_
#define _MEM_PROFILE_H_

#include <stdio.h>

class MemProfile
{
protected:
    double      m_base_ts;
    size_t      m_base_use;
    uint8_t*    m_stack_base;

    double      m_prev_ts;
    size_t      m_prev_use;
    size_t      m_max;

    FILE*       m_file;
    const char* m_funcname;

protected:

    template <typename Type>
    static inline Type absdiff(const Type& a, const Type& b)
    {
        return (a < b) ? (b - a) : (a - b);
    }

    void output(double ts, size_t memcurr)
    {
        fprintf(m_file, "RESULT\tfunc=%s\tts=%f\tmem=%lu\n", m_funcname, ts - m_base_ts, memcurr);
    }

    inline void callback(size_t memcurr)
    {
        size_t use = (memcurr > m_base_use) ? (memcurr - m_base_use) : 0;

        if ((uint8_t*)&use < m_stack_base) // add stack usage
            use += m_stack_base - (uint8_t*)&use;

        double ts = omp_get_wtime();
        if (m_max < use) m_max = use;

        if (ts - m_prev_ts > 0.01 || absdiff(use, m_prev_use) > 16*1024 )
        {
            output(ts, m_max);
            m_max = 0;
            m_prev_ts = ts;
            m_prev_use = use;
        }
    }

    static void static_callback(void* cookie, size_t memcurr)
    {
        return static_cast<MemProfile*>(cookie)->callback(memcurr);
    }

public:

    MemProfile(const char* funcname, const char* filepath)
        : m_funcname(funcname)
    {
        uint8_t stack; m_stack_base = &stack;
        m_file = fopen(filepath, "a");
        malloc_count_set_callback(MemProfile::static_callback, this);
        clear();
    }

    ~MemProfile()
    {
        malloc_count_set_callback(NULL, NULL);
        fclose(m_file);
    }

    void clear()
    {
        m_base_ts = omp_get_wtime();
        m_base_use = malloc_count_current();
        m_prev_ts = 0;
        m_prev_use = 0;
        m_max = 0;
    }

    void finish()
    {
        m_prev_ts = 0;
        m_prev_use = 0;
        callback( malloc_count_current() );
        malloc_count_set_callback(NULL, NULL);
    }
};

#endif // _MEM_PROFILE_H_
