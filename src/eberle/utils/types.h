/******************************************************************************
 * src/eberle/utils/types.h
 *
 * Type definitions used by some src/eberle/ algorithms.
 *
 ******************************************************************************
 * Copyright (C) 2013-2014 Andreas Eberle <email@andreas-eberle.com>
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

#ifndef TYPES_H_
#define TYPES_H_

#include "../../tools/stringtools.h"

namespace types
{

using stringtools::string;
typedef uintptr_t lcp_t;

struct LcpStringPtr
{
    string * strings;
    lcp_t* lcps;

public:
    LcpStringPtr() :
            strings(NULL), lcps(NULL)
    {
    }

    LcpStringPtr(string* strings, lcp_t* lcps) :
            strings(strings), lcps(lcps)
    {
    }

    inline void
    set(string s, lcp_t lcp) const
    {
        *strings = s;
        *lcps = lcp;
    }

    inline void
    setFirst(LcpStringPtr& ptr)
    {
        *strings = *ptr.strings;
        *lcps = *ptr.lcps;
    }

    inline string&
    str() const
    {
        return *strings;
    }

    inline lcp_t&
    lcp() const
    {
        return *lcps;
    }

    inline void
    copyFrom(LcpStringPtr& other, size_t length) const
    {
        memcpy(strings, other.strings, length * sizeof(string));
        memcpy(lcps, other.lcps, length * sizeof(lcp_t));
    }

    inline void
    copyStringsTo(string* destination, size_t length) const
    {
        memcpy(destination, strings, length * sizeof(string));
    }

    inline void
    setLcp(size_t position, lcp_t value) const
    {
        lcps[position] = value;
    }

    // preincrement
    inline LcpStringPtr&
    operator++()
    {
        ++strings;
        ++lcps;
        return *this;
    }

    friend inline LcpStringPtr
    operator+(const LcpStringPtr& ptr, size_t delta);

    friend inline size_t
    operator-(const LcpStringPtr& lhs, const LcpStringPtr& rhs);

    inline bool
    operator<(const LcpStringPtr& rhs)
    {
        return strings < rhs.strings;
    }
};

inline size_t
operator-(const LcpStringPtr& lhs, const LcpStringPtr& rhs)
{
    return lhs.strings - rhs.strings;
}

inline LcpStringPtr
operator+(const LcpStringPtr& ptr, size_t delta)
{
    return LcpStringPtr(ptr.strings + delta, ptr.lcps + delta);
}

} // namespace types

#endif // TYPES_H_
