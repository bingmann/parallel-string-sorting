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
