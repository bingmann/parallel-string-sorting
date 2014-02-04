/******************************************************************************
 * src/tools/stringptr.h
 *
 * StringPtr, StringPtrOut, and NoLcpCalc specializations. Encapsulate string
 * and shadow array pointers.
 *
 * Additionally: LcpStringPtr encapsulates string and lcp arrays, which may be
 * interleaved or separate.
 *
 ******************************************************************************
 * Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
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

#ifndef STRINGPTR_H_
#define STRINGPTR_H_

#include <assert.h>
#include "debug.h"
#include "numa.h"

namespace stringtools {

/// Objectified string array pointer and shadow pointer array for out-of-place
/// swapping of pointers.
template <bool NoLcpCalc>
class StringPtrBase
{
protected:
    //! strings (front) and temporary shadow (back) array
    string      *m_active, *m_shadow;

    //! length of subarray
    size_t      m_size;

    //! false if m_active is original, true if m_shadow is original
    bool        m_flipped;

public:
    /// constructor specifying all attributes
    inline StringPtrBase(string* original, string* shadow = NULL, size_t size = 0, bool flipped = false)
        : m_active(original), m_shadow(shadow), m_size(size), m_flipped(flipped)
    {
    }

    /// true if flipped to back array
    inline bool flipped() const { return m_flipped; }

    /// return currently active array
    inline string* active() const { return m_active; }

    /// return current shadow array
    inline string* shadow() const { return m_shadow; }

    //! return valid length
    inline size_t size() const { return m_size; }

    //! ostream-able
    friend inline std::ostream& operator << (std::ostream& os, const StringPtrBase& sp)
    {
        return os << '(' << sp.active() << '/' << sp.shadow() << '|' << sp.flipped() << ':' << sp.size() << ')';
    }

    /// Advance (both) pointers by given offset, return sub-array
    inline StringPtrBase sub(size_t offset, size_t size) const
    {
        assert(offset + size <= m_size);
        return StringPtrBase(m_active + offset, m_shadow + offset, size, m_flipped);
    }

    /// construct a StringPtrBase object specifying a sub-array with flipping to
    /// other array.
    inline StringPtrBase flip(size_t offset, size_t size) const
    {
        assert(offset + size <= m_size);
        return StringPtrBase(m_shadow + offset, m_active + offset, size, !m_flipped);
    }

    /// Return the original for this StringPtr for LCP calculation
    inline StringPtrBase original() const
    {
        return m_flipped ? flip(0, m_size) : *this;
    }

    /// return subarray pointer to n strings in original array, might copy from
    /// shadow before returning.
    inline StringPtrBase copy_back() const
    {
        if (!m_flipped) {
            return *this;
        }
        else {
            memcpy(m_shadow, m_active, m_size * sizeof(string));
            return flip(0, m_size);
        }
    }

    //! check sorted order of strings
    inline bool check() const
    {
        for (size_t i = 1; i < m_size; ++i)
            assert(scmp(out(i-1), out(i)) <= 0);
        return true;
    }

    /// Return i-th string pointer from m_active
    inline string& str(size_t i) const
    {
        assert(i < m_size);
        return m_active[i];
    }

    /// return reference to the i-th lcp
    inline uintptr_t& lcp(size_t i) const
    {
        if (NoLcpCalc) assert(0);

        assert(!m_flipped);
        assert(i < m_size);
        return ((uintptr_t*)m_shadow)[i];
    }

    /// set the i-th lcp to v and check its value
    inline void set_lcp(size_t i, const uintptr_t& v) const
    {
        if (NoLcpCalc) return;

        assert(i > 0);
        assert(i < m_size);
        assert(v == calc_lcp(out(i-1), out(i)));

        lcp(i) = v;
    }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    inline void fill_lcp(uintptr_t v)
    {
        if (NoLcpCalc) return;

        for (size_t i = 1; i < m_size; ++i)
            set_lcp(i, v);
    }

    /// Return pointer to LCP array
    inline uintptr_t* lcparray() const
    {
        if (NoLcpCalc) assert(0);

        assert(!m_flipped);
        return (uintptr_t*)m_shadow;
    }

    /// Return the output string array
    inline string* output() const
    {
        assert(!m_flipped); // m_active is original/output array
        return m_active;
    }

    /// Return i-th output string pointer from m_active / output()
    inline string& out(size_t i) const
    {
        assert(!m_flipped); // m_active is original/output array
        return str(i);
    }
};

typedef StringPtrBase<false> StringPtr;
typedef StringPtrBase<true> StringPtrNoLcpCalc;

/// Objectified string array pointer and shadow pointer array for out-of-place
/// swapping of pointers. With separate output array. Use class not derived
/// from StringPtrBase!
template <bool NoLcpCalc>
class StringPtrOutBase
{
protected:
    //! encapsuled StringPtrBase type
    typedef StringPtrBase<NoLcpCalc> base_type;

    //! encapsuled StringPtrBase
    base_type   sp;

    //! output string array
    string      *m_output;

public:
    /// constructor specifying all attributes
    inline StringPtrOutBase(string* original, string* shadow = NULL, string* output = NULL,
                            size_t size = 0, bool flipped = false)
        : sp(original, shadow, size, flipped),
          m_output(output)
    { }

    /// true if flipped to back array
    inline bool flipped() const { return sp.flipped(); }

    /// return currently active array
    inline string* active() const { return sp.active(); }

    /// return current shadow array
    inline string* shadow() const { return sp.shadow(); }

    //! return valid length
    inline size_t size() const { return sp.size(); }

    //! ostream-able
    friend inline std::ostream& operator << (std::ostream& os, const StringPtrOutBase& sp)
    {
        return os << '(' << sp.active() << '/' << sp.shadow() << '/' << sp.output()
                  << '|' << sp.flipped() << ':' << sp.size() << ')';
    }

    /// Advance (both) pointers by given offset, return sub-array
    inline StringPtrOutBase sub(size_t offset, size_t size) const
    {
        assert(offset + size <= sp.size());
        return StringPtrOutBase(active() + offset, shadow() + offset, m_output + offset,
                                size, flipped());
    }

    /// construct a StringPtrOutBase object specifying a sub-array with flipping to
    /// other array.
    inline StringPtrOutBase flip(size_t offset, size_t size) const
    {
        assert(offset + size <= sp.size());
        return StringPtrOutBase(shadow() + offset, active() + offset, m_output + offset,
                                size, !flipped());
    }

    /// Return the original for this StringPtrOut for LCP calculation
    inline StringPtrOutBase original() const
    {
        return flipped() ? flip(0, size()) : *this;
    }

    /// return subarray pointer to n strings in original array, might copy from
    /// shadow before returning.
    inline StringPtrOutBase copy_back() const
    {
        memcpy(m_output, active(), size() * sizeof(string));
        return original();
    }

    //! check sorted order of strings
    inline bool check() const
    {
        for (size_t i = 1; i < size(); ++i)
            assert(scmp(out(i-1), out(i)) <= 0);
        return true;
    }

    /// Return i-th string pointer from m_active
    inline string& str(size_t i) const { return sp.str(i); }

    /// return reference to the i-th lcp
    inline uintptr_t& lcp(size_t i) const { return sp.lcp(i); }

    /// set the i-th lcp to v and check its value
    inline void set_lcp(size_t i, const uintptr_t& v) const
    {
        if (NoLcpCalc) return;

        assert(i > 0);
        assert(i < size());
        assert(v == calc_lcp(out(i-1), out(i)));

        lcp(i) = v;
    }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    inline void fill_lcp(uintptr_t v)
    {
        if (NoLcpCalc) return;

        for (size_t i = 1; i < size(); ++i)
            set_lcp(i, v);
    }

    /// Return pointer to LCP array
    inline uintptr_t* lcparray() const { return sp.lcparray(); }

    /// Return the output string array
    inline string* output() const
    {
        return m_output;
    }

    /// Return i-th output string pointer from m_active / output()
    inline string& out(size_t i) const
    {
        assert(i < size());
        return m_output[i];
    }
};

typedef StringPtrOutBase<false> StringPtrOut;
typedef StringPtrOutBase<true> StringPtrOutNoLcpCalc;

//////////////////////////////////////////////////////////////////////

typedef uintptr_t lcp_t;

struct LcpCacheStringPtr
{
public:
    string * strings;
    lcp_t* lcps;
    char* cachedChars;
    size_t size;

public:
    LcpCacheStringPtr()
        : strings(NULL), lcps(NULL), cachedChars(NULL), size(0)
    {
    }

    LcpCacheStringPtr(string* _strings, lcp_t* _lcps, char* _cachedChars, size_t _size)
           : strings(_strings), lcps(_lcps), cachedChars(_cachedChars), size(_size)
    {
    }

    inline bool empty() const
    {
        return (size == 0);
    }

    inline void
    set(string s, lcp_t lcp) const
    {
        assert(size > 0);
        *strings = s;
        *lcps = lcp;
        *cachedChars = s[lcp];
    }

    inline void
    set(string s, lcp_t lcp, char cachedCharacter) const
    {
       assert(size > 0);
       *strings = s;
       *lcps = lcp;
       *cachedChars = cachedCharacter;
    }

    inline void
    setFirst(LcpCacheStringPtr& ptr)
    {
        *strings = *ptr.strings;
        *lcps = *ptr.lcps;
        *cachedChars = *ptr.cachedChars;
    }

    inline string&
    str() const
    {
        assert(size > 0);
        return *strings;
    }

    inline lcp_t&
    lcp() const
    {
        assert(size > 0);
        return *lcps;
    }

    inline char&
    cache() const
    {
        assert(size > 0);
        return *cachedChars;
    }

    inline void
    copyFrom(LcpCacheStringPtr& other, size_t length) const
    {
        memcpy(strings, other.strings, length * sizeof(string));
        memcpy(lcps, other.lcps, length * sizeof(lcp_t));
        memcpy(cachedChars, other.cachedChars, length * sizeof(char));
    }

    inline void
    copyStringsTo(string* destination, size_t length) const
    {
        memcpy(destination, strings, length * sizeof(string));
    }

    inline void
    calculateCache() const
    {
        for(unsigned i = 0; i < size; ++i){
            cachedChars[i] = strings[i][lcps[i]];
        }
    }

    inline void
    allocateNumaMemory(int numaNode, size_t length)
    {
#if 0
        strings     = new string[length];
        lcps        = new lcp_t[length];
        cachedChars = new char[length];

        numa_tonode_memory(strings,     length * sizeof(string), numaNode);
        numa_tonode_memory(lcps,        length * sizeof(lcp_t), numaNode);
        numa_tonode_memory(cachedChars, length * sizeof(char), numaNode);
#else
        strings =     (string*) numa_alloc_onnode(length * sizeof(string), numaNode);
        lcps =        (lcp_t*)  numa_alloc_onnode(length * sizeof(lcp_t),  numaNode);
        cachedChars = (char*)   numa_alloc_onnode(length * sizeof(char),   numaNode);
#endif
        size = length;
    }

    inline void
    freeNumaMemory()
    {
        numa_free(strings, size * sizeof(string));
        numa_free(lcps, size * sizeof(string));
        numa_free(cachedChars, size * sizeof(char));
    }

    // preincrement
    inline LcpCacheStringPtr&
    operator++()
    {
        ++strings;
        ++lcps;
        ++cachedChars;
        --size;
        return *this;
    }

    //! return sub-array of (string,lcp) with offset and size
    inline LcpCacheStringPtr
    sub(size_t offset, size_t n) const
    {
        assert(offset + n <= size);
        return LcpCacheStringPtr(strings + offset, lcps + offset, cachedChars + offset, n);
    }

    //! return empty end array.
    inline LcpCacheStringPtr
    end() const
    {
        return sub(size, 0);
    }

    inline size_t
    operator-(const LcpCacheStringPtr& rhs) const
    {
        return strings - rhs.strings;
    }

    inline bool
    operator<(const LcpCacheStringPtr& rhs) const
    {
        return strings < rhs.strings;
    }
};


struct LcpStringPtr
{
public:
    string * strings;
    lcp_t* lcps;
    size_t size;

public:
    LcpStringPtr()
        : strings(NULL), lcps(NULL), size(0)
    {
    }

    LcpStringPtr(string* _strings, lcp_t* _lcps, size_t _size)
        : strings(_strings), lcps(_lcps), size(_size)
    {
    }

    LcpStringPtr(const LcpCacheStringPtr& ptr)
         : strings(ptr.strings), lcps(ptr.lcps), size(ptr.size)
    {
    }

    LcpStringPtr(const LcpCacheStringPtr& ptr, size_t _size)
            : strings(ptr.strings), lcps(ptr.lcps), size(_size)
    {
    }

    inline bool empty() const
    {
        return (size == 0);
    }

    inline void
    set(string s, lcp_t lcp) const
    {
        assert(size > 0);
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
        assert(size > 0);
        return *strings;
    }

    inline lcp_t&
    lcp() const
    {
        assert(size > 0);
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
        assert(position < size);
        lcps[position] = value;
    }

    // preincrement
    inline LcpStringPtr&
    operator++()
    {
        ++strings;
        ++lcps;
        --size;
        return *this;
    }

    //! return sub-array of (string,lcp) with offset and size
    inline LcpStringPtr
    sub(size_t offset, size_t n) const
    {
        assert(offset + n <= size);
        return LcpStringPtr(strings + offset, lcps + offset, n);
    }

    //! return empty end array.
    inline LcpStringPtr
    end() const
    {
        return sub(size, 0);
    }

    inline size_t
    operator-(const LcpStringPtr& rhs) const
    {
        return strings - rhs.strings;
    }

    inline bool
    operator<(const LcpStringPtr& rhs) const
    {
        return strings < rhs.strings;
    }
};

/// verify LCP array against sorted string array by scanning LCPs
template<bool checkCache>
static inline bool
verify_lcp_cache(string* strings, lcp_t* lcps, char* cache, size_t n, lcp_t expectedFirstLcp)
{
    bool allValid = true;

    if (lcps[0] != expectedFirstLcp)
    {
        std::cout << "lcp[0] = " << lcps[0] << " excepted " << expectedFirstLcp << std::endl;
        allValid = false;
    }
    if(checkCache && *cache != strings[0][lcps[0]]){
        std::cout << "cache[0] = " << cache[0] << " excepted " <<  strings[0][lcps[0]] << std::endl;
        allValid = false;
    }

    for (size_t i = 1; i < n; ++i)
    {
        string s1 = strings[i-1], s2 = strings[i];
        size_t h = calc_lcp(s1, s2);

        if (h != lcps[i])
        {
            std::cout << "lcp[" << i << "] = " << lcps[i] << " excepted " << h << std::endl;
            allValid = false;
        }
        if (checkCache && cache[i] != s2[lcps[i]])
        {
            std::cout << "cache[" << i << "] = " << cache[i] << " excepted " << s2[lcps[i]] << std::endl;
            allValid = false;
        }
    }

    if (allValid)
        std::cout << "All LCPs and cache values valid!" << std::endl;
    else
        std::cout << "Found invalid LCPS and/or cache values!" << std::endl;

    return allValid;
}

static inline bool
verify_lcp(string* strings, lcp_t* lcps, size_t n, lcp_t expectedFirstLcp){
    return verify_lcp_cache<false>(strings, lcps, NULL, n, expectedFirstLcp);
}

static inline bool
verify_lcp_cache(string* strings, lcp_t* lcps, char* cache, size_t n, lcp_t expectedFirstLcp)
{
    return verify_lcp_cache<true>(strings, lcps, cache, n, expectedFirstLcp);
}

} // namespace stringtools

#endif // STRINGPTR_H_
