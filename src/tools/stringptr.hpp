/*******************************************************************************
 * src/tools/stringptr.hpp
 *
 * StringPtr, StringPtrOut, and NoLcpCalc specializations. Encapsulate string
 * and shadow array pointers.
 *
 * Additionally: LcpStringPtr encapsulates string and lcp arrays, which may be
 * interleaved or separate.
 *
 * StringLcpPtr               -> (string,lcp,size)
 * StringLcpCachePtr          -> (string,lcp,charcache,size)
 *
 * StringShadowPtr            -> (string,shadow,size,flip)
 * StringShadowOutPtr         -> (string,shadow,output,size,flip)
 * StringShadowLcpPtr         -> (string,shadow=lcp,size,flip)
 * StringShadowLcpOutPtr      -> (string,shadow=lcp,output,size,flip)
 * StringShadowLcpCacheOutPtr -> (string,shadow=lcp,charcache,output,size,flip)
 *
 *******************************************************************************
 * Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
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
 ******************************************************************************/

#ifndef PSS_SRC_TOOLS_STRINGPTR_HEADER
#define PSS_SRC_TOOLS_STRINGPTR_HEADER

#include <cassert>
#include <stdint.h>
#include "debug.hpp"
#include <numa.h>
#include "stringset.hpp"

namespace stringtools {

typedef uintptr_t lcp_t;

////////////////////////////////////////////////////////////////////////////////

class LcpCacheStringPtr;

class LcpStringPtr
{
public:
    string* strings;
    lcp_t* lcps;
    size_t size;

public:
    LcpStringPtr()
        : strings(NULL), lcps(NULL), size(0)
    { }

    LcpStringPtr(string* _strings, lcp_t* _lcps, size_t _size)
        : strings(_strings), lcps(_lcps), size(_size)
    { }

    bool empty() const
    {
        return (size == 0);
    }

    void setFirst(string s, lcp_t lcp) const
    {
        assert(size > 0);
        *strings = s;
        *lcps = lcp;
    }

    void setFirst(LcpStringPtr& ptr)
    {
        *strings = *ptr.strings;
        *lcps = *ptr.lcps;
    }

    string & firstString() const
    {
        assert(size > 0);
        return *strings;
    }

    lcp_t & firstLcp() const
    {
        assert(size > 0);
        return *lcps;
    }

    void copyFrom(LcpStringPtr& other, size_t length) const
    {
        memcpy(strings, other.strings, length * sizeof(string));
        memcpy(lcps, other.lcps, length * sizeof(lcp_t));
    }

    void copyStringsTo(string* destination, size_t length) const
    {
        memcpy(destination, strings, length * sizeof(string));
    }

    void setLcp(size_t position, lcp_t value) const
    {
        assert(position < size);
        lcps[position] = value;
    }

    // preincrement
    LcpStringPtr& operator ++ ()
    {
        ++strings;
        ++lcps;
        --size;
        return *this;
    }

    //! return sub-array of (string,lcp) with offset and size
    LcpStringPtr sub(size_t offset, size_t n) const
    {
        assert(offset + n <= size);
        return LcpStringPtr(strings + offset, lcps + offset, n);
    }

    //! return empty end array.
    LcpStringPtr end() const
    {
        return sub(size, 0);
    }

    size_t operator - (const LcpStringPtr& rhs) const
    {
        return strings - rhs.strings;
    }

    bool operator < (const LcpStringPtr& rhs) const
    {
        return strings < rhs.strings;
    }
};

////////////////////////////////////////////////////////////////////////////////

class LcpCacheStringPtr
{
public:
    string* strings;
    lcp_t* lcps;
    char_type* cachedChars;
    size_t size;

public:
    LcpCacheStringPtr()
        : strings(NULL), lcps(NULL), cachedChars(NULL), size(0)
    { }

    LcpCacheStringPtr(string* _strings, lcp_t* _lcps, char_type* _cachedChars, size_t _size)
        : strings(_strings), lcps(_lcps), cachedChars(_cachedChars), size(_size)
    { }

    bool empty() const
    {
        return (size == 0);
    }

    void setFirst(string s, lcp_t lcp) const
    {
        assert(size > 0);
        *strings = s;
        *lcps = lcp;
        *cachedChars = s[lcp];
    }

    void setFirst(string s, lcp_t lcp, char cachedCharacter) const
    {
        assert(size > 0);
        *strings = s;
        *lcps = lcp;
        *cachedChars = cachedCharacter;
    }

    void setFirst(LcpCacheStringPtr& ptr)
    {
        *strings = *ptr.strings;
        *lcps = *ptr.lcps;
        *cachedChars = *ptr.cachedChars;
    }

    string & firstString() const
    {
        assert(size > 0);
        return *strings;
    }

    lcp_t & firstLcp() const
    {
        assert(size > 0);
        return *lcps;
    }

    char_type & firstCached() const
    {
        assert(size > 0);
        return *cachedChars;
    }

    void copyFrom(LcpCacheStringPtr& other, size_t length) const
    {
        memcpy(strings, other.strings, length * sizeof(string));
        memcpy(lcps, other.lcps, length * sizeof(lcp_t));
        memcpy(cachedChars, other.cachedChars, length * sizeof(char));
    }

    void copyStringsTo(string* destination, size_t length) const
    {
        memcpy(destination, strings, length * sizeof(string));
    }

    void calculateCache() const
    {
        for (unsigned i = 0; i < size; ++i) {
            cachedChars[i] = strings[i][lcps[i]];
        }
    }

    void allocateNumaMemory(int numaNode, size_t length)
    {
#if 0
        strings = new string[length];
        lcps = new lcp_t[length];
        cachedChars = new char_type[length];

        numa_tonode_memory(strings, length * sizeof(string), numaNode);
        numa_tonode_memory(lcps, length * sizeof(lcp_t), numaNode);
        numa_tonode_memory(cachedChars, length * sizeof(char_type), numaNode);
#else
        strings = (string*)numa_alloc_onnode(length * sizeof(string), numaNode);
        lcps = (lcp_t*)numa_alloc_onnode(length * sizeof(lcp_t), numaNode);
        cachedChars = (char_type*)numa_alloc_onnode(length * sizeof(char_type), numaNode);
#endif
        size = length;
    }

    void freeNumaMemory()
    {
#if 0
        delete[] strings;
        delete[] lcps;
        delete[] cachedChars;
#else
        numa_free(strings, size * sizeof(string));
        numa_free(lcps, size * sizeof(string));
        numa_free(cachedChars, size * sizeof(char));
#endif
    }

    // preincrement
    LcpCacheStringPtr& operator ++ ()
    {
        ++strings;
        ++lcps;
        ++cachedChars;
        --size;
        return *this;
    }

    //! return sub-array of (string,lcp) with offset and size
    LcpCacheStringPtr sub(size_t offset, size_t n) const
    {
        assert(offset + n <= size);
        return LcpCacheStringPtr(strings + offset, lcps + offset, cachedChars + offset, n);
    }

    //! return sub-array of (string,lcp) with offset and size, but remove
    //! cachedChars
    LcpStringPtr subNoCache(size_t offset, size_t n) const
    {
        assert(offset + n <= size);
        return LcpStringPtr(strings + offset, lcps + offset, n);
    }

    //! return empty end array.
    LcpCacheStringPtr end() const
    {
        return sub(size, 0);
    }

    size_t operator - (const LcpCacheStringPtr& rhs) const
    {
        return strings - rhs.strings;
    }

    bool operator < (const LcpCacheStringPtr& rhs) const
    {
        return strings < rhs.strings;
    }

    size_t binarySearch(string searched) const
    {
        size_t idx;
        size_t l = 0;
        size_t r = size - 1;
        size_t hl = 0, hr = 0;

        if (scmp(searched, strings[0], hl) <= 0)
        {
            return 0;
        }
        else if (scmp(searched, strings[r], hr) > 0)
        {
            return size;
        }
        else
        {
            while ((r - l) > 1)
            {
                size_t m = (l + r) / 2;

                size_t h = std::min(hl, hr);

                if (scmp(searched, strings[m], h) <= 0)
                {
                    r = m;
                    hr = h;
                }
                else
                {
                    l = m;
                    hl = h;
                }
            }
            idx = r;
        }

        return idx;
    }
};

////////////////////////////////////////////////////////////////////////////////

//! Objectified string array pointer and shadow pointer array for out-of-place
//! swapping of pointers.
template <bool WithLcp>
class StringShadowPtrBase
{
protected:
    //! strings (front) and temporary shadow (back) array
    string* active_, * shadow_;

    //! length of subarray
    size_t size_;

    //! false if active_ is original, true if shadow_ is original
    bool flipped_;

public:
    //! constructor specifying all attributes
    StringShadowPtrBase(string* original, string* shadow = NULL,
                        size_t size = 0, bool flipped = false)
        : active_(original), shadow_(shadow), size_(size), flipped_(flipped)
    { }

    //! true if flipped to back array
    bool flipped() const { return flipped_; }

    //! return currently active array
    string * active() const { return active_; }

    //! return current shadow array
    string * shadow() const { return shadow_; }

    //! return valid length
    size_t size() const { return size_; }

    //! ostream-able
    friend std::ostream& operator << (std::ostream& os, const StringShadowPtrBase& sp)
    {
        return os << '(' << sp.active() << '/' << sp.shadow()
                  << '|' << sp.flipped() << ':' << sp.size() << ')';
    }

    //! Advance (both) pointers by given offset, return sub-array
    StringShadowPtrBase sub(size_t offset, size_t size) const
    {
        assert(offset + size <= size_);
        return StringShadowPtrBase(active_ + offset, shadow_ + offset, size, flipped_);
    }

    //! construct a StringShadowPtrBase object specifying a sub-array with flipping to
    //! other array.
    StringShadowPtrBase flip(size_t offset, size_t size) const
    {
        assert(offset + size <= size_);
        return StringShadowPtrBase(shadow_ + offset, active_ + offset, size, !flipped_);
    }

    //! Return the original for this StringShadowPtr for LCP calculation
    StringShadowPtrBase original() const
    {
        return flipped_ ? flip(0, size_) : *this;
    }

    //! return subarray pointer to n strings in original array, might copy from
    //! shadow before returning.
    StringShadowPtrBase copy_back() const
    {
        if (!flipped_) {
            return *this;
        }
        else {
            memcpy(shadow_, active_, size_ * sizeof(string));
            return flip(0, size_);
        }
    }

    //! check sorted order of strings
    bool check() const
    {
        for (size_t i = 1; i < size_; ++i)
            assert(scmp(out(i - 1), out(i)) <= 0);
        return true;
    }

    //! Return i-th string pointer from active_
    string & str(size_t i) const
    {
        assert(i < size_);
        return active_[i];
    }

    //! if we want to save the LCPs
    static bool with_lcp()
    {
        return WithLcp;
    }

    //! return reference to the i-th lcp
    uintptr_t & lcp(size_t i) const
    {
        if (!WithLcp) assert(0);

        assert(!flipped_);
        assert(i < size_);
        return ((uintptr_t*)shadow_)[i];
    }

    //! set the i-th lcp to v and check its value
    void set_lcp(size_t i, const uintptr_t& v) const
    {
        if (!WithLcp) return;

        assert(i > 0);
        assert(i < size_);
        assert(v == calc_lcp(out(i - 1), out(i)));

        lcp(i) = v;
    }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    void fill_lcp(uintptr_t v)
    {
        if (!WithLcp) return;

        for (size_t i = 1; i < size_; ++i)
        {
            set_lcp(i, v);
            set_cache(i, 0);
        }
    }

    //! set the i-th distinguishing cache charater to c
    void set_cache(size_t, const char_type&) const
    {
        // no-op
    }

    //! Return pointer to LCP array
    uintptr_t * lcparray() const
    {
        if (!WithLcp) assert(0);

        assert(!flipped_);
        return (uintptr_t*)shadow_;
    }

    //! Return the output string array
    string * output() const
    {
        assert(!flipped_); // active_ is original/output array
        return active_;
    }

    //! Return i-th output string pointer from active_ / output()
    string & out(size_t i) const
    {
        assert(!flipped_); // active_ is original/output array
        return str(i);
    }
};

typedef StringShadowPtrBase<false> StringShadowPtr;
typedef StringShadowPtrBase<true> StringShadowLcpPtr;

////////////////////////////////////////////////////////////////////////////////

//! Objectified string array pointer and shadow pointer array for out-of-place
//! swapping of pointers. With separate output array. Use class, do not derive
//! from StringShadowPtrBase!, since we must adapt all functions using out()!
template <bool WithLcp>
class StringShadowOutPtrBase
{
protected:
    //! encapsuled StringShadowPtrBase type
    typedef StringShadowPtrBase<WithLcp> base_type;

    //! encapsuled StringShadowPtrBase
    base_type sp;

    //! output string array
    string* output_;

public:
    //! constructor specifying all attributes
    StringShadowOutPtrBase(
        string* original, string* shadow = NULL, string* output = NULL,
        size_t size = 0, bool flipped = false)
        : sp(original, shadow, size, flipped),
          output_(output)
    { }

    //! true if flipped to back array
    bool flipped() const { return sp.flipped(); }

    //! return currently active array
    string * active() const { return sp.active(); }

    //! return current shadow array
    string * shadow() const { return sp.shadow(); }

    //! return valid length
    size_t size() const { return sp.size(); }

    //! ostream-able
    friend std::ostream& operator << (std::ostream& os, const StringShadowOutPtrBase& sp)
    {
        return os << '(' << sp.active() << '/' << sp.shadow() << '/' << sp.output()
                  << '|' << sp.flipped() << ':' << sp.size() << ')';
    }

    //! Advance (both) pointers by given offset, return sub-array
    StringShadowOutPtrBase sub(size_t offset, size_t size) const
    {
        assert(offset + size <= sp.size());
        return StringShadowOutPtrBase(
            active() + offset, shadow() + offset, output_ + offset,
            size, flipped());
    }

    //! construct a StringShadowOutPtrBase object specifying a sub-array with flipping to
    //! other array.
    StringShadowOutPtrBase flip(size_t offset, size_t size) const
    {
        assert(offset + size <= sp.size());
        return StringShadowOutPtrBase(
            shadow() + offset, active() + offset, output_ + offset,
            size, !flipped());
    }

    //! Return the original for this StringShadowOutPtr for LCP calculation
    StringShadowOutPtrBase original() const
    {
        return flipped() ? flip(0, size()) : *this;
    }

    //! return subarray pointer to n strings in original array, might copy from
    //! shadow before returning.
    StringShadowOutPtrBase copy_back() const
    {
        memcpy(output_, active(), size() * sizeof(string));
        return original();
    }

    //! check sorted order of strings
    bool check() const
    {
        for (size_t i = 1; i < size(); ++i)
            assert(scmp(out(i - 1), out(i)) <= 0);
        return true;
    }

    //! Return i-th string pointer from active_
    string & str(size_t i) const { return sp.str(i); }

    //! if we want to save the LCPs
    static bool with_lcp() { return base_type::with_lcp(); }

    //! return reference to the i-th lcp
    uintptr_t & lcp(size_t i) const { return sp.lcp(i); }

    //! set the i-th lcp to v and check its value
    void set_lcp(size_t i, const uintptr_t& v) const
    {
        if (!WithLcp) return;

        assert(i > 0);
        assert(i < size());
        assert(v == calc_lcp(out(i - 1), out(i)));

        lcp(i) = v;
    }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    void fill_lcp(uintptr_t v)
    {
        if (!WithLcp) return;

        for (size_t i = 1; i < size(); ++i)
        {
            set_lcp(i, v);
            set_cache(i, 0);
        }
    }

    //! set the i-th distinguishing cache charater to c
    void set_cache(size_t, const char_type&) const
    {
        // no-op
    }

    //! Return pointer to LCP array
    uintptr_t * lcparray() const { return sp.lcparray(); }

    //! Return the output string array
    string * output() const
    {
        return output_;
    }

    //! Return i-th output string pointer from active_ / output()
    string & out(size_t i) const
    {
        assert(i < size());
        return output_[i];
    }
};

typedef StringShadowOutPtrBase<false> StringShadowOutPtr;
typedef StringShadowOutPtrBase<true> StringShadowLcpOutPtr;

////////////////////////////////////////////////////////////////////////////////

//! Objectified string array pointer, shadow pointer array for out-of-place
//! swapping of pointers and lcp output, and character cache of distinguishing
//! characters.
class StringShadowLcpCacheOutPtr : public StringShadowLcpOutPtr
{
protected:
    //! our own type
    typedef StringShadowLcpCacheOutPtr self_type;

    //! encapsuled StringShadowLcpOutPtr type
    typedef StringShadowLcpOutPtr super_type;

    //! character cache array
    char_type* cache_;

public:
    //! constructor specifying all attributes
    StringShadowLcpCacheOutPtr(
        string* original = NULL, string* shadow = NULL, string* output = NULL,
        char_type* cache = NULL,
        size_t size = 0, bool flipped = false)
        : super_type(original, shadow, output, size, flipped),
          cache_(cache)
    { }

    //! return character cache of distinguishing chars
    char_type * cache() const { return cache_; }

    //! Advance (all) pointers by given offset, return sub-array
    self_type sub(size_t offset, size_t size) const
    {
        assert(offset + size <= this->size());
        return self_type(active() + offset, shadow() + offset, output() + offset,
                         cache_ + offset, size, flipped());
    }

    //! construct a StringShadowOutPtrBase object specifying a sub-array with
    //! flipping to other array.
    self_type flip(size_t offset, size_t size) const
    {
        assert(offset + size <= this->size());
        return self_type(shadow() + offset, active() + offset, output_ + offset,
                         cache_ + offset, size, !flipped());
    }

    //! Return the original of this StringPtr for LCP calculation
    self_type original() const
    {
        return flipped() ? flip(0, size()) : *this;
    }

    //! return subarray pointer to n strings in original array, might copy from
    //! shadow before returning.
    self_type copy_back() const
    {
        memcpy(super_type::output_, active(), size() * sizeof(string));
        return original();
    }

    //! set the i-th distinguishing cache charater to c
    void set_cache(size_t i, const char_type& c) const
    {
        assert(i < size());
        cache_[i] = c;
    }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    void fill_lcp(uintptr_t v)
    {
        if (!with_lcp()) return;

        for (size_t i = 1; i < size(); ++i)
        {
            set_lcp(i, v);
            set_cache(i, 0);
        }
    }
};

////////////////////////////////////////////////////////////////////////////////

//! verify LCP array against sorted string array by scanning LCPs
template <bool checkCache, typename StringSet>
static inline bool
verify_lcp_cache(const StringSet& ss,
                 lcp_t* lcps, char_type* cache, lcp_t expectedFirstLcp)
{
    typedef typename StringSet::String String;
    typedef typename StringSet::Iterator Iterator;

    bool allValid = true;
    Iterator begin = ss.begin();

    if (expectedFirstLcp != (lcp_t)-1)
    {
        if (lcps[0] != expectedFirstLcp)
        {
            std::cout << "lcp[0] = " << lcps[0]
                      << " excepted " << expectedFirstLcp << std::endl;
            allValid = false;
        }
        if (checkCache && *cache != ss.get_char(ss[begin], lcps[0]))
        {
            std::cout << "cache[0] = " << cache[0]
                      << " excepted " << ss.get_char(ss[begin], lcps[0])
                      << std::endl;
            allValid = false;
        }
    }

    for (size_t i = 1; i < ss.size(); ++i)
    {
        const String& s1 = ss[begin + i - 1];
        const String& s2 = ss[begin + i];

        size_t h = calc_lcp(ss, s1, s2);

        if (h != lcps[i])
        {
            std::cout << "lcp[" << i << "] = " << lcps[i]
                      << " excepted " << h << std::endl;
            allValid = false;
        }
        if (checkCache && cache[i] != ss.get_char(s2, lcps[i]))
        {
            std::cout << "cache[" << i << "] = " << cache[i]
                      << " excepted " << ss.get_char(s2, lcps[i]) << std::endl;
            allValid = false;
        }
    }

    if (allValid)
        std::cout << "All LCPs and cache values valid!" << std::endl;
    else
        std::cout << "Found invalid LCPS and/or cache values!" << std::endl;

    return allValid;
}

template <typename StringSet>
static inline bool
verify_lcp(const StringSet& ss, lcp_t* lcps, lcp_t expectedFirstLcp)
{
    return verify_lcp_cache<false>(ss, lcps, NULL, expectedFirstLcp);
}

static inline bool
verify_lcp(string* strings, lcp_t* lcps, size_t n, lcp_t expectedFirstLcp)
{
    return verify_lcp_cache<false>(
        parallel_string_sorting::UCharStringSet(strings, strings + n),
        lcps, NULL, expectedFirstLcp);
}

static inline bool
verify_lcp_cache(string* strings, lcp_t* lcps,
                 char_type* cache, size_t n, lcp_t expectedFirstLcp)
{
    return verify_lcp_cache<true>(
        parallel_string_sorting::UCharStringSet(strings, strings + n),
        lcps, cache, expectedFirstLcp);
}

} // namespace stringtools

#endif // !PSS_SRC_TOOLS_STRINGPTR_HEADER

/******************************************************************************/
