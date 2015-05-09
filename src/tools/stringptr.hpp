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

/******************************************************************************/

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

/******************************************************************************/

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

/******************************************************************************/

//! Objectified string array pointer and shadow pointer array for out-of-place
//! swapping of pointers.
template <typename _StringSet>
class StringShadowPtr
{
public:
    typedef _StringSet StringSet;
    typedef typename StringSet::String String;

protected:
    //! strings (front) and temporary shadow (back) array
    StringSet active_, shadow_;

    //! false if active_ is original, true if shadow_ is original
    bool flipped_;

public:
    //! constructor specifying all attributes
    StringShadowPtr(const StringSet& original, const StringSet& shadow,
                    bool flipped = false)
        : active_(original), shadow_(shadow), flipped_(flipped)
    { }

    //! true if flipped to back array
    bool flipped() const { return flipped_; }

    //! return currently active array
    const StringSet & active() const { return active_; }

    //! return current shadow array
    const StringSet & shadow() const { return shadow_; }

    //! return valid length
    size_t size() const { return active_.size(); }

    //! ostream-able
    friend std::ostream& operator << (std::ostream& os,
                                      const StringShadowPtr& sp)
    {
        return os << '(' << sp.active() << '/' << sp.shadow()
                  << '|' << sp.flipped() << ':' << sp.size() << ')';
    }

    //! Advance (both) pointers by given offset, return sub-array
    StringShadowPtr sub(size_t offset, size_t _size) const
    {
        assert(offset + _size <= size());
        return StringShadowPtr(active_.subi(offset, offset + _size),
                               shadow_.subi(offset, offset + _size), flipped_);
    }

    //! construct a StringShadowPtr object specifying a sub-array with flipping
    //! to other array.
    StringShadowPtr flip(size_t offset, size_t _size) const
    {
        assert(offset + _size <= size());
        return StringShadowPtr(shadow_.subi(offset, offset + _size),
                               active_.subi(offset, offset + _size), !flipped_);
    }

    //! Return the original for this StringShadowPtr for LCP calculation
    StringShadowPtr original() const
    {
        return flipped_ ? flip(0, size()) : *this;
    }

    //! return subarray pointer to n strings in original array, might copy from
    //! shadow before returning.
    StringShadowPtr copy_back() const
    {
        if (!flipped_) {
            return *this;
        }
        else {
            std::move(active_.begin(), active_.end(), shadow_.begin());
            return flip(0, size());
        }
    }

    //! check sorted order of strings
    bool check() const
    {
        assert(output().check_order());
        return true;
    }

    //! Return i-th string pointer from active_
    String & str(size_t i) const
    {
        assert(i < size());
        return active_.at(i);
    }

    //! if we want to save the LCPs
    static bool with_lcp()
    {
        return false;
    }

    //! return reference to the i-th lcp
    uintptr_t & lcp(size_t /* i */) const
    {
        abort();
    }

    //! set the i-th lcp to v and check its value
    void set_lcp(size_t /* i */, const uintptr_t& /* v */) const
    { }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    void fill_lcp(uintptr_t /* v */)
    { }

    //! set the i-th distinguishing cache charater to c
    void set_cache(size_t, const char_type&) const
    {
        // no-op
    }

    //! Return pointer to LCP array
    uintptr_t * lcparray() const
    {
        abort();
    }

    //! Return the output string array
    const StringSet & output() const
    {
        assert(!flipped_); // active_ is original/output array
        return active_;
    }

    //! Return i-th output string pointer from active_ / output()
    String & out(size_t i) const
    {
        assert(!flipped_); // active_ is original/output array
        return str(i);
    }
};

/******************************************************************************/

//! Objectified string array pointer and shadow pointer array for out-of-place
//! swapping of pointers. With separate output array. Use class, do not derive
//! from StringShadowPtr!, since we must adapt all functions using out()!
template <typename _StringSet>
class StringShadowOutPtr
{
public:
    typedef _StringSet StringSet;
    typedef typename StringSet::String String;

    //! encapsuled StringShadowPtr type
    typedef StringShadowPtr<StringSet> Super;

protected:
    //! encapsuled StringShadowPtr
    Super sp;

    //! output string array
    StringSet output_;

    static const bool WithLcp = false;

public:
    //! constructor specifying all attributes
    StringShadowOutPtr(
        const StringSet& original, const StringSet& shadow, const StringSet& output,
        bool flipped = false)
        : sp(original, shadow, flipped),
          output_(output)
    { }

    //! true if flipped to back array
    bool flipped() const { return sp.flipped(); }

    //! return currently active array
    const StringSet & active() const { return sp.active(); }

    //! return current shadow array
    const StringSet & shadow() const { return sp.shadow(); }

    //! return valid length
    size_t size() const { return sp.size(); }

    //! ostream-able
    friend std::ostream& operator << (std::ostream& os, const StringShadowOutPtr& sp)
    {
        return os << '(' << sp.active() << '/' << sp.shadow() << '/' << sp.output()
                  << '|' << sp.flipped() << ':' << sp.size() << ')';
    }

    //! Advance (both) pointers by given offset, return sub-array
    StringShadowOutPtr sub(size_t offset, size_t _size) const
    {
        assert(offset + _size <= sp.size());
        return StringShadowOutPtr(
            sp.active().subi(offset, offset + _size),
            sp.shadow().subi(offset, offset + _size),
            output_.subi(offset, offset + _size),
            flipped());
    }

    //! construct a StringShadowOutPtr object specifying a sub-array with flipping to
    //! other array.
    StringShadowOutPtr flip(size_t offset, size_t _size) const
    {
        assert(offset + _size <= sp.size());
        return StringShadowOutPtr(
            sp.shadow().subi(offset, offset + _size),
            sp.active().subi(offset, offset + _size),
            output_.subi(offset, offset + _size),
            !flipped());
    }

    //! Return the original for this StringShadowOutPtr for LCP calculation
    StringShadowOutPtr original() const
    {
        return flipped() ? flip(0, size()) : *this;
    }

    //! return subarray pointer to n strings in original array, might copy from
    //! shadow before returning.
    StringShadowOutPtr copy_back() const
    {
        std::move(sp.active().begin(), sp.active().end(), output_.begin());
        return original();
    }

    //! check sorted order of strings
    bool check() const
    {
        assert(output().check_order());
        return true;
    }

    //! Return i-th string pointer from active_
    string & str(size_t i) const { return sp.str(i); }

    //! if we want to save the LCPs
    static bool with_lcp() { return Super::with_lcp(); }

    //! return reference to the i-th lcp
    uintptr_t & lcp(size_t i) const { return sp.lcp(i); }

    //! set the i-th lcp to v and check its value
    void set_lcp(size_t /* i */, const uintptr_t& /* v */) const
    { }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    void fill_lcp(uintptr_t /* v */)
    { }

    //! set the i-th distinguishing cache charater to c
    void set_cache(size_t, const char_type&) const
    {
        // no-op
    }

    //! Return pointer to LCP array
    uintptr_t * lcparray() const { return sp.lcparray(); }

    //! Return the output string array
    const StringSet & output() const
    {
        return output_;
    }

    //! Return i-th output string pointer from active_ / output()
    String & out(size_t i) const
    {
        assert(i < size());
        return *(output_.begin() + i);
    }
};

/******************************************************************************/

template <typename _StringSet>
class StringShadowLcpPtr : protected StringShadowPtr<_StringSet>
{
public:
    typedef _StringSet StringSet;
    typedef StringShadowPtr<_StringSet> Super;

protected:
    lcp_t* lcps_;

public:
    // *** import protected methods
    using Super::size;
    using Super::flipped;
    using Super::active;
    using Super::shadow;
    using Super::output;
    using Super::out;
    using Super::check;
    using Super::set_cache;

    //! constructor specifying all attributes
    StringShadowLcpPtr(const StringSet& original, const StringSet& shadow,
                       lcp_t* lcps, bool flipped = false)
        : Super(original, shadow, flipped),
          lcps_(lcps)
    { }

    //! Advance (both) pointers by given offset, return sub-array
    StringShadowLcpPtr sub(size_t offset, size_t _size) const
    {
        assert(offset + _size <= size());
        return StringShadowLcpPtr(active().subi(offset, offset + _size),
                                  shadow().subi(offset, offset + _size),
                                  lcps_ + offset, flipped());
    }

    //! construct a StringShadowLcpPtr object specifying a sub-array with
    //! flipping to other array.
    StringShadowLcpPtr flip(size_t offset, size_t _size) const
    {
        assert(offset + _size <= size());
        return StringShadowLcpPtr(shadow().subi(offset, offset + _size),
                                  active().subi(offset, offset + _size),
                                  lcps_ + offset, !flipped());
    }

    //! Return the original for this StringShadowPtr for LCP calculation
    StringShadowLcpPtr original() const
    {
        return flipped() ? flip(0, size()) : *this;
    }

    //! return subarray pointer to n strings in original array, might copy from
    //! shadow before returning.
    StringShadowLcpPtr copy_back() const
    {
        if (!flipped()) {
            return *this;
        }
        else {
            Super::copy_back();
            return flip(0, size());
        }
    }

    //! if we want to save the LCPs
    static bool with_lcp()
    {
        return true;
    }

    //! return reference to the i-th lcp
    uintptr_t & lcp(size_t i) const
    {
        assert(!flipped());
        assert(i < size());
        return lcps_[i];
    }

    //! set the i-th lcp to v and check its value
    void set_lcp(size_t i, const uintptr_t& v) const
    {
        assert(i > 0);
        assert(i < size());
        assert(v == calc_lcp(active(), this->out(i - 1), this->out(i)));

        lcp(i) = v;
    }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    void fill_lcp(uintptr_t v)
    {
        for (size_t i = 1; i < size(); ++i)
        {
            set_lcp(i, v);
            this->set_cache(i, 0);
        }
    }

    //! Return pointer to LCP array
    uintptr_t * lcparray() const
    {
        assert(!flipped());
        return lcps_;
    }
};

/******************************************************************************/

template <typename _StringSet>
class StringShadowLcpOutPtr : protected StringShadowOutPtr<_StringSet>
{
public:
    typedef _StringSet StringSet;
    typedef StringShadowOutPtr<_StringSet> Super;

protected:
    lcp_t* lcps_;

public:
    // *** import protected methods
    using Super::size;
    using Super::flipped;
    using Super::active;
    using Super::shadow;
    using Super::output;
    using Super::out;
    using Super::check;
    using Super::set_cache;

    //! constructor specifying all attributes
    StringShadowLcpOutPtr(
        const StringSet& original, const StringSet& shadow, const StringSet& output,
        lcp_t* lcps, bool flipped = false)
        : Super(original, shadow, output, flipped),
          lcps_(lcps)
    { }

    //! Advance (both) pointers by given offset, return sub-array
    StringShadowLcpOutPtr sub(size_t offset, size_t _size) const
    {
        assert(offset + _size <= size());
        return StringShadowLcpOutPtr(active().subi(offset, offset + _size),
                                     shadow().subi(offset, offset + _size),
                                     output().subi(offset, offset + _size),
                                     lcps_ + offset, flipped());
    }

    //! construct a StringShadowLcpOutPtr object specifying a sub-array with
    //! flipping to other array.
    StringShadowLcpOutPtr flip(size_t offset, size_t _size) const
    {
        assert(offset + _size <= size());
        return StringShadowLcpOutPtr(shadow().subi(offset, offset + _size),
                                     active().subi(offset, offset + _size),
                                     output().subi(offset, offset + _size),
                                     lcps_ + offset, !flipped());
    }

    //! Return the original for this StringShadowPtr for LCP calculation
    StringShadowLcpOutPtr original() const
    {
        return flipped() ? flip(0, size()) : *this;
    }

    //! return subarray pointer to n strings in original array, might copy from
    //! shadow before returning.
    StringShadowLcpOutPtr copy_back() const
    {
        if (!flipped()) {
            return *this;
        }
        else {
            Super::copy_back();
            return flip(0, size());
        }
    }

    //! if we want to save the LCPs
    static bool with_lcp()
    {
        return true;
    }

    //! return reference to the i-th lcp
    uintptr_t & lcp(size_t i) const
    {
        assert(!flipped());
        assert(i < size());
        return lcps_[i];
    }

    //! set the i-th lcp to v and check its value
    void set_lcp(size_t i, const uintptr_t& v) const
    {
        assert(i > 0);
        assert(i < size());
        //assert(v == calc_lcp(active_, out(i - 1), out(i)));

        lcp(i) = v;
    }

    //! Fill whole LCP array with n times the value v, ! excluding the first
    //! LCP[0] position
    void fill_lcp(uintptr_t v)
    {
        for (size_t i = 1; i < size(); ++i)
        {
            set_lcp(i, v);
            this->set_cache(i, 0);
        }
    }

    //! Return pointer to LCP array
    uintptr_t * lcparray() const
    {
        assert(!flipped());
        return lcps_;
    }
};

/******************************************************************************/

//! Objectified string array pointer, shadow pointer array for out-of-place
//! swapping of pointers and lcp output, and character cache of distinguishing
//! characters.
template <typename _StringSet>
class StringShadowLcpCacheOutPtr : public StringShadowLcpOutPtr<_StringSet>
{
public:
    //! our own type
    typedef StringShadowLcpCacheOutPtr Self;

    //! string set type
    typedef _StringSet StringSet;

    //! encapsuled StringShadowLcpOutPtr type
    typedef StringShadowLcpOutPtr<StringSet> Super;

    //! character type
    typedef typename StringSet::Char Char;

protected:
    //! character cache array
    Char* cache_;

    using Super::output_;
    using Super::lcps_;

public:
    // *** import public methods
    using Super::size;
    using Super::flipped;
    using Super::active;
    using Super::shadow;
    using Super::output;
    using Super::with_lcp;
    using Super::set_lcp;

    //! constructor specifying all attributes
    StringShadowLcpCacheOutPtr(
        const StringSet& original, const StringSet& shadow, const StringSet& output,
        lcp_t* lcps, Char* cache = NULL, bool flipped = false)
        : Super(original, shadow, output, lcps, flipped),
          cache_(cache)
    { }

    //! return character cache of distinguishing chars
    Char * cache() const { return cache_; }

    //! Advance (all) pointers by given offset, return sub-array
    Self sub(size_t offset, size_t _size) const
    {
        assert(offset + _size <= size());
        return Self(active().subi(offset, offset + _size),
                    shadow().subi(offset, offset + _size),
                    output().subi(offset, offset + _size),
                    lcps_, cache_ + offset, flipped());
    }

    //! construct a StringShadowOutPtr object specifying a sub-array with
    //! flipping to other array.
    Self flip(size_t offset, size_t _size) const
    {
        assert(offset + _size <= size());
        return Self(shadow().subi(offset, offset + _size),
                    active().subi(offset, offset + _size),
                    output().subi(offset, offset + _size),
                    lcps_, cache_ + offset, !flipped());
    }

    //! Return the original of this StringPtr for LCP calculation
    Self original() const
    {
        return flipped() ? flip(0, size()) : *this;
    }

    //! return subarray pointer to n strings in original array, might copy from
    //! shadow before returning.
    Self copy_back() const
    {
        Super::copy_back();
        return original();
    }

    //! set the i-th distinguishing cache charater to c
    void set_cache(size_t i, const Char& c) const
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

typedef StringShadowLcpCacheOutPtr<parallel_string_sorting::UCharStringSet>
    UCharStringShadowLcpCacheOutPtr;

////////////////////////////////////////////////////////////////////////////////

//! verify LCP array against sorted string array by scanning LCPs
template <bool checkCache = true, typename StringSet>
static inline bool
verify_lcp_cache(
    const StringSet& ss, const lcp_t* lcps, const typename StringSet::Char* cache,
    lcp_t expectedFirstLcp)
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
verify_lcp(const StringSet& ss, const lcp_t* lcps, lcp_t expectedFirstLcp)
{
    return verify_lcp_cache<false, StringSet>(ss, lcps, NULL, expectedFirstLcp);
}

static inline bool
verify_lcp(string* strings, const lcp_t* lcps, size_t n, lcp_t expectedFirstLcp)
{
    typedef parallel_string_sorting::UCharStringSet StringSet;
    return verify_lcp_cache<false, StringSet>(
        StringSet(strings, strings + n), lcps, nullptr, expectedFirstLcp);
}

static inline bool
verify_lcp_cache(string* strings, const lcp_t* lcps,
                 const char_type* cache, size_t n, lcp_t expectedFirstLcp)
{
    typedef parallel_string_sorting::UCharStringSet StringSet;
    return verify_lcp_cache<true, StringSet>(
        StringSet(strings, strings + n), lcps, cache, expectedFirstLcp);
}

} // namespace stringtools

#endif // !PSS_SRC_TOOLS_STRINGPTR_HEADER

/******************************************************************************/
