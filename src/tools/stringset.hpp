/*******************************************************************************
 * src/tools/stringset.hpp
 *
 * Implementations of StringSet concept: UCharStringSet, VectorStringSet,
 * StringSuffixSet.
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
 * Copyright (C) 2015 Timo Bingmann <tb@panthema.net>
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

#ifndef PSS_SRC_TOOLS_STRINGSET_HEADER
#define PSS_SRC_TOOLS_STRINGSET_HEADER

#include <cassert>
#include <stdint.h>
#include <vector>
#include <memory>

#include "debug.hpp"

namespace parallel_string_sorting {

typedef uintptr_t lcp_t;

/******************************************************************************/
// CharIterator -> character group functions

template <typename CharIterator>
inline uint8_t get_char_uint8(CharIterator str, size_t depth)
{
    return (uint8_t)(str[depth]);
}

template <typename CharIterator>
inline uint16_t get_char_uint16(CharIterator str, size_t depth)
{
    uint16_t v = 0;
    if (str[depth] == 0) return v;
    v |= (uint16_t(str[depth]) << 8);
    ++str;
    v |= (uint16_t(str[depth]) << 0);
    return v;
}

template <typename CharIterator>
inline uint32_t get_char_uint32(CharIterator str, size_t depth)
{
#if 0
    uint32_t v = 0;
    if (str[depth] == 0) return v;
    v |= (uint32_t(str[depth]) << 24);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint32_t(str[depth]) << 16);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint32_t(str[depth]) << 8);
    ++str;
    v |= (uint32_t(str[depth]) << 0);
    return v;
#else
    uint32_t v = __builtin_bswap32(*(uint32_t*)(str + depth));
    if ((v & 0xFF000000LU) == 0)
        return 0;
    else if ((v & 0x00FF0000LU) == 0)
        return (v & 0xFFFF0000LU);
    else if ((v & 0x0000FF00LU) == 0)
        return (v & 0xFFFFFF00LU);
    return v;
#endif
}

template <typename CharIterator>
inline uint64_t get_char_uint64_simple(CharIterator str, size_t depth)
{
    uint64_t v = 0;
    if (str[depth] == 0) return v;
    v = (uint64_t(str[depth]) << 56);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 48);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 40);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 32);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 24);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 16);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 8);
    ++str;
    v |= (uint64_t(str[depth]) << 0);
    return v;
}

template <typename CharIterator>
inline uint64_t get_char_uint64_bswap64(CharIterator str, size_t depth)
{
    uint64_t v = __builtin_bswap64(*(uint64_t*)(str + depth));
    if ((v & 0xFF00000000000000LLU) == 0)
        return 0;
    else if ((v & 0x00FF000000000000LLU) == 0)
        return (v & 0xFFFF000000000000LLU);
    else if ((v & 0x0000FF0000000000LLU) == 0)
        return (v & 0xFFFFFF0000000000LLU);
    else if ((v & 0x000000FF00000000LLU) == 0)
        return (v & 0xFFFFFFFF00000000LLU);
    else if ((v & 0x00000000FF000000LLU) == 0)
        return (v & 0xFFFFFFFFFF000000LLU);
    else if ((v & 0x0000000000FF0000LLU) == 0)
        return (v & 0xFFFFFFFFFFFF0000LLU);
    else if ((v & 0x000000000000FF00LLU) == 0)
        return (v & 0xFFFFFFFFFFFFFF00LLU);
    return v;
}

template <typename CharIterator>
inline uint64_t get_char_uint64(CharIterator str, size_t depth);

template <>
inline uint64_t get_char_uint64(unsigned char* str, size_t depth)
{
    return get_char_uint64_bswap64(str, depth);
}

template <>
inline uint64_t get_char_uint64(char* str, size_t depth)
{
    return get_char_uint64_bswap64(str, depth);
}

template <>
inline uint64_t get_char_uint64(std::string::const_iterator str, size_t depth)
{
    return get_char_uint64_bswap64(&(*str), depth);
}

/******************************************************************************/

/*!
 * Base class for common string set functions, included via CRTP.
 */
template <typename StringSet, typename Traits>
class StringSetBase
{
public:
    typename Traits::Char
    get_char(const typename Traits::String& i, size_t depth) const
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);
        return *ss.get_chars(i, depth);
    }

    uint8_t get_uint8(const typename Traits::String& i, size_t depth) const
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);
        return get_char_uint8(ss.get_chars(i, 0), depth);
    }

    uint16_t get_uint16(const typename Traits::String& i, size_t depth) const
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);
        return get_char_uint16(ss.get_chars(i, 0), depth);
    }

    uint32_t get_uint32(const typename Traits::String& i, size_t depth) const
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);
        return get_char_uint32(ss.get_chars(i, 0), depth);
    }

    uint64_t get_uint64(const typename Traits::String& i, size_t depth) const
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);
        return get_char_uint64(ss.get_chars(i, 0), depth);
    }

    //! Subset this string set using begin and size range.
    StringSet subr(const typename Traits::Iterator& begin, size_t size) const
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);
        return ss.sub(begin, begin + size);
    }

    //! Subset this string set using index range.
    StringSet subi(size_t begin, size_t end) const
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);
        return ss.sub(ss.begin() + begin, ss.begin() + end);
    }

    void print() const
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);

        size_t i = 0;
        for (typename Traits::Iterator pi = ss.begin(); pi != ss.end(); ++pi)
        {
            std::cout << "[" << i++ << "] = " << ss[pi]
                      << " = " << ss.get_string(ss[pi], 0) << std::endl;

            // for printing std::unique_ptr<std::string>
            // std::cout << "[" << i++ << "] = " << ss[pi].get();
            // if (ss[pi])
            //     std::cout << " = " << ss.get_string(ss[pi], 0);
            // std::cout << std::endl;
        }
    }

    bool check_order()
    {
        const StringSet& ss = *static_cast<const StringSet*>(this);

        for (typename Traits::Iterator pi = ss.begin() + 1;
             pi != ss.end(); ++pi)
        {
            typename StringSet::CharIterator
            s = ss.get_chars(ss[pi - 1], 0),
            t = ss.get_chars(ss[pi], 0);

            while (*s == *t && *s != 0)
                ++s, ++t;

            if (*s > *t) return false;
        }

        return true;
    }
};

/******************************************************************************/

/*!
 * Traits class implementing StringSet concept for char* and unsigned char* strings.
 */
template <typename CharType>
class GenericCharStringSetTraits
{
public:
    //! exported alias for character type
    typedef CharType Char;

    //! String reference: pointer to first character
    typedef Char* String;

    //! Iterator over string references: pointer over pointers
    typedef String* Iterator;

    //! iterator of characters in a string
    typedef Char* CharIterator;
};

/*!
 * Class implementing StringSet concept for char* and unsigned char* strings.
 */
template <typename CharType>
class GenericCharStringSet
    : public GenericCharStringSetTraits<CharType>,
      public StringSetBase<GenericCharStringSet<CharType>,
                           GenericCharStringSetTraits<CharType> >
{
public:
    typedef GenericCharStringSetTraits<CharType> Traits;

    using typename Traits::Char;
    using typename Traits::String;
    using typename Traits::Iterator;
    using typename Traits::CharIterator;

    //! Construct from begin and end string pointers
    GenericCharStringSet(Iterator begin, Iterator end)
        : begin_(begin), end_(end)
    { }

    //! Return size of string array
    size_t size() const { return end_ - begin_; }
    //! Iterator representing first String position
    Iterator begin() const { return begin_; }
    //! Iterator representing beyond last String position
    Iterator end() const { return end_; }

    //! Iterator-based array access (readable and writable) to String objects.
    String& operator [] (Iterator i) const
    { return *i; }

    //! Return CharIterator for referenced string, which belong to this set.
    CharIterator get_chars(const String& s, size_t depth) const
    { return s + depth; }

    //! Return complete string (for debugging purposes)
    std::string get_string(const String& s, size_t depth = 0) const
    { return std::string(reinterpret_cast<const char*>(s) + depth); }

    //! Subset this string set using iterator range.
    GenericCharStringSet sub(Iterator begin, Iterator end) const
    { return GenericCharStringSet(begin, end); }

protected:
    //! array of string pointers
    Iterator begin_, end_;
};

typedef GenericCharStringSet<char> CharStringSet;
typedef GenericCharStringSet<unsigned char> UCharStringSet;

/******************************************************************************/

/*!
 * Class implementing StringSet concept for a std::vector containing std::string
 * objects.
 */
class VectorStringSetTraits
{
public:
    //! exported alias for assumed string container
    typedef std::vector<std::string> Container;

    //! exported alias for character type
    typedef std::string::value_type Char;

    //! String reference: std::string, which should be reference counted.
    typedef typename Container::value_type String;

    //! Iterator over string references: using std::vector's iterator
    typedef typename Container::iterator Iterator;

    //! iterator of characters in a string
    typedef std::string::const_iterator CharIterator;
};

/*!
 * Class implementing StringSet concept for a std::vector containing std::string
 * objects.
 */
class VectorStringSet
    : public VectorStringSetTraits,
      public StringSetBase<VectorStringSet, VectorStringSetTraits>
{
public:
    //! Construct from begin and end string pointers
    VectorStringSet(const Iterator& begin, const Iterator& end)
        : begin_(begin), end_(end)
    { }

    //! Return size of string array
    size_t size() const { return end_ - begin_; }
    //! Iterator representing first String position
    Iterator begin() const { return begin_; }
    //! Iterator representing beyond last String position
    Iterator end() const { return end_; }

    //! Array access (readable and writable) to String objects.
    String& operator [] (const Iterator& i) const
    { return *i; }

    //! Return CharIterator for referenced string, which belongs to this set.
    CharIterator get_chars(const String& s, size_t depth) const
    { return s.begin() + depth; }

    //! Return complete string (for debugging purposes)
    std::string get_string(const String&& s, size_t depth = 0) const
    { return s.substr(depth); }

    //! Subset this string set using iterator range.
    VectorStringSet sub(Iterator begin, Iterator end) const
    { return VectorStringSet(begin, end); }

protected:
    //! vector of std::string objects
    Iterator begin_, end_;
};

/******************************************************************************/

/*!
 * Class implementing StringSet concept for a std::vector containing
 * std::unique_ptr<std::string> objects, which are non-copyable.
 */
class VectorPtrStringSetTraits
{
public:
    //! exported alias for assumed string container
    typedef std::vector<std::unique_ptr<std::string> > Container;

    //! exported alias for character type
    typedef std::string::value_type Char;

    //! String reference: std::string, which should be reference counted.
    typedef typename Container::value_type String;

    //! Iterator over string references: using std::vector's iterator
    typedef typename Container::iterator Iterator;

    //! iterator of characters in a string
    typedef std::string::const_iterator CharIterator;
};

/*!
 * Class implementing StringSet concept for a std::vector containing std::string
 * objects.
 */
class VectorPtrStringSet
    : public VectorPtrStringSetTraits,
      public StringSetBase<VectorPtrStringSet, VectorPtrStringSetTraits>
{
public:
    //! Construct from begin and end string pointers
    VectorPtrStringSet(const Iterator& begin, const Iterator& end)
        : begin_(begin), end_(end)
    { }

    //! Return size of string array
    size_t size() const { return end_ - begin_; }
    //! Iterator representing first String position
    Iterator begin() const { return begin_; }
    //! Iterator representing beyond last String position
    Iterator end() const { return end_; }

    //! Array access (readable and writable) to String objects.
    String& operator [] (const Iterator& i) const
    { return *i; }

    //! Return CharIterator for referenced string, which belongs to this set.
    CharIterator get_chars(const String& s, size_t depth) const
    { return s->begin() + depth; }

    //! Return complete string (for debugging purposes)
    std::string get_string(const String& s, size_t depth = 0) const
    { return s->substr(depth); }

    //! Subset this string set using iterator range.
    VectorPtrStringSet sub(Iterator begin, Iterator end) const
    { return VectorPtrStringSet(begin, end); }

protected:
    //! vector of std::string objects
    Iterator begin_, end_;
};

/******************************************************************************/

/*!
 * Class implementing StringSet concept for suffix sorting indexes of an
 * unsigned char* text object.
 */
class UCharSuffixSetTraits
{
public:
    //! exported alias for assumed string container
    typedef unsigned char* Container;

    //! exported alias for character type
    typedef unsigned char Char;

    //! String reference: suffix index of the text.
    typedef int String;

    //! Iterator over string references: using std::vector's iterator over
    //! suffix array vector
    typedef String* Iterator;

    //! iterator of characters in a string
    typedef Char* CharIterator;
};

/*!
 * Class implementing StringSet concept for suffix sorting indexes of a
 * std::string text object.
 */
class UCharSuffixSet
    : public UCharSuffixSetTraits,
      public StringSetBase<UCharSuffixSet, UCharSuffixSetTraits>
{
public:
    //! Construct from begin and end string pointers
    UCharSuffixSet(const Container& text,
                   const Iterator& begin, const Iterator& end)
        : text_(text),
          begin_(begin), end_(end)
    { }

    //! Return size of string array
    size_t size() const { return end_ - begin_; }
    //! Iterator representing first String position
    Iterator begin() const { return begin_; }
    //! Iterator representing beyond last String position
    Iterator end() const { return end_; }

    //! Array access (readable and writable) to String objects.
    String& operator [] (const Iterator& i) const
    { return *i; }

    //! Return CharIterator for referenced string, which belongs to this set.
    CharIterator get_chars(const String& s, size_t depth) const
    { return text_ + s + depth; }

    //! Return complete string (for debugging purposes)
    std::string get_string(const String& s, size_t depth = 0) const
    { return std::string(reinterpret_cast<const char*>(text_ + s + depth)); }

    //! Subset this string set using iterator range.
    UCharSuffixSet sub(Iterator begin, Iterator end) const
    { return UCharSuffixSet(text_, begin, end); }

protected:
    //! reference to base text
    const Container& text_;

    //! iterators inside the output suffix array.
    Iterator begin_, end_;
};

/******************************************************************************/

/*!
 * Class implementing StringSet concept for suffix sorting indexes of a
 * std::string text object.
 */
class StringSuffixSetTraits
{
public:
    //! exported alias for assumed string container
    typedef std::string Container;

    //! exported alias for character type
    typedef std::string::value_type Char;

    //! String reference: suffix index of the text.
    typedef typename Container::size_type String;

    //! Iterator over string references: using std::vector's iterator over
    //! suffix array vector
    typedef typename std::vector<String>::iterator Iterator;

    //! iterator of characters in a string
    typedef std::string::const_iterator CharIterator;
};

/*!
 * Class implementing StringSet concept for suffix sorting indexes of a
 * std::string text object.
 */
class StringSuffixSet
    : public StringSuffixSetTraits,
      public StringSetBase<StringSuffixSet, StringSuffixSetTraits>
{
public:
    //! Construct from begin and end string pointers
    StringSuffixSet(const Container& text,
                    const Iterator& begin, const Iterator& end)
        : text_(text),
          begin_(begin), end_(end)
    { }

    //! Initializing constructor which fills output vector sa with indices.
    static StringSuffixSet
    Initialize(const Container& text, std::vector<String>& sa)
    {
        sa.resize(text.size());
        for (size_t i = 0; i < text.size(); ++i)
            sa[i] = i;
        return StringSuffixSet(text, sa.begin(), sa.end());
    }

    //! Return size of string array
    size_t size() const { return end_ - begin_; }
    //! Iterator representing first String position
    Iterator begin() const { return begin_; }
    //! Iterator representing beyond last String position
    Iterator end() const { return end_; }

    //! Array access (readable and writable) to String objects.
    String& operator [] (const Iterator& i) const
    { return *i; }

    //! Return CharIterator for referenced string, which belongs to this set.
    CharIterator get_chars(const String& s, size_t depth) const
    { return text_.begin() + s + depth; }

    //! Return complete string (for debugging purposes)
    std::string get_string(const String& s, size_t depth = 0) const
    { return text_.substr(s + depth); }

    //! Subset this string set using iterator range.
    StringSuffixSet sub(Iterator begin, Iterator end) const
    { return StringSuffixSet(text_, begin, end); }

protected:
    //! reference to base text
    const Container& text_;

    //! iterators inside the output suffix array.
    Iterator begin_, end_;
};

} // namespace parallel_string_sorting

#endif // !PSS_SRC_TOOLS_STRINGSET_HEADER

/******************************************************************************/
