/******************************************************************************
 * src/tools/stringtools.h
 *
 * Some template tools to access strings. Used in sample and radix sorts.
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

namespace stringtools {

/// zero-terminated character strings
typedef unsigned char* string;

/// hacky gcc synthesised 128-bit datatype
typedef unsigned int uint128_t __attribute__((mode(TI)));

/// represent hex octets of large integer datatypes
template <typename Type>
static inline std::string toHex(Type v)
{
    char out[2 * sizeof(v)+1];
    static const char hex[17] = "0123456789ABCDEF";
    for (unsigned int i = 1; i <= sizeof(v); ++i)
    {
        out[2*(sizeof(v)-i)+0] = hex[(v & 0xF0) >> 4];
        out[2*(sizeof(v)-i)+1] = hex[(v & 0x0F) >> 0];
        v >>= 8;
    }
    out[2*sizeof(v)] = 0;
    return out;
}

/// represent binary digits of large integer datatypes
template <typename Type>
static inline std::string toBinary(Type v)
{
    static const int w = (1 << sizeof(v));
    char binstr[w+1];
    binstr[w] = 0;
    for (int i = 0; i < w; i++) {
        binstr[w-i-1] = (v & 1) ? '1' : '0';
        v /= 2;
    }
    return binstr;
}

/// Return traits of key_type
template <typename CharT>
class key_traits {
};

template <>
class key_traits<uint8_t> {
public:
    static const size_t radix = 256;
    static const size_t add_depth = 1;
};


template <>
class key_traits<uint16_t> {
public:
    static const size_t radix = 65536;
    static const size_t add_depth = 2;
};

/// get packed characters from string at certain depth, needed due to
/// zero-termianted strings.
template <typename CharT>
inline CharT get_char(string ptr, size_t depth);

template <>
inline uint8_t get_char<uint8_t>(string str, size_t depth)
{
    return (uint8_t)(str[depth]);
}

template <>
inline uint16_t get_char<uint16_t>(string str, size_t depth)
{
    uint16_t v = 0;
    if (str[depth] == 0) return v;
    v |= (uint16_t(str[depth]) << 8); ++str;
    v |= (uint16_t(str[depth]) << 0);
    return v;
}

template <>
inline uint32_t get_char<uint32_t>(string str, size_t depth)
{
    uint32_t v = 0;
    if (str[depth] == 0) return v;
    v |= (uint32_t(str[depth]) << 24); ++str;
    if (str[depth] == 0) return v;
    v |= (uint32_t(str[depth]) << 16); ++str;
    if (str[depth] == 0) return v;
    v |= (uint32_t(str[depth]) << 8); ++str;
    v |= (uint32_t(str[depth]) << 0);
    return v;
}

template <>
inline uint64_t get_char<uint64_t>(string str, size_t depth)
{
#if 0
    uint64_t v = 0;
    if (str[depth] == 0) return v;
    v = (uint64_t(str[depth]) << 56); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 48); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 40); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 32); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 24); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 16); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 8); ++str;
    v |= (uint64_t(str[depth]) << 0);
    return v;
#else
    uint64_t v = __builtin_bswap64( *(uint64_t*)(str+depth) );
    if ( (v & 0xFF00000000000000LLU) == 0 )
        return 0;
    else if (  (v & 0x00FF000000000000LLU) == 0 )
        return (v & 0xFFFF000000000000LLU);
    else if (  (v & 0x0000FF0000000000LLU) == 0 )
        return (v & 0xFFFFFF0000000000LLU);
    else if (  (v & 0x000000FF00000000LLU) == 0 )
        return (v & 0xFFFFFFFF00000000LLU);
    else if (  (v & 0x00000000FF000000LLU) == 0 )
        return (v & 0xFFFFFFFFFF000000LLU);
    else if (  (v & 0x0000000000FF0000LLU) == 0 )
        return (v & 0xFFFFFFFFFFFF0000LLU);
    else if (  (v & 0x000000000000FF00LLU) == 0 )
        return (v & 0xFFFFFFFFFFFFFF00LLU);
    return v;
#endif
}

template <>
inline uint128_t get_char<uint128_t>(string str, size_t depth)
{
    uint128_t v = 0;
    if (str[depth] == 0) return v;
    v = (uint128_t(str[depth]) << 120); ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 112); ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 104); ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 96); ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 88); ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 80); ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 72); ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 64); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 56); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 48); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 40); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 32); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 24); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 16); ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 8); ++str;
    v |= (uint64_t(str[depth]) << 0);
    return v;
}

/// Templated hardware operation to get the number of zero bits, starting from
/// the most significant.
template <typename T>
inline int count_high_zero_bits(const T& t);

template <>
inline int count_high_zero_bits<uint32_t>(const uint32_t& t)
{
    if (t == 0) return sizeof(t) * 8;
    return __builtin_clz(t);
}

template <>
inline int count_high_zero_bits<uint64_t>(const uint64_t& t)
{
    if (t == 0) return sizeof(t) * 8;
    return __builtin_clzll(t);
}

template <>
inline int count_high_zero_bits<uint128_t>(const uint128_t& t)
{
    if (t == 0) return sizeof(t) * 8;
    uint64_t hi = (t >> 64);
    int b = __builtin_clzll(hi);
    if (b > 0)
        return b;
    else
        return 64 + __builtin_clzll( (uint64_t)t );
}

/// Objectified string array pointer and shadow pointer array for out-of-place
/// swapping of pointers.
class StringPtr
{
public:

    /// strings (front) and temporary shadow (back) array
    string      *m_front, *m_back;

    /// true if back array is active
    bool        m_flip;

    /// constructor specifying all attributes
    inline StringPtr(string* front, string* back = NULL, bool flip = false)
        : m_front(front), m_back(back), m_flip(flip)
    {
    }

    /// true if flipped to back array
    inline bool flipped() const
    {
        return m_flip;
    }

    /// return currently active array
    inline string* active()
    {
        return (m_flip ? m_back : m_front);
    }

    /// return current shadow array
    inline string* shadow()
    {
        return (m_flip ? m_front : m_back);
    }

    /// construct a StringPtr object specifying a subarray with flipping to
    /// other array.
    inline StringPtr flip_ptr(size_t offset) const
    {
        return StringPtr(m_front + offset, m_back + offset, !m_flip);
    }

    /// return subarray pointer to n strings in original array, might copy from
    /// shadow before returning.
    inline string* to_original(size_t n)
    {
        if (m_flip) {
            assert(m_back);
            memcpy(m_front, m_back, n * sizeof(string));
            m_flip = false;
        }
        return m_front;
    }
};

} // namespace stringtools
