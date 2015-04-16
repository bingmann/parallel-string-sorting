/*******************************************************************************
 * src/tools/stringtools.h
 *
 * Some template tools to access strings. Used in sample and radix sorts.
 *
 *******************************************************************************
 * Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
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

#ifndef PSS_SRC_TOOLS_STRINGTOOLS_HEADER
#define PSS_SRC_TOOLS_STRINGTOOLS_HEADER

#include <cassert>
#include "debug.h"

namespace stringtools {

//! zero-terminated character strings
typedef unsigned char* string;

//! type of characters in strings
typedef unsigned char char_type;

/// hacky gcc synthesised 128-bit datatype
typedef unsigned int uint128_t __attribute__ ((mode(TI)));

/// represent hex octets of large integer datatypes
template <typename Type>
static inline std::string toHex(Type v)
{
    char out[2 * sizeof(v) + 1];
    static const char hex[17] = "0123456789ABCDEF";
    for (unsigned int i = 1; i <= sizeof(v); ++i)
    {
        out[2 * (sizeof(v) - i) + 0] = hex[(v & 0xF0) >> 4];
        out[2 * (sizeof(v) - i) + 1] = hex[(v & 0x0F) >> 0];
        v >>= 8;
    }
    out[2 * sizeof(v)] = 0;
    return out;
}

/// represent binary digits of large integer datatypes
template <typename Type>
static inline std::string toBinary(Type v, const int width = (1 << sizeof(Type)))
{
    char binstr[width + 1];
    binstr[width] = 0;
    for (int i = 0; i < width; i++) {
        binstr[width - i - 1] = (v & 1) ? '1' : '0';
        v /= 2;
    }
    return binstr;
}

/// compare strings by scanning
static inline int scmp(const string _s1, const string _s2)
{
    string s1 = _s1, s2 = _s2;

    while (*s1 != 0 && *s1 == *s2)
        s1++, s2++;
    return (*s1 - *s2);
}

// compare strings by scanning. Start at given lcp, which also returns the final lcp.
static inline int
scmp(const string _s1, const string _s2, size_t& lcp)
{
    string s1 = _s1 + lcp, s2 = _s2 + lcp;

    while (*s1 != 0 && *s1 == *s2)
        s1++, s2++, lcp++;
    return (*s1 - *s2);
}

/// calculate lcp by scanning
static inline unsigned int calc_lcp(const string _s1, const string _s2)
{
    string s1 = _s1, s2 = _s2;

    size_t h = 0;
    while (*s1 != 0 && *s1 == *s2)
        ++h, ++s1, ++s2;

    return h;
}

/// Return traits of key_type
template <typename CharT>
class key_traits
{ };

template <>
class key_traits<uint8_t>
{
public:
    static const size_t radix = 256;
    static const size_t add_depth = 1;
};

template <>
class key_traits<uint16_t>
{
public:
    static const size_t radix = 65536;
    static const size_t add_depth = 2;
};

template <>
class key_traits<uint32_t>
{
public:
    static const size_t radix = 4294967296;
    static const size_t add_depth = 4;
};

template <>
class key_traits<uint64_t>
{
public:
    static const size_t add_depth = 8;
};

/// get packed characters from string at certain depth, needed due to
/// zero-termianted strings.
template <typename CharT>
static inline CharT get_char(string ptr, size_t depth)
__attribute__ ((always_inline));

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
    v |= (uint16_t(str[depth]) << 8);
    ++str;
    v |= (uint16_t(str[depth]) << 0);
    return v;
}

template <>
inline uint32_t get_char<uint32_t>(string str, size_t depth)
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

template <>
inline uint64_t get_char<uint64_t>(string str, size_t depth)
{
#if 0
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
#else
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
#endif
}

template <>
inline uint128_t get_char<uint128_t>(string str, size_t depth)
{
    uint128_t v = 0;
    if (str[depth] == 0) return v;
    v = (uint128_t(str[depth]) << 120);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 112);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 104);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 96);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 88);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 80);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 72);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint128_t(str[depth]) << 64);
    ++str;
    if (str[depth] == 0) return v;
    v |= (uint64_t(str[depth]) << 56);
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

/// Templated hardware operation to get the number of zero bits, starting from
/// the most significant.
template <typename T>
static inline int count_high_zero_bits(const T& t)
__attribute__ ((always_inline));

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
        return 64 + __builtin_clzll((uint64_t)t);
}

/// Templated hardware operation to get the number of zero bits, starting from
/// the least significant.
template <typename T>
static inline int count_low_zero_bits(const T& t)
__attribute__ ((always_inline));

template <>
inline int count_low_zero_bits<uint32_t>(const uint32_t& t)
{
    if (t == 0) return sizeof(t) * 8;
    return __builtin_ctz(t);
}

template <>
inline int count_low_zero_bits<uint64_t>(const uint64_t& t)
{
    if (t == 0) return sizeof(t) * 8;
    return __builtin_ctzll(t);
}

//! Class to transform in-order to level-order indexes in a perfect binary tree
template <size_t treebits>
struct TreeCalculations
{
    static const bool   debug = false;

    static const size_t numnodes = (1 << treebits) - 1;

    static inline unsigned int
                        level_to_inorder(unsigned int id)
    {
        assert(id > 0);
        DBG(debug, "index: " << id << " = " << toBinary(id));

        static const int bitmask = numnodes;

        int hi = treebits - 32 + count_high_zero_bits<uint32_t>(id);
        DBG(debug, "high zero: " << hi);

        unsigned int bkt = ((id << (hi + 1)) & bitmask) | (1 << hi);

        DBG(debug, "bkt: " << bkt << " = " << toBinary(bkt));
        return bkt;
    }

    static inline unsigned int
                       in_to_levelorder(unsigned int id)
    {
        assert(id > 0);
        DBG(debug, "index: " << id << " = " << toBinary(id));

        static const int bitmask = numnodes;

        int lo = count_low_zero_bits<uint32_t>(id) + 1;
        DBG(debug, "low zero: " << lo);

        unsigned int bkt = ((id >> lo) & bitmask) | (1 << (treebits - lo));

        DBG(debug, "bkt: " << bkt << " = " << toBinary(bkt));
        return bkt;
    }

    static inline void self_verify()
    {
        for (size_t i = 1; i <= numnodes; ++i)
        {
            std::cout << toBinary(i, treebits) << " -> " << std::endl;

            size_t id = level_to_inorder(i);
            std::cout << toBinary(id, treebits) << " -> " << std::endl;

            id = in_to_levelorder(id);
            std::cout << toBinary(id, treebits) << std::endl;

            std::cout << std::endl;
            assert(id == i);
        }
    }
};

static inline void self_verify_tree_calculations()
{
    TreeCalculations<4>::self_verify();
    TreeCalculations<5>::self_verify();
    TreeCalculations<6>::self_verify();
    TreeCalculations<11>::self_verify();
    TreeCalculations<12>::self_verify();
    TreeCalculations<13>::self_verify();
    TreeCalculations<14>::self_verify();
    TreeCalculations<15>::self_verify();
    TreeCalculations<16>::self_verify();
}

} // namespace stringtools

#include "stringptr.h"

#endif // !PSS_SRC_TOOLS_STRINGTOOLS_HEADER

/******************************************************************************/
