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

/// represent hex octets of larger datatypes
template <typename Type>
std::string toHex(const Type& _v)
{
    Type v = _v;
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
    return __builtin_clz(t);
}

template <>
inline int count_high_zero_bits<uint64_t>(const uint64_t& t)
{
    return __builtin_clzll(t);
}

template <>
inline int count_high_zero_bits<uint128_t>(const uint128_t& t)
{
    uint64_t hi = (t >> 64);
    int b = __builtin_clzll(hi);
    if (b > 0)
        return b;
    else
        return 64 + __builtin_clzll( (uint64_t)t );
}

} // namespace stringtools
