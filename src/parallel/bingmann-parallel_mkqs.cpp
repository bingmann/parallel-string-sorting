/*******************************************************************************
 * src/parallel/bingmann-parallel_mkqs.cpp
 *
 * Parallel multikey-quicksort.
 *
 *******************************************************************************
 * Copyright (C) 2013-2015 Timo Bingmann <tb@panthema.net>
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

#include "bingmann-parallel_mkqs.hpp"

namespace bingmann_parallel_mkqs {

static inline void
bingmann_sequential_mkqs_cache8(unsigned char** strings, size_t n)
{
    return bingmann_sequential_mkqs_cache8(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT(bingmann_sequential_mkqs_cache8,
               "bingmann/sequential_mkqs_cache8",
               "multikey_cache with 8byte cache (non-recursive)")

static inline
void bingmann_parallel_mkqs(unsigned char** strings, size_t n)
{
    return bingmann_parallel_mkqs(
        parallel_string_sorting::UCharStringSet(strings, strings + n), 0);
}

PSS_CONTESTANT_PARALLEL(bingmann_parallel_mkqs,
                        "bingmann/parallel_mkqs",
                        "Parallel MKQS with blocks and cache8")

} // namespace bingmann_parallel_mkqs

/******************************************************************************/
