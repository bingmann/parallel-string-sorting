/*******************************************************************************
 * src/parallel/bingmann-parallel_radix_sort.cpp
 *
 * Parallel radix sort with work-balancing.
 *
 * The set of strings is sorted using a 8- or 16-bit radix sort
 * algorithm. Recursive sorts are processed in parallel using a lock-free job
 * queue and OpenMP threads. Two radix sort implementations are used:
 * sequential in-place and parallelized out-of-place.
 *
 * The sequential radix sort is implemented using an explicit recursion stack,
 * which enables threads to "free up" work from the top of the stack when other
 * threads become idle. This variant uses in-place permuting and a character
 * cache (oracle).
 *
 * To parallelize sorting of large buckets an out-of-place variant is
 * implemented using a sequence of three Jobs: count, distribute and
 * copyback. All threads work on a job-specific context called RadixStepCE,
 * which encapsules all variables of an 8-bit or 16-bit radix sort (templatized
 * with key_type).
 *
 *******************************************************************************
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
 ******************************************************************************/

#include <cstdlib>
#include <cstring>

#include <iostream>
#include <vector>

#include "bingmann-parallel_radix_sort.hpp"

namespace bingmann_parallel_radix_sort {

/******************************************************************************/
// Frontends

static inline void parallel_radix_sort_8bit(string* strings, size_t n)
{
    return parallel_radix_sort_8bit_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n),
        /* depth */ 0);
}

PSS_CONTESTANT_PARALLEL(parallel_radix_sort_8bit,
                        "bingmann/parallel_radix_sort_8bit",
                        "Parallel MSD Radix sort with load balancing, 8-bit BigSorts")

static inline void parallel_radix_sort_16bit(string* strings, size_t n)
{
    return parallel_radix_sort_16bit_generic(
        parallel_string_sorting::UCharStringSet(strings, strings + n),
        /* depth */ 0);
}

PSS_CONTESTANT_PARALLEL(parallel_radix_sort_16bit,
                        "bingmann/parallel_radix_sort_16bit",
                        "Parallel MSD Radix sort with load balancing, 16-bit BigSorts")

} // namespace bingmann_parallel_radix_sort

/******************************************************************************/
