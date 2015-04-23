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
    typedef ParallelMKQS<parallel_string_sorting::UCharStringSet> MKQS;

    MKQS::Context ctx;
    ctx.g_strings = strings;
    ctx.g_totalsize = n;

    MKQS::StrCache* cache = new MKQS::StrCache[n];
    for (size_t i = 0; i < n; ++i) {
        cache[i].str = strings[i];
    }
    JobQueue jobqueue;
    MKQS::SequentialJob<true>(ctx, strings, n, 0, cache).sequential_mkqs(jobqueue);
}

PSS_CONTESTANT(bingmann_sequential_mkqs_cache8,
               "bingmann/sequential_mkqs_cache8",
               "multikey_cache with 8byte cache (non-recursive)")

static inline
void bingmann_parallel_mkqs(unsigned char** strings, size_t n)
{
    typedef ParallelMKQS<parallel_string_sorting::UCharStringSet> MKQS;

    MKQS::Context ctx;
    ctx.g_strings = strings;
    ctx.g_totalsize = n;
    ctx.g_threadnum = omp_get_max_threads();
    ctx.g_sequential_threshold = std::max(g_inssort_threshold,
                                          ctx.g_totalsize / ctx.g_threadnum);

    g_stats >> "block_size" << block_size;

    JobQueue jobqueue;
    new MKQS::ParallelJob<MKQS::BlockSourceInput>(ctx, jobqueue, strings, n, 0);
    jobqueue.loop();
}

PSS_CONTESTANT_PARALLEL(bingmann_parallel_mkqs,
                        "bingmann/parallel_mkqs",
                        "Parallel MKQS with blocks and cache8")

} // namespace bingmann_parallel_mkqs

/******************************************************************************/
