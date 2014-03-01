/******************************************************************************
 * src/sequential/bingmann-sample_sort.h
 *
 * Experiments with sequential Super Scalar String Sample-Sort (S^5).
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

namespace bingmann_sample_sort {

static const bool debug = false;
static const bool debug_splitter = false;
static const bool debug_bucketsize = false;
static const bool debug_recursion = false;
static const bool debug_splitter_tree = false;

using namespace stringtools;

typedef uint64_t key_type;

static const size_t l2cache = 256*1024;

static const size_t g_samplesort_smallsort = 32*1024;

static const size_t oversample_factor = 2;

static size_t g_ss_steps, g_rs_steps;

enum { TM_GENERAL, TM_MAKE_SAMPLE, TM_MAKE_SPLITTER, TM_CLASSIFY, TM_PREFIXSUM, TM_PERMUTE, TM_SMALLSORT };
static TimerArray g_timer(16);

static inline void sample_sort_pre()
{
    g_ss_steps = g_rs_steps = 0;
    g_timer.clear();
}

static inline void sample_sort_post()
{
    g_stats >> "l2cache" << l2cache
            >> "steps_sample_sort" << g_ss_steps
            >> "steps_base_sort" << g_rs_steps;

    g_stats >> "tm_general" << g_timer.get(TM_GENERAL)
            >> "tm_make_sample" << g_timer.get(TM_MAKE_SAMPLE)
            >> "tm_make_splitter" << g_timer.get(TM_MAKE_SPLITTER)
            >> "tm_classify" << g_timer.get(TM_CLASSIFY)
            >> "tm_prefixsum" << g_timer.get(TM_PREFIXSUM)
            >> "tm_permute" << g_timer.get(TM_PERMUTE)
            >> "tm_smallsort" << g_timer.get(TM_SMALLSORT);
}

} // namespace bingmann_sample_sort
