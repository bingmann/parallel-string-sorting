/*******************************************************************************
 * src/parallel/bingmann-parallel_sample_sort.h
 *
 * Parallel Super Scalar String Sample-Sort, many variant via different
 * Classifier templates.
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

#ifndef PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_HEADER
#define PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_HEADER

#include "../tools/stringtools.h"

namespace bingmann_parallel_sample_sort_lcp {

using namespace stringtools;

void parallel_sample_sort_numa(string* strings, size_t n,
                               int numaNode, int numberOfThreads,
                               const LcpCacheStringPtr& output);

void parallel_sample_sort_numa2(const StringShadowLcpCacheOutPtr* input,
                                unsigned numInputs);

} // namespace bingmann_parallel_sample_sort_lcp

#endif // !PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_HEADER

/******************************************************************************/
