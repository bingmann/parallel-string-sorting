/*******************************************************************************
 * src/parallel/bingmann-parallel_sample_sort_lcp.hpp
 *
 * Parallel Super Scalar String Sample-Sort, many variant via different
 * Classifier templates and LCP construction.
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

#ifndef PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_LCP_HEADER
#define PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_LCP_HEADER

// we will reinclude the main PSS header
#ifdef PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_HEADER
#undef PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_HEADER
#endif

#ifdef CALC_LCP
#undef CALC_LCP
#endif

#define CALC_LCP 1

#include "bingmann-parallel_sample_sort.hpp"

#endif // !PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_LCP_HEADER

/******************************************************************************/
