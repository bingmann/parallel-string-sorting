/*******************************************************************************
 * src/tools/globals.cpp
 *
 * Global variables for parallel string sorting algorithms, other than those in
 * psstest.cpp
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

#include "globals.hpp"

stats_writer g_stats;

int g_num_threads = 0;

// number of NUMA nodes (may be faked by cmdline)
size_t g_numa_nodes;

// argument -M, --memory, see tools/input.h
std::string gopt_memory_type;

/******************************************************************************/
