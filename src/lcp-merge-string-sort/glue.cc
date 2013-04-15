/******************************************************************************
 * src/lcp-merge-string-sort/glue.cc
 *
 * Glue to attach Shamsundar's Parallel LCP Merge String Sort code to
 * parallel-string-sort framework.
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

#include "../tools/contest.h"

extern int pss_num_threads;

#define PTHREAD

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#ifdef _WIN32
#ifndef PTHREAD
#include <windows.h>
#endif
#include <io.h>
#else
#include <unistd.h>
#endif
#include <sys/types.h>
#include <fcntl.h>

namespace shamsundar_lcp_merge_string_sort {

#include "lcpsrtMT.c"

AS *sP;

static void prepare(unsigned char **strings, size_t size)
{
    sP = new AS[size];
    for (size_t i = 0; i < size; ++i) {
        sP[i].str = (char*)strings[i];
        sP[i].llcp = 0;
    }

    gSP = new AS*[size];
    for (size_t i = 0; i < size; ++i)
        gSP[i] = &sP[i];
}

static void string_sort(unsigned char **strings, size_t size)
{
    int cnts[4];
    int fcol=1,lcol=1000000000;

    NCPU = pss_num_threads;

    cnts[0]=cnts[1]=cnts[2]=cnts[3]=0;

    gpTMP = (pAS *)calloc((size+11)/2, sizeof(pAS)); 

    merge_sortMT(gSP,size,cnts,gpTMP,fcol,lcol);

    for (size_t i = 0; i < size; ++i) {
        strings[i] = (unsigned char*)gSP[i]->str;
    }

    free(gpTMP);
    delete [] sP;
    delete [] gSP;
}

CONTESTANT_REGISTER_PARALLEL_PREPARE(
    prepare, string_sort,
    "shamsundar/lcp-merge-string-sort",
    "Parallelized LCP Merge sort by N. Shamsundar");

// Note that the assembler version uses i686 instructions, and cannot be used
// on modern x86_64 bit architectures.

} // namespace shamsundar_lcp_merge_string_sort
