/******************************************************************************
 * src/sinha-copy-burstsort/glue.cc
 *
 * Glue to attach Sinha's C-Burstsort code to parallel-string-sort framework.
 *
 ******************************************************************************
 * Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
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

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>

extern const char* g_string_data;
extern size_t g_string_datasize;

namespace sinha_copy_burstsort {

#include "C-burstsort.c"
#include "CP-burstsort.c"
#include "CPL-burstsort.c"
#include "utils.c"
//#include "copy-burstsort.c"

void getsegs(string) {}
void getkeys(string) {}
void listinputs() {}

void sinha_C_burstsort_prepare(string *strings, size_t size)
{
    char CHARCOUNT[256] = { 0 };
    MAXKEYLEN = 0;

    for (size_t i = 0; i < size; ++i)
    {
        size_t j;
        for (j = 0; strings[i][j]; ++j)
        {
            CHARCOUNT[ strings[i][j] ] = 1;
        }
        if (j > MAXKEYLEN) MAXKEYLEN = j;
    }
    LOCHAR = 1; while (!CHARCOUNT[LOCHAR]) ++LOCHAR;
    HICHAR = 255; while (!CHARCOUNT[HICHAR]) --HICHAR;
}

void sinha_C_burstsort(string *strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 0;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS0 = 0;		/* number of low threshold bursts to allow */

    NKEYS = size;
    NBYTES = g_string_datasize;
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Run C-burstsort on the data that was loaded in getsegs() above. */

    /* Create root node of burst trie. */
    NODE1 *rt = GETNODE(NODE1); ++NODES;

    /* Insert a group of keys; seg[] and seglim[] are set up so that
       segments begin and end at key boundaries. */
    sadd((string)g_string_data, (string)g_string_data + g_string_datasize, rt);

    /* Discard any remaining free bursts; see below in straverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    straverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_ssavetrie(rt, strings, (string)g_string_data);

    /* Deallocate the root node. */
    skill(rt);
}

void sinha_fbC_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 0;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS = 100;		/* number of low threshold bursts to allow */

    NKEYS = size;
    NBYTES = g_string_datasize;
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Run C-burstsort on the data that was loaded in getsegs() above. */

    /* Create root node of burst trie. */
    NODE1 *rt = GETNODE(NODE1); ++NODES;

    /* Insert a group of keys; seg[] and seglim[] are set up so that
       segments begin and end at key boundaries. */
    sadd((string)g_string_data, (string)g_string_data + g_string_datasize, rt);

    /* Discard any remaining free bursts; see below in straverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    straverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_ssavetrie(rt, strings, (string)g_string_data);

    /* Deallocate the root node. */
    skill(rt);
}

void sinha_sC_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 4000;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS0 = 0;		/* number of low threshold bursts to allow */

    NKEYS = size;
    NBYTES = g_string_datasize;
    NSEGS = 1;
    string mySEGMENTS = { (string)g_string_data };
    string mySEGLIMITS = { (string)g_string_data + g_string_datasize };

    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Run C-burstsort on the data that was loaded in getsegs() above. */

    /* Create root node of burst trie. */
    NODE1 *rt = GETNODE(NODE1); ++NODES;

    /* Sampling mode prebuilds a trie based on a random sample of keys. */
    if (SAMPLERATE)
    {
        /* Build trie using 1/SAMPLERATE of the keys and bursting any bin
           that fills at BINSIZE0 rather than the usual threshold. */
        ssample(&mySEGMENTS, &mySEGLIMITS, rt, NKEYS/SAMPLERATE);

        /* Clear data but retain node/bin structure of the "skeleton trie". */
        clear(rt);
    }

    /* Insert a group of keys; seg[] and seglim[] are set up so that
       segments begin and end at key boundaries. */
    sadd((string)g_string_data, (string)g_string_data + g_string_datasize, rt);

    /* Discard any remaining free bursts; see below in straverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    straverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_ssavetrie(rt, strings, (string)g_string_data);

    /* Deallocate the root node. */
    skill(rt);
}

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_C_burstsort,
                            "sinha/C_burstsort", "Original C-burstsort");

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_fbC_burstsort,
                            "sinha/fbC_burstsort", "Original fbC-burstsort");

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_sC_burstsort,
                            "sinha/sC_burstsort", "Original sC-burstsort");

void sinha_CP_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 0;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS0 = 0;		/* number of low threshold bursts to allow */

    NKEYS = size;
    NBYTES = g_string_datasize;
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Create root node of burst trie. */
    NODE1 *rt = GETNODE(NODE1); ++NODES;

    /* Read string sequence. */
    radd((string)g_string_data, NKEYS, rt);

    /* Discard any remaining free bursts; see below in rtraverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    rtraverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_savetrie(rt, (string*)strings);

    /* Deallocate the root node. */
    rkill(rt);
}

void sinha_fbCP_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 0;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS = 100;		/* number of low threshold bursts to allow */

    NKEYS = size;
    NBYTES = g_string_datasize;
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Create root node of burst trie. */
    NODE1 *rt = GETNODE(NODE1); ++NODES;

    /* Read string sequence. */
    radd((string)g_string_data, NKEYS, rt);

    /* Discard any remaining free bursts; see below in rtraverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    rtraverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_savetrie(rt, (string*)strings);

    /* Deallocate the root node. */
    rkill(rt);
}

void sinha_sCP_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 4000;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS0 = 0;		/* number of low threshold bursts to allow */

    NKEYS = size;
    NBYTES = g_string_datasize;
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Create root node of burst trie. */
    NODE1 *rt = GETNODE(NODE1); ++NODES;

    if (SAMPLERATE)
    {
        rsample((string)g_string_data, rt, NKEYS / SAMPLERATE);
        clear(rt);
    }

    /* Read string sequence. */
    radd((string)g_string_data, NKEYS, rt);

    /* Discard any remaining free bursts; see below in rtraverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    rtraverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_savetrie(rt, (string*)strings);

    /* Deallocate the root node. */
    rkill(rt);
}

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_CP_burstsort,
                            "sinha/CP_burstsort", "Original CP-burstsort");

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_fbCP_burstsort,
                            "sinha/fbCP_burstsort", "Original fbCP-burstsort");

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_sCP_burstsort,
                            "sinha/sCP_burstsort", "Original sCP-burstsort");

/* Initialization from pbdo() */
void sinha_CPL_burstsort_initialize()
{
    /* The command line arg TAILRATE specifies a
       percentage of the average key length to use as
       the key segment size; we calculate this now. */
    double d = TAILRATE; d /= 100;
    d *= NBYTES; d /= NKEYS;

    /* Command line arg TAILSIZE0 specifies a maximum
       key segment size; we choose the smaller of
       TAILSIZE0 and TAILSIZE calculated above. */
    TAILSIZE = MIN(d, TAILSIZE0);

    /* Both TAILSIZE and TAILSIZE + 4 are used in many
       operations, so we precalculate both; the extra
       four bytes is room for a pointer to the original
       record, which is kept immediately before the
       TAILSIZE bytes of suffix. */
    TAILSIZE4 = TAILSIZE + sizeof(string);

    /* In contrast to C- and CP- burstsort, the key
       segments inserted in bins during CPL-burstsort
       are of known fixed length and we can calculate
       ahead of time how many will fit for each bin
       size. */
    int	lv = 0; int sz = BINSIZE0;
    while (sz < CACHESIZE)
    {
        THR[lv] = sz / TAILSIZE4;
        sz *= 2;
        ++lv;
    }

    /* For the preceding levels, we allowed room for
       TAILSIZE bytes of each key plus a record pointer;
       when the bin reaches cache size, we also need to
       reserve room for the sorting pointerss that will
       be needed during container sorting.*/
    THR[lv] = CACHESIZE / (TAILSIZE + 8);
}

void sinha_CPL_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 0;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS0 = 0;		/* number of low threshold bursts to allow */
    TAILSIZE0 = 12;             /* maximum tail size in CPL-burstsort */
    TAILRATE = 80;		/* used with TAILSIZE0 in CPL-burstsort */

    NKEYS = size;
    NBYTES = g_string_datasize;
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    sinha_CPL_burstsort_initialize();

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Create root node of burst trie. */
    NODE2 *rt = GETNODE(NODE2); ++NODES;

    /* Read string sequence. */
    padd((string)g_string_data, NKEYS, rt);

    /* Discard any remaining free bursts; see below in rtraverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    ptraverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_psavetrie(rt, (string*)strings);

    /* Deallocate the root node. */
    pkill(rt);
}

void sinha_fbCPL_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 0;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS = 100;		/* number of low threshold bursts to allow */
    TAILSIZE0 = 12;             /* maximum tail size in CPL-burstsort */
    TAILRATE = 80;		/* used with TAILSIZE0 in CPL-burstsort */

    NKEYS = size;
    NBYTES = g_string_datasize;
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    sinha_CPL_burstsort_initialize();

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Create root node of burst trie. */
    NODE2 *rt = GETNODE(NODE2); ++NODES;

    /* Read string sequence. */
    padd((string)g_string_data, NKEYS, rt);

    /* Discard any remaining free bursts; see below in rtraverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    ptraverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_psavetrie(rt, (string*)strings);

    /* Deallocate the root node. */
    pkill(rt);
}

void sinha_sCPL_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */
    while (BINSIZE0 < MAXKEYLEN) BINSIZE0 *= 2; /* must fix at least one key */

    SAMPLERATE = 4000;		/* if >0, will sample 1/SAMPLERATE of data */
    FREEBURSTS0 = 0;		/* number of low threshold bursts to allow */
    TAILSIZE0 = 12;             /* maximum tail size in CPL-burstsort */
    TAILRATE = 80;		/* used with TAILSIZE0 in CPL-burstsort */

    NKEYS = size;
    NBYTES = g_string_datasize;
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    sinha_CPL_burstsort_initialize();

    /* Zero counters for burst trie objects. */
    NODES = BINS = NULLBINS = 0;

    /* Create root node of burst trie. */
    NODE2 *rt = GETNODE(NODE2); ++NODES;

    if (SAMPLERATE)
    {
        psample((string)g_string_data, rt, NKEYS / SAMPLERATE);
        pclear(rt);
    }

    /* Read string sequence. */
    padd((string)g_string_data, NKEYS, rt);

    /* Discard any remaining free bursts; see below in rtraverse() for why. */
    FREEBURSTS = 0;

    /* Starting at the root node, traverse the trie and sort each bin. */
    ptraverse(rt);

    /* Write out string characters and string pointers to psstest */
    tb_psavetrie(rt, (string*)strings);

    /* Deallocate the root node. */
    pkill(rt);
}

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_CPL_burstsort,
                            "sinha/CPL_burstsort", "Original CPL-burstsort");

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_fbCPL_burstsort,
                            "sinha/fbCPL_burstsort", "Original fbCPL-burstsort");

CONTESTANT_REGISTER_PREPARE(sinha_C_burstsort_prepare, sinha_sCPL_burstsort,
                            "sinha/sCPL_burstsort", "Original sCPL-burstsort");

}
