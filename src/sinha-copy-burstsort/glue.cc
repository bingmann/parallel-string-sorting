/******************************************************************************
 * src/sinha-copy-burstsort/glue.cc
 *
 * Glue to attach Sinha's CopyBurstsort code to parallel-string-sort framework.
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

void sinha_C_burstsort_prepare(unsigned char **strings, size_t size)
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

    LOCHAR = 32; while (!CHARCOUNT[LOCHAR]) ++LOCHAR;
    HICHAR = 255; while (!CHARCOUNT[HICHAR]) --HICHAR;
}

void sinha_C_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */	

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
    tb_ssavetrie(rt, (string*)strings, (char*)g_string_data);

    /* Deallocate the root node. */
    skill(rt);
}

void sinha_fbC_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */	
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
    tb_ssavetrie(rt, (string*)strings, (char*)g_string_data);

    /* Deallocate the root node. */
    skill(rt);
}

void sinha_sC_burstsort(unsigned char **strings, size_t size)
{
    /* set default values for command line args. */
    CACHESIZE = 1<<19;
    BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */	

    SAMPLERATE = 4000;		/* if >0, will sample 1/SAMPLERATE of data */
    //FREEBURSTS0 = 100;		/* number of low threshold bursts to allow */

    NKEYS = size;
    NBYTES = g_string_datasize;
    NSEGS = 1;
    string mySEGMENTS = { (char*)g_string_data };
    string mySEGLIMITS = { (char*)g_string_data + g_string_datasize };

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
    tb_ssavetrie(rt, (string*)strings, (char*)g_string_data);

    /* Deallocate the root node. */
    skill(rt);
}

CONTESTANT_REGISTER_UCARRAY_PREPARE(sinha_C_burstsort_prepare, sinha_C_burstsort, "Original C-burstsort");

CONTESTANT_REGISTER_UCARRAY_PREPARE(sinha_C_burstsort_prepare, sinha_fbC_burstsort, "Original fbC-burstsort");

CONTESTANT_REGISTER_UCARRAY_PREPARE(sinha_C_burstsort_prepare, sinha_sC_burstsort, "Original sC-burstsort");

}
