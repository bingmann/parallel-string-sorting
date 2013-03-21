/* -*- tab-width: 8 -*- */

#include "copy.h"

/*
 * Input command line, process arguments and initiate testbed program
 *	for demonstrating C-, CP- and CPL-burstsorts.  Displays help screen
 * if any unrecognized argument is encountered.
*/
int main(int argc, char *argv[])
{
	int i, j, nreps, ntoks, set, dset1, srt, v[20], sn[20] = {10, 11, 12, 20, 21, 22, 30, 31, 32, -1};
	string VP = NULL, tok[20];

	string dnam[] = {"set1_nodup", "set2_nodup", "set3_nodup", 
		 "set4_nodup", "set5_nodup", "set6_nodup", "set1_genome", 
		 "set2_genome", "set3_genome", "set4_genome", "set5_genome", 
		 "set6_genome", "set1_random", "set2_random", "set3_random", 
		 "set4_random", "set5_random", "set6_random", "set1_url", 
		 "set2_url", "set3_url", "set4_url", "set5_url"};

	/* set default values for command line args. */
    INPUTDIR = "data/";
    OUTPUTDIR = "sort-output/";
	MAXKEYS = 100000000;	/* maximum keys to load */
	MAXBYTES = 700000000;	/* maximum bytes to load */
	CACHESIZE = 1<<19;

	dset1 = 0;				/* Start at this index in dnam[] above. */
	INPUTORDER = 0;			/* -1/0/1 for sorted/random/reverse order */
	NSEGS0 = 50;			/* number of segments in input buffer */
	SEGSIZE0 = 1<<15;		/* initial segment size; will be adjusted */
	BINSIZE0 = 1<<9;		/* initial size of buffer in terminal node */	
	SAMPLERATE0 = 4000;		/* if >0, will sample 1/SAMPLERATE of data */
	FREEBURSTS0 = 100;		/* number of low threshold bursts to allow */
	TAILSIZE0 = 12;		    /* maximum tail size in CPL-burstsort */
	TAILRATE = 80;			/* used with TAILSIZE0 in CPL-burstsort */
	nreps = 1;				/* desired number of repeats */

	/* Load and parse command line arguments */
	for (i = 1; i < argc; ++i)
	{
		/* break up each white space delimited arg into tokens
			separated by "=" or ",". */
		ntoks = tokenize(argv[i], tok, ',', '=');
		if		  (strcmp(tok[0], "ds") == 0) {sscanf(tok[1], "%d", &dset1);}
		else if (strcmp(tok[0], "nr") == 0) {sscanf(tok[1], "%d", &nreps);}
		else if (strcmp(tok[0], "cs") == 0) {sscanf(tok[1], "%d", &CACHESIZE);}
		else if (strcmp(tok[0], "wr") == 0) {sscanf(tok[1], "%d", &WRITEFILE);}
		else if (strcmp(tok[0], "fb") == 0) {sscanf(tok[1], "%d", &FREEBURSTS0);}
		else if (strcmp(tok[0], "mb") == 0) {sscanf(tok[1], "%d", &MAXBYTES);}
		else if (strcmp(tok[0], "mk") == 0) {sscanf(tok[1], "%d", &MAXKEYS);}
		else if (strcmp(tok[0], "g0") == 0) {sscanf(tok[1], "%d", &NSEGS0);}
		else if (strcmp(tok[0], "or") == 0) {sscanf(tok[1], "%d", &INPUTORDER);}
		else if (strcmp(tok[0], "sr") == 0) {sscanf(tok[1], "%d", &SAMPLERATE0);}
		else if (strcmp(tok[0], "s0") == 0) {sscanf(tok[1], "%d", &BINSIZE0);}
		else if (strcmp(tok[0], "tf") == 0) {sscanf(tok[1], "%d", &TAILRATE);}
		else if (strcmp(tok[0], "tl") == 0) {sscanf(tok[1], "%d", &TAILSIZE0);}
		else if (strcmp(tok[0], "sn") == 0)
			/* load list of sort numbers following arg "sn" */
			{
				 for (srt = 0; srt < ntoks - 1; ++srt)
					  sscanf(tok[srt + 1], "%d", &sn[srt]);
				 sn[srt] = -1;
			}
			
		else if (strcmp(tok[0], "vl") == 0)
			/* load list of alternative values for arg introduced
				by arg "vl"; the first token after "vl" is the tag 
				for the arg, and the remaining tokens are the values. */
			{
				 VP = tok[1];
				 for(j = 2; j < ntoks; ++j)
					  sscanf(tok[j], "%d", &v[j-2]); v[ntoks-2] = -1;
			}
		else {showhelp(); return(0);}
	}

	/* Print row of titles for output fields. */
	sayln("DS SORTNAME CACHESIZE FREEBURSTS . NKEYS . Ti Tg Tb Tt . Tmed Tmin Tnorm ");

	/* Run the specified tests for each data set in turn. */
	for (set = dset1; set <= dset1; ++set)
	{
		DATANAME = dnam[set]; 
		
		/* For each data set, run the specified sort variants. */
		for (srt = 0; sn[srt] >= 0; ++srt)
		{
			/* Initialize the sort variant and do the specified number
				of repeats; if VP is not null, the sort will be repeated
				using a list of alternative values for the arg specified
				in VP. */
			if (VP == NULL) {setsort(sn[srt]); doreps(nreps);}
				else
				{
					/* for (i = 0; v[i] >= 0; ++i) */
					/* { */
						if (strcmp(VP, "cs") == 0) CACHESIZE = v[i];
						else if (strcmp(VP, "fb") == 0) FREEBURSTS0 = v[i];
						else if (strcmp(VP, "mb") == 0) MAXBYTES = v[i];
						else if (strcmp(VP, "mk") == 0) MAXKEYS = v[i];
						else if (strcmp(VP, "g0") == 0) NSEGS0 = v[i];
						else if (strcmp(VP, "or") == 0) INPUTORDER = v[i];
						else if (strcmp(VP, "sr") == 0) SAMPLERATE0 = v[i];
						else if (strcmp(VP, "s0") == 0) BINSIZE0 = v[i];
						else if (strcmp(VP, "tf") == 0) TAILRATE = v[i];
						else if (strcmp(VP, "tl") == 0) TAILSIZE0 = v[i];
						setsort(sn[srt]); 
						doreps(nreps);
					/* } */
				}
		} 
		NBYTES = MAXALLOCATED = ALLOCATED = 0;
	} 
	return(0);

} /* End of main function */

/*
 * Attempts to explain itself.
 *
 */
void showhelp()
{
	sayln(" ??\n");
	sayln("Available arguments (use in any order) are:\n");
	sayln("	nr		# of repeats to run			ds	  first data set to run");
	sayln("	cs		cache size(bytes)			s0	  initial container size");
	sayln("	fb		# of free bursts			g0	  # of buffer segments");
	sayln("	mk		maximum keys				mb	  maximum bytes");
	sayln("	or		order of input data			sr	  sampling rate");
	sayln("	tl		tail length for CPL-bs		tf	  tail factor for CPL-bs");
	sayln("	wr		output sorted strings		sn    sort variant(s) to run");

	sayln("Arguments are separated by white space and may not contain white space. Each argument must");
	sayln("start with a tag followed by '=' or ', ' and one or more numbers or additional tags.");
	sayln("Tag 'sn' takes a list of sort numbers.  Tag 'vl' takes a second tag, followed by a list");
	sayln("of up to 20 values for the argument specified by the second tag.\n");

	sayln("Available sorts are:\n");
	sayln("	 10	C-bs		Copying burstsort");
	sayln("	 11	Cf-bs		Same + free bursts set by arg FREEBURSTS");
	sayln("	 12	Cs-bs		Same + sampling set by arg SAMPLERATE");
	sayln("	 20	CP-bs		Record burstsort");
	sayln("	 21	CPf-bs		Same + free bursts set by arg FREEBURSTS");
	sayln("	 22	CPs-bs		Same + sampling set by arg SAMPLERATE");
	sayln("	 30	CPL-bs		Paging burstsort");
	sayln("	 31	CPLf-bs		Same + free bursts set by arg FREEBURSTS");
	sayln("	 32	CPLs-bs		Same + sampling set by arg SAMPLERATE");

	sayln("Available data sets are:\n");

	sayln("	 0	Set 1 of No Duplicates (100, 000 strings)");
	sayln("	 1	Set 2 of No Duplicates (316, 230 strings)");
	sayln("	 2	Set 3 of No Duplicates (1, 000, 000 strings)");
	sayln("	 3	Set 4 of No Duplicates (3, 162, 300 strings)");
	sayln("	 4	Set 5 of No Duplicates (10, 000, 000 strings)");
	sayln("	 5	Set 6 of No Duplicates (31, 623, 000 strings)\n");

	sayln("	 6	Set 1 of Genome (100, 000 strings)");
	sayln("	 7	Set 2 of Genome (316, 230 strings)");
	sayln("	 8	Set 3 of Genome (1, 000, 000 strings)");
	sayln("	 9	Set 4 of Genome (3, 162, 300 strings)");
	sayln("	 10 Set 5 of Genome (10, 000, 000 strings)");
	sayln("	 11 Set 6 of Genome (31, 623, 000 strings)\n");

	sayln("	 12	 Set 1 of Random (100, 000 strings)");
	sayln("	 13	 Set 2 of Random (316, 230 strings)");
	sayln("	 14	 Set 3 of Random (1, 000, 000 strings)");
	sayln("	 15	 Set 4 of Random (3, 162, 300 strings)");
	sayln("	 16	 Set 5 of Random (10, 000, 000 strings)");
	sayln("	 17	 Set 6 of Random (31, 623, 000 strings)\n");

	sayln("	 18	 Set 1 of URL (100, 000 strings)");
	sayln("	 19	 Set 2 of URL (316, 230 strings)");
	sayln("	 20	 Set 3 of URL (1, 000, 000 strings)");
	sayln("	 21	 Set 4 of URL (3, 162, 300 strings)");
	sayln("	 22	 Set 5 of URL (10, 000, 000 strings)\n");

	sayln("	 23	Set 3 of Random Strings, each of length 100 characters\n");

	sayln("Usage: sn=10, 11, 20 ds=5 nr=1 will run sorts 10, 11 and 20 once on dataset 5.\n\n");
} /* End of showhelp function */

/*
 * Takes a sort number specified in the command line args, and selects
 * and initializes the specified sort variant.
 *
 */
void setsort(int sn)
{
	FREEBURSTS = SAMPLERATE = 0; 
	
	switch(sn)
	{
		case 10: TEST = sbdo; SORTNAME = "C-burstsort"; 
					break;
		case 11: TEST = sbdo; SORTNAME = "fbC-burstsort"; 
					FREEBURSTS = FREEBURSTS0; break;
		case 12: TEST = sbdo; SORTNAME = "sC-burstsort"; 
					SAMPLERATE = SAMPLERATE0; break;

		case 20: TEST = rbdo; SORTNAME = "CP-burstsort"; 
					break;
		case 21: TEST = rbdo; SORTNAME = "fbCP-burstsort"; 
					FREEBURSTS = FREEBURSTS0; break;
		case 22: TEST = rbdo; SORTNAME = "sCP-burstsort"; 
					SAMPLERATE = SAMPLERATE0; break;

		case 30: TEST = pbdo; SORTNAME = "CPL-burstsort"; 
					break;
		case 31: TEST = pbdo; SORTNAME = "fbCPL-burstsort"; 
					FREEBURSTS = FREEBURSTS0; break;
		case 32: TEST = pbdo; SORTNAME = "sCPL-burstsort"; 
					SAMPLERATE = SAMPLERATE0; break;
	}
}

/*
 * Output parameters set for a particular sort run. 
 *
 */
void listinputs()
{
	 if (!REP)
	 {
		  say(DATANAME); say(SORTNAME); isay(CACHESIZE); 
		  isay(FREEBURSTS); dot();
	 }
}

/*
 * Do the specified number of repeats for a run involving a data set, 
 * a sort variant and a set of parameter values (in global variables
 * set in main() or setsort()).
 *
 */
void doreps(int nr)
{
	double rt, pt, st, it, gt, bt, ct, ht , kt, wt, ist, xst, mo, mbs;

	/* Reset timing system. */
	treset(MAXTMRS, nr);
	rt = pt = st = it = gt = bt = ct = ht = 
	kt = wt = ist = xst = mo = mbs = 0;

	/* Based on arg INPUTORDER, set FILETYPE to load an unsorted, 
		presorted or reverse sorted input data file. */
	if (INPUTORDER == 0)
		 FILETYPE = "dat";
	else if (INPUTORDER == 1)
		 FILETYPE = "srt";
	else
		 FILETYPE = "rev";
	
	/* Run and time the repetitions. */
	TEST(nr, fp(INPUTDIR, DATANAME, FILETYPE));

	isay(NKEYS); dot();

	/* Collect median times for each timing phase. */
	rt = tmed(tr); pt = tmed(tx); st = tmed(ts); it = tmed(ti); gt = tmed(tg);
	bt = tmed(tb); ct = tmed(tc); ht = tmed(th); kt = tmed(tk); wt = tmed(tw);
	ist = pt+st+it+gt+bt+ct; xst = ist+rt+wt;

	/* Calculate "memory overhead" = peak allocated bytes during sorting
		divided by bytes sorted. */
	mo = DR(MAXALLOCATED, NBYTES); 
	
	/* Calculate megabytes sorted per second. */
	mbs = DR(NBYTES, DR(ist, 1000));

	/* Find shortest overall time. */
	tmin(nr);

	/* Output times for insertion phase, growth phase, burst phase, 
		container sorting phase, total internal sorting time, shortest 
		time among iterations, and time normalized for n log n. */
	dsay(it, 0); dsay(gt, 0); dsay(bt, 0); dsay(ct, 0); dot();
	dsay(ist, 0); tminsay(); tnormsay(); cr();
}

/*
* Load input data file as a number of segments, each of which can be freed
* once its keys have been inserted in the burst trie; used by C-burstsort,
* which sorts only keys; CP- and CPL-burstsort, which sort records based on
* keys, use the following function getkeys().
*
*/
void getsegs(string filepath)
{
    FILE *f; string s, t; int c, delim, seg, seglen, nread, n, keylen, fbytes;

    if ((f = fopen(filepath, "rb")) == NULL)
    {
        say("couldn't open"); brp(filepath);
    }
    
    /* First check the length of the file. */
    fseek(f, 0, SEEK_END); fbytes = ftell(f); rewind(f);

    /* If the requested number of bytes is more than the file length, or is not
       set, adjust MAXBYTES to the file length. */
    
    if (MAXBYTES > fbytes || MAXBYTES == 0) MAXBYTES = fbytes;

    /* Set a reasonable segment size.  First, divide the number of bytes that will
       be read and sorted by the requested number of segments set in argument
       NSEGS0.  Adjust this to an even multiple of the optimal system I/O buffer
       size BUFSIZ.  Finally, make sure the segments are not too small.  (Small
       data files with large requested segment counts could otherwise create
       segments shorter than a long string key.) */
    
    SEGSIZE = MAXBYTES / NSEGS0;
    SEGSIZE -= (SEGSIZE % BUFSIZ);

    while (SEGSIZE < BINSIZE0) SEGSIZE *= 2;

    /* Now calculate how many segments we will actually end up with. */
    NSEGS = MAXBYTES / SEGSIZE; if (MAXBYTES % SEGSIZE) ++NSEGS;
    
    /* Set up pointers to the beginning and end of each segment. */
    SEGMENTS = (string*) CALLOC(NSEGS, sizeof(string));
    SEGLIMITS = (string*) CALLOC(NSEGS, sizeof(string));
   
    /* Clear array used to count char usage. */
    memset((char*) CHARCOUNT, 0, 256 * sizeof(int));

    /* Clear global variables that will record the number of bytes and keys
       loaded and the length of the longest key loaded. */
    NBYTES = NKEYS = MAXKEYLEN = 0;
    
    /* Create the segments and load the data. */
    for (seg = 0; seg < NSEGS; ++seg) {
        s = SEGMENTS[seg] = (char*) MALLOC(SEGSIZE);
        
        /* Read SEGSIZE bytes unless that would exceed MAXBYTES. */
        nread=MIN(SEGSIZE, MAXBYTES - NBYTES);
        nread=fread(SEGMENTS[seg], 1, nread, f);
        
        /* If this is the first segment, we need to determine the EOL char in
           use; the following loop scans forward until it encounters '\n' or '\r'
           and assumes that char to be the delimiter for the rest of the file. */
        if (seg == 0)
        {
            t = s; while ((c = *t) != '\n' && c != '\r') ++t; delim = c;
        }
        
        /* keylen is the length of the current key; seglen is the length of the
           current segment, and is incremented only in whole key lengths so that
           segments always end at a key boundary. */
        n = keylen = seglen = 0;

        /* Now scan keys until the end of the read, or until the byte or key
           limit is hit. */
        while (n < nread) {
            if ((c = s[n]) == delim) {
                /* This is a key boundary, so update the segment length. */
                s[n++] = 0; seglen = n; ++NKEYS;

                /* If this is the longest key so far, update MAXKEYLEN. */
                if (++keylen > MAXKEYLEN) MAXKEYLEN = keylen; keylen = 0;
            } 
            else 
            {
                /* Replace any nonprinting chars with '.'. */
                if (c < 32 || c == 127) s[n] = c = 46;
                CHARCOUNT[c] = 1; ++keylen; ++n;
            }                       
        }
            
        SEGLIMITS[seg] = s += seglen; NBYTES += seglen;
        
        /* If the last read did not end at a key boundary, this seek moves back
           to the first byte of the incomplete key. */
        fseek(f, NBYTES, SEEK_SET);
    }

    /* MAXKEYLEN will be used to determine when a bin (terminal node) in a burst
       trie does not have enough room to add another key; here we set LIMSIZE0,
       the offset of the limit pointer in a new bin created at BINSIZE0, to
       BINSIZE0 - MAXKEYLEN. */
    LIMSIZE0 = BINSIZE0 - MAXKEYLEN;

    /* Scan the counts of char usage to find the lowest and highest chars
       actually encountered; only this range need be scanned in traversing
       the burst trie. */
    LOCHAR = 32; while (!CHARCOUNT[LOCHAR]) ++LOCHAR;
    HICHAR = 255; while (!CHARCOUNT[HICHAR]) --HICHAR;
    
} /* End of getsegs function */


/*
 * Load input data file for use in CP- or CPL-burstsort, which sort records 
 * based on string keys; since these sorts need to keep the original record
 * data in memory until all keys are sorted, the segmented approach used in 
 * getsegs() is not applicable. 
 *
 */
void getkeys(string filepath)
{
	FILE *f; string s, lim; int c, delim, keylen;

	if ((f = fopen(filepath, "rb")) == NULL) 
		{say("couldn't open"); brp(filepath);} 
	
	/* The two extra bytes leave room for terminating chars before and after
		the data; see below for why. */
	INBUF = (char *) MALLOC(MAXBYTES + 2); 
	*INBUF = 0; ++INBUF;

	/* If file is shorter than MAXBYTES, adjust MAXBYTES. */
	MAXBYTES=fread(INBUF, 1, MAXBYTES, f);
	fclose(f);

	/* Determine what char is being used as EOL delimiter. */
	s = INBUF;
	while ((c = *s) != '\n' && c != '\r') ++s; 
	delim = c; 
	
	/* Place a delimiter so that the last key will be truncated at MAXBYTES;
		this behavior is inconsistent with getsegs(), which discards any final
		partial key instead of truncating it. */
	/* lim = INBUF + MAXBYTES; *(lim - 1) = delim;  */
	
	/* To get behavior consistent with getsegs(), replace the line above with: */
		lim = INBUF + MAXBYTES; while (*(lim - 1) != delim) --lim;
	
	/* Clear array used to count char usage. */
	memset((char*) CHARCOUNT, 0, 256*sizeof(int));

	/* Clear global variables that will record the number of bytes and keys
		loaded and the length of the longest key loaded. */
	NBYTES = NKEYS = MAXKEYLEN = 0;

	/* comments */
	s = INBUF;
	while (s<lim && NKEYS<MAXKEYS)
	{
		keylen = 0;
		while ((c = *s) != delim)
		{
			 /* Replace any nonprinting chars with '.'. */
			 if (c < 32 || c == 127) *s = c = 46; 
			 CHARCOUNT[c] = 1; ++keylen; ++s;
		}

		*s++ = 0; ++NKEYS;
		
		/* If this is the longest key so far, update MAXKEYLEN. */
		if (++keylen > MAXKEYLEN) MAXKEYLEN = keylen;
	}

	/* Change the delimiter following the input buffer to a null terminator. */
	*s = 0; 
	
	NBYTES = s - INBUF; 
	
	/* MAXKEYLEN will be used to determine when a bin (terminal node) in a burst
		trie does not have enough room to add another key; here we set LIMSIZE0, 
		the offset of the limit pointer in a new bin created at BINSIZE0, to
		BINSIZE0 - MAXKEYLEN. */
	LIMSIZE0 = BINSIZE0 - MAXKEYLEN; 
	
	/* Scan the char counts to find the lowest and highest chars encountered;
		only this range need be scanned in traversing the burst trie. */
	LOCHAR = 32; while (!CHARCOUNT[LOCHAR]) ++LOCHAR; 
	HICHAR = 255; while (!CHARCOUNT[HICHAR]) --HICHAR;
} /* End of getkeys function */
