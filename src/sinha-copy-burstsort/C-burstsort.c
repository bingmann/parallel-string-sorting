#include "copy.h"

/*
 * Run a specified number of iterations of C-burstsort.
 * This version of copying burstsort sorts string keys 
 * that are not part of larger records.  The 's' prefix
 * in the function names is from its original nickname
 * "suffix" burstsort (dropped to avoid any confusion
 * with suffix trees).
 *
 */
void sbdo(int nreps, string inputfilepath)
{
	for (REP=0; REP<nreps; ++REP)
	{
		/* Set timer phase to "read". */
		ton(tr);

		/* Zero current & maximum allocated memory counters. */
		MAXALLOCATED = ALLOCATED = 0;

		/* Load data from source file into segmented buffer. */
		getsegs(inputfilepath);
		
		/* Set timer phase to "prep". */
		ton(tx);

		/* Zero counters for burst trie objects. */
		NODES = BINS = NULLBINS = 0;

		/* List input values and settings for this run. */
		listinputs();

		/* Set timer phase to "insert". */
		ton(ti);

		/* Run C-burstsort on the data that was loaded in getsegs() above. */
		sbs(SEGMENTS, SEGLIMITS);
	}
	ton(0);
}

/*
 * Sort a buffer of string keys that has been divided into arbitrary seg.
 *
 */
void sbs(string *seg, string *seglim)
{
	int i; NODE1 *rt;

	/* Create root node of burst trie. */
	rt = GETNODE(NODE1); ++NODES;

	/* Sampling mode prebuilds a trie based on a random sample of keys. */
	if (SAMPLERATE)
	{
		 /* Set timer phase to "sample". */
		 ton(ts);

		 /* Build trie using 1/SAMPLERATE of the keys and bursting any bin
			  that fills at BINSIZE0 rather than the usual threshold. */
		 ssample(seg, seglim, rt, NKEYS/SAMPLERATE);

		 /* Clear data but retain node/bin structure of the "skeleton trie". */
		 clear(rt);

		 /* Set timer phase back to "insert". */
		 ton(ti);
	}

	/* Insert the keys, freeing each segment when its keys have been added. */
	for (i = 0;i < NSEGS; ++i)
	{
		 /* Insert a group of keys; seg[] and seglim[] are set up so that
			 segments begin and end at key boundaries. */
		 sadd(seg[i], seglim[i], rt);

		 /* If allocated memory has reached a new high, save it; then, free
			 the segment of keys that was just inserted. */
		 UPDATEMEM; FREE(seg[i], SEGSIZE);
	}

	/* When all segments are inserted, free the segmented buffer. */
	FREE(seg, NSEGS*sizeof(string)); FREE(seglim, NSEGS*sizeof(string));

	/* Discard any remaining free bursts; see below in straverse() for why. */
	FREEBURSTS = 0;

	/* Set timer phase to "container sort". */
	ton(tc);

	/* Starting at the root node, traverse the trie and sort each bin. */
	straverse(rt);

	UPDATEMEM;

	if (WRITEFILE)
	{
		 /* Set timer phase to "write" and save sorted data to disk file. */
		 ton(tw); ssavetrie(rt, OUTPUTDIR, DATANAME);
	}

	/* Set timer phase to "kill" and deallocate the root node. */
	ton(tk); skill(rt);
}

/*
 * Build a model burst trie using n keys randomly chosen from the current
 * data set; during sampling, bins do not grow -- all bins that become full
 * at their original size are burst.  Sampling and "free bursts" (see below)
 * are alternate strategies for reducing the copying costs incurred when
 * large bins are burst.
 */
void ssample(string *seg, string *seglim, NODE1 *rt, int n)
{
	NODE1 *nd; NODE1 *bn;
	int ct, i;
	char c; string b, s, lim;

	srand(clock());
	while (n--)
	{
		 /* Choose a random segment of keys */
		 i = rand()%NSEGS; lim = seglim[i];

		 /* and a random position within that segment; */
		 b = seg[i] + rand() % (lim-seg[i]);

		 /* then scan downward to the beginning of a key. */
		 while (*b != 0 && b > seg[i]) --b;

		 /* If we hit a null, the key begins at the next byte;
			 if we reached the start of the segment, 
			 a key begins on that exact byte. */
		 if (b>seg[i]) ++b;

		 /* Use the first char of the key to index from the root node to the
			 for keys that start with that char. */
		 bn = rt + (c = *b++);

		 /* A count field set to -1 indicates that the bin has previously burst
			 and is now a branched node indexing a subtrie. */
		 while (bn != NULL && bn->ct == -1)
		 {
			 /* The burst bin's base pointer now actually points to the base node
				 of the subtrie, so cast it as such */
			 nd = (NODE1*) bn->bp;

			 /* and use the next char of the key to index to the appropriate
				 bin of the subtrie, and so on ... until a terminal bin is
				 found (and the count field has a value other than -1). */
			 bn = nd + (c = *b++);
		 }

		 /* If c is a NULL, the key is exhausted; when this happens in ssample(), 
			 we simply ignore that key, since it does not affect bursting and thus
			 does not influence the trie structure formed during sampling */
		 if (c == 0) continue;

		 /* For any other char, we increment the count of keys in the bin. */
		 ct = ++bn->ct;

		 /* If the count is 1, we need to create a new bin. */
		 if (ct == 1)
		 {
			 /* Allocate a string buffer of size BINSIZE0 at the base pointer
				 (bp) of the bin, and update the tally of bin objects. */
			 s = bn->bp = (char*) MALLOC(BINSIZE0); ++BINS;

			 /* Set the bin's limit pointer to its base + LIMSIZE0, where
				 LIMSIZE0 = BINSIZE0 - MAXKEYLEN; MAXKEYLEN is guaranteed to be >= to
				 the longest possible key length. */
			 bn->lp = s + LIMSIZE0;

			 /* Copy the remainder of the key into the bin's string buffer. */
			 while ((*s++ = *b++) != 0) ;

			 /* Set the bin's active pointer to point at the next free byte. */
			 bn->ap = s;

		 } else {

			  /* Retrieve the active pointer of an existing bin. */
			  s = bn->ap;

			  /* Copy the string tail and update the active pointer. */
			  while ((*s++ = *b++) != 0) ; bn->ap = s;

			  /* If the active pointer passes the limit pointer, fewer than
				  MAXKEYLEN bytes remain, and the next key tail may not fit;
				  since bins do not grow in ssample(), the bin must be burst. */
			  if (s>bn->lp) samburst(bn);
		 }
	}
}

/*
 * Burst a bin that fills (at its original size BINSIZE0) during
 * sampling by function ssample() above.
 *
 */
void samburst(NODE1 *bn)
{
	string b0, b, s; NODE1 *rt; int c, ct, n;

	/* Create a new branched node that will hold the bin being burst. */
	rt = GETNODE(NODE1); ++NODES; 

	/* Set the input pointer to the base of the bin being burst. */
	b = b0 = bn->bp; 
	
	/* Set the burst bin's base pointer to point to the new subtrie. */
	bn->bp = (char*) rt;
	
	/* Get the count of string tails from the burst bin. */
	n = bn->ct; 
	
	/* Set the count to -1 to flag that the bin has been burst. */
	bn->ct = -1;

	/* Insert the string tails from the burst bin into the new subtrie;
		see identical code in ssample() for comments. */
	while (n--)
	{
		bn = rt +(c = *b++);
		if (c == 0) continue; ct = ++bn->ct;

		if (ct == 1)
		{
			++BINS; s = bn->bp = (char*) MALLOC(BINSIZE0); 
			bn->lp = s+LIMSIZE0;
			while ((*s++ = *b++) != 0) ; bn->ap = s;
		}
		else
		{
			s = bn->ap;
			while ((*s++ = *b++) != 0) ; bn->ap = s;
			if (s>bn->lp) samburst(bn);
		}
	 } 
	 
	 /* When the tails have been moved, free the old bin's container. */
	 UPDATEMEM; ALLOCATED -= BINSIZE0; free(b0); --BINS;
}

/*
 * After sampling, discard the sampled data, but retain the trie skeleton, 
 * which will then be filled using the full input data.
 *
 */
void clear(NODE1 *nd)
{
	int c, ct; NODE1 *bn;

	/* Starting at the root node, scan the range of char values set
		when the data was read. */
	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		
		/* Recurse to clear burst bins flagged by count = -1. */ 
		if ((ct = bn->ct) == -1) clear((NODE1*) bn->bp);
		
		/* For other non-empty bins... */
		else if (ct)
		{
			 /* If the bin did not actually burst during sampling, 
				 but it is projected to burst based on its current
				 content and the known sampling rate, burst it 
				 preemptively. */
			if (((ct<<2) + bn->ap-bn->bp) * SAMPLERATE>CACHESIZE)
			{
				samburst(bn); clear((NODE1*) bn->bp);
			}
			
			/* Otherwise, reset its active pointer and key tail count. */
			else
			{
				bn->ap = bn->bp; bn->ct = 0;
			}
		}
	}
}

/*
 * Main routine called by sbs() to insert string keys into the burst
 * trie starting at rt.	 If SAMPLERATE>0, rt will represent a skeleton
 * trie structure built during sampling; otherwise, it will be a single
 * node that expands by bursting as keys are added.  The code in sadd()
 * is similar to the insertion routines in ssample() and samburst(), but
 * covers a few more possibilities.	 Only the new parts are commented.
 *
 */
void sadd(string b, string blim, NODE1 *rt)
{
	NODE1 *bn, *nd; 
	string s, t, bp, ap, lp; 
	int c, ct, sz, sz2;

	/* sbs() ensures that b and blim coincide with key boundaries. */
	while (b < blim)
	{
		/* Traverse the trie until a terminal node (bin) is reached. */
		bn = rt + (c = *b++);
		while (bn != NULL && bn->ct == -1) {
			nd = (NODE1*) bn->bp; 
			bn = nd + (c = *b++);
		} 
		ct = ++bn->ct;

		/* We ignored exhausted keys during sampling (since they do not lead to 
			bursting or bin formation), but now we need to count each exhausted 
			key.	Such keys have no tail left to be copied to a bin, and are 
			recorded only as counts, but they can be reconstructed from the chars 
			traversed to reach a given bin.	I.e., the exhausted keys counted in 
			a bin are identical to the prefix traversed to reach it. */
		if (c == 0) {
			if (ct == 1) ++NULLBINS; continue;
		}
		
		/* Create a new bin if needed; see ssample() for comments. */
		if ((bp = bn->bp) == NULL)
		{
			++BINS;
			s = bn->bp = (string) MALLOC(BINSIZE0);
			bn->lp = s + LIMSIZE0;
			
			while ((*s++ = *b++) != 0) ; 
			bn->ap = s;
		
		/* Or add key to existing bin... */
		} else {
			s = bn->ap;
			while ((*s++ = *b++) != 0) ; 
			bn->ap = s;

			/* When an existing bin is full (next byte > limit pointer lp) it 
				must either burst or grow. */
			if (s > (lp = bn->lp))
			{
				sz = lp + MAXKEYLEN - bp;
			
				/* If FREEBURSTS > 0, burst bin regardless of size; if bin size == 
				CACHESIZE, burst bin regardless of FREEBURSTS. */ 
				if (FREEBURSTS > 0 || sz == CACHESIZE)
				{
					
					/* Set timer phase to "burst". */
					ton(tb); 
					
					/* Burst the bin into a subtrie. */
					sburst(bn);
				}

				/* if FREEBURSTS == 0, and bin size < CACHESIZE, copy the bin's 
					string tails to a larger container. */
				else
				{
					/* Set timer phase to "grow". */
					ton(tg); 
					
					/* Double the container size; since the old container is freed, 
						net allocated memory increases by the old size. */
					sz2 = sz<<1; ALLOCATED += sz;
					
					/* Allocate the new container and set its limit pointer to
						leave room for the largest possible key. */
					ap = bn->bp = (string)malloc(sz2);
					bn->lp = ap + sz2 - MAXKEYLEN;

					/* Copy string tails from the old container to the new. */
					t = bp; while (t < s) *ap++ = *t++;
					
					/* Free the old container and record the active pointer. */
					free(bp); bn->ap = ap;
				} 
				
				/* Set timer back to "insert" phase. */
				ton(ti);
			}
		}
	}
}

void sburst(NODE1 *bg) {
	NODE1 *rt, *t; string b, b0, s, p, q, r; 
	int c, ct, n, nc, nn, sz, sz0, sz2;
	
	nc=BINS; nn=NULLBINS; do {
		b=b0=bg->bp; sz0=bg->lp+MAXKEYLEN-b; 
		n=bg->ct; bg->ct=-1; 
		rt=GETNODE(NODE1); ++NODES; bg->bp=(char*)rt;
		while (n--) {
			bg=rt+(c=*b++); 
			while (bg->ct==-1) {
				t=(NODE1*)bg->bp; bg=t+(c=*b++);
			} 
			ct=++bg->ct; 
			if (c==0) {if (ct==1) ++NULLBINS; continue;}
			if ((p=bg->bp)==NULL) {
				++BINS; s=bg->bp=(string)MALLOC(BINSIZE0); 
				bg->lp=s+LIMSIZE0; 
				while ((*s++=*b++)!=0) ; bg->ap=s;
			} else {
				s=bg->ap; 
				while ((*s++=*b++)!=0) ; bg->ap=s; 
				if (s>(r=bg->lp)) {
					sz=r+MAXKEYLEN-p; sz2=sz<<1; 
					r=p; q=s; 
					s=bg->bp=(string)malloc(sz2); 
					ALLOCATED+=sz;
					bg->lp=s+sz2-MAXKEYLEN; 
					while (p<q) *s++=*p++; 
					free(r); bg->ap=s; 
				}
			}
		} UPDATEMEM; FREE(b0, sz0); --BINS;
	} while (BINS == nc && NULLBINS == nn); --FREEBURSTS;
}


/*
 * Traverse the trie and sort each container of string tails with an
 * appropriate helper sort.
 *
 */
void straverse(NODE1 *nd)
{
	int c, ct, tpsz;
	NODE1 *bn;
	string s, bp, ap, *tp, *tp0;

	/* Scan only the range of char values actually observed during key
		insertion. */
	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		
		/* If a given node is a subtrie, process it recursively. */
		if (bn != NULL && (ct = bn->ct) == -1) straverse((NODE1*) bn->bp);

		/* If the node was empty or a subtrie, skip to next char value. */
		if (ct < 1) continue;

		/* Calculate room needed for pointers to each string tail. */
		tpsz = ct<<2; 
		
		/* Get the base pointer and active pointer (next free byte). */
		s = bp = bn->bp; ap = bn->ap;

		/* If the string tails plus the sorting pointers add up to more than
			cache size, burst the bin and process it recursively, just as if
			it had burst during insertion. */
		if (ap + tpsz > bp + CACHESIZE)
		{
			 ton(tb); sburst(bn);
			 ton(tc); straverse((NODE1*) bn->bp);
			 continue;
		}

		/* Create a buffer to hold sorting pointers for the string tails. */
		tp = tp0 = (string*) MALLOC(tpsz);
		
		/* The bin's active pointer was copied to local variable ap above; now
			we redeploy field ap to hold the address of the sorting pointer
			array. */
		bn->ap = (string) tp0;
		
		/* Scan thru the string tails to build the sorting pointers. */
		while (s < ap) {*tp++ = s; while (*s++) ;}

		/* If the bin contains more than one string tail, sort it; multikey
			quicksort (defaulting to insertion sort for small numbers of keys)
			is used here, but forward radix sort is also a good choice; in 
			practice, multikey quicksort tends to be faster on data sets that
			generate tries with less branching, and forward radix on sets with 
			high branching; faster sorting can be achieved by adaptively choosing
			mkq when the ratio of terminal bins to branched nodes < 25 and radix
			when it is greater. */
		if (ct > 1) mkqsort(tp0, ct, 0);
	}
}

/*
 * De-allocate a trie or subtrie.
 *
 */
void skill(NODE1 *nd)
{
	int c, ct;
	NODE1 *bn;
	string p;

	/* Scan the range of observed char values. */
	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		
		/* Recursively kill any subtries. */
		if (bn != NULL && (ct = bn->ct) == -1) skill((NODE1*) bn->bp);
		
		else if (ct)
		{
			p = bn->bp; --BINS;
		
			/* Free buffer holding string tails. */
			FREE(p, MAXKEYLEN + (bn->lp-p));
			
			/* Free buffer holding sorted pointers. */
			FREE(bn->ap, bn->ct<<2);
		}
	}
	
	/* Decrement the null bin tally. */
	if (nd->ct) --NULLBINS;
	
	/* Free the node and decrement node tally. */
	FREE(nd, 256 * sizeof(NODE1)); --NODES;
}

/*
 * Output sorted trie to a disk file of sorted string keys.
 *
 */
void ssavetrie(NODE1 *nd, string path, string dset)
{
	FILE *f;
	string s;
	char pfx[400];

	/* Open output file. */
	f = fopen(s = fp(path, dset, "c-burstsort"), "w+");
	if (f == NULL)
	{
		say("could not open "); brp(s);
	}
	else
	{
		/* Write trie ot subtrie to file. */
		sputtrie(nd, f, pfx, 0);
		fflush(f);
		fclose(f);
	}
}

/*
 * Append sorted keys from a (sub)trie to a disk file.
 *
 */
void sputtrie(NODE1 *nd, FILE *f, string pfx, int d)
{
	int ct, c; string *tp; NODE1 *bn;

	/* pfx holds any prefix chars accumulated at char depth 0 .. d - 1. */
	pfx[d + 1] = 0; 
	
	/* For the count of exhausted keys in this node ... */
	ct = nd->ct;
	while (ct--)
	{
		/* Output only the shared prefix. */
		fputs(pfx, f);
		fputc('\n', f);
	}

	/* For the keys with tails ... */
	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		/* Build shared prefix from traversed chars. */
		bn = nd + c; pfx[d] = c;
		
		/* Recursively output subtries. */
		if (bn != NULL && (ct = bn->ct) == -1)
			 sputtrie((NODE1*) bn->bp, f, pfx, d + 1);

		/* For non-empty bins ... */
		else if (ct)
		{
			/* Use sorted pointers to access string tails in order. */
			tp = (string*) bn->ap;
			while (ct--)
			{
				/* Output shared prefix. */
				fputs(pfx, f);
				
				/* Append tail and delimiter. */
				fputs(*tp++, f);
				fputc('\n', f);
			}
		}
	}
	pfx[d] = 0;
}

/*
 * Insertion sort - among the most efficient sorting alternatives for short
 * lists; a good cutoff is n <= 12 for mkqsort or n <= 23 for forward radix.
 *
 */
void inssort(string *tp, int n, int d)
{
	string *rp, *lp, ls, rs, rd, ld;

	for (rp = tp + 1; --n > 0; ++rp)
	{
		for (rs = *rp, lp = rp; lp > tp; --lp)
		{
			 ls = *(lp - 1);
			 rd = rs + d;
			 ld = ls + d;

			 while (*ld == *rd && *rd != 0)
			 {
				 ++ld; ++rd;
			 }

			 if (*rd < *ld)
				  *lp = ls;
			 else
				  break;
		 }
		*lp = rs;
	}
}

/*
 * Based on multikey quicksort of Bentley and Sedgewick.	 Hands off short
 * sublists to insertion sort. 
 *
 */
void mkqsort(string *tp, int n, int d)
{
	 int k, nl, ne, nr, pv;
	 char **ol, **il, **ir, **or, **l, **m, **r, *t;

	 if (n < 13)
	 {
		 inssort(tp, n, d);
		 return;
	 }
	 l = tp; m = tp + (n>>1); r = tp + n - 1;

	 if (n > 31)
	 {
		 k = n>>3;
		 l = med3(l, l + k, l + k * 2, d);
		 m=med3(m - k, m, m + k, d);
		 r=med3(r - k * 2, r - k, r, d);
	 }

	 m=med3(l, m, r, d);
	 PSWAP(tp, m, t);
	 pv=P2C(tp, d);
	 r = tp + n; ol = il = tp + 1; ir = or = r - 1;

	 while (il <= ir && P2C(il, d) == pv)
	 {
		 ++ol; ++il;
	 }

	 while (il <= ir && P2C(ir, d) == pv)
	 {
		 --or; --ir;
	 }

	 for (;;)
	 {
		 while (il <= ir && (k = P2C(il, d) - pv) <= 0)
		 {
			 if (k == 0)
			 {
				 PSWAP(ol, il, t); ++ol;
			 }
			 ++il;
		 }

		 while (il <= ir && (k = P2C(ir, d) - pv) >= 0)
		 {
			 if (k == 0)
			 {
				 PSWAP(ir, or, t); --or;
			 }
			 --ir;
		 }

		 if (il > ir)
			 break;

		 PSWAP(il, ir, t);
		 ++il; --ir;
	 }

	 nl = il - ol; nr = or - ir; ne = n - (nl + nr);

	 k = MIN(ol - tp, nl); vswap(tp, il - k, k);
	 k = MIN(nr, r - or - 1); vswap(r - k, il, k);

	 if (ne > 1 && pv > 0)
		 mkqsort(tp + nl, ne, d + 1);

	 if (nl > 1)
		 mkqsort(tp, nl, d);

	 if (nr > 1)
		 mkqsort(r - nr, nr, d);
}
