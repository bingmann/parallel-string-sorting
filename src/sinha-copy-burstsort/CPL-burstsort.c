#include "copy.h"

/*
 * Run a specified number of interations of CPL-burstsort.
 * Like CP-, CPL-burstsort stably sorts records by their
 * string keys, but CPL- copies only a fixed length 
 * segment of the key tails to its bins; in cases where
 * more chars are needed for sorting, additional fixed
 * length segments are "paged" in from the original 
 * records (hence the 'p' prefix for "paging" burstsort).
 * Only code specific to CPL-burstsort is commented here;
 * shared code is commented under C-burstsort and/or
 * CP-burstsort, q.v.
 *
 */
void pbdo(int nr, string fp)
{
	int lv, sz;
	double d;

	ton(tr); ALLOCATED = 0;
	
	getkeys(fp);
	
	for (REP = 0; REP < nr; ++REP)
	{
		ton(tx);
		
		/* The command line arg TAILRATE specifies a
			percentage of the average key length to use as
			the key segment size; we calculate this now. */
		d = TAILRATE; d /= 100; 
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
		lv = 0; sz = BINSIZE0;
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
		
		NODES = BINS = NULLBINS = 0;
		MAXALLOCATED = ALLOCATED;
		listinputs();

		/* Run CPL-burstsort. */
		ton(ti); pbs(INBUF);

	}
	--INBUF; FREE(INBUF, MAXBYTES + 1);
	ton(0);
}

/*
 * Sort buffer ib, containing records indexed by string
 * keys.
 *
 */
void pbs(string ib)
{
	NODE2 *rt;

	rt = GETNODE(NODE2); ++NODES;

	if (SAMPLERATE)
	{
		ton(ts);
		psample(ib, rt, NKEYS / SAMPLERATE);
		pclear(rt);
		ton(ti);
	}
	
	padd(ib, NKEYS, rt);

	UPDATEMEM; FREEBURSTS = 0; ton(tc); 
	
	ptraverse(rt); UPDATEMEM;

	if (WRITEFILE)
	{
		ton(tw); psavetrie(rt, OUTPUTDIR, DATANAME);
	}
	
	ton(tk); pkill(rt);
}

/*
 * Build model burst trie using sample of n random keys;
 * except for types and the samburst() routine called, 
 * the code is identical to rsample().
 *
 */
void psample(string b0, NODE2 *rt, int n)
{

	NODE2 *t; NODE2 *bn;
	int ct;
	char c;
	string b, s, lim;

	srand(clock());
	lim = b0 + NBYTES;
	while (n--)
	{
		b = b0 + rand() % NBYTES;
		while (*b) --b; ++b;

		bn = rt + (c = *b++);
		while (bn != NULL && bn->lv == -1)
		{
			t = (NODE2*) bn->bp;
			bn = t + (c = *b++);
		}

		if (!c) continue;
		
		ct = ++bn->ct;

		if (ct == 1)
		{
			++BINS;
			s = bn->bp = (char*) MALLOC(BINSIZE0);
			while ((*s++ = *b++) != 0) ;
			bn->ap = s;
		}
		else
		{
			s=bn->ap;
			while ((*s++ = *b++) !=0 ) ;
			bn->ap = s;
			if (s + MAXKEYLEN > bn->bp + BINSIZE0)
			psamburst(bn);
		}
	}
}

/*
 * Burst a bin that fills during sampling by function
 * psample().	The only reason we can't use samburst()
 * is that CPL-burstsort uses NODE2 objects rather than
 * NODE1; otherwise, the code is almost identical.
 *
 */
void psamburst(NODE2 *bn)
{
	string b0, b, s; NODE2 *rt; int c, ct, n;

	rt = GETNODE(NODE2); ++NODES;
	
	b = b0 = bn->bp;
	
	bn->bp = (char*) rt;

	n = bn->ct;

	/* NODE2 objects flag burst nodes by setting their 
		"lv" field to -1 rather than their "ct" field. */
	bn->lv = -1;

	while (n--)
	{
		bn = rt + (c = *b++);
		if (c == 0) continue; ct = ++bn->ct;

		if (ct == 1)
		{
			++BINS; s = bn->bp = (char*) MALLOC (BINSIZE0);
			while ((*s++ = *b++) != 0) ; bn->ap = s;
		}
		else
		{
			s = bn->ap;
			while ((*s++ = *b++) != 0) ;
				bn->ap = s;
				if (s + MAXKEYLEN > bn->bp + BINSIZE0)
					psamburst(bn);
		}
	}
	
	UPDATEMEM; ALLOCATED -= BINSIZE0; free(b0); --BINS;
}

/*
 * Retain the trie structure built by sampling, but clear
 * the data in it; a separate version is needed in CPL-
 * burstsort because of the difference in NODE objects.
 *
 */
void pclear(NODE2 *nd)
{
	int c;
	NODE2 *bn;

	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		
		if (bn != NULL && bn->lv == -1) pclear((NODE2*) bn->bp);

		else if (bn->ct)
		{
			bn->ap = bn->bp; bn->ct = 0;
		}
	}
}

/*
 * Main routine called by pbs() to insert string keys into the burst
 * trie starting at rt.  Like CP-burstsort, CPL-burstsort stores a
 * record pointer before each string tail in a bin.  Unlike C- and CP-,
 * CPL- copies only a fixed number of tail bytes from each key to its
 * bin.  If more bytes are needed to determine a key's sorted order,
 * they are fetched via its record pointer (at the probable cost of an
 * out of cache memory access), but if not, CPL- comes out ahead,
 * particularly if many keys contain significant terminal sequences
 * that do not affect their sorting.
 *
 */
void padd(string b, int n, NODE2 *rt) {
	NODE2 *bn, *nd; string rp, s, bp, ap, t; 
	int c, ct, lv, sz;
	
	while (n--) 
	{
		rp = b;
		bn = rt + (c = *b++); 
		
		/* In CPL-burstsort, the 'lv' field instead of the
			'ct' field is used to flag burst nodes. */
		while (bn != NULL && bn->lv == -1) {
			nd = (NODE2*) bn->bp; 
			bn = nd + (c = *b++);
		} 
		ct=++bn->ct;
		
		/* Exhausted keys */
		if (c == 0)
		{
			if (ct == 1)
			{
				++NULLBINS; 
				s = bn->bp = (string) MALLOC(BINSIZE0); 
				RP(s) = rp; bn->ap = s + 4;
			} 
			else 
			{
				RP(bn->ap) = rp; 
				bn->ap += 4; 
				sz = BINSIZE0<<bn->lv;

				/* Since CPL-burstsort copies only fixed length 
					segments of string tail data, it does not need 
					the limit pointer approach used to detect full 
					bins in C- and CP-burstsort. Instead, we just 
					compare bin size with key count and the fixed 
					size per key (size of record pointer + TAILSIZE 
					bytes for a key or 0 for an exhausted key). */
				if (ct<<2 == sz)
				{
					bp = bn->bp = realloc(bn->bp, sz<<1); 
					ALLOCATED += sz; ++bn->lv; bn->ap = bp + sz;
				}
			}
		} 
		else
		
		/* Keys with tails. */
		{
			if ((bp = bn->bp)==NULL) {
				++BINS; 
				s = bn->bp = (string) MALLOC(BINSIZE0); 

			  /* Instead of copying up to a null terminator, we
			  	  increment the active pointer by TAILSIZE + 4
			  	  and copy up to it. */
				ap = bn->ap = s + TAILSIZE4; 
				RP(s) = b; s += 4; 
				do {*s++ = c = *b++;} while (c > 0 && s < ap);
			} 
			else 
			{
				s = bn->ap; 
				ap = bn->ap += TAILSIZE4; 
				RP(s) = b; s += 4; 
				do {*s++ = c = *b++;} while (c > 0 && s < ap);
				
				if (ct == THR[lv = bn->lv])
				{
					sz = BINSIZE0<<lv; 
					if (FREEBURSTS > 0 || sz == CACHESIZE) 
					{
						ton(tb); pburst(bn);
					} else {
						ton(tg); 
						s = bn->bp = (string) malloc(sz<<1); 
						ALLOCATED += sz; ++bn->lv; 
						t = bp; while (t < ap) *s++ = *t++; 
						free(bp); bn->ap=s; 						
					} 
					ton(ti);
				}
			}
		  
		  /* If the segment copied did not finish the tail, 
		  		we need to scan thru the remaining bytes to 
		     the beginning of the next key. */
			if (c) while (*b++) ;
		}
	}
}

/*
 * Convert a full terminal bin to a branched node / subtrie;
 *
 */
void pburst(NODE2 *bn) {
	NODE2 *rt; string b, b0, b1, rp, s, bp, ap, t; 
	int c, ct, lv, n, nc, nn, sz, sz0;
	
	nc=BINS; nn=NULLBINS; 
	do 
	{
		b1 = b0 = bn->bp; 
		sz0 = BINSIZE0<<bn->lv; n = bn->ct; 

		/* CPL-burstsort flags 'lv' rather than 'ct'. */
		bn->lv = -1; 
		
		rt = GETNODE(NODE2); ++NODES; 
		bn->bp = (string) rt;
		
		while (n-->0) 
		{
			/* Copy the record pointer from the source bin. */
			rp = RP(b1);
			
			/* Skip over the record pointer to get the first tail byte. */
			c = *(b1 + 4);
			
			/* Leave b1 pointing to the next record pointer. */
			b1 += TAILSIZE4;

			if (c == 0)
			{
				/* Since we won't be retrieving any more segments
					from an exhausted key, allof() resets the record 
					pointer to the first byte of the key. */
				rp = allof(rp);
				ct = ++rt->ct;

				if (ct == 1)
				{
					++NULLBINS; 
					s = rt->bp = (string) MALLOC(BINSIZE0); 
					RP(s) = rp; 
					rt->ap = s + 4;
				}
				else
				{
					RP(rt->ap) = rp; 
					rt->ap += 4; 
					sz = BINSIZE0<<rt->lv;
					
					/* An exhausted key needs four bytes for its
						record pointer; since bin size is always a 
						multiple of 4, the bin is full when ct * 4
						equals bin size. */
					if (ct<<2 == sz) 
					{
						rt->ap = rt->bp = realloc(rt->bp, sz<<1); 
						ALLOCATED += sz; ++rt->lv; rt->ap += sz;
					}
				}
			}
			else
			{
				/* This simple version of cburst() goes back to 
					the original key for input bytes, most likely 
					incurring an out-of-cache reference; even so, 
					CPL-burstsort beats CP-burstsort on datasets 
					where keys are significantly longer than the 
					byte depth needed to sort them; a faster (but 
					much trickier) approach is to get input bytes 
					from the tail segment in cache until it runs 
					out, and only then to fetch the next segment. */
				b = ++rp;

				bn = rt + c; ct = ++bn->ct;

				if ((bp = bn->bp) == NULL) 
				{
					++BINS; 
					s = bn->bp = (string) MALLOC(BINSIZE0); 
					ap = bn->ap = s + TAILSIZE4; 
					RP(s) = b; s += 4; 
					do {*s++ = c = *b++;} while (c > 0 && s < ap);
				}
				else
				{
					s = bn->ap; 
					ap = bn->ap += TAILSIZE4; 
					RP(s) = b; s += 4; 
					do {*s++ = c = *b++;} while (c > 0 && s < ap);
					
					if (ct == THR[lv = bn->lv]) 
					{
						sz = BINSIZE0<<lv; 
						s = bn->bp = (string) malloc(sz<<1); 
						ALLOCATED += sz; ++bn->lv; 
						t = bp; while (t < ap) *s++ = *t++; 
						free(bp); bn->ap = s; 
					}
				}
			}
		}
		UPDATEMEM; FREE(b0, sz0); --BINS;
	} while (BINS == nc && NULLBINS == nn); --FREEBURSTS;
}

/*
 * Traverse the trie and sort each container with a helper sort
 * designed for the segmented data access approach used in CP-
 * burstsort.
 *
 */
void ptraverse(NODE2 *nd)
{
	int c, ct, tpsz;
	NODE2 *bn;
	string bp, ap, s, *tp, *tp0;

	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		
		if (bn != NULL && bn->lv == -1)
		{
			ptraverse((NODE2*) bn->bp);
			continue;
		}
		
		if ((bp = bn->bp) == NULL) continue;

		ct = bn->ct;
		tpsz = ct<<2;
		s = bp;
		ap = bn->ap;

		tp = tp0 = (string*) MALLOC(tpsz);
		bn->ap = (string) tp0; s += 4;
		
		while (s < ap) {*tp++ = s; s += TAILSIZE4;}

		if (ct > 1) mkqsortp(tp0, ct, 0);
		
		tp=tp0;
		while (ct--)
		{
			s = *tp - 4;
			*tp++ = allof(RP(s));
		}

		UPDATEMEM; ALLOCATED -= BINSIZE0<<bn->lv;
		free(bp); --BINS;
	}
}

/*
 * De-allocate trie or subtrie; same as rkill() except
 * for differences tied to NODE2 object.
 *
 */
void pkill(NODE2 *nd)
{
	int c, ct;
	NODE2 *bn;

	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		if (bn != NULL && bn->lv == -1)
			pkill((NODE2*) bn->bp);
		else if ((ct = bn->ct) > 0)
		{
			FREE(bn->ap, ct<<2);
		}
	}
	if (nd->ct)
	{
		--NULLBINS;
		FREE(nd->bp, BINSIZE0<<nd->lv);
	}
	FREE(nd, 256*sizeof(NODE2));
	--NODES;
}

/*
 * Output a sorted file of records. 
 *
 */
void psavetrie(NODE2 *nd, string path, string dset)
{
	FILE *f;
	string s;

	if ((f = fopen(s = fp(path, dset, "cpl-burstsort"), "w+")) == NULL)
	{
		say("could not open "); brp(s);
	}
	else
	{
		pputtrie(nd, f);
		fflush(f);
		fclose(f);
	}
}

/*
 * Append sorted records to a file.
 *
 */
void pputtrie(NODE2 *nd, FILE *f)
{
	int ct, c; string *tp; NODE2 *bn;

	ct = nd->ct;
	tp = (string*) nd->bp;
	
	while (ct--)
	{
		fputs(*tp++, f);
		fputc('\n', f);
	}

	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		if (bn != NULL && bn->lv == -1)
			pputtrie((NODE2*) bn->bp, f);
		else if ((ct = bn->ct) > 0)
		{
			tp = (string*) bn->ap;
			while (ct--)
			{
				fputs(*tp++, f);
				fputc('\n', f);
			}
		}
	}
}

/*
 * A version of stabilizing inssorts() adapted to the segmented
 * data access of CPL-burstsort.
 *
 */
void inssortp(string *tp, int n, int d)
{
	string *rp, *lp, ls, rs, rd, ld;
	int o;

	for (rp = tp + 1; --n > 0; ++rp)
	{
		for (rs = *rp, lp = rp; lp > tp; --lp)
		{

			ls = *(lp - 1);
			rd = rs + d;
			ld = ls + d;
			o = d;

			while (*ld == *rd && *rd != 0)
			{
				if (++o == TAILSIZE)
				{
					ld = *(string*) (ls - 4) + TAILSIZE;
					rd = *(string*) (rs - 4) + TAILSIZE;
				}
				else
				{
					++ld; ++rd;
				}
			}
			if (*rd < *ld || (*rd == 0 && rs < ls))
				*lp = ls;
			else
				break;
		}
		*lp = rs;
	}
}

/*
 * A version of stabilizing mkqsorts() adapted to the segmented
 * data access of CPL-burstsort.
 *
 */
void mkqsortp(string *tp, int n, int d)
{
	int k, nl, ne, nr, pv;
	char **ol, **il, **ir, **or, **l, **m, **r, *t, *s, *rp, **q;

	if (d == TAILSIZE)
	{
		d = 0;
		for (k = 0; k < n; ++k)
		{
			s = tp[k];
			q = (string*) (s - 4);
			rp = *q += TAILSIZE;
			t = s + TAILSIZE;
			do
			{
				*s++ = *rp++;
			} while (s < t);}
	}

	if (n < 13)
	{
		inssortp(tp, n, d);
		return;
	}
	l = tp;
	m = tp + (n>>1);
	r = tp + n - 1;

	if (n > 31)
	{
		k = n>>3;
		l = med3(l, l + k, l + k * 2, d);
		m = med3(m - k, m, m + k, d);
		r = med3(r - k * 2, r - k, r, d);
	}

	m = med3(l, m, r, d);
	PSWAP(tp, m, t);
	pv = P2C(tp, d);
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
	nl = il - ol;
	nr = or - ir;
	ne = n - (nl + nr);

	k = MIN(ol - tp, nl); vswap(tp, il - k, k);
	k = MIN(nr, r - or - 1); vswap(r - k, il, k);

	if (ne > 1)
	{
		if (pv == 0)
			stabilize(tp + nl, ne);
		else
			mkqsortp(tp + nl, ne, d + 1);
	}

	if (nl > 1)
		mkqsortp(tp, nl, d);
	if (nr > 1)
		mkqsortp(r - nr, nr, d);
}
