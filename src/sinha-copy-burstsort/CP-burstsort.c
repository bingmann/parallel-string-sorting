/* -*- tab-width: 8 -*- */

#include "copy.h"

/*
 * Run a specified number of interations of CP-burstsort.
 * This version of copying burstsort stably sorts records
 * indexed by string keys. The 'r' prefix refers to nick-
 * name "record" burstsort. Only CP-specific code details
 * are commented here; refer to C-burstsort for comments 
 * on shared code.
 *
 */
void rbdo(int nr, string fp)
{
	ton(tr); ALLOCATED = 0; 
	
	/* Load data from source file into char array INBUF. */
	getkeys(fp);

	for (REP = 0; REP < nr; ++REP)
	{
		ton(tx);

		TAILSIZE = NODES = BINS = NULLBINS = 0;
		MAXALLOCATED = ALLOCATED;
		listinputs();

		ton(ti); 
		
		/* Run CP-burstsort. */
		rbs(INBUF);
	}
	--INBUF; FREE(INBUF, MAXBYTES + 1); ton(0);
}

/*
 * Sort buffer ib, containing records indexed by string
 * keys.
 *
 */
void rbs(string ib)
{
	NODE1 *rt;
	
	rt = GETNODE(NODE1); ++NODES;

	if (SAMPLERATE)
	{
		ton(ts); rsample(ib, rt, NKEYS / SAMPLERATE);
		clear(rt);
		ton(ti);
	}

	radd(ib, NKEYS, rt);

	UPDATEMEM; FREEBURSTS = 0; ton(tc); 
	
	rtraverse(rt); UPDATEMEM;

	if (WRITEFILE)
	{
		ton(tw); savetrie(rt, OUTPUTDIR, DATANAME);
	}
	
	ton(tk); rkill(rt);
}

/*
 * Build model burst trie using sample of n random keys;
 * except that we don't used a segmented buffer, sampling
 * is almost identical to that in C-burstsort; since the
 * sample trie will be cleared before use, we don't
 * bother to keep record pointers until the main insertion
 * in radd() below.
 * 
 */
void rsample(string b0, NODE1 *rt, int n)
{
	NODE1 *t; NODE1 *bn; 
	int ct; 
	unsigned char c; 
	string b, s, lim;

	srand(clock()); 
	
	/* Set upper limit of input buffer to the total number
		of bytes read. */
	lim = b0 + NBYTES;
	while (n--)
	{
		/* Choose a random offset into the input buffer. */
		b = b0 + rand() % (NBYTES);

		/* Back up until a null terminator is found; then
			move one byte forward to the first byte of a
			key.	The NULL placed before the input buffer in
			getbytes() ensures that this loop will terminate
			when b == b0 - 1, and allows a simpler test. */
		while (*b) --b; ++b;
		
		/* Without the sentinel value, a double test
			would be needed:
		
		while (*b && (b > b0)) --b; if (b > b0) ++b;
		
		*/

		bn = rt + (c = *b++);
		while (bn != NULL && bn->ct == -1)
		{
			t = (NODE1*) bn->bp;
			bn = t + (c = *b++);
		}

		if (c == 0) continue;
		
		ct = ++bn->ct;

		if (ct == 1)
		{
			++BINS;
			s = bn->bp = (string) MALLOC(BINSIZE0);
			bn->lp = s + LIMSIZE0;

			while ((*s++ = *b++) != 0) ;
			bn->ap = s;
		}
		else
		{
			s = bn->ap;
			while ((*s++ = *b++) != 0) ;
			bn->ap = s;
			
			/* The same samburst() function used in C-
				burstsort works fine in CP-burstsort;
				similarly, there is no need for a separate
				version of clear(). */
			if (s > bn->lp) samburst(bn);
		}
	}
}

/*
 * Main routine called by rbs() to insert string keys into the burst
 * trie starting at rt.	 The main differences between radd() and sadd()
 * involve record pointers.  While C-burstsort copies each string tail
 * in a bin immediately after the last, CP-burstsort puts a four byte
 * record pointer before each tail.	 This pointer is used to access or 
 * output records in order once their keys have been sorted.
 *
 */
void radd(string b, int n, NODE1 *rt)
{
	NODE1 *bn, *nd;
	string rp, s, t, bp, ap, lp;
	int c, ct, sz, sz2;

	while (n--)
	{
		/* Set record pointer rp to the first byte of the string key 
			(which CP-burstsort assumes to be the first byte of the record).
			If the key is not the first field of the record, a modified
			version of CP- using information on the key offset may be
			needed. */			
		rp = b;
		
		bn = rt + (c = *b++);
		while (bn != NULL && bn->ct == -1)
		{
			nd = (NODE1*) bn->bp; 
			bn = nd + (c = *b++);
		}
		ct = ++bn->ct;

		/* In C-burstsort, a NULLBIN is just a count of exhausted strings, 
			and has no container (char array) since there are no string 
			tails to put in it.	In contrast, CP-burstsort needs to store
			a record pointer for each exhausted string, so a real container
			is needed, even though there are no tails to put in it. */
		if (c == 0)
		{
			if (ct == 1)
			{
				++NULLBINS;
				s = bn->bp = (string) MALLOC(BINSIZE0);
				
				/* This NULLBIN will fill with four byte record pointers
					only, and BINSIZE0 should always be a multiple of the
					system pointer size, so we can set the limit at the
					full size of the bin's container. */
					
				bn->lp = s + BINSIZE0;
				/* Macro RP() inserts the record pointer and castes it as
					string bytes preceding the string tail (which in this
					does not exist, since we're handling exhausted keys). */
				RP(s) = rp; 
				
				/* Having inserted the record pointer, we then set the bin's
				active pointer to the following byte. */
				bn->ap = s + RPSIZE;
			}
			else
			{
				RP(bn->ap) = rp; 
				s = bn->ap += RPSIZE;
				
				if (s == bn->lp)
				{
					sz = s - bn->bp; 
					sz2 = sz<<1; ALLOCATED += sz;
					bp = bn->bp = (string)realloc(bn->bp, sz2);
					bn->lp = bp + sz2;
					bn->ap = bp + sz;
				}
			}
		}
		else
		{
			if ((bp = bn->bp) == NULL)
			{
				++BINS;
				s = bn->bp = (string) MALLOC(BINSIZE0);
				bn->lp = s + LIMSIZE0;

				RP(s) = rp; s += RPSIZE;
				while ((*s++ = *b++) != 0) ; 
				bn->ap = s;
			}
			else
			{
				s = bn->ap; 
				RP(s) = rp; s += RPSIZE;
				while ((*s++ = *b++) != 0) ; 
				bn->ap = s;

				/* The code handling full bins here is identical
					to that in C-burstsort. */
				if (s >= (lp = bn->lp))
				{
					sz = lp + MAXKEYLEN - bp;
					if (FREEBURSTS > 0 || sz == CACHESIZE)
					{
						ton(tb); rburst(bn);
					}
					else
					{
						ton(tg); 
						sz2 = sz<<1; ALLOCATED += sz;
						ap = bn->bp = (string) malloc(sz2);
						bn->lp = ap + sz2 - MAXKEYLEN;
						t = bp; while (t < s) *ap++ = *t++;
						free(bp); bn->ap = ap;
					}
					ton(ti);
				}
			}
		}
	}
}

/*
 * Convert a full terminal bin to a branched node / subtrie.
 *
 */
void rburst(NODE1 *bn)
{
	NODE1 *rt, *nd;
	string b, b0, s, t, bp, ap, lp, rp;
	int c, ct, n, nc, nn, sz0, sz, sz2;

	nc = BINS; nn = NULLBINS;
	do
	{
		b = b0 = bn->bp;
		sz0 = bn->lp + MAXKEYLEN - b; n = bn->ct;
		
		bn->ct = -1;
		rt = GETNODE(NODE1); ++NODES;
		bn->bp = (string) rt;

		while (n--)
		{
			/* During a burst, we are reading from a previous bin
				rather than the original record, so (a) we need to
				copy the existing record pointer and (b) we need
				to skip over its bytes and use only the remaining
				tail bytes in building the subtrie structure. */
			rp = RP(b); b += RPSIZE;
			
			bn = rt + (c = *b++);
			while (bn != NULL && bn->ct == -1)
			{
				nd = (NODE1*) bn->bp;
				bn = nd + (c = *b++);
			}
			ct = ++bn->ct;

			if (c == 0)
			{
				if (ct == 1)
				{
					++NULLBINS;
					s = bn->bp = (string) MALLOC(BINSIZE0);
					bn->lp = s + BINSIZE0;
					RP(s) = rp; bn->ap = s + RPSIZE;
				}
				else
				{
					s = bn->ap; 
					RP(s) = rp; s += RPSIZE; 
					bn->ap = s;
					
					if (s == bn->lp)
					{
						sz = s - bn->bp; 
						sz2 = sz<<1; ALLOCATED += sz;
						bp = bn->bp = (string)realloc(bn->bp, sz2);
						bn->lp = bp + sz2; bn->ap = bp + sz; 
					}
				}
			}
			else
			{
				if ((bp = bn->bp) == NULL)
				{
					++BINS;
					s = bn->bp = (string) MALLOC(BINSIZE0);
					bn->lp = s + LIMSIZE0;
					RP(s) = rp; s += RPSIZE;

					while ((*s++ = *b++) != 0) ;
					bn->ap = s;
				}
				else
				{
					s = bn->ap; 
					RP(s) = rp; s += RPSIZE;
					while ((*s++ = *b++) != 0) ;
					bn->ap = s;

					if (s >= (lp = bn->lp))
					{
						sz = lp + MAXKEYLEN - bp; 
						sz2 = sz<<1; ALLOCATED += sz;
						ap = bn->bp = (string) malloc(sz2);
						bn->lp = ap + sz2 - MAXKEYLEN;
						t = bp; while (t < s) *ap++ = *t++;
						free(bp); bn->ap = ap;
					}
				}
			}
		}
	
		UPDATEMEM; FREE(b0, sz0); --BINS;
	
	} while (BINS == nc && NULLBINS == nn); 
	
	--FREEBURSTS;
}

/*
 * Traverse the trie and sort each container of string tails with an
 * appropriate helper sort.
 *
 */
void rtraverse(NODE1 *nd)
{
	int c, ct, tpsz; 
	NODE1 *bn;
	string s, bp, ap, *tp, *tp0;

	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		if (bn != NULL && (ct = bn->ct) == -1) rtraverse((NODE1*) bn->bp);
		if (ct < 1) continue;

		tpsz = ct * sizeof(string*);
		s = bp = bn->bp; ap = bn->ap;

		if (ap + tpsz > bp + CACHESIZE)
		{
			ton(tb); rburst(bn);
			ton(tc); rtraverse((NODE1*) bn->bp);
			continue;
		}

		tp = tp0 = (string*) MALLOC(tpsz);
		bn->ap = (string) tp0;
		
		/* We set the sorting pointers to point to the tails 
			rather than the record pointers that precede them. */
		while (s < ap)
		{
			/* Skip over the record pointer. */
			s += RPSIZE; 
			/* Set the tail pointer. */
			*tp++ = s;
			/* Scan thru the rest of the tail. */
			while (*s++) ;
		}

		/* mkqsorts() (q.v.) is a stable version of mkqsort(). */
		if (ct > 1) mkqsorts(tp0, ct, 0);
		
		/* After sorting, we replace each tail pointer with the
			record pointer of the corresponding tail; this puts
			the records pointers in order in a minimal contiguous
			block, and we can then release the char array that held
			the interspersed record pointers and tails. */
		for (tp = tp0; ct-- > 0; ++tp) *tp = RP(*tp - RPSIZE);

		UPDATEMEM; ALLOCATED -= (MAXKEYLEN + bn->lp - bp);
		free(bp); --BINS;
	}
}

/*
 * De-allocate trie or subtrie.
 *
 */
void rkill(NODE1 *nd)
{
	int c, ct;
	NODE1 *bn;

	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		if (bn != NULL && (ct = bn->ct) == -1) rkill((NODE1*) bn->bp);
		else if (ct)
		{
			/* At this point in CP-burstsort, we've already freed
				the bin's original container at bp, but we need to
				free the sorted record pointers at ap. */
                    FREE(bn->ap, bn->ct * sizeof(string*));
		}
	}
	if (nd->ct)
	{
		--NULLBINS;
		FREE(nd->bp, nd->lp-nd->bp);
	}
	FREE(nd, 256*sizeof(NODE1)); --NODES;
}

/*
 * Output a sorted file of records.
 *
 */
void savetrie(NODE1 *nd, const char* path, const char* dset)
{
	FILE *f;
	char* s;

	if ((f = fopen(s = fp(path, dset, "cp-burstsort"), "w+")) == NULL)
	{
		say("could not open "); brp(s);
	}
	else
	{
		puttrie(nd, f);
		fflush(f);
		fclose(f);
	}
}

/*
 * Append sorted records to a file.
 *
 */
void puttrie(NODE1 *nd, FILE *f)
{
	int ct, c; string *tp; NODE1 *bn;

	ct = nd->ct;
	
	/* While sputtrie() had to reconstruct a prefix from
		traversed chars, puttrie() can simply use a record
		pointer to retrieve the original key.	This version
		just outputs the key; to output the whole record, 
		puttrie() would need to know about the fields of 
		the record and how they are delimited. */
	tp = (string*) nd->bp;
	
	/* Output any exhausted keys. */
	while (ct--)
	{
                fputs((char*)*tp++, f);
		fputc('\n', f);
	}

	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		if (bn != NULL && (ct = bn->ct) == -1)
			puttrie((NODE1*) bn->bp, f);
		else if (ct)
		{
			tp = (string*) bn->ap;
			while (ct--)
			{
				fputs((char*)*tp++, f);
				fputc('\n', f);
			}
		}
	}
}

/** Added by Timo Bingmann to output the trie to a string pointer array */

string* puttrie_strout;        /* string pointer output iterator */

/*
 * Append sorted records to a file.
 *
 */
void tb_puttrie(NODE1 *nd)
{
	int ct, c; string *tp; NODE1 *bn;

	ct = nd->ct;
	
	/* While sputtrie() had to reconstruct a prefix from
		traversed chars, puttrie() can simply use a record
		pointer to retrieve the original key.	This version
		just outputs the key; to output the whole record, 
		puttrie() would need to know about the fields of 
		the record and how they are delimited. */
	tp = (string*) nd->bp;
	
	/* Output any exhausted keys. */
	while (ct--)
	{
                *puttrie_strout++ = *tp++;
	}

	for (c = LOCHAR; c <= HICHAR; ++c)
	{
		bn = nd + c;
		if (bn != NULL && (ct = bn->ct) == -1)
			tb_puttrie((NODE1*) bn->bp);
		else if (ct)
		{
			tp = (string*) bn->ap;
			while (ct--)
			{
                                *puttrie_strout++ = *tp++;
			}
		}
	}
}

/*
 * Output a sorted file of records.
 */
void tb_savetrie(NODE1 *nd, string* strout)
{
        puttrie_strout = strout;
        tb_puttrie(nd);
}

/*
 * A version of insertion sort that can restore stability to
 * data after partial sorting with an unstable helper sort;
 * 
 *
 */
void inssorts(string *tp, int n, int d)
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

				/* The re-stabilization occurs here; if tails are 
					equal, we just ensure that final order matches 
					the order of the tails in the bin, which reflects 
					their original order of insertion.*/
			  if (*rd < *ld || (*rd == 0 && rs < ls)) 
				  *lp = ls;
			  else
				  break;
		}
		*lp = rs;
	}
}

/*
 * Stabilize sorting by otherwise unstable helper sorts such
 * as multikey quicksort or displacing forward radix sort; this
 * routine is used with CP- and CPL-burstsort and works by
 * ordering a run of equal key tails according to the position
 * of those tails within the bin (which preserves the original
 * order in which they were inserted.
 *
 */
void stabilize(string *tp, int n)
{
	int n1, n2, n3;
	string *l, *m, *r, pv, t;

	if (n < 13)
	{
		 for (r = tp + 1; --n > 0; ++r)
		 {
			 pv = *r;
			 for (l = r; l > tp; --l)
			 {
				 t = *(l - 1);
				 if (pv < t) *l = t;
				 else
					 break;
			 }
			 *l = pv;
		 }
		 return;
	}

	n1 = n>>1; n2 = n1>>1; n3 = n2>>1;
	l = tp; m = tp + n1; r = tp + n - 1;

	if (n > 31)
	{
		l = med3s(l, l + n3, l + n2);
		m = med3s(m - n3, m, m + n3);
		r = med3s(r - n2, r - n3, r);
	}

	m = med3s(l, m, r);
	PSWAP(tp, m, t);
	pv = *tp;
	l = tp + 1; r = tp + n - 1;

	for (;;)
	{
		while (l <= r && *l < pv) ++l;
		while (l <= r && *r > pv) --r;
		if (l > r)
			break;
		PSWAP(l, r, t);
		++l; --r;
	}

	PSWAP(tp, r, t);
	--r; n1 = l - tp; n2 = n - n1;
	if (n1 > 1)
		stabilize(tp, n1);
	if (n2 > 1)
		stabilize(tp + n1, n2);
}

/*
 * This version of mkqsort restores stability by calling 
 * stabilize() to handle runs of equal keys and by passing
 * short runs to inssorts();
 *
 */
void mkqsorts(string *tp, int n, int d)
{
	int k, nl, ne, nr, pv;
	string *ol, *il, *ir, *_or, *l, *m, *r, t;

	if (n < 13)
	{
		inssorts(tp, n, d);
		return;
	}
	l = tp; m = tp + (n>>1); r = tp + n - 1;

	if (n > 31)
	{
		k = n>>3;
		l = med3(l, l + k, l + k * 2, d);
		m = med3(m - k, m, m + k, d);
		r = med3(r - k * 2, r - k, r, d);
	}

	m = med3(l, m, r, d);
	PSWAP(tp, m, t);
	pv=P2C(tp, d);
	r = tp + n; ol = il =tp + 1; ir = _or = r - 1;

	while (il <= ir && P2C(il, d) == pv)
	{
		++ol; ++il;
	}
	while (il <= ir && P2C(ir, d) == pv)
	{
		--_or; --ir;
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
				PSWAP(ir, _or, t); --_or;
			}
			--ir;
		}

		if (il > ir)
			break;
		PSWAP(il, ir, t);
		++il; --ir;
	}
	nl = il - ol; nr = _or - ir; ne = n - (nl + nr);

	k = MIN(ol - tp, nl);
	vswap(tp, il - k, k);
	k = MIN(nr, r - _or - 1);
	vswap(r - k, il, k);

	if (ne > 1)
	{
		if (pv == 0)
			stabilize(tp + nl, ne);
		else
			 mkqsorts(tp + nl, ne, d + 1);
	}

	if (nl > 1)
		mkqsorts(tp, nl, d);
	if (nr > 1)
		mkqsorts(r - nr, nr, d);
}
