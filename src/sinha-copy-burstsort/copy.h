/* -*- tab-width: 8 -*- */
#ifndef COPY_H
#define COPY_H
	#include <time.h>
	#include <stdio.h>
	#include <string.h>
	#include <stdlib.h>
	#include <math.h>
	#include <float.h>
	#include <limits.h>

/* miscellaneous utility macros */

	/* returns minimum of two numbers */
	#define MIN(X, Y) ((X) <= (Y) ? (X) : (Y))

	/* returns maximum of two numbers */
	#define MAX(X, Y) ((X) >= (Y) ? (X) : (Y))

	/* returns ratio of two numbers as a double */
	#define DR(X, Y) ((double) (X) / (double) (Y))

	/* returns product of two numbers as a double */
	#define DP(X, Y) ((double) (X) * (double) (Y))

	/* returns (X * Y) / Z as a double */
	#define DPR(X, Y, Z) (((double) (X) * (double) (Y)) / (double) (Z))

	/* swaps the contents of two pointers */
	#define PSWAP(X, Y, Z) {Z = *(X); *(X) = *(Y); *(Y) = Z;}

	/* retrieve char at depth D in string tail pointed to by TP */
	#define P2C(TP, D) (*(*(TP) + D))

	/* copy N chars from S to T; (N, S, T are all changed) */
	#define SC(T, S, N) while (N-- > 0) *T++ = *S++

/* burstsort macros */

	/* returns an empty trie node of the specified type */
	#define GETNODE(typ) (typ*) CALLOC(256, sizeof(typ))

	/* in CP-burstsort, returns a pointer to the record matching a key tail */
	#define RP(x) *((string*) (x))

        /*-tb the size of the record pointer */
        #define RPSIZE  sizeof(string*)

/* memory allocation macros */ /* tb: stripped of counting */

	/* compare current/maximum allocated memory; update maximum if needed */
        /*#define UPDATEMEM if (ALLOCATED > MAXALLOCATED) MAXALLOCATED = ALLOCATED*/
        #define UPDATEMEM

	/* return x bytes from malloc(), update running total of allocated memory */
	/*#define MALLOC(x) malloc(x); ALLOCATED += (x)*/
        #define MALLOC(x) malloc(x);

	/* return mem from calloc() and update running total */
        /*#define CALLOC(x, y) calloc(x, y); ALLOCATED += (x * y)*/
        #define CALLOC(x, y) calloc(x, y);

	/* free y bytes of memory at pointer x and adjust running total */
        /*#define FREE(x, y) ALLOCATED -= y; free(x)*/
        #define FREE(x, y) free(x)

/* timer macros */

	/* maximum number of repetitions of a sorting run */
	#define MAXREPS 12

	/* maximum number of timer phases recording different program activities */
	#define MAXTMRS 12

	/* factor converting clock ticks to milliseconds */
	#define MSEC_PER_CLOCK ((double) 1000.0 / (double) CLOCKS_PER_SEC)

	/* timer phases */
	#define tr 1 /* read; input data from disk file */
	#define tx 2 /* prep; miscellaneous stuff not fitting elsewhere */
	#define ts 3 /* sample; analyse sample of data to set sort parameters */
	#define ti 4 /* insert; add keys to burst trie */
	#define tb 5 /* burst; expand a terminal node into a subtrie */
	#define tg 6 /* grow; copy terminal node string tails to a larger array */
	#define tc 7 /* traverse; step thru the burst trie sorting each bin in turn */
	#define th 8 /* check; test for correct ordering and stability */
	#define tk 9 /* kill; deallocate trie objects and buffers  */
	#define tw 10 /* write; output data to disk file */

#define DUMP {\
	isaye(TAILSIZE,"TAILSIZE"); cr();\
	isaye(NKEYS,"NKEYS"); isaye(NBYTES,"NBYTES"); cr();\
	isaye(LOCHAR,"LOCHAR"); isaye(HICHAR,"HICHAR"); cr();\
	isaye(NCHARS,"NCHARS"); isaye(MAXKEYLEN,"MAXKEYLEN"); cr();\
	isaye(NSEGS,"NSEGS"); isaye(SEGSIZE,"SEGSIZE"); cr();\
	isaye(BINSIZE0,"BINSIZE0"); isaye(LIMSIZE0,"LIMSIZE0"); cr();\
	isaye(FREEBURSTS,"FREEBURSTS"); isaye(SAMPLERATE,"SAMPLERATE"); cr();\
	isaye(INPUTORDER,"INPUTORDER"); isaye(WRITEFILE,"WRITEFILE"); cr();\
	isaye(CACHESIZE,"CACHESIZE"); isaye(NODES,"NODES"); cr();\
	isaye(BINS,"BINS"); isaye(NULLBINS,"NULLBINS"); cr();\
	isaye(ALLOCATED,"ALLOCATED"); isaye(MAXALLOCATED,"MAXALLOCATED"); br();\
}

	typedef unsigned char *string;

	/* object used by helper sorts to stack sort runs */
	typedef struct {
		string *a;	/* pointer to strings being sorted */
		int n;		/* count of strings in run */
		int d;		/* depth at which to compare chars */
	} run;

	/* object used to time program phases */
	typedef struct {int n; double t[MAXREPS];} timer;

	/* node object used in burstsort A */
	typedef struct trierec {char **ntp; char lv[256]; int ct[256]; char **p[256];} trie0;

	/* node object used in C- and CP- burstsort */
	typedef struct {
		int ct;	  /* count of string keys in a terminal node */
		string bp; /* pointer to base of array holding string tails */
		string ap; /* pointer to active position where next tail will be added */
		string lp; /* limit pointer signaling end of array */
	} NODE1;

	/* node object used in CPL-burstsort */
	typedef struct {
		int lv, 		/* level used to find size of bin */
		ct;			/* count of string keys in bin */
		string bp, 	/* base pointer */
		ap;			/* active position pointer */
	} NODE2;
	
	typedef void (*sort1) (string *a, int n, int d);
	typedef void (*test) (int nr, string fp);

	/* globals */
	size_t	
                MAXKEYS, 		/* maximum keys to read and sort (arg) */
	        MAXBYTES, 		/* maximum bytes to read and sort (arg) */
		NKEYS, 			/* number of keys sorted */
		NBYTES, 			/* number of bytes sorted */
                CHARCOUNT[256], 	/* array used to determine char usage */
		LOCHAR, 			/* lowest char value found in data */
		HICHAR, 			/* highest char value found in data */
                NCHARS, 			/* equals 1 + HICHAR - LOCHAR */
		
                NSEGS0, 			/* default segment count of input buffer (arg) */
		NSEGS, 			/* adjusted segment count of input buffer */
                SEGSIZE0, 		/* default size of input buffer segments */
		SEGSIZE, 		/* adjusted size of input buffer segments */
		CACHESIZE, 		/* size of cache memory in bytes (arg) */
		BINSIZE0, 		/* initial container size (arg) */
		MAXKEYLEN, 		/* should be >= longest expected key */
		LIMSIZE0, 		/* equals BINSIZE0 - MAXKEYLEN */
                FREEBURSTS0, 	/* number of free bursts to allow (arg) */
	 	FREEBURSTS, 	/* number of free bursts available */
                SAMPLERATE0, 	/* sampling rate (arg) */
		SAMPLERATE,		/* sampling rate in use */
		TAILSIZE0, 		/* default tail size for CPL-burstsort (arg) */
		TAILRATE, 		/* used in setting TAILSIZE (arg) */
		TAILSIZE, 		/* length of tail segments used in CPL-burstsort */
		TAILSIZE4, 		/* TAILSIZE+RPSIZE */
		THR[20], 		/* thresholds for CPL-burstsort */
		
                INPUTORDER, 	/* load sorted/random/reversed input data (arg) */
		WRITEFILE, 		/* toggles output of sorted strings (arg) */		
		TIMERPHASE, 	/* index of currently active timer */
		REP, 				/* # of current rep */

		NODES, 			/* number of branched nodes created */
		BINS, 			/* number of terminal nodes created */
		NULLBINS, 		/* number of exhausted string counts created */
		OUTOFORDER, 	/* number of keys out of sorted order */
		UNSTABLE, 		/* number of equal keys not in original order */
		UNIQUE, 			/* number of unique keys */
		SORTEDBYTES, 	/* number of bytes checked after sorting */
		EXHAUSTEDKEYS, /* number of exhausted keys counted */
		ALLOCATED, 		/* running total of allocated memory bytes */
		MAXALLOCATED;	/* peak allocated memory bytes during a sort */
	double
		TMIN, 			/* lowest time required by a repeat of a sort */
		TNORM;			/* equals TMIN/(NKEYS * log10(NKEYS)) */
        char
		*INPUTDIR, 		/* directory containing input files */
                *DATANAME, 		/* name of file to be sorted */
		*FILETYPE, 		/* file extension */
		*OUTPUTDIR, 		/* directory in which to create output files */
		*SORTNAME, 		/* name of sort algorithm in use */
                *HELPERNAME;	        /* name of helper sort in use */

        string
		INBUF, 			/* input buffer used in CP- and CPL-burstsort */
		*SEGMENTS, 		/* segmented input buffer used in C-burstsort */
		*SEGLIMITS;		/* limits for the segments of SEGMENTS	 */
	FILE
		*ERRORLOG;		/* error log */
	test
		TEST;				/* set by setsort() to be run in doreps() */
	timer
		*TMR, 			/* timer for current timing phase */
		TR[MAXTMRS];	/* array of timers for each timing phase */

	/***********************************************************************************************/

	int main(int argc, char *argv[]);
	void showhelp();
	void setsort(int sn);
	void listinputs();
	void doreps(int nr);
	void getsegs(string fp);
	void getkeys(string fp);

	/* C-burstsort */
	void sbdo(int nr, string fp);
        void sbs(string *seg, string *seglim);
	void ssample(string *seg, string *seglim, NODE1 *rt, int n);
	void samburst(NODE1 *bn);
	void clear(NODE1 *t);
	void sadd(string b, string bn, NODE1 *rt);
	void sburst(NODE1 *bn);
	void straverse(NODE1 *t);
	void scheck(NODE1 *t, int d);
	void skill(NODE1 *t);
	void ssavetrie(NODE1 *t, const char* path, const char* dset);
	void sputtrie(NODE1 *t, FILE *f, string pfx, int d);
	void inssort(string *a, int n, int d);
	void mkqsort(string *a, int n, int d);

	/* CP-burstsort */
	void rbdo(int nr, string fp);
	void rbs(string b);
	void rsample(string b, NODE1 *rt, int n);
	void radd(string b, int n, NODE1 *rt);
	void rburst(NODE1 *bn0);
	void rtraverse(NODE1 *t);
	void rcheck(NODE1 *t, int d);
	void rkill(NODE1 *t);
	void savetrie(NODE1 *t, const char* path, const char* dset);
	void puttrie(NODE1 *t, FILE *f);
	void inssorts(string *a, int n, int d);
	void stabilize(string *a, int n);
	void mkqsorts(string *a, int n, int d);

	/* CPL-burstsort */
	void pbdo(int nr, string fp);
	void pbs(string b);
	void psample(string b, NODE2 *rt, int n);
	void psamburst(NODE2 *bn);
	void pclear(NODE2 *t);
	void padd(string b, int n, NODE2 *rt);
	void pburst(NODE2 *bn0);
	void ptraverse(NODE2 *t);
	void pcheck(NODE2 *t, int d);
	void pkill(NODE2 *t);
	void psavetrie(NODE2 *t, const char* path, const char* dset);
	void pputtrie(NODE2 *t, FILE *f);
	void inssortp(string *a, int n, int d);
	void mkqsortp(string *a, int n, int d);

	/* utility routines */
	void vswap(string *a, string *b, int n);
	string *med3(string *a, string *b, string *c, int d);
	string *med3s(string *a, string *b, string *c);
	
	void cr();
	void dot();
	void say(const char* s);
	void sayln(const char* s);
	void isay(int i);	  void isaye(int i, const char* s);
	void dsay(double d, int n);
	void esay(const char* s);
	void br();
	void brp(const char* s);
	
	string allof(string s);
	int tokenize(string s, string *tok, char d, char dd);
	
	char* fp(const char* dir, const char* fn, const char* ft);
	
	void treset(int nt, int nr);
        static inline void ton(int x) {}
	timer *tsort(int i);
	double tmean(int i);
	double tmed(int i);
	void tmin(int nr);
	void tminsay();
	void tnormsay();
	void isaye(int i,string s);	
#endif
