/*
   This is an implementation of Adaptive radixsort, a simple
   modification of the standard MSD radixsort using distributive
   partitioning. A more complete discussion of the algorithm,
   the implementation and a comparison with other well known sorting
   algorithms, both radix sorting and comparison-based methods, can
   be found in:

      S. Nilsson. Radix Sorting & Searching. PhD thesis, Department
      of Computer Science, Lund University, 1996.

   The thesis can be fetched from my homepage:

      http://www.cs.hut.fi/~sni/

   The code presented in this file has been tested with care but is
   not guaranteed for any purpose. The writer does not offer any
   warranties nor does he accept any liabilities with respect to
   the code.

   Stefan Nilsson, 26 oct 1996.

   Laboratory of Information Processing Science
   Helsinki University of Technology
   Stefan.Nilsson@hut.fi
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include "../tools/contest.hpp"

#define CHARS 256
typedef int character;
typedef unsigned char char_t;
typedef char_t *string;

typedef int boolean;
#define TRUE 1
#define FALSE 0

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define MAXBLOCKS 100

/******************* A simple memory allocator ********************/

typedef struct {
	void *block[MAXBLOCKS];
	int allocnr;
	int nr;
	int blocksize;
	void *current, *first, *last;
} memory;

/*
   Initialize a memory allocator m that can allocate at most
   blocksize * MAXBLOCKS memory slots each consisting of
   elemsize bytes.
*/
static inline void initmem(memory *m, int elemsize, int blocksize)
{
   m->blocksize = MAX(blocksize, 1000) * elemsize;
   m->block[0] = malloc(m->blocksize);
   m->allocnr = 0;
   m->nr = -1;
   m->current = m->first = m->last = NULL;
}

static inline void *allocmem(memory *m, int elemsize)
{
   if (m->current < m->last)
      m->current = (void *) ((char *) (m->current) + elemsize);
   else {
      if (++m->nr >= MAXBLOCKS) return NULL;
      if (m->nr > m->allocnr) {
         m->block[m->nr] = malloc(m->blocksize);
         m->allocnr++;
      }
      m->current = m->first = m->block[m->nr];
      m->last = (void *)
         ((char *) (m->first) + m->blocksize - elemsize);
   }
   return m->current;
}

static inline void *deallocmem(memory *m, int elemsize)
{
   if (m->current > m->first)
      m->current = (void *) ((char *) (m->current) - elemsize);
   else {
      m->nr--;
      if (m->nr >= 0) {
         m->first = m->block[m->nr];
         m->current = m->last = (void *)
            ((char *) (m->first) + m->blocksize - elemsize);
      } else
         m->current = m->first = m->last = NULL;
   }
   return m->current;
}

static inline void resetmem(memory *m)
{
   m->nr = -1;
   m->current = m->first = m->last = NULL;
}

static inline void freemem(memory *m)
{
   int i;

   for (i = 0; i <= m->allocnr; i++)
      free(m->block[i]);
   m->nr = m->allocnr = -1;
   m->current = m->first = m->last = NULL;
}

/********************* Insertion sort ************************/

typedef struct listrec *list;
struct listrec {
   string str;
   list next;
   int length;
};

static inline int scmp(unsigned char *s1, unsigned char *s2)
{
    while( *s1 != '\0' && *s1 == *s2 )
        s1++, s2++;
    return( *s1-*s2 );
}

/* Insertion sort for linked lists of character strings.
   The strings all have a common prefix of length p. */
static inline list Insertsort(list r, list *tail, int p)
{
   list fi, la, t;

   for (fi = la = r, r = r->next; r; r = la->next)
      if (scmp(r->str+p, la->str+p) >= 0) /* add to tail */
         la = r;
      else if (scmp(r->str+p, fi->str+p) <= 0) { /* add to head */
         la->next = r->next;
         r->next = fi;
         fi = r;
      } else { /* insert into middle */
         for (t = fi; scmp(r->str+p, t->next->str+p) >= 0; )
            t = t->next;
         la->next = r->next;
         r->next = t->next;
         t->next = r;
      }
   *tail = la;
   return fi;
}

/************ Adaptive Radixsort (8,16-bit version) ***************/

#define INSERT_BREAK 12
#define BYTE_BREAK 1500

#define BUCKETS CHARS*CHARS
#define CHAR(s, p) s[p]
#define SHORT(s, p) s[p] << 8 | (s[p] ? s[p+1] : 0)
#define IS_ENDMARK(ch) ((ch & 255) == 0)

#define HIGH(ch) ch >> 8
#define LOW(ch) ch & 255

typedef struct bucketrec {
   list head, tail;
   int size; /* size of list, 0 if finished */
} bucket;

typedef struct stackrec {
   list head, tail;
   int size; /* size of list, 0 if finished */
   int pos;  /* current position in string */
} stack;

static memory stackmem[1];
static stack *stackp;

void push(list head, list tail, int size, int pos)
{   
   stackp = (stack *) allocmem(stackmem, sizeof(struct stackrec));
   stackp->head = head;
   stackp->tail = tail;
   stackp->size = size;
   stackp->pos = pos;
}

stack *pop()
{
   stack *temp;

   temp = stackp;
   stackp = (stack *) deallocmem(stackmem, sizeof(struct stackrec));
   return temp;
}

stack *top()
{
   return stackp;
}

boolean stackempty()
{

   return (!stackp);
}

void intobucket1(bucket *b, list h, list t, int size,
                 character ch, character *chmin, character *chmax)
{
   if (!b->head) {
      b->head = h;
      b->tail = t;
      b->size = size;
      if (ch != '\0'){
         if (ch < *chmin) *chmin = ch;
         if (ch > *chmax) *chmax = ch;
      }
   } else {
      b->tail->next = h;
      b->tail = t;
      b->size += size;
   }
}

void intobucket2(bucket *b, list h, list t, int size,
                 character ch, int *used1, int *used2)
{
   if (!b->head) {
      b->head = h;
      b->tail = t;
      b->size = size;
      used1[HIGH(ch)] = used2[LOW(ch)] = TRUE;
   } else {
      b->tail->next = h;
      b->tail = t;
      b->size += size;
   }
}

void ontostack(bucket *b, int pos)
{
   b->tail->next = NULL;
   if (b->size <= INSERT_BREAK) {
      if (b->size > 1)
         b->head = Insertsort(b->head, &b->tail, pos);
      b->size = 0; /* finished */
   }
   if (!b->size && !stackempty() && !top()->size) {
      top()->tail->next = b->head;
      top()->tail = b->tail;
   }
   else {
      push(b->head, b->tail, b->size, pos);
      b->size = 0;
   }
   b->head = NULL;
}

void onebyte(list a, int pos)
{
   static bucket b[CHARS];
   bucket *bp;
   character ch, prevch;
   character chmin = CHARS-1, chmax = 0;
   list t = a, tn;
   int size = 1;

   prevch = CHAR(t->str, pos);           /* into buckets */
   for ( ; (tn = t->next); t = tn) {
      ch = CHAR(tn->str, pos); size++;
      if (ch == prevch) continue;
      intobucket1(b+prevch,a, t, size-1, prevch, &chmin, &chmax);
      a = tn;
      prevch = ch; size = 1;
   }
   intobucket1(b+prevch, a, t, size, prevch, &chmin, &chmax);

   if (b->head) {    /* ch = '\0', end of string */
      b->size = 0;   /* finished */
      ontostack(b, pos);
   }
   for (bp = b + chmin; bp <= b + chmax; bp++)
      if (bp->head) ontostack(bp, pos+1);
}

void twobytes(list a, int pos)
{
   static bucket b[BUCKETS]; /* buckets */
   character ch, prevch;
   list t = a, tn;
   int size = 1;
   int used1[CHARS]; /* What buckets are used? */
   int used2[CHARS];
   int buckets1 = 0, buckets2 = 0;
   character ch1, ch2, high;

   for (ch = 0; ch < CHARS; ch++)
      used1[ch] = used2[ch] = FALSE;

   prevch = SHORT(t->str, pos);              /* into buckets */
   for ( ; (tn = t->next); t = tn) {
      ch = SHORT(tn->str, pos); size++;
      if (ch == prevch) continue;
      intobucket2(b+prevch,a, t, size-1, prevch, used1, used2);
      a = tn;
      prevch = ch; size = 1;
   }
   intobucket2(b+prevch, a, t, size, prevch, used1, used2);

   for (ch = 0; ch < CHARS; ch++) {
      if (used1[ch]) used1[buckets1++] = ch;
      if (used2[ch]) used2[buckets2++] = ch;
   }

   for (ch1 = 0; ch1 < buckets1; ch1++) {  /* put onto stack */
      high = used1[ch1] << 8;
      for (ch2 = 0; ch2 < buckets2; ch2++) {
         ch = high | used2[ch2];
         if (b[ch].head) {
            if IS_ENDMARK(ch) b[ch].size = 0;    /* finished */
            ontostack(b+ch, pos+2);
         }
      }
   }
}

list MSD(list a, int n)
{
   list res = NULL;
   stack *s;

   if (n < 2) return a;
   initmem(stackmem, sizeof(struct stackrec), n/50);
   push(a, NULL, n, 0);

   while (!stackempty()) {
      s = pop();
      if (!s->size) { /* finished */
         s->tail->next = res;
         res = s->head;
         continue;
      }
      if (s->size <= BYTE_BREAK)
         onebyte(s->head, s->pos);
      else
         twobytes(s->head, s->pos);
   }

   freemem(stackmem);
   return res;
}

/****************************** Glue by Rantala *************************/

void arssort(string strings[], size_t scnt)
{
    list listnodes;
    size_t i;

    /* allocate memory based on the number of strings in the array */
    listnodes = (list ) calloc(scnt, sizeof(struct listrec));

    /* point the linked list nodes to the strings in the array */
    for( i=0; i<scnt; i++)
    {
        listnodes[i].str = strings[i];
        if (i<(scnt-1))
            listnodes[i].next = &listnodes[i+1]; 
        else
            listnodes[i].next = NULL;
    }

    /* sort */
    list sortednodes = MSD(listnodes, scnt);

    /* write the strings back into the array */
    for (i = 0;  i < scnt ; i++, sortednodes=sortednodes->next)
        strings[i] = sortednodes->str;

    free(listnodes);
}

void nilsson_adaptive_msd(unsigned char **strings, size_t n)
{
    return arssort((char_t**)strings, n);
}

PSS_CONTESTANT(nilsson_adaptive_msd,
               "nilsson/adaptive_msd",
               "Adaptive MSD Radix Sort by Stefan Nilsson")

/******************************************************************************/
