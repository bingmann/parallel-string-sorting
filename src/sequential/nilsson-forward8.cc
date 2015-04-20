/*
   Forward radixsort with a fixed sized alphabet. The algorithm
   inspects one character at a time. This code will work well for
   alphabets of small size (8 bits). Larger alphabets (16 bits or
   more) may, however, require some heuristic to avoid inspecting
   empty buckets.

   S. Nilsson. Radix Sorting and Searching. PhD thesis, Department
   of Computer Science, Lund University, 1990.

   The code presented in this file has been tested with care but is
   not guaranteed for any purpose. The writer does not offer any
   warranties nor does he accept any liabilities with respect to
   the code.

   Stefan Nilsson, 8 jan 1997.

   Laboratory of Information Processing Science
   Helsinki University of Technology
   Stefan.Nilsson@hut.fi
*/

#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "../tools/contest.h"

#define CHARS 256
#define INSERTBREAK 20
typedef unsigned char* string;

#define MAXBLOCKS 100
#define TRUE 1
#define FALSE 0
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef int boolean;
typedef int character;

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

/**************** Forward radix sort (bitstrings) *******************/

#define IS_ENDMARK(ch) (ch == '\0')
#define CHAR(s, p) s[p]

typedef struct grouprec *group;
typedef struct bucketrec *bucket;

struct grouprec {
   list head, tail; /* a list of elements */
   group next;      /* the next group */
   group nextunf;   /* the next unfinished group */
   group insp;      /* insertion point */
   boolean finis;   /* is the group finished? */
};
/* The group structure member insp is used to make splitting of
   groups possible during the phase where elements are moved from
   buckets back into their previous groups. The group structure
   member finis indicates if the elements in the group are sorted;
   this information makes it easy to skip finished groups during a
   traversal of the group data structure */

struct bucketrec {
   list head, tail; /* a list of elements */
   int size;        /* list length */
   group tag;       /* group tag */
   bucket next;     /* next bucket item */
};

static memory groupmem[1];
static memory bucketmem[1];

/* Put a list of elements into a bucket. We distinguish between two
   cases. If the first bucket item has the same tag as the list to
   be inserted the list is just appended, otherwise a new bucket
   is created. */
static void intobucket(bucket *b, list head, list tail,
                       int size, group g)
{
   bucket btemp = *b, newb;

   if (!btemp || btemp->tag != g) {   /* create new tag */
      newb = (bucket) allocmem(bucketmem, sizeof(struct bucketrec));
      newb->next = btemp;
      newb->head = head;
      newb->size = size;
      newb->tag = g;
      *b = btemp = newb;
   } else {   /* append */
      btemp->tail->next = head;
      btemp->size += size;
   }
   tail->next = NULL;
   btemp->tail = tail;
}

/* Travers the groups and put the elements into buckets.
   The parameter pos indicates the current position in the string.
   To be able to skip groups that are already sorted we keep track
   of the previous group. Also, the previously read character is
   recorded. In this way it is possible to move the elements in
   blocks consisting of strings that have a common character in
   position pos. Furthermore, a group that is not split during this
   phase is left behind and not put into a bucket. */
static void intobuckets(group g, bucket b[], int pos)
{
   group prevg;
   character ch, prevch;
   boolean split;
   list tail, tailn;
   int size;

   resetmem(bucketmem);
   for (prevg = g, g = g->nextunf ; g; g = g->nextunf) {
      if (g->finis)
         {prevg->nextunf = g->nextunf; continue;}
      tail = g->head; split = FALSE;
      prevch = CHAR(tail->str, pos); size = 1;
      for ( ; (tailn = tail->next); tail = tailn) {
         ch = CHAR(tailn->str, pos); size++;
         if (ch == prevch) continue;
         intobucket(b+prevch, g->head, tail, size-1, g);
         g->head = tailn; split = TRUE;
         prevch = ch; size = 1;
      }
      if (split) {
         intobucket(b+prevch, g->head, tail, size, g);
         g->head = NULL;
         prevg = g;
      } else if (IS_ENDMARK(prevch))
         prevg->nextunf = g->nextunf;
      else
         prevg = g;
   }
}

/* Put a list into group g and, at the same time, split g.
   If two consecutive groups are both finished, there is no need
   to perform any splitting. */
static void intogroup(group g, list head, list tail, boolean finis)
{
   group newg;

   if (!g->head) {   /* back into old group */
      g->head = head;
      g->tail = tail;
      g->finis = finis;
      g->insp = g;
   } else if (finis && g->insp->finis) {  /* don't split if both */
      g->insp->tail->next = head;         /* groups are finished */
      g->insp->tail = tail;
   }
   else {   /* split */
      newg = (group) allocmem(groupmem, sizeof(struct grouprec));
      newg->head = head;
      newg->tail = tail;
      newg->next = g->insp->next;
      newg->nextunf = g->insp->nextunf;
      newg->finis = finis;
      g->insp = g->insp->nextunf = g->insp->next = newg;
   }
}

/* Traverse the buckets and put the elements back into their groups.
   Split the groups and mark all finished groups.
   The elements are moved in blocks. */
static void intogroups(bucket b[], int pos)
{
   character ch;
   bucket s;
   boolean finis;

   for (ch = 0; ch < CHARS; ch++) {
      if (!b[ch]) continue;
      for (s = b[ch]; s; s = s->next) {
         finis = IS_ENDMARK(ch);
         if (s->size < INSERTBREAK && !finis) {
            if (s->size > 1)
               s->head = Insertsort(s->head, &s->tail, pos);
            finis = TRUE;
         }
         intogroup(s->tag, s->head, s->tail, finis);
      }
      b[ch] = NULL;
   }
}

/* Travers the groups and return the elements in sorted order. */
static list collect(group g)
{
   list head, tail;

   g = g->next;
   head = g->head;
   tail = g->tail;
   for (g = g->next; g; g = g->next) {
      tail->next = g->head;
      tail = g->tail;
   }
   return head;
}

static inline list forward1(list t, int n)
{
   static bucket b[CHARS];   /* buckets */
   group g, g2;              /* groups */
   int pos = 0;              /* pos in string */

   if (n<2) return t;

   initmem(groupmem, sizeof(struct grouprec), n/15);
   initmem(bucketmem, sizeof(struct bucketrec), n/5);

   /* We use a dummy group g as the header of the group data
      structure. It does not contain any elements, but only a
      pointer to the first unfinished group. */
   g = (group) allocmem(groupmem, sizeof(struct grouprec));
   g2 = (group) allocmem(groupmem, sizeof(struct grouprec));
   g->next = g->nextunf = g2;
   g2->head = t;
   g2->next = g2->nextunf = NULL; 
   g2->finis = FALSE;

   intobuckets(g, b, pos);
   while (g->nextunf) {
      pos++;
      intogroups(b, pos);
      intobuckets(g, b, pos);
   }
   t = collect(g);

   freemem(bucketmem);
   freemem(groupmem);

   return t;
}

/****************************** Glue by Rantala *************************/

void frssort1(string strings[], size_t scnt)
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
    list sortednodes = forward1(listnodes, scnt);

    /* write the strings back into the array */
    for (i = 0;  i < scnt ; i++, sortednodes=sortednodes->next)
        strings[i] = sortednodes->str;

    free(listnodes);
}

void nilsson_forward8(unsigned char **strings, size_t n)
{
    return frssort1(strings, n);
}

PSS_CONTESTANT(nilsson_forward8,
               "nilsson/forward8",
               "Forward Radix Sort 8-bit by Stefan Nilsson")
