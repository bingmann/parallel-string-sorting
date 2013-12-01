/*
   Multikey quicksort, a radix sort algorithm for arrays of character
   strings by Bentley and Sedgewick.

   J. Bentley and R. Sedgewick. Fast algorithms for sorting and
   searching strings. In Proceedings of 8th Annual ACM-SIAM Symposium
   on Discrete Algorithms, 1997.

   http://www.CS.Princeton.EDU/~rs/strings/index.html

   The code presented in this file has been tested with care but is
   not guaranteed for any purpose. The writer does not offer any
   warranties nor does he accept any liabilities with respect to
   the code.

   Stefan Nilsson, 8 jan 1997.

   Laboratory of Information Processing Science
   Helsinki University of Technology
   Stefan.Nilsson@hut.fi

   Updated with std::swap and std::min by Timo Bingmann in 2012.
*/

#include "inssort.h"

namespace bs_mkqs {

typedef unsigned char* string;

#define i2c(i) x[i][depth]

static void vecswap(int i, int j, int n, string x[])
{
    while (n-- > 0) {
        std::swap(x[i], x[j]);
        i++;
        j++;
    }
}

static inline void ssort1(string x[], int n, int depth)
{
    int    a, b, c, d, r, v;
    if (n <= 1)
        return;
    a = rand() % n;
    std::swap(x[0], x[a]);
    v = i2c(0);
    a = b = 1;
    c = d = n-1;
    for (;;) {
        while (b <= c && (r = i2c(b)-v) <= 0) {
            if (r == 0) { std::swap(x[a], x[b]); a++; }
            b++;
        }
        while (b <= c && (r = i2c(c)-v) >= 0) {
            if (r == 0) { std::swap(x[c], x[d]); d--; }
            c--;
        }
        if (b > c) break;
        std::swap(x[b], x[c]);
        b++;
        c--;
    }
    r = std::min(a, b-a);     vecswap(0, b-r, r, x);
    r = std::min(d-c, n-d-1); vecswap(b, n-r, r, x);
    r = b-a; ssort1(x, r, depth);
    if (i2c(r) != 0)
        ssort1(x + r, a + n-d-1, depth+1);
    r = d-c; ssort1(x + n-r, r, depth);
}

static inline void multikey1(string x[], int n)
{ ssort1(x, n, 0); }


/* ssort2 -- Faster Version of Multikey Quicksort */

static inline void vecswap2(string *a, string *b, int n)
{
    while (n-- > 0) {
        string t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

#define ptr2char(i) (*(*(i) + depth))

static inline string *med3func(string *a, string *b, string *c, int depth)
{
    int va, vb, vc;
    if ((va=ptr2char(a)) == (vb=ptr2char(b)))
        return a;
    if ((vc=ptr2char(c)) == va || vc == vb)
        return c;       
    return va < vb
        ? (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}

static inline void ssort2(string a[], size_t n, int depth)
{
    int d, r, partval;
    string *pa, *pb, *pc, *pd, *pl, *pm, *pn;
    if (n < 64) {
        return inssort::inssort(a, n, depth);
    }
    pl = a;
    pm = a + (n/2);
    pn = a + (n-1);
    if (n > 30) { /* On big arrays, pseudomedian of 9 */
        d = (n/8);
        pl = med3func(pl, pl+d, pl+2*d, depth);
        pm = med3func(pm-d, pm, pm+d,   depth);
        pn = med3func(pn-2*d, pn-d, pn, depth);
    }
    pm = med3func(pl, pm, pn, depth);
    std::swap(*a, *pm);
    partval = ptr2char(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;
    for (;;) {
        while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
            if (r == 0) std::swap(*pa++, *pb);
            pb++;
        }
        while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
            if (r == 0) std::swap(*pc, *pd--);
            pc--;
        }
        if (pb > pc) break;
        std::swap(*pb++, *pc--);
    }
    pn = a + n;
    r = std::min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
    if ((r = pb-pa) > 1)
        ssort2(a, r, depth);
    if (ptr2char(a + r) != 0)
        ssort2(a + r, pa-a + pn-pd-1, depth+1);
    if ((r = pd-pc) > 1)
        ssort2(a + n-r, r, depth);
}

static inline void multikey2(string a[], size_t n)
{ ssort2(a, n, 0); }

static void bs_mkqsort(unsigned char **strings, size_t n)
{
    return multikey2(strings, n);
}

CONTESTANT_REGISTER(bs_mkqsort,
                    "bs/mkqsort",
                    "bs_mkqs Original Multikey-Quicksort")

#undef i2c
#undef ptr2char

} // namespace bs_mkqs

// global procedure for base sorting
static inline void mkqsort(unsigned char **strings, size_t n, int depth)
{
    return bs_mkqs::ssort2(strings, n, depth);
}
