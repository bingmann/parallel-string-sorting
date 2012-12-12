/* lzw.h -- define the lzw functions.

   Copyright (C) 1992-1993 Jean-loup Gailly.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */

#ifndef BITS
#  define BITS 16
#endif
#define INIT_BITS 9              /* Initial number of bits per code */

#define	LZW_MAGIC  "\037\235"   /* Magic header for lzw files, 1F 9D */

#define BIT_MASK    0x1f /* Mask for 'number of compression bits' */
/* Mask 0x20 is reserved to mean a fourth header byte, and 0x40 is free.
 * It's a pity that old uncompress does not check bit 0x20. That makes
 * extension of the format actually undesirable because old compress
 * would just crash on the new format instead of giving a meaningful
 * error message. It does check the number of bits, but it's more
 * helpful to say "unsupported format, get a new version" than
 * "can only handle 16 bits".
 */

#define BLOCK_MODE  0x80
/* Block compression: if table is full and compression rate is dropping,
 * clear the dictionary.
 */

#define LZW_RESERVED 0x60 /* reserved bits */

#define	CLEAR  256       /* flush the dictionary */
#define FIRST  (CLEAR+1) /* first free entry */

#include <stdint.h>
#include <ucontext.h>
#include "zioP.h"

#if defined _REENTRANT || defined _THREAD_SAFE
# include <pthread.h>
weak_symbol(pthread_sigmask);
#endif

__extension__
static __inline__ int sigucmask(int how, const sigset_t * __restrict set, sigset_t * __restrict oset)
{
#if defined _REENTRANT || defined _THREAD_SAFE
    if (&pthread_sigmask)
	return pthread_sigmask(how, set, oset);
    else
#endif
    return sigprocmask(how, set, oset);
}

#if defined(__ia64__) && defined(uc_sigmask)	/* Broken header sys/ucontext.h -> bits/sigcontext.h */
__extension__
static __inline__ unsigned long int sig_ia64_mask(const sigset_t set)
{
    unsigned long int mask = 0;
    int cnt = (8 * sizeof(unsigned long int));
    if (cnt > NSIG) cnt = NSIG;
    while (--cnt >= 0) {
	if (!sigismember(&set, cnt))
	    continue;
	mask |= (1 << (cnt - 1));		/* sigmask() macro is is not usable for BSD way */
    }
    return mask;
}
#endif

typedef struct _LZW_s {
    uint8_t  *inbuf;
    uint8_t  *outbuf;
    uint16_t *d_buf;
    uint8_t  *tab_suffix;
    uint16_t *tab_prefix;
    uint8_t  *transfer;
    off_t     bytes_in;
    off_t     bytes_out;
    size_t    insize;
    size_t    inptr;
    size_t    outpos;
    size_t    tcount;
    ssize_t   tsize;
    ssize_t   rsize;
    struct unlzw_s {
	uint8_t *stackp;
	int32_t  code;
	int      finchar;
	int32_t  oldcode;
	int32_t  incode;
	uint32_t inbits;
	uint32_t posbits;
	uint32_t bitmask;
	int32_t  free_ent;
	int32_t  maxcode;
	int32_t  maxmaxcode;
	off_t    newdif;
	int      n_bits;
	int      block_mode;
	int      maxbits;
    } n;
    int ifd;
    char *ifname;
    ucontext_t *uc;
    uint8_t *stack;
} LZW_t;

extern LZW_t *openlzw(const char * __restrict, const char * __restrict);
extern LZW_t *dopenlzw(int fildes, const char * __restrict);
extern ssize_t readlzw(LZW_t * __restrict, char * __restrict, const size_t);
extern void closelzw(LZW_t * __restrict);
