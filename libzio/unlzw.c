/* unlzw.c -- decompress files in LZW format.
 * The code in this file is directly derived from the public domain 'compress'
 * written by Spencer Thomas, Joe Orost, James Woods, Jim McKie, Steve Davies,
 * Ken Turkowski, Dave Mack and Peter Jannesen.
 *
 * This is a temporary version which will be rewritten in some future version
 * to accommodate in-memory decompression.
 *
 * Tue Dec 12 17:54:07 CET 2006 - werner@suse.de
 * Be able to emulate a zlib-like behaviour: open, read, and close .Z files
 * in memory. I'm using context switchting and a global allocated structure
 * to be able to read during the main loop in unlzw() does its work. For this
 * nearly _all_ variables affected by the context switch are forward to this
 * structure, even the stack and the context type its self.
 * The oringal source was adopted from the gzip version 1.3.7.
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#define DIST_BUFSIZE	0x2000
#define INBUFSIZ	0x2000
#define WSIZE		0x2000
#define INBUF_EXTRA	0x40
#define OUTBUFSIZ	0x2000
#define OUTBUF_EXTRA	0x800
#define STACK_SIZE	0x1000

#include "lzw.h"

static void unlzw(LZW_t *in);

/* ===========================================================================
 * Fill the input buffer. This is called only when the buffer is empty.
 * Adopted from gzip version 1.3.7 util.c
 */
__extension__
static __inline__ int fill_inbuf(LZW_t *in)
{
    /* Read as much as possible */
    in->insize = 0;
    do {
	ssize_t len = read(in->ifd, in->inbuf + in->insize, INBUFSIZ - in->insize);
	if (len < 0) {
	    if (errno == EINTR)
		continue;
	    perror(__FUNCTION__);
	    break;
	}
	if (len == 0)
	    break;
	in->insize += len;
    } while (in->insize < INBUFSIZ);

    if (in->insize == 0)
	return EOF;

    in->bytes_in += (off_t)in->insize;
    in->inptr = 1;
    return (in->inbuf)[0];
}

/* ===========================================================================
 * Does the same as write(), but also handles partial pipe writes and checks
 * for error return.
 * Adopted from gzip version 1.3.7 util.c
 * Note that this version uses context switching, switch back to old context.
 */
__extension__
static __inline__ void write_buf(LZW_t *in, const unsigned char* buf, size_t cnt)
{
    do {
	if ((in->tsize = (in->tcount > cnt) ? cnt : in->tcount)) {
	    (void)memcpy(in->transfer, buf, in->tsize);
	    buf += in->tsize;
	    cnt -= in->tsize;
	}
	swapcontext(&in->uc[1], &in->uc[0]);
    } while (cnt);
}

#define get_byte(in)	((in)->inptr < (in)->insize ? (in)->inbuf[(in)->inptr++] : fill_inbuf((in)))
#define memzero(s,n)	memset ((void*)(s), 0, (n))

#define MAXCODE(n)	(1L << (n))

#ifndef	BYTEORDER
#	define	BYTEORDER	0000
#endif

#ifndef	NOALLIGN
#	define	NOALLIGN	0
#endif


union	bytes {
    long  word;
    struct {
#if BYTEORDER == 4321
	uint8_t	b1;
	uint8_t	b2;
	uint8_t	b3;
	uint8_t	b4;
#else
#if BYTEORDER == 1234
	uint8_t	b4;
	uint8_t	b3;
	uint8_t	b2;
	uint8_t	b1;
#else
#	undef	BYTEORDER
	int  dummy;
#endif
#endif
    } bytes;
};

#if BYTEORDER == 4321 && NOALLIGN == 1
#  define input(b,o,c,n,m){ \
     (c) = (*(uint32_t *)(&(b)[(o)>>3])>>((o)&0x7))&(m); \
     (o) += (n); \
   }
#else
#  define input(b,o,c,n,m){ \
     uint8_t *p = &(b)[(o)>>3]; \
     (c) = ((((long)(p[0]))|((long)(p[1])<<8)| \
     ((long)(p[2])<<16))>>((o)&0x7))&(m); \
     (o) += (n); \
   }
#endif

#define tab_prefixof(i)		in->tab_prefix[i]
#define clear_tab_prefixof(in)	memzero((in)->tab_prefix, sizeof(unsigned short)*(1<<(BITS)));
#define de_stack		((uint8_t *)(&in->d_buf[DIST_BUFSIZE-1]))
#define tab_suffixof(i)		in->tab_suffix[i]

/* ============================================================================
 * Decompress in to out.  This routine adapts to the codes in the
 * file building the "string" table on-the-fly; requiring no table to
 * be stored in the compressed file.
 * IN assertions: the buffer inbuf contains already the beginning of
 *   the compressed data, from offsets iptr to insize-1 included.
 *   The magic header has already been checked and skipped.
 *   bytes_in and bytes_out have been initialized.
 *
 * Adopted from gzip version 1.3.7 unlzw.c
 * This is mainly the head of the old unlzw() before its main loop.
 */
LZW_t *openlzw(const char * path, const char *mode)
{
    int fildes;
    if (!mode || *mode != 'r')
	goto out;
    if ((fildes = open(path, O_RDONLY)) < 0)
        goto out;
    return dopenlzw(fildes, mode);
out:
    return (LZW_t*)0;
}

LZW_t *dopenlzw(int fildes, const char *mode)
{
    LZW_t *in = (LZW_t*)0;
    uint8_t magic[2];
    sigset_t sigmask, oldmask;

    if (!mode || *mode != 'r')
	goto out;

    if ((in = (LZW_t*)malloc(sizeof(LZW_t))) == (LZW_t*)0)
	goto err;
    memset(in, 0, sizeof(LZW_t));

    if ((in->inbuf = (uint8_t*)malloc(sizeof(uint8_t)*(INBUFSIZ))) == (uint8_t*)0)
	goto err;
    if ((in->outbuf = (uint8_t*)malloc(sizeof(uint8_t)*(OUTBUFSIZ+OUTBUF_EXTRA))) == (uint8_t*)0)
	goto err;
    if ((in->d_buf = (uint16_t*)malloc(sizeof(uint16_t)*DIST_BUFSIZE)) == (uint16_t*)0)
	goto err;
    if ((in->tab_suffix = (uint8_t*) malloc(sizeof(uint8_t )*(2L*WSIZE))) == (uint8_t*)0)
	goto err;
    if ((in->tab_prefix = (uint16_t*)malloc(sizeof(uint16_t)*(1<<(BITS)))) == (uint16_t*)0)
	goto err;
    if ((in->ifd = fildes) < 0)
	goto err;

    if ((in->stack = (uint8_t*)malloc(STACK_SIZE)) == (uint8_t*)0)
	goto err;
    if ((in->uc = (ucontext_t*)malloc(2*sizeof(ucontext_t))) == (ucontext_t*)0)
	goto err;
    if (getcontext(&in->uc[1]) < 0)
	goto err;
    in->uc[1].uc_link = &in->uc[0];
    in->uc[1].uc_stack.ss_sp = in->stack;
    in->uc[1].uc_stack.ss_size = STACK_SIZE;
    if (sigucmask(SIG_SETMASK, (sigset_t*)0, &sigmask) < 0)
	goto err;
    if (sigaddset(&sigmask, SIGINT) < 0)
	goto err;
    if (sigaddset(&sigmask, SIGQUIT) < 0)
	goto err;
#if defined(__ia64__) && defined(uc_sigmask)	/* On ia64 the type of uc_sigmask is ulong not sigset_t */
    in->uc[1].uc_sigmask = sig_ia64_mask(sigmask);
#else
    in->uc[1].uc_sigmask = sigmask;
#endif
    makecontext(&in->uc[1], (void(*)(void))unlzw, 1, in);

    sigucmask(SIG_SETMASK, &sigmask, &oldmask);
    magic[0] = get_byte(in);
    magic[1] = get_byte(in);
    sigucmask(SIG_SETMASK, &oldmask, &sigmask);

    if (memcmp(magic, LZW_MAGIC, sizeof(magic)))
	goto err;

    in->n.block_mode = BLOCK_MODE;  /* block compress mode -C compatible with 2.0 */
    in->rsize = in->insize;

    in->n.maxbits = get_byte(in);
    in->n.block_mode = in->n.maxbits & BLOCK_MODE;

    if ((in->n.maxbits & LZW_RESERVED) != 0) {
	fprintf(stderr, "%s: warning, unknown flags 0x%x\n",
		__FUNCTION__, in->n.maxbits & LZW_RESERVED);
    }
    in->n.maxbits &= BIT_MASK;
    in->n.maxmaxcode = MAXCODE(in->n.maxbits);

    if (in->n.maxbits > BITS) {
	fprintf(stderr, "%s: compressed with %d bits, can only handle %d bits\n",
		__FUNCTION__, in->n.maxbits, BITS);
	goto err;
    }

    in->n.maxcode = MAXCODE(in->n.n_bits = INIT_BITS)-1;
    in->n.bitmask = (1<<in->n.n_bits)-1;
    in->n.oldcode = -1;
    in->n.finchar = 0;
    in->n.posbits = in->inptr<<3;

    in->n.free_ent = ((in->n.block_mode) ? FIRST : 256);

    clear_tab_prefixof(in);	/* Initialize the first 256 entries in the table. */

    for (in->n.code = 255 ; in->n.code >= 0 ; --in->n.code) {
	tab_suffixof(in->n.code) = (uint8_t)in->n.code;
    }
out:
    return in;
err:
    closelzw(in);
    return (LZW_t*)0;
}

/*
 * New function, simply to free all allocated objects in
 * reverse order and close the input file.
 */
void closelzw(LZW_t * in)
{
    if (in == (LZW_t*)0)
	return;
    if (in->uc) free(in->uc);
    if (in->stack) free(in->stack);
    if (in->ifd >= 0) close(in->ifd);
    if (in->tab_prefix) free(in->tab_prefix);
    if (in->tab_suffix) free(in->tab_suffix);
    if (in->d_buf)  free(in->d_buf);
    if (in->outbuf) free(in->outbuf);
    if (in->inbuf)  free(in->inbuf);
    free(in);
    in = (LZW_t*)0;
}

/*
 * Adopted from gzip version 1.3.7 unlzw.c
 * This is mainly the body of the old unlzw() which is its main loop.
 */
static void unlzw(LZW_t *in)
{
    do {
	int i;
	int e;
	int o;

    resetbuf:
	e = in->insize - (o = (in->n.posbits>>3));

	for (i = 0 ; i < e ; ++i) {
	    in->inbuf[i] = in->inbuf[i+o];
	}
	in->insize = e;
	in->n.posbits = 0;

	if (in->insize < INBUF_EXTRA) {
	    do {
		in->rsize = read(in->ifd, in->inbuf + in->insize, INBUFSIZ - in->insize);
		if (in->rsize < 0) {
		    if (errno == EINTR)
			continue;
		    perror(__FUNCTION__);
		    break;
		}
		if (in->rsize == 0)
		    break;
		in->insize += in->rsize;
	    } while (in->insize < INBUFSIZ);
	    in->bytes_in += (off_t)in->insize;
	}
	in->n.inbits = ((in->rsize != 0) ? ((long)in->insize - in->insize%in->n.n_bits)<<3 :
			((long)in->insize<<3)-(in->n.n_bits-1));

	while (in->n.inbits > in->n.posbits) {
	    if (in->n.free_ent > in->n.maxcode) {
		in->n.posbits = ((in->n.posbits-1) +
			   ((in->n.n_bits<<3)-(in->n.posbits-1+(in->n.n_bits<<3))%(in->n.n_bits<<3)));
		++in->n.n_bits;
		if (in->n.n_bits == in->n.maxbits) {
		    in->n.maxcode = in->n.maxmaxcode;
		} else {
		    in->n.maxcode = MAXCODE(in->n.n_bits)-1;
		}
		in->n.bitmask = (1<<in->n.n_bits)-1;
		goto resetbuf;
	    }
	    input(in->inbuf,in->n.posbits,in->n.code,in->n.n_bits,in->n.bitmask);

	    if (in->n.oldcode == -1) {
		if (256 <= in->n.code)
		    fprintf(stderr, "%s: corrupt input.\n", __FUNCTION__);
		in->outbuf[in->outpos++] = (uint8_t)(in->n.finchar = (int)(in->n.oldcode=in->n.code));
		continue;
	    }
	    if (in->n.code == CLEAR && in->n.block_mode) {
		clear_tab_prefixof(in);
		in->n.free_ent = FIRST - 1;
		in->n.posbits = ((in->n.posbits-1) +
			   ((in->n.n_bits<<3)-(in->n.posbits-1+(in->n.n_bits<<3))%(in->n.n_bits<<3)));
		in->n.maxcode = MAXCODE(in->n.n_bits = INIT_BITS)-1;
		in->n.bitmask = (1<<in->n.n_bits)-1;
		goto resetbuf;
	    }
	    in->n.incode = in->n.code;
	    in->n.stackp = de_stack;

	    if (in->n.code >= in->n.free_ent) { /* Special case for KwKwK string. */
		if (in->n.code > in->n.free_ent) {
#ifdef DEBUG
		    uint8_t *p;
		    in->n.posbits -= in->n.n_bits;
		    p = &in->inbuf[in->n.posbits>>3];
		    fprintf(stderr,
			    "code:%ld free_ent:%ld n_bits:%d insize:%lu\n",
			    in->n.code, in->n.free_ent, in->n.n_bits, in->insize);
		    fprintf(stderr,
			    "posbits:%ld inbuf:%02X %02X %02X %02X %02X\n",
			    in->n.posbits, p[-1],p[0],p[1],p[2],p[3]);
#endif
		    if (in->outpos > 0) {
			write_buf(in, in->outbuf, in->outpos);
			in->bytes_out += (off_t)in->outpos;
			in->outpos = 0;
		    }
		    fprintf(stderr, "%s: corrupt input.\n", __FUNCTION__);
		}
		*--in->n.stackp = (uint8_t)in->n.finchar;
		in->n.code = in->n.oldcode;
	    }

	    /* Generate output characters in reverse order */
	    while ((uint32_t)in->n.code >= (uint32_t)256) {
		*--in->n.stackp = tab_suffixof(in->n.code);
		in->n.code = tab_prefixof(in->n.code);
	    }
	    *--in->n.stackp =	(uint8_t)(in->n.finchar = tab_suffixof(in->n.code));

	    /* And put them out in forward order */
	    if (in->outpos + (in->n.newdif = (de_stack - in->n.stackp)) >= OUTBUFSIZ) {
		do {
		    if (in->n.newdif > OUTBUFSIZ - in->outpos)
			in->n.newdif = OUTBUFSIZ - in->outpos;

		    if (in->n.newdif > 0) {
			memcpy(in->outbuf + in->outpos, in->n.stackp, in->n.newdif);
			in->outpos += in->n.newdif;
		    }
		    if (in->outpos >= OUTBUFSIZ) {
			write_buf(in, in->outbuf, in->outpos);
			in->bytes_out += (off_t)in->outpos;
			in->outpos = 0;
		    }
		    in->n.stackp+= in->n.newdif;
		} while ((in->n.newdif = (de_stack - in->n.stackp)) > 0);
	    } else {
		memcpy(in->outbuf + in->outpos, in->n.stackp, in->n.newdif);
		in->outpos += in->n.newdif;
	    }

	    if ((in->n.code = in->n.free_ent) < in->n.maxmaxcode) { /* Generate the new entry. */

		tab_prefixof(in->n.code) = (uint16_t)in->n.oldcode;
		tab_suffixof(in->n.code) = (uint8_t)in->n.finchar;
		in->n.free_ent = in->n.code+1;
	    }
	    in->n.oldcode = in->n.incode;	/* Remember previous code.	*/
	}
    } while (in->rsize != 0);

    if (in->outpos > 0) {
	write_buf(in, in->outbuf, in->outpos);
	in->bytes_out += (off_t)in->outpos;
	in->tsize = EOF;
	in->outpos = 0;
    }
}

/*
 * New function, simply to read from the output buffer of unlzw().
 * We do this by switching into the context of unlzw() and back
 * to our old context if the provided buffer is filled.
 */
ssize_t readlzw(LZW_t * in, char* buffer, const size_t size)
{
    in->transfer = (uint8_t*)buffer;
    in->tcount = size;
    in->tsize = 0;
    if (in->uc == (ucontext_t*)0)
	return 0;			/* For (f)lex scanner ... */
    swapcontext(&in->uc[0], &in->uc[1]);
    if (in->tsize < 0) {
	free(in->uc);			/* ... do not enter next */
	in->uc = (ucontext_t*)0;
	free(in->stack);
	in->stack = (uint8_t*)0;
	return 0;
    }
    return in->tsize;
}

#ifdef TEST
int main()
{
    ssize_t len;
    char buffer[1024];

    LZW_t *lzw = openlzw("man.1.Z", "r");
    if (!lzw)
	return -1;

    do {
	len = readlzw(lzw, &buffer[0], sizeof(buffer));
	write(1, &buffer[0], len);
    } while (len != 0);

    closelzw(lzw);
    return 0;
}
#endif
