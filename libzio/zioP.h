/*
 * zioP.h	Internal header for libzio, including required
 *		standard glibc header, zlib.h, and bzlib.h.
 *		Making the used libz and bzlib functions weak symbols.
 *
 * Copyright 2004 Werner Fink, 2004 SuSE LINUX AG, Germany.
 * Copyright 2006 Werner Fink, 2006 SuSE Products GmbH, Germany.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Author:      Werner Fink <werner@suse.de>
 */

#ifndef _ZIO_P_H
#define _ZIO_P_H

#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#ifdef HAVE_LIBIO_H
# include <libio.h>
#endif
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>

#ifndef  unused
# define unused		__attribute__((__unused__))
#endif
#ifndef  nonnull
# define nonnull(parm)	__attribute__((__nonnull__ parm))
#endif
#ifndef  wur
# define wur		__attribute__((__warn_unused_result__))
#endif
#define alignof(type)	(sizeof(type)+(sizeof(type)%sizeof(void*)))
#define strsize(str)	((strlen(str)+1)*sizeof(char))

#if !defined(HAVE_FOPENCOOKIE) && !defined(HAVE_FUNOPEN)
# error Requires fopencookie(3GNU) or funopen(3BSD)
#endif

#if defined(HAVE_LIBIO_H) || defined(HAVE_FOPENCOOKIE)
# if defined __GLIBC__ && __GLIBC__ > 1
#  undef  LIBIO_IS_FIXED
#  if __GLIBC__ > 2 || (__GLIBC__ >= 2 && __GLIBC_MINOR__ > 0)
#   define LIBIO_IS_FIXED
#  endif
# else
#  error The libzio requires the GLIBC
# endif
#endif

#if defined __GNUC__
#  if defined __USE_ISOC99
#    define _cat_pragma(exp)	_Pragma(#exp)
#    define _weak_pragma(exp)	_cat_pragma(weak name)
#  else
#    define _weak_pragma(exp)
#  endif
#  define _declare(name)	__extension__ extern __typeof__(name) name
#  define weak_symbol(name)	_weak_pragma(name) _declare(name) __attribute__((weak))
#else
#  error The libzio requires the GCC
#endif

#if defined(HAS_ZLIB_H)
# include <zlib.h>
# ifndef NO_WEAK
weak_symbol(gzopen);
weak_symbol(gzdopen);
weak_symbol(gzread);
weak_symbol(gzwrite);
weak_symbol(gzseek);
weak_symbol(gzflush);
weak_symbol(gzclose);
# endif
#endif

#if defined(HAS_BZLIB_H)
# include <bzlib.h>
# ifndef NO_WEAK
weak_symbol(BZ2_bzopen);
weak_symbol(BZ2_bzdopen);
weak_symbol(BZ2_bzread);
weak_symbol(BZ2_bzwrite);
/* no BZ2_bzseek */
weak_symbol(BZ2_bzflush);
weak_symbol(BZ2_bzclose);
# endif
#endif

#if defined(HAS_LZMA_H)
# include <stdint.h>
# include <lzma.h>
# ifndef NO_WEAK
weak_symbol(lzma_easy_encoder);
weak_symbol(lzma_lzma_preset);
weak_symbol(lzma_alone_encoder);
weak_symbol(lzma_auto_decoder);
weak_symbol(lzma_code);
weak_symbol(lzma_end);
# endif
#else /* !HAS_LZMA_H */
# if defined(HAS_LZMADEC_H)
#  include <stdint.h>
#  include <lzmadec.h>
#  ifndef NO_WEAK
weak_symbol(lzmadec_open);
weak_symbol(lzmadec_dopen);
weak_symbol(lzmadec_read);
/* no lzmadec_write() */
weak_symbol(lzmadec_seek);
weak_symbol(lzmadec_close);
/* no lzmadec_flush() */
#  endif
# endif
#endif /* !HAS_LZMA_H */

#if defined(HAVE_FOPENCOOKIE)
# undef HAVE_FUNOPEN
__extension__ typedef off_t   zio_off_t;
__extension__ typedef int     zio_int_t;
# if !defined(LIBIO_IS_FIXED)
__extension__ typedef _IO_cookie_io_functions_t cookie_io_functions_t;
__extension__ typedef ssize_t cookie_read_function_t  __P ((void *, char *, size_t));
__extension__ typedef ssize_t cookie_write_function_t __P ((void *, const char *, size_t));
__extension__ typedef int     cookie_seek_function_t  __P ((void *, off_t, int));
__extension__ typedef int     cookie_close_function_t __P ((void *));
# endif
#endif
#if defined(HAVE_FUNOPEN)
__extension__ typedef size_t zio_off_t;
__extension__ typedef fpos_t zio_int_t;
__extension__ typedef int    cookie_read_function_t  __P ((void *, char *, int));
__extension__ typedef int    cookie_write_function_t __P ((void *, const char *, int));
__extension__ typedef fpos_t cookie_seek_function_t  __P ((void *, fpos_t, int));
__extension__ typedef int    cookie_close_function_t __P ((void *));
__extension__ typedef struct
{
    cookie_read_function_t  *read;
    cookie_write_function_t *write;
    cookie_seek_function_t  *seek;
    cookie_close_function_t *close;
} cookie_io_functions_t;
static __inline__ FILE *fopencookie(void *__restrict,
				    const char *__restrict,
				    cookie_io_functions_t) nonnull((1,2)) wur;
static __inline__ FILE *fopencookie(void *__restrict cookie,
				    const char *__restrict mode unused,
				    cookie_io_functions_t io_funcs)
{
    return funopen(cookie, io_funcs.read, io_funcs.write, io_funcs.seek, io_funcs.close);
}
#endif
#endif /* _ZIO_P_H */
