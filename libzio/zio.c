/*
 * zio.c	Provide an streamable interface to gziped/bzip2ed/LZW/LZMA files
 *
 * Copyright 2004 Werner Fink, 2004 SuSE LINUX AG, Germany.
 * Copyright 2006 Werner Fink, 2006 SuSE Products GmbH, Germany.
 * Copyright 2009 Werner Fink, 2009 SuSE Products GmbH, Germany.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Author:	Werner Fink <werner@suse.de>
 */

#include "zioP.h"
#include "zio.h"
#include "lzw.h"

#if defined(HAS_ZLIB_H)
static ssize_t   zread(void *cookie, char *buf, size_t count)
{
    if (!cookie) {
	errno = EINVAL;
	return -1;
    }
    return (ssize_t)gzread((gzFile)cookie, (voidp)buf, count);
}

static ssize_t   zwrite(void *cookie, const char *buf, size_t count)
{
    if (!cookie) {
	errno = EINVAL;
	return -1;
    }
    return (ssize_t)gzwrite((gzFile)cookie, (const voidp)buf, count);
}

static zio_int_t zseek(void *cookie, zio_off_t *poffset, int whence)
{
    if (!cookie) {
	errno = EINVAL;
	return -1;
    }
    return (zio_int_t)gzseek((gzFile)cookie, (z_off_t)(*poffset), whence);
}

static int       zclose(void *cookie)
{
    int status;
    if (!cookie) {
	errno = EINVAL;
	return -1;
    }
    (void)gzflush((gzFile)cookie, Z_FINISH);
    status = gzclose((gzFile)cookie);
    return (status >= 0) ? 0 : EOF;
}

__extension__
static cookie_io_functions_t ioz = {
    .read  = (cookie_read_function_t*) &zread,
    .write = (cookie_write_function_t*)&zwrite,
    .seek  = (cookie_seek_function_t*) &zseek,
    .close = (cookie_close_function_t*)&zclose,
};
#else /* !HAS_ZLIB_H */
# error No support for `.gz' nor `.z' format
#endif /* !HAS_ZLIB_H */

#if defined(HAS_BZLIB_H)
# ifndef MIN
#  define MIN(x,y) ((x) < (y) ? (x) : (y))
# endif

typedef struct bzfile_s {
    size_t total_out;
    BZFILE *file;
    char *mode;
    char *path;
    int fd;
} bzfile_t;

static ssize_t   bzread(void *cookie, char *buf, size_t count)
{
    bzfile_t *bzf = (bzfile_t*)cookie;
    ssize_t len = -1;
    if (!bzf)
	goto out;
    if (bzf->file)
	len = (ssize_t)BZ2_bzread(bzf->file, (void*)buf, count);
    if (len > 0)
	bzf->total_out += len;
out:
    if (len < 0)
	errno = EINVAL;
    return len;
}

static ssize_t   bzwrite(void *cookie, const char *buf, size_t count)
{
    bzfile_t *bzf = (bzfile_t*)cookie;
    ssize_t len = -1;
    if (!bzf)
	goto out;
    if (bzf->file)
	len = (ssize_t)BZ2_bzread(bzf->file, (void*)buf, count);
    if (len > 0)
	bzf->total_out += len;
out:
    if (len < 0)
	errno = EINVAL;
    return len;
}

static zio_int_t bzseek(void *cookie, zio_off_t *poffset, int whence)
{
    bzfile_t *bzf = (bzfile_t*)cookie;
    off_t offset = (off_t)*poffset;
    off_t oldpos, newpos;
    if (!bzf) {
	errno = EINVAL;
	return -1;
    }

    oldpos = (off_t)bzf->total_out;
    switch (whence) {
    case SEEK_SET:
	if (offset < 0)
	    return -1;
	newpos = offset;
	break;
    case SEEK_CUR:
	if ((offset < 0 && (off_t)(-1 * offset) > oldpos) || (offset > 0 && (offset+oldpos) < oldpos))
	    return -1;
	newpos = (off_t)bzf->total_out + offset;
	break;
    case SEEK_END:
	newpos = -1;
	break;
    default:
	errno = EINVAL;
	return -1;
    }

    if (whence != SEEK_END && newpos < oldpos) {
	int status = BZ2_bzflush(bzf->file);
	BZ2_bzclose(bzf->file);
	if (status < 0) {
	    errno = EINVAL;
	    return -1;
	}
	if (bzf->fd >= 0) {
	    lseek(bzf->fd, 0, SEEK_SET);
	    bzf->file = BZ2_bzdopen(bzf->fd, bzf->mode);
	} else if (bzf->path) {
	    bzf->file = BZ2_bzopen(bzf->path, bzf->mode);
	} else {
	    errno = EINVAL;
	    return -1;
	}
	if (bzf->file == (BZFILE*)0) {
	    errno = EINVAL;
	    return -1;
	}
	bzf->total_out = 0;
    }
    if (newpos == oldpos)
	return oldpos;
    else {
	char buf[1<<12];
	while (newpos > oldpos || newpos == -1) {
	    size_t  req_size = MIN(sizeof(buf), newpos - oldpos);
	    ssize_t got_size = BZ2_bzread(bzf->file, buf, req_size);
	    if (got_size != (ssize_t)(req_size)) {
		if (got_size < 0)
		    return -1;
		else {
		    newpos = oldpos + got_size;
		    break;
		}
	    }
	    oldpos += got_size;
	}
	return newpos;
    }
}

static int       bzclose(void *cookie)
{
    bzfile_t *bzf = (bzfile_t*)cookie;
    int status = -1;
    if (!bzf) {
	errno = EINVAL;
	goto out;
    }
    if (bzf->file) {
	status = BZ2_bzflush(bzf->file);
	BZ2_bzclose(bzf->file);
    }
    free(cookie);
out:
    return (status >= 0) ? 0 : EOF;
}

__extension__
static cookie_io_functions_t iobz = {
    .read  = (cookie_read_function_t*) &bzread,
    .write = (cookie_write_function_t*)&bzwrite,
    .seek  = (cookie_seek_function_t*) &bzseek,
    .close = (cookie_close_function_t*)&bzclose,
};
#else /* !HAS_BZLIB_H */
# warning No support for .bz2 format
#endif /* !HAS_BZLIB_H */

#if defined(HAS_LZMA_H)
# ifndef MIN
#  define MIN(x,y) ((x) < (y) ? (x) : (y))
# endif

typedef struct lzfile_s {
    uint8_t buf[1<<12];
    lzma_stream strm;
    FILE *file;
    int encoding;
    int level;
    int what;
    int eof;
} lzfile_t;

static lzma_ret lzmaopen(lzma_stream *__restrict strm, const char mode, const char what, int level)
{
    lzma_ret ret;
    if (mode == 'w') {
	if (what == 'x')
	    ret = lzma_easy_encoder(strm, level, LZMA_CHECK_CRC32);
	else {
	    lzma_options_lzma opt;
	    lzma_lzma_preset(&opt, level);
	    ret = lzma_alone_encoder(strm, &opt);
	}
   } else
	ret = lzma_auto_decoder(strm, 100<<20, 0);
   return ret;
}

static ssize_t   lzmaread(void *cookie, char *buf, size_t count)
{
    lzfile_t *lzma = (lzfile_t*)cookie;
    lzma_stream *strm;
    ssize_t eof = 0;
    if (!lzma || lzma->encoding)
	return -1;
    if (lzma->eof)
	return 0;
    strm = &lzma->strm;

    strm->next_out = (uint8_t*)buf;
    strm->avail_out = count;
    while (1) {
	lzma_ret lret;
	if (!strm->avail_in) {
	    strm->next_in = lzma->buf;
	    strm->avail_in = fread(lzma->buf, 1, sizeof(lzma->buf), lzma->file);
	    if (strm->avail_in == 0)
		eof = 1;
	}
	lret = lzma_code(strm, LZMA_RUN);
	if (lret == LZMA_STREAM_END) {
	    lzma->eof = 1;
	    return count - strm->avail_out;
	}
	if (lret != LZMA_OK)
	    return -1;
	if (strm->avail_out == 0)
	    return count;
	if (eof) {
	    eof = feof(lzma->file);
	    clearerr(lzma->file);
	    if (eof)
		break;
	    return -1;
	}
    }
    return 0;
}

static ssize_t   lzmawrite(void *cookie, const char *buf, size_t count)
{
    lzfile_t *lzma = (lzfile_t*)cookie;
    lzma_stream *strm;
    if (!lzma || !lzma->encoding)
	return -1;
    if (!count)
	return 0;
    strm = &lzma->strm;

    strm->next_in = (uint8_t*)buf;
    strm->avail_in = count;
    while (1) {
	lzma_ret lret;
	size_t len;
	strm->next_out = lzma->buf;
	strm->avail_out = sizeof(lzma->buf);
	lret = lzma_code(strm, LZMA_RUN);
	if (lret != LZMA_OK)
	    break;
	len = sizeof(lzma->buf) - strm->avail_out;
	if (len && fwrite(lzma->buf, 1, len, lzma->file) != len)
	    break;
	if (strm->avail_in == 0)
	    return len;
    }
    return -1;
}

static zio_int_t lzmaseek(void *cookie, zio_off_t *poffset, int whence)
{
    lzfile_t *lzma = (lzfile_t*)cookie;
    off_t offset = (off_t)*poffset;
    lzma_stream *strm;
    off_t oldpos, newpos;
    if (!lzma)
	return -1;
    strm = &lzma->strm;

    oldpos = (off_t)strm->total_out;
    switch (whence) {
    case SEEK_SET:
	if (offset < 0)
	    return -1;
	newpos = offset;
	break;
    case SEEK_CUR:
	if ((offset < 0 && (off_t)(-1 * offset) > oldpos) || (offset > 0 && (offset+oldpos) < oldpos))
	    return -1;
	newpos = (off_t)strm->total_out + offset;
	break;
    case SEEK_END:
	newpos = -1;
	break;
    default:
	errno = EINVAL;
	return -1;
    }

    if (whence != SEEK_END && newpos < oldpos) {
	lzma_ret ret;
	lzma_end(strm);
	rewind(lzma->file);
	ret = lzmaopen(strm, lzma->encoding ? 'w' : 'r', lzma->what, lzma->level);
	if (ret != LZMA_OK || strm->total_out != 0) {
	    fclose(lzma->file);
	    errno = EINVAL;
	    return -1;
	}
    }
    if (newpos == oldpos)
	return oldpos;
    else {
	char buf[sizeof(lzma->buf)];
	while (newpos > oldpos || newpos == -1) {
	    size_t  req_size = MIN(sizeof(buf), newpos - oldpos);
	    ssize_t got_size = lzmaread(cookie, buf, req_size);
	    if (got_size != (ssize_t)(req_size)) {
		if (got_size < 0)
		    return -1;
		else {
		    newpos = oldpos + got_size;
		    break;
		}
	    }
	    oldpos += got_size;
	}
	return newpos;
    }
}

static int       lzmaclose(void *cookie)
{
    lzfile_t *lzma = (lzfile_t*)cookie;
    lzma_ret lret = LZMA_STREAM_END;
    lzma_stream *strm;
    int fret = -1;
    if (!lzma)
	return -1;
    if (!lzma->encoding)
	goto out;
    strm = &lzma->strm;
    while (1) {
	size_t len;
	strm->avail_out = sizeof(lzma->buf);
	strm->next_out = (uint8_t*)lzma->buf;
	lret = lzma_code(strm, LZMA_FINISH);
	if (lret != LZMA_OK && lret != LZMA_STREAM_END)
	    goto out;
	len = sizeof(lzma->buf) - strm->avail_out;
	if (len && fwrite(lzma->buf, 1, len, lzma->file) != len)
	    goto out;
	if (lret == LZMA_STREAM_END)
	    break;
    }
    lzma_end(strm);
out:
    fret = fclose(lzma->file);
    free(lzma);
    if (lret != LZMA_STREAM_END)
	fret = -1;
    return fret;
}

__extension__
static cookie_io_functions_t iolzma = {
    .read  = (cookie_read_function_t*) &lzmaread,
    .write = (cookie_write_function_t*)&lzmawrite,
    .seek  = (cookie_seek_function_t*) &lzmaseek,
    .close = (cookie_close_function_t*)&lzmaclose,
};
#else /* !HAS_LZMA_H */
# if defined(HAS_LZMADEC_H)
static ssize_t   lzmaread(void *cookie, char *buf, size_t count)
{
    if (!cookie) {
	errno = EINVAL;
	return -1;
    }
    return lzmadec_read((lzmadec_FILE*)cookie, (uint8_t*)buf, count);
}

static ssize_t   lzmawrite(void *cookie, const char *buf, size_t count)
{
    errno = ENOTSUP;
    return -1;
}

static zio_int_t lzmaseek(void *cookie, zio_off_t *poffset, int whence)
{
    if (!cookie) {
	errno = EINVAL;
	return -1;
    }
    return (zio_int_t)lzmadec_seek((lzmadec_FILE*)cookie, (off_t)(*poffset), whence);
}

static int       lzmaclose(void *cookie)
{
    if (!cookie) {
	errno = EINVAL;
	return -1;
    }
    int_fast8_t status = lzmadec_close((lzmadec_FILE*)cookie);
    return (status >= 0) ? 0 : EOF;
}

__extension__
static cookie_io_functions_t iolzma = {
    .read  = (cookie_read_function_t*) &lzmaread,
    .write = (cookie_write_function_t*)&lzmawrite,
    .seek  = (cookie_seek_function_t*) &lzmaseek,
    .close = (cookie_close_function_t*)&lzmaclose,
};
# else /* !HAS_LZMADEC_H */
#  warning No support for .lzma format
# endif /* !HAS_LZMADEC_H */
#endif /* !HAS_LZMA_H */

typedef struct lzwfile_s {
    size_t total_out;
    LZW_t *file;
    char *mode;
    char *path;
    int fd;
} lzwfile_t;

static ssize_t   lzwread(void *cookie, char *buf, size_t count)
{
    lzwfile_t *lzw = (lzwfile_t*)cookie;
    ssize_t len = -1;
    if (!lzw)
	goto out;
    if (lzw->file)
	len = readlzw(lzw->file, buf, count);
    if (len > 0)
	lzw->total_out += len;
out:
    if (len < 0)
	errno = EINVAL;
    return len;
}

static ssize_t   lzwwrite(void *cookie, const char *buf, size_t count)
{
    errno = ENOTSUP;
    return -1;
}

static zio_int_t lzwseek(void *cookie, zio_off_t *poffset, int whence)
{
    lzwfile_t *lzw = (lzwfile_t*)cookie;
    off_t offset = (off_t)*poffset;
    off_t oldpos, newpos;
    if (!lzw) {
	errno = EINVAL;
	return -1;
    }

    oldpos = (off_t)lzw->total_out;
    switch (whence) {
    case SEEK_SET:
	if (offset < 0)
	    return -1;
	newpos = offset;
	break;
    case SEEK_CUR:
	if ((offset < 0 && (off_t)(-1 * offset) > oldpos) || (offset > 0 && (offset+oldpos) < oldpos))
	    return -1;
	newpos = (off_t)lzw->total_out + offset;
	break;
    case SEEK_END:
	newpos = -1;
	break;
    default:
	errno = EINVAL;
	return -1;
    }

    if (whence != SEEK_END && newpos < oldpos) {
	closelzw(lzw->file);
	if (lzw->fd >= 0) {
	    lseek(lzw->fd, 0, SEEK_SET);
	    lzw->file = dopenlzw(lzw->fd, lzw->mode);
	} else if (lzw->path) {
	    lzw->file = openlzw(lzw->path, lzw->mode);
	} else {
	    errno = EINVAL;
	    return -1;
	}
	if (lzw->file == (LZW_t*)0) {
	    errno = EINVAL;
	    return -1;
	}
	lzw->total_out = 0;
    }
    if (newpos == oldpos)
	return oldpos;
    else {
	char buf[1<<12];
	while (newpos > oldpos || newpos == -1) {
	    size_t  req_size = MIN(sizeof(buf), newpos - oldpos);
	    ssize_t got_size = readlzw(lzw->file, buf, req_size);
	    if (got_size != (ssize_t)(req_size)) {
		if (got_size < 0)
		    return -1;
		else {
		    newpos = oldpos + got_size;
		    break;
		}
	    }
	    oldpos += got_size;
	}
	return newpos;
    }
}

static int       lzwclose(void *cookie)
{
    lzwfile_t *lzw = (lzwfile_t*)cookie;
    if (!lzw) {
	errno = EINVAL;
	return -1;
    }
    if (lzw->file)
	closelzw(lzw->file);
    free(cookie);
    return 0;
}

__extension__
static cookie_io_functions_t iolzw = {
    .read  = (cookie_read_function_t*) &lzwread,
    .write = (cookie_write_function_t*)&lzwwrite,
    .seek  = (cookie_seek_function_t*) &lzwseek,
    .close = (cookie_close_function_t*)&lzwclose,
};

static inline char autodetect(char **__restrict path, const char *__restrict check)
{
    const size_t len = strlen(*path);
    char *suff = strrchr(*path, '.');
    char *ext = *path;
    char what = 'n';

    if (suff) {
	suff++;
	if      (strcmp(suff, "z"   ) == 0)
	    what = 'z';
	else if (strcmp(suff, "gz"  ) == 0)
	    what = 'g';
	else if (strcmp(suff, "Z"   ) == 0)
	    what = 'Z';
	else if (strcmp(suff, "bz2" ) == 0)
	    what = 'b';
	else if (strcmp(suff, "lzma") == 0)
	    what = 'l';
	else if (strcmp(suff, "xz"  ) == 0)
	    what = 'x';
    }

    if (what == 'n' && *check == 'r') {
	int olderr, fd;
	struct stat st;
	char m[5];
	ext = malloc(sizeof(char)*(len + 5 + 1));
	if (!ext)
	    goto out;
	strcpy(ext, *path);
	suff = (ext+len);

	olderr = errno;
	if (stat(strcat(ext, ".gz"),  &st) == 0) {
	    what = 'g';
	    goto skip;
	}
	*suff = '\0';
	if (stat(strcat(ext, ".bz2"), &st) == 0) {
	    what = 'b';
	    goto skip;
	}
	*suff = '\0';
	if (stat(strcat(ext, ".z"),   &st) == 0) {
	    what = 'z';
	    goto skip;
	}
	*suff = '\0';
	if (stat(strcat(ext, ".Z"),   &st) == 0) {
	    what = 'Z';
	    goto skip;
	}
	*suff = '\0';
	if (stat(strcat(ext, ".lzma"), &st) == 0) {
	    what = 'l';
	    goto skip;
	}
	*suff = '\0';
	if (stat(strcat(ext, ".xz"), &st) == 0) {
	    what = 'x';
	    goto skip;
	}
	*suff = '\0';

	if ((fd = open(ext, O_RDONLY|O_NOCTTY)) < 0)
	    goto skip;
	if (read(fd, m, sizeof(m)) == sizeof(m)) {
	    if (m[0] == '\037' && m[1] == '\213')
		what = 'g';
	    if (m[0] == '\037' && m[1] == '\235')
		what = 'Z';
	    if (m[0] == '\037' && m[1] == '\236')
		what = 'z';
	    else if (m[0] == 'B' && m[1] == 'Z' && m[2] == 'h')
		what = 'b';
	    else if (m[0] == ']' && m[1] == '\0' && m[2] == '\0' && m[3] == '\200') /* weak!! */
		what = 'l';
	    else if (m[0] == '\377' && m[1] == 'L' && m[2] == 'Z' && m[3] == 'M' && m[4] == 'A')
		what = 'l';
	    else if (m[0] == '\375' && m[1] == '7' && m[2] == 'z' && m[3] == 'X' && m[4] == 'Z')
		what = 'x';
	}
	close(fd);
    skip:
	errno = olderr;
    }
out:
    *path = ext;
    return what;
}

FILE * fzopen(const char * path, const char * mode)
{
    FILE * ret = (FILE *)0;
    char * check = (char*)0, * ext = (char*)0;
    size_t n = 0, len;
    unsigned int i;
    char what = 'n';

    if (!mode || !(n = strlen(mode))) {
	errno = EINVAL;
	goto out;
    }

    if (!(check = (char*)malloc(n*sizeof(char))))
	goto out;

    /* No append mode possible */
    switch (*mode) {
	case 'r': check[0] = 'r'; break;
	case 'w': check[0] = 'w'; break;
	default:  errno = EINVAL; goto out;
    }

    for (i = 1; i < n; i++) {
	/* We can only open for reading OR writing but NOT for both */
	switch (mode[i]) {
	    case '\0': break;
	    case '+': errno = EINVAL; goto out;
	    case 'b': case 'x': check[i] = mode[i]; continue;
	    /* Ingore switches for gzopen() */
	    case 'f': case 'h': check[i] = '\0';    continue;
	    default:		check[i] = '\0';    continue;
	}
	break;
    }

    if (!path || !(len = strlen(path))) {
	errno = EINVAL;
	goto out;
    }

    ext = (char *)path;
    what = autodetect(&ext, check);

    switch (what) {
    case 'g':
    case 'z':		/* Is this correct? Old gzip magic */
#if defined(HAS_ZLIB_H)
	{
	    gzFile cookie;

	    if (&gzopen == NULL) {
		errno = ENOSYS;
		goto out;
	    }

	    if (!(cookie = gzopen(ext, mode))) {
		if (!errno)
		    errno = ENOMEM;
		goto out;
	    }

	    if (!(ret = fopencookie((void*)cookie, check, ioz))) {
		gzclose(cookie);
		errno = EINVAL;
		goto out;
	    }
# ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
# endif
	}
#else /* !HAS_ZLIB_H */
	errno = ENOTSUP;
#endif /* !HAS_ZLIB_H */
	break;
    case 'Z':
	{
	    lzwfile_t *__restrict cookie;
	    if (*mode != 'r') {
		errno = ENOTSUP;
		goto out;
	    }

	    if (posix_memalign((void*)&cookie, sizeof(void*), alignof(lzwfile_t)+strsize(ext)+strsize(mode)) != 0)
		goto out;
	    memset(cookie, 0, alignof(lzwfile_t)+strsize(ext)+strsize(mode));

	    cookie->fd = -1;
	    cookie->mode = ((char*)cookie)+alignof(lzwfile_t);
	    cookie->path = cookie->mode + strsize(mode);
	    strcpy(cookie->mode, mode);
	    strcpy(cookie->path, ext);

	    if (!(cookie->file = openlzw(ext, mode))) {
		free(cookie);
		if (!errno)
		    errno = ENOMEM;
		goto out;
	    }

	    if (!(ret = fopencookie((void*)cookie, check, iolzw))) {
		closelzw(cookie->file);
		free(cookie);
		errno = EINVAL;
		goto out;
	    }
# ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
# endif
	}
	break;
    case 'b':
#if defined(HAS_BZLIB_H)
	{
	    bzfile_t *__restrict cookie;
	    int level = 5;

	    if (&BZ2_bzopen == NULL) {
		errno = ENOSYS;
		goto out;
	    }

	    if (posix_memalign((void*)&cookie, sizeof(void*), alignof(bzfile_t)+strsize(ext)+strsize(mode)) != 0)
		goto out;
	    memset(cookie, 0, alignof(bzfile_t)+strsize(ext)+strsize(mode));

	    for (i = 1; i < n; i++) {
		if (mode[i] >= '0' && mode[i] <= '9') {
		    level = (int)mode[i];
		    break;
		}
	    }

	    cookie->fd = -1;
	    cookie->mode = ((char*)cookie)+alignof(bzfile_t);
	    cookie->path = cookie->mode+strsize(mode);
	    strcpy(cookie->mode, mode);
	    strcpy(cookie->path, ext);

	    if (!(cookie->file = BZ2_bzopen(ext, mode))) {
		free(cookie);
		if (!errno)
		    errno = ENOMEM;
		goto out;
	    }

	    if (!(ret = fopencookie((void*)cookie, check, iobz))) {
		BZ2_bzclose(cookie->file);
		free(cookie);
		errno = EINVAL;
		goto out;
	    }
# ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
# endif
	}
#else /* !HAS_BZLIB_H */
	errno = ENOTSUP;
#endif /* !HAS_BZLIB_H */
	break;
#if defined(HAS_LZMA_H)
    case 'l':
    case 'x':
	{
	    int level = LZMA_PRESET_DEFAULT;
	    lzfile_t *__restrict cookie;
	    lzma_ret lret;

	    if (&lzma_auto_decoder == NULL) {
		errno = ENOSYS;
		goto out;
	    }

	    if (posix_memalign((void*)&cookie, sizeof(void*), alignof(lzfile_t)) != 0)
		goto out;
	    memset(cookie, 0, alignof(lzfile_t));

	    if ((cookie->file = fopen(ext, check)) == NULL) {
		free(cookie);
		goto out;
	    }

	    for (i = 1; i < n; i++) {
		if (mode[i] >= '0' && mode[i] <= '9') {
		    level = (int)mode[i];
		    break;
		}
	    }

	    cookie->eof = 0;
	    cookie->what = what;
	    cookie->strm = (lzma_stream)LZMA_STREAM_INIT;
	    cookie->level = level;
	    cookie->encoding = (check[0] == 'w') ? 1 : 0;
	    lret = lzmaopen(&cookie->strm, check[0], what, level);

	    if (lret != LZMA_OK) {
		fclose(cookie->file);
		free(cookie);
		errno = EINVAL;
		goto out;
	    }
	    if (!(ret = fopencookie((void*)cookie, check, iolzma))) {
		lzma_end(&cookie->strm);
		fclose(cookie->file);
		free(cookie);
		errno = EINVAL;
		goto out;
	    }
#  ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
#  endif
	}
	break;
#else /* !HAS_LZMA_H */
    case 'x':
	errno = ENOTSUP;
	break;
    case 'l':
# if defined(HAS_LZMADEC_H)
	{
	    lzmadec_FILE* cookie;

	    if (*mode != 'r') {
		errno = ENOTSUP;
		goto out;
	    }

	    if (&lzmadec_open == NULL) {
		errno = ENOSYS;
		goto out;
	    }

	    if (!(cookie = lzmadec_open(ext))) {
		if (!errno)
		    errno = ENOMEM;
		goto out;
	    } 

	    if (!(ret = fopencookie((void*)cookie, check, iolzma))) {
		lzmadec_close(cookie);
		errno = EINVAL;
		goto out;
	    }
#  ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
#  endif
	}
# else /* !HAS_LZMADEC_H */
	errno = ENOTSUP;
# endif /* !HAS_LZMADEC_H */
#endif /* !HAS_LZMA_H */
	break;
    default:
	ret = fopen(ext, mode);
	break;
    }

out:
    if (check)
	free(check);
    if (ext && ext != path)
	free(ext);

    return ret;
}

FILE * fdzopen(int fildes, const char * mode, const char *what)
{
    FILE * ret = (FILE *)0;
    char * check = (char*)0;
    size_t n = 0;
    unsigned int i;

    if (!mode || !(n = strlen(mode))) {
	errno = EINVAL;
	goto out;
    }

    if (!(check = (char*)malloc(n*sizeof(char))))
	goto out;

    /* No append mode possible */
    switch (*mode) {
	case 'r': check[0] = 'r'; break;
	case 'w': check[0] = 'w'; break;
	default:  errno = EINVAL; goto out;
    }

    for (i = 1; i < n; i++) {
	/* We can only open for reading OR writing but NOT for both */
	switch (mode[i]) {
	    case '\0': break;
	    case '+': errno = EINVAL; goto out;
	    case 'b': case 'x': check[i] = mode[i]; continue;
	    /* Ingore switches for gzopen() */
	    case 'f': case 'h': check[i] = '\0';    continue;
	    default:		check[i] = '\0';    continue;
	}
	break;
    }

    switch (*what) {
    case 'g':
    case 'z':		/* Is this correct? Old gzip magic */
#if defined(HAS_ZLIB_H)
	{
	    gzFile cookie;

	    if (&gzdopen == NULL) {
		errno = ENOSYS;
		goto out;
	    }

	    if (!(cookie = gzdopen(fildes, mode))) {
		if (!errno)
		    errno = ENOMEM;
		goto out;
	    }

	    if (!(ret = fopencookie((void*)cookie, check, ioz))) {
		gzclose(cookie);
		errno = EINVAL;
		goto out;
	    }
# ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
# endif
	}
#else /* !HAS_ZLIB_H */
	errno = ENOTSUP;
#endif /* !HAS_ZLIB_H */
	break;
    case 'Z':
	{
	    lzwfile_t *__restrict cookie;
	    if (*mode != 'r') {
		errno = ENOTSUP;
		goto out;
	    }

	    if (posix_memalign((void*)&cookie, sizeof(void*), alignof(lzwfile_t)+strsize(mode)) != 0)
		goto out;
	    memset(cookie, 0, alignof(lzwfile_t)+strsize(mode));

	    cookie->fd = fildes;
	    cookie->mode = ((char*)cookie)+alignof(lzwfile_t);
	    strcpy(cookie->mode, mode);

	    if (!(cookie->file = dopenlzw(fildes, mode))) {
		free(cookie);
		if (!errno)
		    errno = ENOMEM;
		goto out;
	    }

	    if (!(ret = fopencookie((void*)cookie, check, iolzw))) {
		closelzw(cookie->file);
		free(cookie);
		errno = EINVAL;
		goto out;
	    }
# ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
# endif
	}
	break;
    case 'b':
#if defined(HAS_BZLIB_H)
	{
	    bzfile_t *__restrict cookie;
	    int level = 5;

	    if (&BZ2_bzdopen == NULL) {
		errno = ENOSYS;
		goto out;
	    }

	    if (posix_memalign((void*)&cookie, sizeof(void*), alignof(bzfile_t)+strsize(mode)) != 0)
		goto out;
	    memset(cookie, 0, alignof(bzfile_t)+strsize(mode));

	    for (i = 1; i < n; i++) {
		if (mode[i] >= '0' && mode[i] <= '9') {
		    level = (int)mode[i];
		    break;
		}
	    }

	    cookie->fd = fildes;
	    cookie->mode = ((char*)cookie)+alignof(bzfile_t);
	    strcpy(cookie->mode, mode);

	    if (cookie->mode == (char*)0) {
		free(cookie);
		goto out;
	    }
	    if (!(cookie->file = BZ2_bzdopen(fildes, mode))) {
		free(cookie);
		if (!errno)
		    errno = ENOMEM;
		goto out;
	    }

	    if (!(ret = fopencookie((void*)cookie, check, iobz))) {
		BZ2_bzclose(cookie->file);
		free(cookie);
		errno = EINVAL;
		goto out;
	    }
# ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
# endif
	}
#else /* !HAS_BZLIB_H */
	errno = ENOTSUP;
#endif /* !HAS_BZLIB_H */
	break;
#if defined(HAS_LZMA_H)
    case 'l':
    case 'x':
	{
	    int level = LZMA_PRESET_DEFAULT;
	    lzfile_t *__restrict cookie;
	    lzma_ret lret;

	    if (&lzma_auto_decoder == NULL) {
		errno = ENOSYS;
		goto out;
	    }

	    if (posix_memalign((void*)&cookie, sizeof(void*), alignof(lzfile_t)) != 0)
		goto out;
	    memset(cookie, 0, alignof(lzfile_t));

	    if ((cookie->file = fdopen(fildes, check)) == NULL) {
		free(cookie);
		goto out;
	    }

	    for (i = 1; i < n; i++) {
		if (mode[i] >= '0' && mode[i] <= '9') {
		    level = (int)mode[i];
		    break;
		}
	    }

	    cookie->eof = 0;
	    cookie->what = *what;
	    cookie->strm = (lzma_stream)LZMA_STREAM_INIT;
	    cookie->level = level;
	    cookie->encoding = (check[0] == 'w') ? 1 : 0;
	    lret = lzmaopen(&cookie->strm, check[0], *what, level);
	    if (lret != LZMA_OK) {
		fclose(cookie->file);
		free(cookie);
		errno = EINVAL;
		goto out;
	    }
	    if (!(ret = fopencookie((void*)cookie, check, iolzma))) {
		lzma_end(&cookie->strm);
		fclose(cookie->file);
		free(cookie);
		errno = EINVAL;
		goto out;
	    }
#  ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
#  endif
	}
	break;
#else /* !HAS_LZMA_H */
    case 'x':
	errno = ENOTSUP;
	break;
    case 'l':
# if defined(HAS_LZMADEC_H)
	{
	    lzmadec_FILE* cookie;

	    if (*mode != 'r') {
		errno = ENOTSUP;
		goto out;
	    }

	    if (&lzmadec_open == NULL) {
		errno = ENOSYS;
		goto out;
	    }

	    if (!(cookie = lzmadec_dopen(fildes))) {
		if (!errno)
		    errno = ENOMEM;
		goto out;
	    } 

	    if (!(ret = fopencookie((void*)cookie, check, iolzma))) {
		lzmadec_close(cookie);
		errno = EINVAL;
		goto out;
	    }
#  ifndef LIBIO_IS_FIXED
	    if (ret->_fileno < 0)
		ret->_fileno = 0;
#  endif
	}
# else /* !HAS_LZMADEC_H */
	errno = ENOTSUP;
# endif /* !HAS_LZMADEC_H */
#endif /* !HAS_LZMA_H */
	break;
    default:
	ret = fdopen(fildes, mode);
	break;
    }

out:
    if (check)
	free(check);
    return ret;
}
