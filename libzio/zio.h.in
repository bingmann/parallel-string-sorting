/*
 * zio.h	Provide an streamable interface to gziped/bzip2ed files
 *
 * Copyright 2004 Werner Fink, 2004 SuSE LINUX AG, Germany.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Author:	Werner Fink <werner@suse.de>
 */

#ifndef _ZIO_H
#define _ZIO_H

#include <stdio.h>
#define ZIO_VERSION @@VERSION@@
#ifdef __cplusplus
extern "C" {
#endif
/*
 * This function can open a gziped file for reading OR writing, but
 * NOT both (not `+' possible) NOR can it open for appending (no `a').
 */
extern FILE *fzopen __P((__const char *__restrict, __const char *__restrict));
extern FILE *fdzopen __P((int, __const char *__restrict mode, __const char *__restrict));

#ifdef __cplusplus
}
#endif

#endif /* _ZIO_H */
