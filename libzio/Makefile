#
# Makefile for compiling libzio
#
# Author: Werner Fink, <werner@suse.de>
#

LARGE	= $(shell getconf LFS_CFLAGS)
CFLAGS	= $(RPM_OPT_FLAGS) -pipe -Wall -D_GNU_SOURCE -D_REENTRANT $(LARGE)
CC	= gcc
MAJOR	= 1
MINOR	= 00
VERSION	= $(MAJOR).$(MINOR)
SONAME	= libzio.so.$(MAJOR)
LDMAP	= -Wl,--version-script=zio.map

prefix	= /usr
libdir	= $(prefix)/lib
datadir	= $(prefix)/share
mandir	= $(datadir)/man
incdir	= $(prefix)/include

FILES	= README	\
	  COPYING	\
	  Makefile	\
	  zio.c		\
	  zioP.h	\
	  zio.h.in	\
	  testt.c	\
	  tests.c	\
	  lzw.h		\
	  unlzw.c	\
	  zio.map	\
	  fzopen.3.in

### Includes and Defines (add further entries here):

cc-include = $(shell $(CC) $(INCLUDES) -include $(1) -S -o /dev/null -xc /dev/null > /dev/null 2>&1 && echo "-D$(2)")
cc-library = $(shell echo 'int main () { return 0; }' |$(CC) -l$(1:lib%=%) -o /dev/null -xc - > /dev/null 2>&1 && echo yes)
cc-function = $(shell echo 'extern void $(1)(void); int main () { $(1)(); return 0; }' |$(CC) -o /dev/null -xc - > /dev/null 2>&1 && echo "-D$(2)")

CFLAGS += $(call cc-include,libio.h,HAVE_LIBIO_H)
CFLAGS += $(call cc-function,fopencookie,HAVE_FOPENCOOKIE)
CFLAGS += $(call cc-function,funopen,HAVE_FUNOPEN)

CFLAGS += $(call cc-include,zlib.h,HAS_ZLIB_H)
CFLAGS += $(call cc-include,bzlib.h,HAS_BZLIB_H)
CFLAGS += $(call cc-include,lzmadec.h,HAS_LZMADEC_H)
CFLAGS += $(call cc-include,lzma.h,HAS_LZMA_H)

LIBS  = -lz
ifeq ($(call cc-library,libbz2),yes)
 LIBS += -lbz2
endif
ifeq ($(call cc-library,liblzma),yes)
 LIBS += -llzma
else
ifeq ($(call cc-library,lzmadec),yes)
 LIBS += -llzmadec
endif
endif

all: shared static
noweak: CFLAGS += -DNO_WEAK
noweak: LINK += $(LIBS)
noweak: all
shared:	libzio.so.$(VERSION) zio.map
static:	libzio.a

obj/zio.o: zio.c zioP.h zio.h
	test -d obj/ || mkdir obj/
	$(CC) $(CFLAGS) -o $@ -c $<

obs/zio.o: zio.c zioP.h zio.h
	test -d obs/ || mkdir obs/
	$(CC) $(CFLAGS) -fPIC -o $@ -c $<

obj/unlzw.o: unlzw.c lzw.h
	test -d obj/ || mkdir obj/
	$(CC) $(CFLAGS) -funroll-loops -o $@ -c $<

obs/unlzw.o: unlzw.c lzw.h
	test -d obs/ || mkdir obs/
	$(CC) $(CFLAGS) -fPIC -o $@ -c $<

libzio.a: obj/zio.o obj/unlzw.o
	ar -rv $@ $^
	ranlib $@

libzio.so.$(VERSION): obs/zio.o obs/unlzw.o
	gcc -shared -Wl,-soname,$(SONAME),-stats $(LDMAP) -o $@ $^ $(LINK)

zioP.h: /usr/include/bzlib.h /usr/include/zlib.h
zio.h:  zio.h.in /usr/include/stdio.h
	sed 's/@@VERSION@@/$(VERSION)/' < $< > $@
fzopen.3: fzopen.3.in
	sed 's/@@VERSION@@/$(VERSION)/' < $< > $@

unlzw.c: lzw.h

install: install-shared install-static install-data

install-shared: libzio.so.$(VERSION)
	mkdir -p $(DESTDIR)$(libdir)
	install -m 0755 libzio.so.$(VERSION) $(DESTDIR)$(libdir)/
	ln -sf libzio.so.$(VERSION)	     $(DESTDIR)$(libdir)/libzio.so.$(MAJOR)
	ln -sf libzio.so.$(MAJOR)	     $(DESTDIR)$(libdir)/libzio.so

install-static: libzio.a
	mkdir -p $(DESTDIR)$(libdir)
	install -m 0644 libzio.a             $(DESTDIR)$(libdir)/

install-data: zio.h fzopen.3
	mkdir -p $(DESTDIR)$(incdir)
	mkdir -p $(DESTDIR)$(mandir)/man3
	install -m 0644 zio.h                $(DESTDIR)$(incdir)/
	install -m 0644 fzopen.3             $(DESTDIR)$(mandir)/man3/

clean:
	rm -f *.a *.so* testt tests zio.h
	rm -rf obj/ obs/
	rm -f libzio-$(VERSION).tar.gz

dest:   clean
	mkdir libzio-$(VERSION)
	cp -p $(FILES) libzio-$(VERSION)/
	tar czf libzio-$(VERSION).tar.gz libzio-$(VERSION)/
	rm -rf libzio-$(VERSION)/

testt:  testt.c libzio.a
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

tests:  tests.c libzio.a
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
