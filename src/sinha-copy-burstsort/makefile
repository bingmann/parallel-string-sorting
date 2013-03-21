# The name of our program:
PROGRAM= copy-burstsort

# These are our C source files:
SRCS= copy-burstsort.c C-burstsort.c CP-burstsort.c utils.c CPL-burstsort.c 

# The name of the C-compiler we want to use:
CC= gcc

# The flags (options) that the C-compiler needs:
#CFLAGS = -O3            # optimization
CFLAGS = -g -Wall -ansi  -pedantic    # debugging, warnings

# Names of libraries we need to link with our program:
LDLIBS= -lm

# Object files are derived from the sources, with .c replaced by .o
OBJECTS = ${SRCS:%.c=%.o}

# The first target, the name of our program:
${PROGRAM}: $(OBJECTS)
	$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

# To remove generated object files:
clean:
	-rm -f ${OBJECTS}
