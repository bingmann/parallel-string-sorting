#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "zio.h"

int main(int argc, char *argv[])
{
    FILE *file;
    char line[1024];
    size_t len;

    if (!(file = fdzopen(fileno(stdin), "r", argc > 1 ? argv[1] : "g"))) {
	fprintf(stderr, "%s\n", strerror(errno));
	return 1;
    }

    while ((len = fread(line, sizeof(char), sizeof (line), file))) {
	size_t ret = fwrite(line, sizeof(char), len, stdout);
	if ((ret != len) && ferror(stdout)) {
	    clearerr(stdout);
	}
    }

    fclose(file);

    return 0;
}
