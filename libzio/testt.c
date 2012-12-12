#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "zio.h"

int main(int argc, char *argv[])
{
    while (argc > 1) {
	FILE *file;
	char line[1024];
	size_t len;

	argv++;
	argc--;

	if (!(file = fzopen(*argv, "r"))) {
	    fprintf(stderr, "%s\n", strerror(errno));
	    continue;
	}

	while ((len = fread(line, sizeof(char), sizeof (line), file))) {
	    size_t ret = fwrite(line, sizeof(char), len, stdout);
	    if ((ret != len) && ferror(stdout)) {
		clearerr(stdout);
	    }
	}

	fclose(file);
    }

    return 0;
}
