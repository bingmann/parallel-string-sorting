/*
 * Simple C program to count the character frequency in an input stream.
 * 2013-03-28 Timo Bingmann <tb@panthema.net>
 *
 * Compile: gcc -O3 -W -Wall charcount.c -o charcount
 * 
 * Usage: charcount < input.txt | tee freq.txt
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>

int main()
{
    int c;
    unsigned long long count = 0;
    unsigned int i, used = 0;

    unsigned long long mark[256];
    memset(mark, 0, sizeof(mark));

    while ( (c = getc(stdin)) != EOF )
        ++count, ++mark[c];

    for (i = 0; i < 256; ++i)
    {
        if (mark[i] == 0) continue;

        if (isprint(i)) {
            printf("char['%c'] = %llu\n", (char)i, mark[i]);
        }
        else {
            printf("char[%u] = %llu\n", i, mark[i]);
        }
        ++used;
    }

    printf("Total: %u chars used in %llu bytes.\n", used, count);

    return 0;
}
