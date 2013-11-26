/*
 * Simple C program to count the line lengths in an input stream.
 * 2013-04-13 Timo Bingmann <tb@panthema.net>
 *
 * Compile: gcc -O3 -W -Wall linecount.c -o linecount
 * 
 * Usage: linecount < input.txt | tee lfreq.txt
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

int main()
{
    typedef unsigned long long ull_type;

    int c;
    ull_type line_count = 0, char_count = 0, length_sum = 0;
    unsigned int i, char_used;
    ull_type ln, length;

    ull_type charmark[256];

    const unsigned int max_linemark = 1024*1024;

    // clear charmark
    memset(charmark, 0, sizeof(charmark));

    // allocate and clear linemark
    ull_type* linemark = malloc( (max_linemark+1) * sizeof(ull_type) );
    memset(linemark, 0, (max_linemark+1) * sizeof(ull_type));

    ln = 0; /* line start */
    while ( (c = getc(stdin)) != EOF )
    {
        if (c == '\n') {
            length = char_count - ln;
            length_sum += length;
            ++line_count;
            
            if (length < max_linemark)
                linemark[length]++;
            else
                linemark[max_linemark]++;

            ln = char_count+1; /* next line start */
        }
        ++char_count;
        ++charmark[c];
    }
    /* count last line */
    if (ln <= char_count) {
        length = char_count - ln;
        length_sum += length;
        ++line_count;
            
        if (length < max_linemark)
            linemark[length]++;
        else
            linemark[max_linemark]++;
    }

    /* count used characters */
    char_used = 0;
    for (i = 0; i < 256; ++i)
    {
        if (charmark[i] == 0) continue;
        ++char_used;
    }

    printf("Total: %llu lines in %llu bytes, alphabet %u, average line length %.6f.\n",
           line_count, char_count, char_used, (length_sum + line_count) / (double)line_count);

    printf("Excluding newline: %llu characters, alphabet %u, average line length %.6f.\n",
           length_sum, char_used - (charmark['\n'] ? 1 : 0), length_sum / (double)line_count);

    for (i = 0; i < 256; ++i)
    {
        if (charmark[i] == 0) continue;

        if (isprint(i)) {
            printf("char['%c'] = %llu\n", (char)i, charmark[i]);
        }
        else {
            printf("char[%u] = %llu\n", i, charmark[i]);
        }
    }

    printf("Line Length Distribution:\n");

    for (i = 0; i < max_linemark; ++i)
    {
        if (linemark[i] == 0) continue;
        printf("%u\t%llu\n", i, linemark[i]);
    }

    return 0;
}
