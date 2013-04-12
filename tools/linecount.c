/*
 * Simple C program to count the line lengths in an input stream.
 * 2013-04-12 Timo Bingmann <tb@panthema.net>
 *
 * Compile: gcc -O3 -W -Wall linecount.c -o linecount
 * 
 * Usage: linecount < input.txt | tee lfreq.txt
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main()
{
    typedef unsigned long long ull_type;

    int c;
    ull_type count_lines = 0, count_chars = 0, length_sum = 0;
    unsigned int i;
    ull_type ln, length;

    const unsigned int maxmark = 1024*1024;

    ull_type* mark = malloc( (maxmark+1) * sizeof(ull_type) );
    memset(mark, 0, sizeof(mark));

    ln = 0; /* line start */
    while ( (c = getc(stdin)) != EOF )
    {
        if (c == '\n') {
            length = count_chars - ln;
            length_sum += length;
            ++count_lines;
            
            if (length < maxmark)
                mark[length]++;
            else
                mark[maxmark]++;

            ln = count_chars+1; /* next line start */
        }
        ++count_chars;
    }
    /* count last line */
    if (ln <= count_chars) {
        length = count_chars - ln;
        length_sum += length;
        ++count_lines;
            
        if (length < maxmark)
            mark[length]++;
        else
            mark[maxmark]++;
    }

    printf("Total: %llu lines in in %llu bytes, average length %.6f.\n",
           count_lines, count_chars, length_sum / (double)count_lines);

    for (i = 0; i < maxmark; ++i)
    {
        if (mark[i] == 0) continue;
        printf("%u\t%llu\n", i, mark[i]);
    }

    return 0;
}
