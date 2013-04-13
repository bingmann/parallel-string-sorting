/*
 * Simple C program to calculate lcp and distingishing prefix of a _sorted_
 * input stream.
 * 2013-04-13 Timo Bingmann <tb@panthema.net>
 *
 * Compile: gcc -O3 -W -Wall lcp-dprefix.c -o lcp-dprefix
 * 
 * Usage: lcp-dprefix < input.txt
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main()
{
    typedef unsigned long long ull_type;

    ull_type line_count = 0, char_count = 0;
    ull_type depth = 0, depth_prev = 0;

    ull_type lcptotal = 0, dprefix = 0;

    const int maxline = 1024*1024;

    char* line = (char*)malloc(maxline);
    char* line_prev = (char*)malloc(maxline);
    char* s;

    while ( fgets(line, maxline, stdin) )
    {
        /* trim off newline if exists */
        if (line[ strlen(line)-1 ] == '\n') {
            line[ strlen(line)-1 ] = 0;
        }

        /* check ascending order of strings */
        if (line_count != 0 &&
            strcmp(line, line_prev) < 0)
        {
            printf("Input violates sorted order in lines %llu and %llu\n",
                   line_count, line_count+1);
            return 0;
        }

        /* calculate depth == LCP of prev and this */
        depth = 0;
        while ( line_prev[depth] == line[depth] &&
                line_prev[depth] != 0 ) ++depth;

        lcptotal += depth;

        /* distingishing length of prev and this */
        ++depth;

        if (depth_prev < depth) {
            dprefix += depth - depth_prev; /* add extra distinguishing characters */
        }
        dprefix += depth;
        depth_prev = depth;

        /* swap pointers */
        s = line_prev, line_prev = line, line = s;

        line_count++;
        char_count += strlen(line); /* excludes newlines! */
    }

    printf("line_count = %llu, char_count = %llu, average LCP %.3f\n",
           line_count, char_count, lcptotal / (double)line_count);
    
    printf("distingishing prefix = %llu, percentage of char_count = %.3f\n",
           dprefix, dprefix * 100.0 / (double)char_count);

    return 0;
}
