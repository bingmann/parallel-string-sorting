#ifndef EBERLE_LCP_MERGESORT_H_
#define EBERLE_LCP_MERGESORT_H_

#include <iostream>

#include "../utils/types.h"

#include "../../tools/stringtools.h"
#include "../../tools/statsfile.h"

namespace eberle_mergesort_lcp
{

using namespace types;

using namespace stringtools;

static inline
void
eberle_lcp_merge(string* input1, unsigned* lcps1, size_t length1, string* input2, unsigned* lcps2, size_t length2, string* output,
        unsigned* outputLcps)
{
    const string* end1 = input1 + length1;
    const string* end2 = input2 + length2;

    unsigned lcp1 = *lcps1;
    unsigned lcp2 = *lcps2;

    //do the merge
    while (input1 < end1 && input2 < end2)
    {
        if (lcp1 == lcp2)
        { // CASE 1 lcps are equal => do string comparision starting at lcp+1st position
            string s1 = *input1 + lcp1;
            string s2 = *input2 + lcp1;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != '\0' && *s1 == *s2)
                s1++, s2++;

            const unsigned lcp = s1 - *input1;

            if (*s1 <= *s2)
            { 	// CASE 1.1: lcp1 <= lcp2
                *output = *input1;
                *outputLcps = lcp1;
                ++input1;
                ++lcps1;
                lcp1 = *lcps1;
                lcp2 = lcp;
            }
            else
            { 	// CASE 1.2: lcp1 > lcp2
                *output = *input2;
                *outputLcps = lcp2;
                ++input2;
                ++lcps2;
                lcp1 = lcp;
                lcp2 = *lcps2;
            }
        }

        else if (lcp1 < lcp2)
        {   // CASE 2: lcp1 > lcp2
            *output = *input2;
            *outputLcps = lcp2;
            ++input2;
            ++lcps2;
            lcp2 = *lcps2;
        }

        else
        {   // CASE 3: lcp1 < lcp2
            *output = *input1;
            *outputLcps = lcp1;
            ++input1;
            ++lcps1;
            lcp1 = *lcps1;
        }

        ++output;
        ++outputLcps;
    }

    if (input1 < end1)
    {   // if there are remaining elements in stream1, copy them to the end
        memcpy(output, input1, (end1 - input1) * sizeof(string));
        memcpy(outputLcps, lcps1, (end1 - input1) * sizeof(unsigned));
        *outputLcps = lcp1;
    }
    else
    {
        memcpy(output, input2, (end2 - input2) * sizeof(string));
        memcpy(outputLcps, lcps2, (end2 - input2) * sizeof(unsigned));
        *outputLcps = lcp2;
    }
}

static inline
void
eberle_lcp_mergesort(string *strings, string* tmp, unsigned* tmpLcps, string* output, unsigned* outputLcps, size_t length)
{

    if (length <= 1)
    {
        *output = *strings;
        *outputLcps = 0;
        return;
    }

    const size_t length1 = length / 2;
    const size_t length2 = length - length1;

    eberle_lcp_mergesort(strings, output, outputLcps, tmp, tmpLcps, length1);
    eberle_lcp_mergesort(strings + length1, output + length1, outputLcps + length1, tmp + length1, tmpLcps + length1, length2);

    eberle_lcp_merge(tmp, tmpLcps, length1, tmp + length1, tmpLcps + length1, length2, output, outputLcps);
}

void
eberle_lcp_mergesort_seperate(string *strings, size_t n)
{
    // Allocate memory for LCPs and temporary string array

    unsigned* outputLcps = new unsigned[n];
    string* tmpStrings = new string[n];
    unsigned* tmpLcps = new unsigned[n];

    eberle_lcp_mergesort(strings, tmpStrings, tmpLcps, strings, outputLcps, n);

    delete outputLcps;
    delete tmpStrings;
    delete tmpLcps;
}

CONTESTANT_REGISTER(eberle_lcp_mergesort_seperate, "eberle/mergesort_lcp_binary_seperate", "Binary Mergesort with LCP-usage by Andreas Eberle")

// implementation follows
static inline
void
eberle_lcp_merge(const LcpStringPtr& input1, size_t length1, const LcpStringPtr& input2, size_t length2, const LcpStringPtr& output)
{
    eberle_lcp_merge(input1.strings, input1.lcps, length1, input2.strings, input2.lcps, length2, output.strings, output.lcps);
}

static inline
void
eberle_lcp_mergesort(string *strings, const LcpStringPtr& tmp, const LcpStringPtr& output, size_t length)
{
    if (length <= 1)
    {
        output.set(*strings, 0);
        return;
    }

    const size_t length1 = length / 2;
    const size_t length2 = length - length1;
    const LcpStringPtr tmp2 = tmp + length1;

    eberle_lcp_mergesort(strings, output, tmp, length1);
    eberle_lcp_mergesort(strings + length1, output + length1, tmp2, length2);

    eberle_lcp_merge(tmp, length1, tmp2, length2, output);
}

void
eberle_lcp_mergesort_ptr(string *strings, size_t n)
{
    // Allocate memory for LCPs and temporary string array

    unsigned* outputLcps = new unsigned[n];
    string* tmpStrings = new string[n];
    unsigned* tmpLcps = new unsigned[n];

    LcpStringPtr output(strings, outputLcps);
    LcpStringPtr tmp(tmpStrings, tmpLcps);

    eberle_lcp_mergesort(strings, tmp, output, n);

    delete outputLcps;
    delete tmpStrings;
    delete tmpLcps;
}

CONTESTANT_REGISTER(eberle_lcp_mergesort_ptr, "eberle/mergesort_lcp_binary_ptr", "Binary Mergesort with LCP-usage by Andreas Eberle")

static inline
void eberle_lcp_merge(AS* input1, size_t length1, AS* input2, size_t length2,
        AS* output)
{
    const AS* end1 = input1 + length1;
    const AS* end2 = input2 + length2;
    AS* outStream = output;

    //do the merge
    while (input1 < end1 && input2 < end2)
    {
        AS* a = input1;
        AS* b = input2;

        if (a->lcp == b->lcp)
        { // CASE 1 lcps are equal
            string s1 = a->text + a->lcp;
            string s2 = b->text + a->lcp;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != '\0' && *s1 == *s2)
            s1++, s2++;

            const UINT lcp = s1 - a->text;

            if (*s1 <= *s2)
            { 	// CASE 1.1: a <= b
                *outStream = *a;
                input1++;
                b->lcp = lcp;

            }
            else
            { 			// CASE 1.2: a > b
                *outStream = *b;
                input2++;
                a->lcp = lcp;
            }

        }
        else if (a->lcp < b->lcp)
        { // CASE 2: a > b
            *outStream = *b;
            input2++;

        }
        else
        { // CASE 3: a < b
            *outStream = *a;
            input1++;
        }

        outStream++;
    }

    if (input1 < end1)
    { // if there are remaining elements in stream1, copy them to the end
        memcpy(outStream, input1, (end1 - input1) * sizeof(AS));
    }
    else
    {
        memcpy(outStream, input2, (end2 - input2) * sizeof(AS));
    }
}

static inline
void
eberle_lcp_mergesort(string *strings, AS *tmp, AS *output, size_t length)
{

    if (length <= 1)
    {
        output->lcp = 0;
        output->text = *strings;
        return;
    }

    const size_t length1 = length / 2;
    const size_t length2 = length - length1;
    AS * tmp2 = tmp + length1;

    eberle_lcp_mergesort(strings, output, tmp, length1);
    eberle_lcp_mergesort(strings + length1, output + length1, tmp2, length2);

    eberle_lcp_merge(tmp, length1, tmp2, length2, output);
}

void
eberle_lcp_mergesort(string *strings, size_t n)
{
    //allocate memory for annotated strings
    AS *tmp = static_cast<AS *>(malloc(n * sizeof(AS)));
    AS *output = static_cast<AS *>(malloc(n * sizeof(AS)));

    eberle_lcp_mergesort(strings, tmp, output, n);

    for (size_t i = 0; i < n; i++)
    {
        strings[i] = output[i].text;
    }

    free(tmp);
    free(output);
}

CONTESTANT_REGISTER(eberle_lcp_mergesort, "eberle/mergesort_lcp_binary", "Binary Mergesort with LCP-usage by Andreas Eberle")

static inline
void
eberle_lcp_merge(AS* input1, size_t length1, AS* input2, size_t length2, string* output)
{
    const AS* end1 = input1 + length1;
    const AS* end2 = input2 + length2;

    //do the merge
    while (input1 < end1 && input2 < end2)
    {
        AS* a = input1;
        AS* b = input2;

        if (a->lcp == b->lcp)
        { // CASE 1 lcps are equal
            string s1 = a->text + a->lcp;
            string s2 = b->text + a->lcp;

            // check the strings starting after lcp and calculate new lcp
            while (*s1 != '\0' && *s1 == *s2)
            s1++, s2++;

            const UINT lcp = s1 - a->text;

            if (*s1 <= *s2)
            {   // CASE 1.1: a <= b
                *output = a->text;
                input1++;
                b->lcp = lcp;

            }
            else
            {           // CASE 1.2: a > b
                *output = b->text;
                input2++;
                a->lcp = lcp;
            }

        }
        else if (a->lcp < b->lcp)
        { // CASE 2: a > b
            *output = b->text;
            input2++;

        }
        else
        { // CASE 3: a < b
            *output = a->text;
            input1++;
        }

        output++;
    }

    if (input1 < end1)
    { // if there are remaining elements in stream1, copy them to the end
      //memcpy(output, input1, (end1 - input1) * sizeof(AS));
        for (; input1 < end1; input1++, output++)
        {
            *output = input1->text;
        }
    }
    else
    {
        //memcpy(output, input2, (end2 - input2) * sizeof(AS));
        for (; input2 < end2; input2++, output++)
        {
            *output = input2->text;
        }
    }
}

}
// namespace eberle_lcp_mergesort

#endif // EBERLE_LCP_MERGESORT_H_

