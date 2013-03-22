/* This source code is from the following article:
 *
 * Ng, Waihong, and Katsuhiko Kakehi. "Merging String Sequences by Longest
 * Common Prefixes." IPSJ Digital Courier 4 (2008): 69–78.
 *
 * Appendix: Source Code of LCP Merge sort
 */

namespace ng_lcpmergesort {

/*
  Instead of just attaching the LLCPs to the strings to reduce memory usage,
  this code allocates n * ((size of key pointer) + (size of integer)) bytes of
  memory to store the annotated strings for clarity and simplicity.
*/

typedef unsigned int UINT;

typedef struct _AS {
    UINT llcp; char* str;
} AS, *PAS, *ASS; // annotated string sequence

void lcpm(PAS pX, PAS pXE, PAS pY, PAS pYE, ASS pZ)
{
    while ((pX<pXE) && (pY<pYE)) {
        if (*((UINT*)pX)!=*((UINT*)pY)) {
            if (*((UINT*)pX)>*((UINT*)pY)) {
                *pZ=*pX; pX++; pZ++;
            }
            else {
                *pZ=*pY; pY++; pZ++;
            }
        }
        else {
            char *x=pX->str+*((UINT*)pX);
            char *y=pY->str+*((UINT*)pX);
            while ((*x==*y) && (*x!=0)) { x++; y++; }
            if (*x>*y) {
                *(UINT*)pX=(x-(pX->str));
                *pZ=*pY; pY++; pZ++;
            }
            else {
                *(UINT*)pY=(y-(pY->str));
                *pZ=*pX; pX++; pZ++;
            }
        }
    }
    if (pX<pXE) {
        memcpy((char*)pZ, (char*)pX,
               (char*)pXE-(char*)pX);
    }
}

// for the last merge
void lcpm_f(PAS pX, PAS pXE, PAS pY, PAS pYE, char** pZ)
{
    while ((pX<pXE) && (pY<pYE)) {
        if (*(UINT*)pX!= *(UINT*)pY) {
            if (*(UINT*)pX > *(UINT*)pY) {

                *pZ=pX->str; pX++; pZ++;
            }
            else {
                *pZ=pY->str; pY++; pZ++;
            }
        }
        else {
            char *x=pX->str+*(UINT*)pX;
            char *y=pY->str+*(UINT*)pX;
            while ((*x==*y) && (*x!=0)) { x++; y++; }
            if (*x>*y) {
                *(UINT*)pX=x-pX->str;
                *pZ=pY->str; pY++; pZ++;
            }
            else {
                *(UINT*)pY=y-pY->str;
                *pZ=pX->str; pX++; pZ++;
            }
        }
    }
    while (pX<pXE) { *pZ=pX->str; pZ++; pX++; }
    while (pY<pYE) { *pZ=pY->str; pZ++; pY++; }
}

void lcpms_i(ASS pX, UINT n, char** str, ASS pZ)
{
    if (n>1) {
        UINT l=n/2;
        lcpms_i(pX+0, l-0, str+0, pZ);
        lcpms_i(pX+l, n-l, str+l, pZ);
        memcpy(pZ, pX, l*sizeof(AS));
        lcpm(pZ+0, pZ+l, pX+l, pX+n, pX);
    }
    else { pX->llcp=0; pX->str=*str; }
}

void lcpms(char** pStr, UINT n)
{
    ASS pX, pT;
    pX=pT=(ASS)malloc(sizeof(AS)*n);
    UINT l=n/2;
    lcpms_i(pX+0, l-0, pStr+0, (ASS)pStr);
    lcpms_i(pX+l, n-l, pStr+l, (ASS)pStr);
    lcpm_f(pX+0, pX+l, pX+l, pX+n, pStr);
    free(pT);
}

void ng_lcpmergesort(unsigned char **strings, size_t n)
{
    return lcpms((char**)strings, n);
}

CONTESTANT_REGISTER_UCARRAY(ng_lcpmergesort, "LCP-Mergesort Original by Waihong Ng and Katsuhiko Kakehi")


} // namespace ng_lcpmergesort
