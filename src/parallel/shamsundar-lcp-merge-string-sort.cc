/* Recursive Posix Thread OR Windows Threads version of
// Longest Common Prefix Length based Merge sort, as described in:

http://dspace.wul.waseda.ac.jp/dspace/bitstream/2065/28672/4/Honbun-4624_01.pdf

// A special feature of the merge: the first part is put into a temporary
// array. The output is into the full array. Elements in the second part
// will be moved down as needed. Any remnants in the second part are already
// in place.

// Written by N. Shamsundar, University of Houston, May 2009.

// To build,

   (i) Configure by setting values in the source code for

          NCPU (number of cores/CPUs),
          SIZ  (maximum size of buffer to hold input text),
          NLIN (maximum number of lines in input text)
          fcol (first column of key field, column numbers starting with 1)
          lcol (last  column of key field)

   (ii-a) Build using the provided assembler routines (for speed)

                                       OR

   (ii-b) Comment out the prototype declaration of 'merge' and two calls to 'merge'.
           Uncomment the two calls to 'Merge'.

    ml /c /coff mrgc.asm
    cl /DPTHREAD lcpsrtmt.c mrgc.obj pthreadVC2.lib (for Posix Threads in Windows)
    cl /MT lcpsrtmt.c mrgc.obj (for Windows Threads)

    gcc -pthread -O lcpsrtmt.c mrgc.s -lpthread -s -o lcpsrtmt (Cygwin, Linux, Unix)

    (iii) Run

    lcpsrtmt <input file> <output file>

    * Added glue to attach Shamsundar's Parallel LCP Merge String Sort code to
    * parallel-string-sort framework in 2013 Timo Bingmann <tb@panthema.net>.
*/

#include "../tools/contest.h"

extern int pss_num_threads;

#define PTHREAD

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#ifdef _WIN32
#ifndef PTHREAD
#include <windows.h>
#endif
#include <io.h>
#else
#include <unistd.h>
#endif
#include <sys/types.h>
#include <fcntl.h>

namespace shamsundar_lcp_merge_string_sort {

/* Written by Shamsundar as a poor man's replacement for
   the vararglist function error(status,errnum,fmt,...) */
int error_message_count=0;
#define ERR0(status,errnum,fmt) \
   {fprintf(stderr,"Program %s, File %s, Line %d: ", \
    argv[0],__FILE__,__LINE__); \
   fprintf(stderr,fmt); \
   if(errnum)fprintf(stderr,": Err.Num. %d\n",errnum); \
   if(status)exit(status); \
   error_message_count++; \
   }
#define ERR1(status,errnum,fmt,p1) \
   {fprintf(stderr,"Program %s, File %s, Line %d: ", \
    argv[0],__FILE__,__LINE__); \
   fprintf(stderr,fmt,p1); \
   if(errnum)fprintf(stderr,": Err.Num. %d\n",errnum); \
   if(status)exit(status); \
   error_message_count++; \
   }

#define FT '\0'
//#define NCPU 2
#define MAXCPU 65
int NCPU = 1;
#define MIN(x,y) (x) < (y) ? (x) : (y)

typedef unsigned int UINT;
typedef struct _AS {
    UINT llcp; char* str;
} AS, *pAS, *aAS; // annotated string sequence

pAS *gpTMP,*gSP;

// The assembler version uses i686 instructions, and cannot be used on modern
// x86_64 bit architectures. -tb

extern void  merge( pAS a[], int la, pAS B[], int lb, pAS C[], int lcol, int *ncomp ); //ASM version

int DEBUG=0;

void

#ifdef _WIN32
   _inline
#else
   inline
#endif

Merge( aAS a[], int m, aAS b[], int n, aAS c[], int lcol, int *ncomp )
{
/* Merge a[0..m-1] and b[0..n-1] into c[0..m+n-1] */

    int i,j,k;
    if(DEBUG){
        printf("List A\n"); for(i=0; i<m; i++)printf("%3d %2d %s\n",i,a[i]->llcp,a[i]->str);
        printf("\nList B\n"); for(i=0; i<n; i++)printf("%3d %2d %s\n",i,b[i]->llcp,b[i]->str);
    }

    if(c+m-b){fprintf(stderr,"Improper overlap\n"); exit(1);}
    i = 0; j = 0; k = 0;
    while (i < m && j < n){
        if(a[i]->llcp - b[j]->llcp){
            if(a[i]->llcp > b[j]->llcp)c[k++]=a[i++]; else c[k++]=b[j++];
        }
        else{
            char *x=a[i]->str+a[i]->llcp,
                *y=b[j]->str+a[i]->llcp,
                *xe=a[i]->str+lcol-1;
            (*ncomp)++;

            while ((*x==*y) && (x < xe) && (*x-FT)) {
                x++; y++;
            }
            if (*x > *y) {
                a[i]->llcp=x-a[i]->str;        // update prefix
                c[k++]=b[j++];
            }
            else {
                b[j]->llcp=y-b[j]->str;        // update prefix
                c[k++]=a[i++];
            }
        }
    }
    while ( i < m ) c[k++] = a[i++];           // output remnants in A
    // Any remnants in b are already in place!
    if(DEBUG){
        printf("\nList C\n"); for(i=0; i<m+n; i++)printf("%3d %2d %s\n",i,c[i]->llcp,c[i]->str);
        printf("\n");
    }
}

void merge_sort( pAS pX[], int n, int cnts[], pAS *lPTMP, int fcol, int lcol )       // Recursive!
{
    UINT i;
// cnts[0..3] contain updated counts of comparisons, merges, copies and skips
    if (n>2) {
        UINT l=n/2;
        char *sX,*sY,*p,*q,*XE;
        if(l > 1)merge_sort(pX+0, l-0,cnts, lPTMP, fcol, lcol);   /* sort left half into lower half of pX */
        if(n-l > 1)merge_sort(pX+l, n-l,cnts, lPTMP, fcol, lcol); /* sort right half into upper half of pX */
        sX=pX[l-1]->str; sY=pX[l]->str; cnts[0]++;
        XE=pX[l-1]->str+lcol;
        p=sX+fcol-1; q=sY+fcol-1;
        while(*p-FT && *p==*q && p < XE){p++; q++;}
        if(*p==FT || *p <= *q){
            pX[l]->llcp=p-sX;      // update prefix
            return;                          /* no need to merge if sorted */
        }
        else{
            if(lPTMP-pX){
                memcpy(lPTMP, pX, l*sizeof(pAS));  /* copy lower half of pX into pTMP */
                cnts[2]+=l;
            }
            cnts[1]++;
            if(DEBUG){
                printf("List A\n"); for(i=0; i<l; i++)printf("%3d %2d %s\n",i,pX[i]->llcp,pX[i]->str);
                printf("\nList B\n"); for(i=l; i<n; i++)printf("%3d %2d %s\n",i,pX[i]->llcp,pX[i]->str);
            }
            Merge(lPTMP+0, l, pX+l, n-l, pX, lcol, cnts);       /* merge pX_lower and pX_upper, putting the result into pX */
//            merge(lPTMP+0, l, pX+l, n-l, pX, lcol, cnts);         // Assembler version
            if(DEBUG){
                printf("\nList C\n"); for(i=0; i<n; i++)printf("%3d %2d %s\n",i,pX[i]->llcp,pX[i]->str);
                printf("\n");
            }
        }
        return;
    }
    else{                                    // n = 2; n = 1 should not happen
        char *sX,*sY,*p,*q,*XE;
        if(n==1)return; cnts[3]++;
        sX=pX[0]->str; sY=pX[1]->str;    // merging not needed, but update LLCP
        XE=sX+lcol-1;
        cnts[0]++;
        p=sX+fcol-1; q=sY+fcol-1;
        while(*p-FT && *p==*q){
            if(p==XE)break;
            p++; q++;
        }
        if(*p > *q){
            pX[0]->str=sY; pX[1]->str=sX;
        }
        if(p-sX)pX[1]->llcp=p-sX;
    }
}

typedef struct {
   int cnts[4];
   pAS *pX;
   int n;
   pAS *lPTMP;
   int fcol; int lcol;
   } msortargs;
msortargs margs[MAXCPU];

// Wrapper to fold merge_sort into the pThread worker function protocol
#ifdef PTHREAD
void *worker( void *warg){
#else
DWORD WINAPI worker( LPVOID warg){
#endif
    msortargs *pStrct=(msortargs *)warg;
    pAS *pX=pStrct->pX; int n=pStrct->n; pAS *lpTMP=pStrct->lPTMP;
    int *lcnts=pStrct->cnts; int fcol=pStrct->fcol; int lcol=pStrct->lcol;
    merge_sort(pX,n,lcnts,lpTMP,fcol,lcol);
#ifdef PTHREAD
    return NULL;
#else
    return 0;
#endif
}

void merge_sortMT( pAS pX[], int n, int cnts[], pAS *lPTMP, int fcol, int lcol ) // NOT Recursive!
{
    int i,blksiz,blen;
#ifdef PTHREAD
    pthread_t thrd[MAXCPU];
#else
    DWORD   dwThreadIdArray[MAXCPU]; HANDLE  hThreadArray[MAXCPU];
#endif
    /* Divide the sort array into NCPU blocks. Hand off blocks 0,1,...NCPU-1 to other threads */
    blksiz=n/(2.0*NCPU)+0.5; blksiz+=blksiz;
    if(blksiz%2){
        fprintf(stderr,"Blksiz is not even with n=%d\n",n);
        exit(1);
    }
    for(i=0; i<NCPU-1; i++){
        margs[i].pX=pX+i*blksiz; margs[i].n=blksiz;
        margs[i].cnts[0]=margs[i].cnts[1]=margs[i].cnts[2]=margs[i].cnts[3]=0;
        margs[i].lPTMP=lPTMP+i*blksiz/2; margs[i].fcol=fcol; margs[i].lcol=lcol;
    }
    for(i=NCPU-2; i>=0; i--){
#ifdef PTHREAD
        if(pthread_create(thrd+i,NULL,worker,(void *)(margs+i))){
            fprintf(stderr,"Error creating thread %d\n",i);
            exit(1);
        }
#else
        hThreadArray[i] =
            CreateThread(
                NULL,                   // default security attributes
                0,                      // use default stack size
                worker,       // thread function name
                &margs[i],          // argument to thread function
                0,                      // use default creation flags
                &dwThreadIdArray[i]);   // returns the thread identifier
#endif
        if(DEBUG)
            printf("Worker %d, %ld to %ld\n",i,margs[i].pX-gSP,margs[i].pX-gSP+blksiz-1);
        // This line only when PThreads are not used :
        // worker(margs+i);
    }
    // The threads have all been started. Sort the high-end block in the manager thread
    // Even the boss does his bit.
    //
    //
    cnts[0]=cnts[1]=cnts[2]=cnts[3]=0; blen=MIN(n-(NCPU-1)*blksiz,blksiz);
    if(DEBUG)printf("Boss 0, %ld to %ld\n",pX+(NCPU-1)*blksiz-gSP,pX-gSP+(NCPU-1)*blksiz+blen-1);
    merge_sort(pX+(NCPU-1)*blksiz,blen,cnts,lPTMP+(NCPU-1)*blksiz/2,fcol,lcol);

    // Wait for the threads to complete; as each completes, merge its results with the part that
    // has already been sorted.

    for(i=NCPU-2; i>=0; i--){

#ifdef PTHREAD
        if(pthread_join(thrd[i],NULL)){
            fprintf(stderr,"Error joining thread %d\n",i);
            exit(1);
        }
#else
        WaitForSingleObject(hThreadArray[i],INFINITE);
        CloseHandle(hThreadArray[i]);
#endif
        cnts[0]+=margs[i].cnts[0]; cnts[1]+=margs[i].cnts[1];
        cnts[2]+=margs[i].cnts[2]; cnts[3]+=margs[i].cnts[3];

        if(DEBUG)printf("Merge    %ld to %ld\n",margs[i].pX-gSP,margs[i].pX-gSP+blksiz-1);
        cnts[0]++;
        if(strcmp(pX[(i+1)*blksiz-1]->str,pX[(i+1)*blksiz]->str) > 0){  // merge is needed
            memcpy(lPTMP+i*blksiz/2,margs[i].pX,blksiz*sizeof(pAS)); cnts[2]+=blksiz;
            Merge(lPTMP+i*blksiz/2, blksiz, pX+(i+1)*blksiz, blen, pX+i*blksiz, lcol, cnts);
//         merge(lPTMP+i*blksiz/2, blksiz, pX+(i+1)*blksiz, blen, pX+i*blksiz, lcol, cnts);         // Assembler version
            cnts[1]++;
        }
        else cnts[3]++;  /* Skip count */
        blen+=blksiz;
    }
    return;
}

#if 0 // -tb
char *bint2str(long l, char *s){
char *p;
int m,t,u;
u=l % 1000; l/=1000; t=l % 1000; l/=1000; m=l % 1000;
sprintf(s,"%03d,%03d,%03d",m,t,u);
p=s; while(*p=='0' || *p==',')*p++=' ';
if(*p=='\0')*(--p)='0';
return s;
}

#define SIZ  0x4000000
#define NLIN 0x400000

char buf[SIZ];
AS slin[NLIN], *sP[NLIN];
int brks[NLIN];

int main(int argc,char *argv[]){
    int fd,blen,nlin=0,i; char *p=buf; FILE *fil;
    int c1,c2,c3,c4,cnts[4]; char tstrs[5][13]; int fcol=17,lcol=24;
    c1=clock(); gSP=sP;
    if(argc-3)ERR0(1,0,"usage: aflgsrt <ifil> <ofil>")
                  fd=open(argv[1],O_RDONLY);
    if(fd < 0)ERR1(1,0,"Could not open %s for reading",argv[1])
                  blen=read(fd,buf,SIZ);
    close(fd);
    if(blen==SIZ)ERR0(1,0,"buffer not big enough");
/* Read file, replace EOL by null */
    while(p-buf < blen){
        sP[nlin]=slin+nlin;
        slin[nlin].llcp=fcol-1; slin[nlin++].str=p; p=strchr(p,'\n');
        if(p[-1]=='\r')p[-1]='\0'; *(p++)='\0';
        if(nlin==NLIN)ERR0(1,0,"NLIN not big enough")
    }

    gpTMP=(pAS *)calloc((nlin+11)/2,sizeof(pAS));
    cnts[0]=cnts[1]=cnts[2]=cnts[3]=0;
    c2=clock();
    merge_sortMT(sP,nlin,cnts,gpTMP,fcol,lcol);
    c3=clock();
    free(gpTMP);

    if(cnts[1]){
/* Write output file */
        fil=fopen(argv[2],"wb");
        if(fil==NULL)ERR1(1,0,"Could not open %s for writing",argv[2])
        for(i=0; i<nlin; i++){
            fputs(sP[i]->str,fil); fputc('\n',fil);
        }
        fclose(fil);
        c4=clock();
    }
    else{
        ERR0(0,0,"Input file is already sorted\n");
        c4=c3;
    }
    ERR1(0,0,"Read input           time : %6.3f\n",(double)(c2-c1)/CLOCKS_PER_SEC);
    ERR1(0,0,"MergeSort            time : %6.3f\n",(double)(c3-c2)/CLOCKS_PER_SEC);
    if(cnts[1])ERR1(0,0,"Write output         time : %6.3f\n",(double)(c4-c3)/CLOCKS_PER_SEC);
    ERR1(0,0,"Total task           time : %6.3f\n",(double)(c4-c1)/CLOCKS_PER_SEC);
    printf("\n%s records, %s merges, %s skipped, %s comparisons, %s copies\n",
           bint2str(nlin,tstrs[0]),
           bint2str(cnts[1],tstrs[1]),
           bint2str(cnts[3],tstrs[2]),
           bint2str(cnts[0],tstrs[3]),
           bint2str(cnts[2],tstrs[4]));
    return 0;
}
#endif

// Glue to attach to string sorting test framework

AS *sP;

static void prepare(unsigned char **strings, size_t size)
{
    sP = new AS[size];
    for (size_t i = 0; i < size; ++i) {
        sP[i].str = (char*)strings[i];
        sP[i].llcp = 0;
    }

    gSP = new AS*[size];
    for (size_t i = 0; i < size; ++i)
        gSP[i] = &sP[i];
}

static void string_sort(unsigned char **strings, size_t size)
{
    int cnts[4];
    int fcol=1,lcol=1000000000;

    NCPU = pss_num_threads;

    cnts[0]=cnts[1]=cnts[2]=cnts[3]=0;

    gpTMP = (pAS *)calloc((size+11)/2, sizeof(pAS)); 

    merge_sortMT(gSP,size,cnts,gpTMP,fcol,lcol);

    for (size_t i = 0; i < size; ++i) {
        strings[i] = (unsigned char*)gSP[i]->str;
    }

    free(gpTMP);
    delete [] sP;
    delete [] gSP;
}

CONTESTANT_REGISTER_PARALLEL_PREPARE(
    prepare, string_sort,
    "shamsundar/lcp-merge-string-sort",
    "Parallelized LCP Merge sort by N. Shamsundar");

} // namespace shamsundar_lcp_merge_string_sort
