# README to parallel-string-sorting

This source code archive contains our test framework for parallel and
sequential string sorting algorithms, as described in our technical report and
paper "Engineering Parallel Sorting". See the project web page
<http://panthema.net/2013/parallel-string-sorting/> for more information.

## Source Code Overview

The main algorithm implementations are in the directory "`src/*/`", with only
"`src/tools/`" containing auxiliary headers.

`src/sequential/` is a collection of old and new sequential string sorting
algorithms, including the original multikey-quicksort from Bentley and
Sedgewick, various older radix sorts from McIlroy, Bostic, and McIlroy, as well
as Stefan Nilsson. The CRadix and LCPMergesort algorithms by Waihong Ng and
Ranjan Sinha's original array and list burstsorts are also included.

In the `src/sequential/` directory we also included sequential versions of
super scalar string sample sort (S^5), which were developed to test the
algorithms before parallelizing them.

Tommi Rantala's huge library of original string sorting implementation is
located in `src/rantala/`, including extra tools.

Ranjan Sinha's copy-burstsort variants are located in
`src/sinha-copy-burstsort`. These are heavily modified and connected to the
framework via the glue.cc source file.

All parallel string sorting implementations, including our fastest one,
**pS^5**, is located in `src/parallel/`. The new algorithm with highest speedups
"pS^5-Unroll" is located in `src/parallel/bingmann-parallel_sample_sort.cc`. For
NUMA systems we developed multiway LCP-mergesort, which is located in
`src/parallel/eberle-ps5-parallel-toplevel-merge.h` and uses the LCP-losertree
in `src/tools/eberle-lcp-losertree.h`.

Only one binary program `psstest` is built when compiling with default options,
which contains all implementations in the collection. The main source code file
`psstest.cc` includes most simple algorithms via header files. Further
implementations are separated either via C++ namespaces or individual
compilation units (`.o`s).

The framework combines all algorithms into a "Contest" using the classes in
`src/tools/contest.h`. An implementation registers a sorting function via
macros. All registered implementation can be selected using the "`-a`" command
line switch of psstest. Psstest is used to run speed tests on the various input
instances and automatically checks the calculated result.

## Using psstest

Psstest is the main program used to run, check and measure performance of the
string sorting algorithms. Running the program without arguments will print a
command line help.

The algorithms implementation are selected using the "`-a`" switch. Run `-a
list` for a (long) list of string sorting algorithms. For example `-a mkqs` will
select all algorithms *containing* "mkqs". You can use `-a` multiple times to
select different sets (multiple `-a` are logically inclusive). To select an
algorithm name fully matching a string, use the "`-A`" switch.

The input is selected by giving a file name or artificial random source
name. Available random inputs are "`randomASCII`", "`random10`", "`random4`"
and "`random255`", where the number specifies the alphabet size and ASCII is
described in our paper.

The program will automatically decompress files ending in "`.gz`", "`.bz2`",
"`.xz`" and "`.lzo`" by spawning the appropriate decompressors as a child
program. However, the *decompressed file size* must be added to the file name
for this to be efficient: you must rename a compressed input file like
`test.gz` to `test.12345.gz`, where 12345 is the decompressed file size!

One can *limit* the input size (number of bytes) using `-s <size>`, where size
can be expressed with suffixes like "512mb".

The input strings can be saved using `--input` and the sorted output is saved
with `--output` (otherwise it is discarded).

For parallel string sorters, the number of threads is specified by one of the
arguments: `--threads`, `--some-threads`, `--all-threads` or `--thread-list
<#>`. Check the command line help for details.

To isolate multiple algorithms runs, use the `--fork` argument, which will load
data only once, or the `--datafork` parameter to reload the input in each
forked child program. These can be combined with `--timeout` to abort a child
program after the specified time.

## Building psstest

To build the program a recent gcc C++ compiler, "cmake" version 2.8 or higher,
and the GNU GMP library are required.

Having "parallel-string-sorting-X.Y.Z.tar" unpacked in the current directory,
the following commands will compile and run the psstest program:

```
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
$ cd src/
$ ./psstest
```

If the cmake configuration routine does not finish successfully, you may need
to install some development packages via your Linux distribution's package
manager.

Here is a simple example to sort random number strings (10 MiB total chars)
using pS^5-Unroll, and write them unsorted and sorted to files:

```
$ ./psstest -a parallel_sample_sortBTCU2 -s 10mb -i unsorted.txt -o sorted.txt random10
```

## Further Information

We have collected source code, papers, some test instances, and further
information at

<http://panthema.net/2013/parallel-string-sorting/>

## Exits

Written 2013-05-02 by Timo Bingmann, updated 2014-03-09

## Appendix: List of String Sorting Implementations

```
Running parallel-string-sorting 0a2c4fa89c1c2148dfd6a0903bd225141a158c5b on i10host
Called as ./psstest -a list
Available string sorting algorithms:
akiba/parallel_radix_sort                               Parallel MSD radix sort by Takuya Akiba
bingmann/lcp_insertion_sort                             LCP-aware insertion sort
bingmann/lcp_insertion_sort_cache                       LCP-aware insertion sort (with distinguishing character cache)
bingmann/lcp_insertion_sort_nolcp                       LCP-aware insertion sort (without LCP output)
bingmann/lcp_insertion_sort_pseudocode                  LCP-aware insertion sort close to pseudo-code, with checking
bingmann/lcp_mergesort_128way                           128-way LCP-Mergesort by Andreas Eberle and Timo Bingmann
bingmann/lcp_mergesort_16way                            16-way LCP-Mergesort by Andreas Eberle and Timo Bingmann
bingmann/lcp_mergesort_4way                             4-way LCP-Mergesort by Andreas Eberle and Timo Bingmann
bingmann/lcp_mergesort_8way                             8-way LCP-Mergesort by Andreas Eberle and Timo Bingmann
bingmann/lcp_mergesort_binary                           Binary Mergesort with LCP-merge by Andreas Eberle and Timo Bingmann
bingmann/msd_CE                                         bingmann/msd_CE (rantala CE original)
bingmann/msd_CE2                                        bingmann/msd_CE2 (CE with reused prefix sum)
bingmann/msd_CE3                                        bingmann/msd_CE3 (CE2 with iterators)
bingmann/msd_CE_nr                                      bingmann/msd_CE_nr (CE non-recursive)
bingmann/msd_CE_nr2                                     bingmann/msd_CE_nr2 (CE non-recursive)
bingmann/msd_CI                                         bingmann/msd_CI (rantala CI original with oracle)
bingmann/msd_CI2                                        bingmann/msd_CI2 (CI without oracle)
bingmann/msd_CI3                                        bingmann/msd_CI3 (CI2 with swap operations)
bingmann/msd_CI4                                        bingmann/msd_CI4 (CI3 with swap cache)
bingmann/msd_CI5                                        bingmann/msd_CI5 (CI4 with charcache)
bingmann/msd_CI5_16bit                                  bingmann/msd_CI5_16bit (CI5 with 16-bit radix)
bingmann/msd_CI_nr                                      bingmann/msd_CI_nr (CI non-recursive)
bingmann/msd_CI_nr2                                     bingmann/msd_CI_nr2 (CI non-recursive)
bingmann/msd_CI_nr3                                     bingmann/msd_CI_nr3 (CI non-recursive, charcache)
bingmann/parallel_mkqs                                  Parallel MKQS with blocks and cache8
bingmann/parallel_radix_sort_16bit                      Parallel MSD Radix sort with load balancing, 16-bit BigSorts
bingmann/parallel_radix_sort_8bit                       Parallel MSD Radix sort with load balancing, 8-bit BigSorts
bingmann/parallel_sample_sortBSC                        pS5: binary search, bktcache
bingmann/parallel_sample_sortBSC_lcp                    pS5: binary search, bktcache_lcp
bingmann/parallel_sample_sortBTC                        pS5: binary tree, bktcache
bingmann/parallel_sample_sortBTCE                       pS5: binary tree, equality, bktcache
bingmann/parallel_sample_sortBTCEU1                     pS5: binary tree, equality, bktcache, unroll tree
bingmann/parallel_sample_sortBTCEU1_lcp                 pS5: binary tree, equality, bktcache, unroll tree_lcp
bingmann/parallel_sample_sortBTCE_lcp                   pS5: binary tree, equality, bktcache_lcp
bingmann/parallel_sample_sortBTCT                       pS5: binary tree, bktcache, tree calc
bingmann/parallel_sample_sortBTCTU1                     pS5: binary tree, bktcache, unroll tree, tree calc
bingmann/parallel_sample_sortBTCTU1_lcp                 pS5: binary tree, bktcache, unroll tree, tree calc_lcp
bingmann/parallel_sample_sortBTCTU2                     pS5: binary tree, bktcache, unroll tree and strings, tree calc
bingmann/parallel_sample_sortBTCTU2_lcp                 pS5: binary tree, bktcache, unroll tree and strings, tree calc_lcp
bingmann/parallel_sample_sortBTCT_lcp                   pS5: binary tree, bktcache, tree calc_lcp
bingmann/parallel_sample_sortBTCU1                      pS5: binary tree, bktcache, unroll tree
bingmann/parallel_sample_sortBTCU1_lcp                  pS5: binary tree, bktcache, unroll tree_lcp
bingmann/parallel_sample_sortBTCU2                      pS5: binary tree, bktcache, unroll tree and strings
bingmann/parallel_sample_sortBTCU2_lcp                  pS5: binary tree, bktcache, unroll tree and strings_lcp
bingmann/parallel_sample_sortBTCU2_out                  pS5: binary tree, bktcache, unroll tree and strings, separate output
bingmann/parallel_sample_sortBTCU2_out_lcp              pS5: binary tree, bktcache, unroll tree and strings, separate output_lcp
bingmann/parallel_sample_sortBTC_lcp                    pS5: binary tree, bktcache_lcp
bingmann/qsort1                                         Run stdlib qsort with string comparsions (bytewise)
bingmann/qsort4                                         Run stdlib qsort with string comparsions (4 bytewise)
bingmann/qsort8                                         Run stdlib qsort with string comparsions (8 bytewise)
bingmann/qsort_strcmp                                   Run stdlib qsort with strcmp comparsion
bingmann/sample_sortBS                                  bingmann/sample_sortBS (binary search, no cache)
bingmann/sample_sortBSC                                 bingmann/sample_sortBSC (binary search, bkt cache)
bingmann/sample_sortBSCA                                bingmann/sample_sortBSCA (binary search, assembler CMOV, bkt cache)
bingmann/sample_sortBSCE                                bingmann/sample_sortBSCE (binary search equal, bkt cache)
bingmann/sample_sortBSCEA                               bingmann/sample_sortBSCEA (binary search equal, assembler CMOV, bkt cache)
bingmann/sample_sortBT                                  bingmann/sample_sortBT (binary tree, no cache)
bingmann/sample_sortBTC                                 bingmann/sample_sortBTC (binary tree, bkt cache)
bingmann/sample_sortBTCA                                bingmann/sample_sortBTCA (binary tree, asm CMOV, bkt cache)
bingmann/sample_sortBTCE2                               bingmann/sample_sortBTCE2 (binary tree equal, bkt cache)
bingmann/sample_sortBTCE2A                              bingmann/sample_sortBTCE2A (binary tree equal, asm CMOV, bkt cache)
bingmann/sample_sortBTCE2U                              bingmann/sample_sortBTCE2U (binary tree equal unroll, asm CMOV, bkt cache)
bingmann/sample_sortBTCE3                               bingmann/sample_sortBTCE3 (adapt binary tree equal, bkt cache)
bingmann/sample_sortBTCE3A                              bingmann/sample_sortBTCE3A (adapt binary tree equal, asm CMOV, bkt cache)
bingmann/sample_sortBTCE3U                              bingmann/sample_sortBTCE3U (adapt binary tree equal unroll, asm CMOV, bkt cache)
bingmann/sample_sortBTCU                                bingmann/sample_sortBTCU (binary tree, unrolled, bkt cache)
bingmann/sample_sortRBTCE                               bingmann/sample_sortRBTCE (adapt binary tree equal, bkt cache)
bingmann/sample_sortRBTCEA                              bingmann/sample_sortRBTCEA (adapt binary tree equal, asm CMOV, bkt cache)
bingmann/sequential_mkqs_cache8                         multikey_cache with 8byte cache (non-recursive)
bingmann/stdsort1                                       Run std::sort with string comparsions (bytewise)
bingmann/stdsort4                                       Run std::sort with string comparsions (4 bytewise)
bingmann/stdsort8                                       Run std::sort with string comparsions (8 bytewise)
bs/mkqsort                                              bs_mkqs Original Multikey-Quicksort
eberle/lcp_insertion_sort                               LCP aware inssertion sort by Andreas Eberle
eberle/lcp_insertion_sort_cache                         LCP aware insertion sort with cached characters calculation by Andreas Eberle
eberle/mergesort_lcp_binary                             Binary Mergesort with LCP-usage by Andreas Eberle
eberle/mergesort_lcp_losertree_16way                    Mergesort with lcp aware Losertree by Andreas Eberle
eberle/mergesort_lcp_losertree_32way                    Mergesort with lcp aware Losertree by Andreas Eberle
eberle/mergesort_lcp_losertree_4way                     Mergesort with lcp aware Losertree by Andreas Eberle
eberle/mergesort_lcp_losertree_64way                    Mergesort with lcp aware Losertree by Andreas Eberle
eberle/parallel-lcp-mergesort                           parallel LCP aware mergesort by Andreas Eberle
eberle/ps5-parallel-toplevel-merge                      NUMA aware sorting algorithm running pS5 on local memory and then doing a parallel merge by Andreas Eberle
eberle/ps5-parallel-toplevel-merge-assisting            pS5-LCP-Merge with JobQueue assisting each other by Andreas Eberle and Timo Bingmann
insertion_sort                                          String Insertion-Sort
mbm/radixsort                                           MSD Radix Sort by P. M. McIlroy, K. Bostic, and M. D. McIlroy
ng/cradix                                               CRadix Original by Waihong Ng and Katsuhiko Kakehi
ng/lcpmergesort                                         LCP-Mergesort Original by Waihong Ng and Katsuhiko Kakehi
ng/rantala_cradix                                       CRadix by Waihong Ng and Katsuhiko Kakehi modified by Rantala
nilsson/adaptive_msd                                    Adaptive MSD Radix Sort by Stefan Nilsson
nilsson/forward16                                       Forward Radix Sort 16-bit by Stefan Nilsson
nilsson/forward8                                        Forward Radix Sort 8-bit by Stefan Nilsson
rantala/burstsort2_bagwell                              burstsort2 with vector_bagwell bucket type
rantala/burstsort2_brodnik                              burstsort2 with vector_brodnik bucket type
rantala/burstsort2_sampling_bagwell                     burstsort2 sampling with vector_bagwell bucket type
rantala/burstsort2_sampling_brodnik                     burstsort2 sampling with vector_brodnik bucket type
rantala/burstsort2_sampling_superalphabet_bagwell       burstsort2 sampling superalphabet with vector_bagwell bucket type
rantala/burstsort2_sampling_superalphabet_brodnik       burstsort2 sampling superalphabet with vector_brodnik bucket type
rantala/burstsort2_sampling_superalphabet_vector        burstsort2 sampling superalphabet with std::vector bucket type
rantala/burstsort2_sampling_superalphabet_vector_block  burstsort2 sampling superalphabet with vector_block bucket type
rantala/burstsort2_sampling_vector                      burstsort2 sampling with std::vector bucket type
rantala/burstsort2_sampling_vector_block                burstsort2 sampling with vector_block bucket type
rantala/burstsort2_superalphabet_bagwell                burstsort2 superalphabet with vector_bagwell bucket type
rantala/burstsort2_superalphabet_brodnik                burstsort2 superalphabet with vector_brodnik bucket type
rantala/burstsort2_superalphabet_vector                 burstsort2 superalphabet with std::vector bucket type
rantala/burstsort2_superalphabet_vector_block           burstsort2 superalphabet with vector_block bucket type
rantala/burstsort2_vector                               burstsort2 with std::vector bucket type
rantala/burstsort2_vector_block                         burstsort2 with vector_block bucket type
rantala/burstsort_bagwell                               burstsort with vector_bagwell bucket type
rantala/burstsort_brodnik                               burstsort with vector_brodnik bucket type
rantala/burstsort_mkq_recursiveburst_1                  burstsort_mkq 1byte alphabet with recursiveburst
rantala/burstsort_mkq_recursiveburst_2                  burstsort_mkq 2byte alphabet with recursiveburst
rantala/burstsort_mkq_recursiveburst_4                  burstsort_mkq 4byte alphabet with recursiveburst
rantala/burstsort_mkq_simpleburst_1                     burstsort_mkq 1byte alphabet with simpleburst
rantala/burstsort_mkq_simpleburst_2                     burstsort_mkq 2byte alphabet with simpleburst
rantala/burstsort_mkq_simpleburst_4                     burstsort_mkq 4byte alphabet with simpleburst
rantala/burstsort_sampling_bagwell                      burstsort sampling with vector_bagwell bucket type
rantala/burstsort_sampling_brodnik                      burstsort sampling with vector_brodnik bucket type
rantala/burstsort_sampling_superalphabet_bagwell        burstsort sampling superalphabet with vector_bagwell bucket type
rantala/burstsort_sampling_superalphabet_brodnik        burstsort sampling superalphabet with vector_brodnik bucket type
rantala/burstsort_sampling_superalphabet_vector         burstsort sampling superalphabet with std::vector bucket type
rantala/burstsort_sampling_superalphabet_vector_block   burstsort sampling superalphabet with vector_block bucket type
rantala/burstsort_sampling_vector                       burstsort sampling with std::vector bucket type
rantala/burstsort_sampling_vector_block                 burstsort sampling with vector_block bucket type
rantala/burstsort_superalphabet_bagwell                 burstsort superalphabet with vector_bagwell bucket type
rantala/burstsort_superalphabet_brodnik                 burstsort superalphabet with vector_brodnik bucket type
rantala/burstsort_superalphabet_vector                  burstsort superalphabet with std::vector bucket type
rantala/burstsort_superalphabet_vector_block            burstsort superalphabet with vector_block bucket type
rantala/burstsort_vector                                burstsort with std::vector bucket type
rantala/burstsort_vector_block                          burstsort with vector_block bucket type
rantala/funnelsort_128way_bfs                           funnelsort 128way bfs
rantala/funnelsort_128way_dfs                           funnelsort 128way dfs
rantala/funnelsort_16way_bfs                            funnelsort 16way bfs
rantala/funnelsort_16way_dfs                            funnelsort 16way dfs
rantala/funnelsort_32way_bfs                            funnelsort 32way bfs
rantala/funnelsort_32way_dfs                            funnelsort 32way dfs
rantala/funnelsort_64way_bfs                            funnelsort 64way bfs
rantala/funnelsort_64way_dfs                            funnelsort 64way dfs
rantala/funnelsort_8way_bfs                             funnelsort 8way bfs
rantala/funnelsort_8way_dfs                             funnelsort 8way dfs
rantala/mergesort_2way                                  mergesort with 2way merger
rantala/mergesort_2way_parallel                         mergesort_parallel with 2way merger
rantala/mergesort_2way_unstable                         mergesort 2way unstable
rantala/mergesort_3way                                  mergesort with 3way merger
rantala/mergesort_3way_parallel                         mergesort_parallel with 3way merger
rantala/mergesort_3way_unstable                         mergesort 3way unstable
rantala/mergesort_4way                                  mergesort with 4way merger
rantala/mergesort_4way_parallel                         mergesort_parallel with 4way merger
rantala/mergesort_4way_unstable                         mergesort 4way unstable
rantala/mergesort_cache1_lcp_2way                       mergesort LCP with 2way merger and 1byte cache
rantala/mergesort_cache1_lcp_2way_parallel              mergesort_parallel LCP with 2way merger and 1byte cache
rantala/mergesort_cache2_lcp_2way                       mergesort LCP with 2way merger and 2byte cache
rantala/mergesort_cache2_lcp_2way_parallel              mergesort_parallel LCP with 2way merger and 2byte cache
rantala/mergesort_cache4_lcp_2way                       mergesort LCP with 2way merger and 4byte cache
rantala/mergesort_cache4_lcp_2way_parallel              mergesort_parallel LCP with 2way merger and 4byte cache
rantala/mergesort_lcp_2way                              mergesort LCP with 2way merger
rantala/mergesort_lcp_2way_parallel                     mergesort_parallel LCP with 2way merger
rantala/mergesort_lcp_2way_unstable                     mergesort Unstable LCP with 2way merger
rantala/mergesort_lcp_2way_unstable_parallel            mergesort_parallel unstable LCP with 2way merger
rantala/mergesort_losertree_1024way                     mergesort 1024way loser tree based
rantala/mergesort_losertree_1024way_parallel            mergesort parallel 1024way loser tree based
rantala/mergesort_losertree_128way                      mergesort 128way loser tree based
rantala/mergesort_losertree_128way_parallel             mergesort parallel 128way loser tree based
rantala/mergesort_losertree_256way                      mergesort 256way loser tree based
rantala/mergesort_losertree_256way_parallel             mergesort parallel 256way loser tree based
rantala/mergesort_losertree_512way                      mergesort 512way loser tree based
rantala/mergesort_losertree_512way_parallel             mergesort parallel 512way loser tree based
rantala/mergesort_losertree_64way                       mergesort 64way loser tree based
rantala/mergesort_losertree_64way_parallel              mergesort parallel 64way loser tree based
rantala/msd_A                                           msd_A
rantala/msd_A2                                          msd_A2
rantala/msd_A2_adaptive                                 msd_A2_adaptive
rantala/msd_A_adaptive                                  msd_A_adaptive
rantala/msd_CE0                                         msd_CE0: baseline
rantala/msd_CE1                                         msd_CE1: oracle
rantala/msd_CE2                                         msd_CE2: oracle+loop fission
rantala/msd_CE3                                         msd_CE3: oracle+loop fission+adaptive
rantala/msd_CE4                                         msd_CE4: oracle+loop fission+adaptive+16bit counter
rantala/msd_CE5                                         msd_CE5: oracle+loop fission+adaptive+16bit counter+prealloc
rantala/msd_CE6                                         msd_CE6: oracle+loop fission+adaptive+16bit counter+prealloc+unroll
rantala/msd_CE7                                         msd_CE7: oracle+loop fission+adaptive+16bit counter+prealloc+unroll+sortedness
rantala/msd_CI                                          msd_CI
rantala/msd_CI_adaptive                                 msd_CI: adaptive
rantala/msd_DB                                          msd_DB
rantala/msd_D_std_deque                                 msd_D_std_deque
rantala/msd_D_std_deque_adaptive                        msd_D_std_deque_adaptive
rantala/msd_D_std_list                                  msd_D_std_list
rantala/msd_D_std_list_adaptive                         msd_D_std_list_adaptive
rantala/msd_D_std_vector                                msd_D_std_vector
rantala/msd_D_std_vector_adaptive                       msd_D_std_vector_adaptive
rantala/msd_D_vector_bagwell                            msd_D_vector_bagwell
rantala/msd_D_vector_bagwell_adaptive                   msd_D_vector_bagwell_adaptive
rantala/msd_D_vector_block                              msd_D_vector_block
rantala/msd_D_vector_block_adaptive                     msd_D_vector_block_adaptive
rantala/msd_D_vector_brodnik                            msd_D_vector_brodnik
rantala/msd_D_vector_brodnik_adaptive                   msd_D_vector_brodnik_adaptive
rantala/msd_D_vector_malloc                             msd_D_vector_malloc
rantala/msd_D_vector_malloc_adaptive                    msd_D_vector_malloc_adaptive
rantala/msd_D_vector_malloc_counter_clear               msd_D_vector_malloc_counter_clear
rantala/msd_D_vector_malloc_counter_clear_adaptive      msd_D_vector_malloc_counter_clear_adaptive
rantala/msd_D_vector_realloc                            msd_D_vector_realloc
rantala/msd_D_vector_realloc_adaptive                   msd_D_vector_realloc_adaptive
rantala/msd_D_vector_realloc_counter_clear              msd_D_vector_realloc_counter_clear
rantala/msd_D_vector_realloc_counter_clear_adaptive     msd_D_vector_realloc_counter_clear_adaptive
rantala/msd_D_vector_realloc_shrink_clear               msd_D_vector_realloc_shrink_clear
rantala/msd_D_vector_realloc_shrink_clear_adaptive      msd_D_vector_realloc_shrink_clear_adaptive
rantala/multikey_block1                                 multikey_block with 1byte alphabet
rantala/multikey_block2                                 multikey_block with 2byte alphabet
rantala/multikey_block4                                 multikey_block with 4byte alphabet
rantala/multikey_cache4                                 multikey_cache with 4byte cache
rantala/multikey_cache8                                 multikey_cache with 8byte cache
rantala/multikey_dynamic_bagwell1                       multikey_dynamic with vector_bagwell bucket type and 1byte alphabet
rantala/multikey_dynamic_bagwell2                       multikey_dynamic with vector_bagwell bucket type and 2byte alphabet
rantala/multikey_dynamic_bagwell4                       multikey_dynamic with vector_bagwell bucket type and 4byte alphabet
rantala/multikey_dynamic_brodnik1                       multikey_dynamic with vector_brodnik bucket type and 1byte alphabet
rantala/multikey_dynamic_brodnik2                       multikey_dynamic with vector_brodnik bucket type and 2byte alphabet
rantala/multikey_dynamic_brodnik4                       multikey_dynamic with vector_brodnik bucket type and 4byte alphabet
rantala/multikey_dynamic_vector1                        multikey_dynamic with std::vector bucket type and 1byte alphabet
rantala/multikey_dynamic_vector2                        multikey_dynamic with std::vector bucket type and 2byte alphabet
rantala/multikey_dynamic_vector4                        multikey_dynamic with std::vector bucket type and 4byte alphabet
rantala/multikey_dynamic_vector_block1                  multikey_dynamic with vector_block bucket type and 1byte alphabet
rantala/multikey_dynamic_vector_block2                  multikey_dynamic with vector_block bucket type and 2byte alphabet
rantala/multikey_dynamic_vector_block4                  multikey_dynamic with vector_block bucket type and 4byte alphabet
rantala/multikey_multipivot_brute_simd1                 multikey_multipivot brute_simd with 1byte alphabet
rantala/multikey_multipivot_brute_simd2                 multikey_multipivot brute_simd with 2byte alphabet
rantala/multikey_multipivot_brute_simd4                 multikey_multipivot brute_simd with 4byte alphabet
rantala/multikey_simd1                                  multikey_simd with 1byte alphabet
rantala/multikey_simd2                                  multikey_simd with 2byte alphabet
rantala/multikey_simd4                                  multikey_simd with 4byte alphabet
rantala/multikey_simd_parallel1                         multikey_simd_parallel with 1byte alphabet
rantala/multikey_simd_parallel2                         multikey_simd_parallel with 2byte alphabet
rantala/multikey_simd_parallel4                         multikey_simd_parallel with 4byte alphabet
shamsundar/lcp-merge-string-sort                        Parallelized LCP Merge sort by N. Shamsundar
sinha/CPL_burstsort                                     Original CPL-burstsort
sinha/CP_burstsort                                      Original CP-burstsort
sinha/C_burstsort                                       Original C-burstsort
sinha/burstsortA                                        burstsortA Original Burstsort with arrays
sinha/burstsortL                                        burstsortL Original Burstsort with linked-lists
sinha/fbCPL_burstsort                                   Original fbCPL-burstsort
sinha/fbCP_burstsort                                    Original fbCP-burstsort
sinha/fbC_burstsort                                     Original fbC-burstsort
sinha/sCPL_burstsort                                    Original sCPL-burstsort
sinha/sCP_burstsort                                     Original sCP-burstsort
sinha/sC_burstsort                                      Original sC-burstsort
```
