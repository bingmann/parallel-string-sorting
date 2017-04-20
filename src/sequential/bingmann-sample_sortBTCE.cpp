/*******************************************************************************
 * src/sequential/bingmann-sample_sortBTCE.cpp
 *
 * Experiments with sequential Super Scalar String Sample-Sort (S^5).
 *
 * Binary tree search with equality branch and bucket cache.
 *
 *******************************************************************************
 * Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "bingmann-sample_sort.hpp"
#include "bingmann-sample_sort_common.hpp"

namespace bingmann_sample_sortBTCE {

using namespace bingmann_sample_sort;

#if 0
static const size_t numsplitters2 = 16;
static const size_t treebits = tlx::Log2Floor<numsplitters2>::value;
static const size_t numsplitters = (1 << treebits) - 1;
#else
//static const size_t l2cache = 256*1024;

// bounding equations:
// splitters            + bktsize
// n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

//static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));
//static const size_t numsplitters2 = l2cache / sizeof(key_type);
static const size_t numsplitters2 = (l2cache - sizeof(size_t)) / (2 * sizeof(size_t));

//static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / ( sizeof(key_type) );

static const size_t treebits = tlx::Log2Floor<numsplitters2>::value;
static const size_t numsplitters = (1 << treebits) - 1;
#endif

// ----------------------------------------------------------------------------

// Variant of string sample-sort: use super-scalar binary search on splitters
// with equality check and index caching.

template <size_t treebits>
class ClassifySimple
{
public:
    static const size_t numsplitters = (1 << treebits) - 1;

    //! search in splitter tree for bucket number
    static inline unsigned int
    find_bkt(const key_type& key, const key_type* splitter_tree0)
    {
        // binary tree traversal without left branch

        const key_type* splitter_tree = splitter_tree0 - 1;
        unsigned int i = 1;

        while (i <= numsplitters)
        {
            if (key == splitter_tree[i])
                return 2 * TreeCalculations<treebits>::level_to_preorder(i) - 1;
            else if (key < splitter_tree[i])
                i = 2 * i + 0;
            else    // (key > splitter_tree[i])
                i = 2 * i + 1;
        }

        i -= numsplitters + 1;

        return 2 * i; // < or > bucket
    }
};

template <size_t treebits>
class ClassifyAssembler
{
public:
    static const size_t numsplitters = (1 << treebits) - 1;

    // binary search on splitter array for bucket number
    static inline
    unsigned int find_bkt(const key_type& key, const key_type* splitter_tree0)
    {
        const key_type* splitter_tree = splitter_tree0 - 1;
        unsigned int i;

        // hand-coded assembler binary tree traversal with equality
        asm ("mov    $1, %%rax \n"             // rax = i
             // body of while loop
             "1: \n"
             "cmpq	(%[splitter_tree],%%rax,8), %[key] \n"
             "je     2f \n"
             "lea    (%%rax,%%rax), %%rax \n"
             "lea    1(%%rax), %%rcx \n"
             "cmova  %%rcx, %%rax \n"            // CMOV rax = 2 * i + 1
             "cmp    %[numsplitters1], %%rax \n" // i < numsplitters+1
             "jb     1b \n"
             "sub    %[numsplitters1], %%rax \n" // i -= numsplitters+1;
             "lea    (%%rax,%%rax), %%rax \n"    // i = i*2
             "jmp    3f \n"
             "2: \n"
             "bsr    %%rax, %%rdx \n"            // dx = bit number of highest one
             "mov    %[treebits], %%rcx \n"
             "sub    %%rdx, %%rcx \n"            // cx = treebits - highest
             "shl    %%cl, %%rax \n"             // shift ax to left
             "and    %[numsplitters], %%rax \n"  // mask off other bits
             "lea    -1(%%rcx), %%rcx \n"
             "mov    $1, %%rdx \n"               // dx = (1 << (hi-1))
             "shl    %%cl, %%rdx \n"             //
             "or     %%rdx, %%rax \n"            // ax = OR of both
             "lea    -1(%%rax,%%rax), %%rax \n"  // i = i * 2 - 1
             "3: \n"
             : "=&a" (i)
             :[key] "r" (key), [splitter_tree] "r" (splitter_tree),
             [numsplitters1] "g" (numsplitters + 1),
             [treebits] "g" (treebits),
             [numsplitters] "g" (numsplitters)
             : "rcx", "rdx");

        return i;
    }
};

template <size_t treebits>
struct ClassifyUnroll
{
#define SPLITTER_TREE_STEP                         \
    /* inside the tree */                          \
    "cmpq   (%[splitter_tree],%%rax,8), %[key] \n" \
    "je     2f \n"                                 \
    "lea    (%%rax,%%rax), %%rax \n"               \
    "lea    1(%%rax), %%rcx \n"                    \
    "cmova  %%rcx, %%rax \n"             /* CMOV rax = 2 * i + 1 */

#define SPLITTER_TREE_END                                                     \
    /* at leaf level */                                                       \
    "sub    %[numsplitters1], %%rax \n"  /* i -= numsplitters+1; */           \
    "lea    (%%rax,%%rax), %%rax \n"     /* i = i*2 */                        \
    "jmp    3f \n"                                                            \
    "2: \n"                                                                   \
    "bsr    %%rax, %%rdx \n"             /* dx = bit number of highest one */ \
    "mov    %[treebits], %%rcx \n"                                            \
    "sub    %%rdx, %%rcx \n"             /* cx = treebits - highest */        \
    "shl    %%cl, %%rax \n"              /* shift ax to left */               \
    "and    %[numsplitters], %%rax \n"   /* mask off other bits */            \
    "lea    -1(%%rcx), %%rcx \n"                                              \
    "mov    $1, %%rdx \n"                /* dx = (1 << (hi-1)) */             \
    "shl    %%cl, %%rdx \n"                                                   \
    "or     %%rdx, %%rax \n"             /* ax = OR of both */                \
    "lea    -1(%%rax,%%rax), %%rax \n"   /* i = i * 2 - 1 */                  \
    "3: \n"
};

template <>
struct ClassifyUnroll<11>
{
    static const size_t treebits = 11;
    static const size_t numsplitters = (1 << treebits) - 1;

    //! binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt(const key_type& key, const key_type* splitter_tree0)
    {
        const key_type* splitter_tree = splitter_tree0 - 1;
        unsigned int i;

        // hand-coded assembler binary tree traversal with equality, using CMOV
        asm ("mov    $1, %%rax \n"             // rax = i
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP

             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP

             SPLITTER_TREE_STEP

             SPLITTER_TREE_END

             : "=&a" (i)
             :[key] "r" (key), [splitter_tree] "r" (splitter_tree),
             [numsplitters1] "g" (numsplitters + 1),
             [treebits] "g" (treebits),
             [numsplitters] "g" (numsplitters)
             : "rcx", "rdx");

        return i;
    }
};

template <>
struct ClassifyUnroll<12>
{
    static const size_t treebits = 12;
    static const size_t numsplitters = (1 << treebits) - 1;

    //! binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt(const key_type& key, const key_type* splitter_tree0)
    {
        const key_type* splitter_tree = splitter_tree0 - 1;
        unsigned int i;

        // hand-coded assembler binary tree traversal with equality, using CMOV
        asm ("mov    $1, %%rax \n"             // rax = i
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP

             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP

             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP

             SPLITTER_TREE_END

             : "=&a" (i)
             :[key] "r" (key), [splitter_tree] "r" (splitter_tree),
             [numsplitters1] "g" (numsplitters + 1),
             [treebits] "g" (treebits),
             [numsplitters] "g" (numsplitters)
             : "rcx", "rdx");

        return i;
    }
};

template <>
struct ClassifyUnroll<13>
{
    static const size_t treebits = 13;
    static const size_t numsplitters = (1 << treebits) - 1;

    //! binary search on splitter array for bucket number
    static inline unsigned int
    find_bkt(const key_type& key, const key_type* splitter_tree0)
    {
        const key_type* splitter_tree = splitter_tree0 - 1;
        unsigned int i;

        // hand-coded assembler binary tree traversal with equality, using CMOV
        asm ("mov    $1, %%rax \n"             // rax = i
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP

             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP

             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP
             SPLITTER_TREE_STEP

             SPLITTER_TREE_END

             : "=&a" (i)
             :[key] "r" (key), [splitter_tree] "r" (splitter_tree),
             [numsplitters1] "g" (numsplitters + 1),
             [treebits] "g" (treebits),
             [numsplitters] "g" (numsplitters)
             : "rcx", "rdx");

        return i;
    }
};

// *****************************************************************************
// Variant of string sample-sort: use super-scalar binary search on splitters,
// with index caching.

template <template <size_t> class Classify>
void sample_sortBTCE2(string* strings, size_t n, size_t depth)
{
    if (n < g_samplesort_smallsort)
    {
        g_rs_steps++;
        g_timer.change(TM_SMALLSORT);
        sample_sort_small_sort(strings, n, depth);
        g_timer.change(TM_GENERAL);
        return;
    }
    g_ss_steps++;

    // step 1: select splitters with oversampling
    g_timer.change(TM_MAKE_SAMPLE);

    const size_t samplesize = oversample_factor * numsplitters;

    key_type samples[samplesize];

    LCGRandom rng(&strings);

    for (unsigned int i = 0; i < samplesize; ++i)
    {
        samples[i] = get_char<key_type>(strings[rng() % n], depth);
    }

    std::sort(samples, samples + samplesize);

    g_timer.change(TM_MAKE_SPLITTER);

    key_type* splitter_tree = new key_type[numsplitters];
    unsigned char* splitter_lcp = new unsigned char[numsplitters];

    {
        unsigned int treelvl[treebits]; // treebits == treelevels

        // fill in initial level indexes
        int idx = 1;
        for (unsigned int i = treebits; i > 0; ) {
            treelvl[--i] = idx;
            idx <<= 1;
        }

        key_type* splitter_tree1 = splitter_tree - 1;
        key_type prevsplitter = 0;

        LOGC(debug_splitter) << "splitter:";
        splitter_lcp[0] = 0; // sentinel for first < everything bucket
        for (size_t i = 0, j = oversample_factor / 2; i < numsplitters; ++i)
        {
            const key_type& splitter = samples[j];

            int l = __builtin_ctz(i + 1);
            LOGC(debug_splitter)
                << "splitter[" << i << "] on level " << l
                << " = tree[" << treelvl[l] << "] = key "
                << tlx::hexdump_type(splitter);
            splitter_tree1[treelvl[l]++] = samples[j];

            if (i != 0) {
                key_type xorSplit = prevsplitter ^ splitter;

                LOGC(debug_splitter)
                    << "    XOR -> " << tlx::hexdump_type(xorSplit) << " - "
                    << count_high_zero_bits(xorSplit)
                    << " bits = " << count_high_zero_bits(xorSplit) / 8
                    << " chars lcp";

                splitter_lcp[i] =
                    (count_high_zero_bits(xorSplit) / 8) |
                    // marker for done splitters
                    ((splitter & 0xFF) ? 0 : 0x80);
            }

            prevsplitter = splitter;
            j += oversample_factor;
        }
    }

    if (debug_splitter_tree)
    {
        DBG1(1, "splitter_tree: ");
        for (size_t i = 0; i < numsplitters; ++i)
        {
            DBG2(1, splitter_tree[i] << " ");
        }
        DBG3(1, "");
    }

    // step 2.2: classify all strings and count bucket sizes
    g_timer.change(TM_CLASSIFY);

    uint16_t* bktcache = new uint16_t[n];

    static const size_t bktnum = 2 * numsplitters + 1;

    for (size_t si = 0; si < n; ++si)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(strings[si], depth);

        unsigned int b = Classify<treebits>::find_bkt(key, splitter_tree);

        assert(b < bktnum);

        bktcache[si] = b;
    }

    delete[] splitter_tree;

    size_t* bktsize = new size_t[bktnum];
    memset(bktsize, 0, bktnum * sizeof(size_t));

    for (size_t si = 0; si < n; ++si)
        ++bktsize[bktcache[si]];

    if (debug_bucketsize)
    {
        DBG1(1, "bktsize: ");
        for (size_t i = 0; i < bktnum; ++i)
        {
            DBG2(1, bktsize[i] << " ");
        }
        DBG3(1, "");
    }

    // step 3: prefix sum
    g_timer.change(TM_PREFIXSUM);

    size_t bktindex[bktnum];
    bktindex[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned int i = 1; i < bktnum; ++i) {
        bktindex[i] = bktindex[i - 1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }
    assert(bktindex[bktnum - 1] == n);

    // step 4: premute in-place
    g_timer.change(TM_PERMUTE);

    for (size_t i = 0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        uint16_t permbkt = bktcache[i];

        while ((j = --bktindex[permbkt]) > i)
        {
            std::swap(perm, strings[j]);
            std::swap(permbkt, bktcache[j]);
        }

        strings[i] = perm;
        i += bktsize[permbkt];
    }

    delete[] bktcache;

    // step 5: recursion
    g_timer.change(TM_GENERAL);

    size_t i = 0, bsum = 0;
    while (i < bktnum - 1)
    {
        // i is even -> bkt[i] is less-than bucket
        if (bktsize[i] > 1)
        {
            LOGC(debug_recursion)
                << "Recurse[" << depth << "]: < bkt " << bsum
                << " size " << bktsize[i]
                << " lcp " << int(splitter_lcp[i / 2] & 0x7F);

            if (!g_toplevel_only)
                sample_sortBTCE2<Classify>(
                    strings + bsum, bktsize[i],
                    depth + (splitter_lcp[i / 2] & 0x7F));
        }
        bsum += bktsize[i++];

        // i is odd -> bkt[i] is equal bucket
        if (bktsize[i] > 1)
        {
            if (splitter_lcp[i / 2] & 0x80) {
                // equal-bucket has NULL-terminated key, done.
                LOGC(debug_recursion)
                    << "Recurse[" << depth << "]: = bkt " << bsum
                    << " size " << bktsize[i] << " is done!";
            }
            else {
                LOGC(debug_recursion)
                    << "Recurse[" << depth << "]: = bkt " << bsum
                    << " size " << bktsize[i] << " lcp keydepth!";

                if (!g_toplevel_only)
                    sample_sortBTCE2<Classify>(strings + bsum, bktsize[i],
                                               depth + sizeof(key_type));
            }
        }
        bsum += bktsize[i++];
    }
    if (bktsize[i] > 0)
    {
        LOGC(debug_recursion)
            << "Recurse[" << depth << "]: > bkt " << bsum
            << " size " << bktsize[i] << " no lcp";

        if (!g_toplevel_only)
            sample_sortBTCE2<Classify>(strings + bsum, bktsize[i], depth);
    }
    bsum += bktsize[i++];
    assert(i == bktnum && bsum == n);

    delete[] splitter_lcp;
    delete[] bktsize;
}

void bingmann_sample_sortBTCE2(string* strings, size_t n)
{
    sample_sort_pre();
    g_stats >> "numsplitters" << numsplitters
        >> "splitter_treebits" << treebits;
    sample_sortBTCE2<ClassifySimple>(strings, n, 0);
    sample_sort_post();
}

PSS_CONTESTANT(bingmann_sample_sortBTCE2, "bingmann/sample_sortBTCE2",
               "bingmann/sample_sortBTCE2 (binary tree equal, bkt cache)")

void bingmann_sample_sortBTCE2A(string* strings, size_t n)
{
    sample_sort_pre();
    g_stats >> "numsplitters" << numsplitters
        >> "splitter_treebits" << treebits;
    sample_sortBTCE2<ClassifyAssembler>(strings, n, 0);
    sample_sort_post();
}

PSS_CONTESTANT(bingmann_sample_sortBTCE2A, "bingmann/sample_sortBTCE2A",
               "bingmann/sample_sortBTCE2A (binary tree equal, asm CMOV, bkt cache)")

void bingmann_sample_sortBTCE2U(string* strings, size_t n)
{
    sample_sort_pre();
    g_stats >> "numsplitters" << numsplitters
        >> "splitter_treebits" << treebits;
    sample_sortBTCE2<ClassifyUnroll>(strings, n, 0);
    sample_sort_post();
}

PSS_CONTESTANT(bingmann_sample_sortBTCE2U, "bingmann/sample_sortBTCE2U",
               "bingmann/sample_sortBTCE2U (binary tree equal unroll, asm CMOV, bkt cache)")

// ----------------------------------------------------------------------------

// Variant of string sample-sort: use recursive binary function on sample array
// to generate splitters.

class SampleSortBTCE3
{
public:
    //! Variant of string sample-sort: use super-scalar binary search on
    //! splitters, with index caching.
    template <template <size_t> class Classify>
    static void sort(string* strings, size_t n, size_t depth)
    {
        if (n < g_samplesort_smallsort)
        {
            g_rs_steps++;
            g_timer.change(TM_SMALLSORT);
            sample_sort_small_sort(strings, n, depth);
            g_timer.change(TM_GENERAL);
            return;
        }
        g_ss_steps++;

        // step 1: select splitters with oversampling
        g_timer.change(TM_MAKE_SAMPLE);

        const size_t samplesize = oversample_factor * numsplitters;

        key_type samples[samplesize];

        LCGRandom rng(&strings);

        for (unsigned int i = 0; i < samplesize; ++i)
        {
            samples[i] = get_char<key_type>(strings[rng() % n], depth);
        }

        std::sort(samples, samples + samplesize);

        g_timer.change(TM_MAKE_SPLITTER);

        key_type splitter_tree[numsplitters + 1];
        unsigned char splitter_lcp[numsplitters + 1];

        TreeBuilderFullDescentTreeCalc<numsplitters>(
            splitter_tree, splitter_lcp, samples, samplesize);

#if 0
        key_type cmp_splitter_tree[numsplitters + 1];
        unsigned char cmp_splitter_lcp[numsplitters];

        {
            unsigned int treelvl[treebits]; // treebits == treelevels

            // fill in initial level indexes
            int idx = 1;
            for (unsigned int i = treebits; i > 0; ) {
                treelvl[--i] = idx;
                idx <<= 1;
            }

            key_type prevsplitter;

            LOGC(debug_splitter) << "splitter:";
            cmp_splitter_lcp[0] = 0; // sentinel for first < everything bucket
            for (size_t i = 0, j = oversample_factor / 2; i < numsplitters; ++i)
            {
                const key_type& splitter = samples[j];

                int l = __builtin_ctz(i + 1);
                LOGC(debug_splitter)
                    << "splitter[" << i << "] on level " << l
                    << " = tree[" << treelvl[l] << "] = key "
                    << tlx::hexdump_type(splitter);
                cmp_splitter_tree[treelvl[l]++] = splitter;

                if (i != 0) {
                    key_type xorSplit = prevsplitter ^ splitter;

                    LOGC(debug_splitter)
                        << "    XOR -> " << tlx::hexdump_type(xorSplit) << " - ";

                    LOGC(debug_splitter)
                        << count_high_zero_bits(xorSplit) << " bits = "
                        << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

                    cmp_splitter_lcp[i] =
                        (count_high_zero_bits(xorSplit) / 8) |
                        // marker for done splitters
                        ((splitter & 0xFF) ? 0 : 0x80);
                }

                prevsplitter = splitter;
                j += oversample_factor;
            }
        }

        for (size_t i = 0; i < numsplitters; ++i)
        {
            LOG1 << "splitter_tree[" << i + 1 << "] = "
                 << tree.splitter_tree[i + 1] << " -? "
                 << cmp_splitter_tree[i + 1];
            assert(tree.splitter_tree[i + 1] == cmp_splitter_tree[i + 1]);

            LOG1 << "splitter_lcp[" << i << "] = "
                 << int(tree.splitter_lcp[i])
                 << " -? " << int(cmp_splitter_lcp[i]);
            assert(tree.splitter_lcp[i] == cmp_splitter_lcp[i]);
        }
#endif

        // step 2.2: classify all strings and count bucket sizes
        g_timer.change(TM_CLASSIFY);

        uint16_t* bktcache = new uint16_t[n];

        static const size_t bktnum = 2 * numsplitters + 1;

        for (size_t si = 0; si < n; ++si)
        {
            // binary search in splitter with equal check
            key_type key = get_char<key_type>(strings[si], depth);

            unsigned int b = Classify<treebits>::find_bkt(key, splitter_tree + 1);

            assert(b < bktnum);

            bktcache[si] = b;
        }

        size_t* bktsize = new size_t[bktnum];
        memset(bktsize, 0, bktnum * sizeof(size_t));

        for (size_t si = 0; si < n; ++si)
            ++bktsize[bktcache[si]];

        if (debug_bucketsize)
        {
            DBG1(1, "bktsize: ");
            for (size_t i = 0; i < bktnum; ++i)
            {
                DBG2(1, bktsize[i] << " ");
            }
            DBG3(1, "");
        }

        // step 3: prefix sum
        g_timer.change(TM_PREFIXSUM);

        size_t bktindex[bktnum];
        bktindex[0] = bktsize[0];
        size_t last_bkt_size = bktsize[0];
        for (unsigned int i = 1; i < bktnum; ++i) {
            bktindex[i] = bktindex[i - 1] + bktsize[i];
            if (bktsize[i]) last_bkt_size = bktsize[i];
        }
        assert(bktindex[bktnum - 1] == n);

        // step 4: premute in-place
        g_timer.change(TM_PERMUTE);

        for (size_t i = 0, j; i < n - last_bkt_size; )
        {
            string perm = strings[i];
            uint16_t permbkt = bktcache[i];

            while ((j = --bktindex[permbkt]) > i)
            {
                std::swap(perm, strings[j]);
                std::swap(permbkt, bktcache[j]);
            }

            strings[i] = perm;
            i += bktsize[permbkt];
        }

        delete[] bktcache;

        // step 5: recursion
        g_timer.change(TM_GENERAL);

        size_t i = 0, bsum = 0;
        while (i < bktnum - 1)
        {
            // i is even -> bkt[i] is less-than bucket
            if (bktsize[i] > 1)
            {
                LOGC(debug_recursion)
                    << "Recurse[" << depth << "]: < bkt " << bsum
                    << " size " << bktsize[i]
                    << " lcp " << int(splitter_lcp[i / 2] & 0x7F);

                sort<Classify>(strings + bsum, bktsize[i],
                               depth + (splitter_lcp[i / 2] & 0x7F));
            }
            bsum += bktsize[i++];

            // i is odd -> bkt[i] is equal bucket
            if (bktsize[i] > 1)
            {
                if (splitter_lcp[i / 2] & 0x80) { // equal-bucket has NULL-terminated key, done.
                    LOGC(debug_recursion)
                        << "Recurse[" << depth << "]: = bkt " << bsum
                        << " size " << bktsize[i] << " is done!";
                }
                else {
                    LOGC(debug_recursion)
                        << "Recurse[" << depth << "]: = bkt " << bsum
                        << " size " << bktsize[i] << " lcp keydepth!";

                    sort<Classify>(strings + bsum, bktsize[i],
                                   depth + sizeof(key_type));
                }
            }
            bsum += bktsize[i++];
        }
        if (bktsize[i] > 0)
        {
            LOGC(debug_recursion)
                << "Recurse[" << depth << "]: > bkt " << bsum
                << " size " << bktsize[i] << " no lcp";

            sort<Classify>(strings + bsum, bktsize[i], depth);
        }
        bsum += bktsize[i++];
        assert(i == bktnum && bsum == n);

        delete[] bktsize;
    }
};

void bingmann_sample_sortBTCE3(string* strings, size_t n)
{
    sample_sort_pre();
    g_stats >> "numsplitters" << numsplitters
        >> "splitter_treebits" << treebits;
    SampleSortBTCE3::sort<ClassifySimple>(strings, n, 0);
    sample_sort_post();
}

PSS_CONTESTANT(bingmann_sample_sortBTCE3, "bingmann/sample_sortBTCE3",
               "bingmann/sample_sortBTCE3 (adapt binary tree equal, bkt cache)")

void bingmann_sample_sortBTCE3A(string* strings, size_t n)
{
    sample_sort_pre();
    g_stats >> "numsplitters" << numsplitters
        >> "splitter_treebits" << treebits;
    SampleSortBTCE3::sort<ClassifyAssembler>(strings, n, 0);
    sample_sort_post();
}

PSS_CONTESTANT(bingmann_sample_sortBTCE3A, "bingmann/sample_sortBTCE3A",
               "bingmann/sample_sortBTCE3A (adapt binary tree equal, asm CMOV, bkt cache)")

void bingmann_sample_sortBTCE3U(string* strings, size_t n)
{
    sample_sort_pre();
    g_stats >> "numsplitters" << numsplitters
        >> "splitter_treebits" << treebits;
    SampleSortBTCE3::sort<ClassifyUnroll>(strings, n, 0);
    sample_sort_post();
}

PSS_CONTESTANT(bingmann_sample_sortBTCE3U, "bingmann/sample_sortBTCE3U",
               "bingmann/sample_sortBTCE3U (adapt binary tree equal unroll, asm CMOV, bkt cache)")

// ----------------------------------------------------------------------------

} // namespace bingmann_sample_sortBTCE

/******************************************************************************/
