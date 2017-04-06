/*******************************************************************************
 * src/parallel/bingmann-parallel_sample_sort_equal.hpp
 *
 * Parallel Super Scalar String Sample-Sort, Classifier variant with equal
 * checks at each node.
 *
 * This file is included by bingmann_parallel_sample_sort.cc
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

#ifndef PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_EQUAL_HEADER
#define PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_EQUAL_HEADER

template <size_t numsplitters>
class TreeBuilderEqual
{
public:
    key_type* m_tree;
    unsigned char* m_lcp_iter;
    key_type* m_samples;

    TreeBuilderEqual(key_type* splitter_tree, unsigned char* splitter_lcp,
                     key_type* samples, size_t samplesize)
        : m_tree(splitter_tree),
          m_lcp_iter(splitter_lcp),
          m_samples(samples)
    {
        key_type sentinel = 0;
        recurse(samples, samples + samplesize, 1, sentinel);

        assert(m_lcp_iter == splitter_lcp + numsplitters);
        splitter_lcp[0] &= 0x80;        // overwrite sentinel lcp for first < everything bucket
        splitter_lcp[numsplitters] = 0; // sentinel for > everything bucket
    }

    ptrdiff_t snum(key_type* s) const
    {
        return (ptrdiff_t)(s - m_samples);
    }

    key_type recurse(key_type* lo, key_type* hi, unsigned int treeidx,
                     key_type& rec_prevkey)
    {
        DBG(debug_splitter, "rec_buildtree(" << snum(lo) << "," << snum(hi)
                                             << ", treeidx=" << treeidx << ")");

        // pick middle element as splitter
        key_type* mid = lo + (ptrdiff_t)(hi - lo) / 2;

        DBG(debug_splitter, "tree[" << treeidx << "] = samples[" << snum(mid) << "] = "
                                    << tlx::hexdump_type(*mid));

        key_type mykey = m_tree[treeidx] = *mid;
#if 1
        key_type* midlo = mid;
        while (lo < midlo && *(midlo - 1) == mykey) midlo--;

        key_type* midhi = mid;
        while (midhi + 1 < hi && *midhi == mykey) midhi++;

        if (midhi - midlo > 1)
            DBG(0, "key range = [" << snum(midlo) << "," << snum(midhi) << ")");
#else
        key_type* midlo = mid, * midhi = mid + 1;
#endif
        if (2 * treeidx < numsplitters)
        {
            key_type prevkey = recurse(lo, midlo, 2 * treeidx + 0, rec_prevkey);

            key_type xorSplit = prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " <<
                tlx::hexdump_type(prevkey) << " XOR " << tlx::hexdump_type(mykey) << " = " <<
                tlx::hexdump_type(xorSplit) << " - " << count_high_zero_bits(xorSplit) <<
                " bits = " << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            // marker for done splitters
            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                            | ((mykey & 0xFF) ? 0 : 0x80);

            return recurse(midhi, hi, 2 * treeidx + 1, mykey);
        }
        else
        {
            key_type xorSplit = rec_prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " <<
                tlx::hexdump_type(rec_prevkey) << " XOR " << tlx::hexdump_type(mykey) << " = " <<
                tlx::hexdump_type(xorSplit) << " - " << count_high_zero_bits(xorSplit) <<
                " bits = " << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            // marker for done splitters
            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                            | ((mykey & 0xFF) ? 0 : 0x80);

            return mykey;
        }
    }
};

template <size_t treebits>
class ClassifyEqual
{
public:
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters + 1];

    /// binary search on splitter array for bucket number
    unsigned int find_bkt_tree(const key_type& key) const
    {
        unsigned int i;

#if 1
        // hand-coded assembler binary tree traversal with equality, using CMOV
        asm ("mov    $1, %%rax \n"             // rax = i
             // body of while loop
             "1: \n"
             "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
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
#else
        // hand-coded assembler binary tree traversal with equality, using SETA
        // this version is slightly slower than the CMOV above.
        asm ("mov    $1, %%rax \n"               // rax = i
             // body of while loop
             "1: \n"
             "cmpq   (%[splitter_tree],%%rax,8), %[key] \n"
             "seta   %%cl \n"                    // cl = 1 if key > splitter_tree
             "movzb  %%cl, %%rcx \n"             // pad cl with zeros
             "je     2f \n"
             "lea    (%%rcx,%%rax,2), %%rax \n"  // rax += 2 * i + 0/1
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
#endif
        return i;
    }

    /// classify all strings in area by walking tree and saving bucket id
    template <typename StringSet>
    void classify(
        const StringSet& strset,
        typename StringSet::Iterator begin, typename StringSet::Iterator end,
        uint16_t* bktout, size_t depth) const
    {
        while (begin != end)
        {
            key_type key = strset.get_uint64(*begin++, depth);
            *bktout++ = find_bkt_tree(key);
        }
    }

    //! return a splitter
    key_type get_splitter(unsigned int i) const
    {
        return splitter_tree[TreeCalculations < treebits > ::pre_to_levelorder(i)];
    }

    /// build tree and splitter array from sample
    void build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderEqual<numsplitters>(splitter_tree, splitter_lcp,
                                       samples, samplesize);
    }
};

template <size_t TreeBits>
class ClassifyEqualUnrollTree
{
public:
    static const size_t treebits = TreeBits;
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters + 1];

    /// specialized implementation of this find_bkt_tree are below
    unsigned int
    find_bkt_tree(const key_type& key) const;

    /// classify all strings in area by walking tree and saving bucket id
    template <typename StringSet>
    void classify(
        const StringSet& strset,
        typename StringSet::Iterator begin, typename StringSet::Iterator end,
        uint16_t* bktout, size_t depth) const
    {
        while (begin != end)
        {
            key_type key = strset.get_uint64(*begin++, depth);
            *bktout++ = find_bkt_tree(key);
        }
    }

    //! return a splitter
    key_type get_splitter(unsigned int i) const
    {
        return splitter_tree[TreeCalculations < treebits > ::pre_to_levelorder(i)];
    }

    /// build tree and splitter array from sample
    void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderEqual<numsplitters>(splitter_tree, splitter_lcp,
                                       samples, samplesize);
    }

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
inline unsigned int
ClassifyEqualUnrollTree<3>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm ("mov    $1, %%rax \n"             // rax = i
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<4>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm ("mov    $1, %%rax \n"             // rax = i
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<5>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm ("mov    $1, %%rax \n"             // rax = i
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<6>::find_bkt_tree(const key_type& key) const
{
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality, using CMOV
    asm ("mov    $1, %%rax \n"             // rax = i
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<7>::find_bkt_tree(const key_type& key) const
{
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

         SPLITTER_TREE_END

         : "=&a" (i)
         :[key] "r" (key), [splitter_tree] "r" (splitter_tree),
         [numsplitters1] "g" (numsplitters + 1),
         [treebits] "g" (treebits),
         [numsplitters] "g" (numsplitters)
         : "rcx", "rdx");

    return i;
}

template <>
inline unsigned int
ClassifyEqualUnrollTree<8>::find_bkt_tree(const key_type& key) const
{
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

         SPLITTER_TREE_END

         : "=&a" (i)
         :[key] "r" (key), [splitter_tree] "r" (splitter_tree),
         [numsplitters1] "g" (numsplitters + 1),
         [treebits] "g" (treebits),
         [numsplitters] "g" (numsplitters)
         : "rcx", "rdx");

    return i;
}

template <>
inline unsigned int
ClassifyEqualUnrollTree<9>::find_bkt_tree(const key_type& key) const
{
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

         SPLITTER_TREE_END

         : "=&a" (i)
         :[key] "r" (key), [splitter_tree] "r" (splitter_tree),
         [numsplitters1] "g" (numsplitters + 1),
         [treebits] "g" (treebits),
         [numsplitters] "g" (numsplitters)
         : "rcx", "rdx");

    return i;
}

template <>
inline unsigned int
ClassifyEqualUnrollTree<10>::find_bkt_tree(const key_type& key) const
{
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

         SPLITTER_TREE_END

         : "=&a" (i)
         :[key] "r" (key), [splitter_tree] "r" (splitter_tree),
         [numsplitters1] "g" (numsplitters + 1),
         [treebits] "g" (treebits),
         [numsplitters] "g" (numsplitters)
         : "rcx", "rdx");

    return i;
}

template <>
inline unsigned int
ClassifyEqualUnrollTree<11>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<12>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<13>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<14>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<15>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<16>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<17>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<18>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<19>::find_bkt_tree(const key_type& key) const
{
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

template <>
inline unsigned int
ClassifyEqualUnrollTree<20>::find_bkt_tree(const key_type& key) const
{
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

#undef SPLITTER_TREE_STEP
#undef SPLITTER_TREE_END

static inline void
parallel_sample_sortBTCE(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyEqual>(strings, n, 0);
}

PSS_CONTESTANT_PARALLEL_LCP(
    parallel_sample_sortBTCE,
    "bingmann/parallel_sample_sortBTCE",
    "pS5: binary tree, equality, bktcache")

static inline void
parallel_sample_sortBTCEU1(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyEqualUnrollTree>(strings, n, 0);
}

PSS_CONTESTANT_PARALLEL_LCP(
    parallel_sample_sortBTCEU1,
    "bingmann/parallel_sample_sortBTCEU1",
    "pS5: binary tree, equality, bktcache, unroll tree")

#endif // !PSS_SRC_PARALLEL_BINGMANN_PARALLEL_SAMPLE_SORT_EQUAL_HEADER

/******************************************************************************/
