/*******************************************************************************
 * src/sequential/bingmann-sample_sortBTCE.hpp
 *
 * (Parallel) Super Scalar String Sample-Sort, Classifier variant with equal
 * checks at each node.
 *
 * This file is included by bingmann_parallel_sample_sort.cc
 *
 *******************************************************************************
 * Copyright (C) 2013-2017 Timo Bingmann <tb@panthema.net>
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

#ifndef PSS_SRC_SEQUENTIAL_BINGMANN_SAMPLE_SORTBTCE_HEADER
#define PSS_SRC_SEQUENTIAL_BINGMANN_SAMPLE_SORTBTCE_HEADER

#include "bingmann-sample_sort.hpp"

namespace bingmann_sample_sort {

template <size_t TreeBits>
class ClassifyEqual
{
public:
    static const size_t treebits = TreeBits;
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters + 1];

    /// binary search on splitter array for bucket number
    unsigned int find_bkt(const key_type& key) const
    {
        // binary tree traversal without left branch

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

    //! classify all strings in area by walking tree and saving bucket id
    void classify(string* strB, string* strE, uint16_t* bktout,
                  size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);
            *bktout++ = find_bkt(key);
        }
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
            *bktout++ = find_bkt(key);
        }
    }

    //! return a splitter
    key_type get_splitter(unsigned int i) const
    {
        return splitter_tree[
            TreeCalculations < treebits > ::pre_to_levelorder(i)];
    }

    /// build tree and splitter array from sample
    void build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderLevelOrder<numsplitters>(
            splitter_tree, splitter_lcp, samples, samplesize);
    }
};

template <size_t TreeBits>
class ClassifyEqualAssembler
{
public:
    static const size_t treebits = TreeBits;
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters + 1];

    /// binary search on splitter array for bucket number
    unsigned int find_bkt(const key_type& key) const
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

    //! classify all strings in area by walking tree and saving bucket id
    void classify(string* strB, string* strE, uint16_t* bktout,
                  size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);
            *bktout++ = find_bkt(key);
        }
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
            *bktout++ = find_bkt(key);
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
        TreeBuilderLevelOrder<numsplitters>(
            splitter_tree, splitter_lcp, samples, samplesize);
    }
};

template <size_t TreeBits>
class ClassifyEqualUnrollTree
{
public:
    static const size_t treebits = TreeBits;
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters + 1];

    /// specialized implementation of this find_bkt are below
    unsigned int
    find_bkt(const key_type& key) const;

    //! classify all strings in area by walking tree and saving bucket id
    void classify(string* strB, string* strE, uint16_t* bktout,
                  size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);
            *bktout++ = find_bkt(key);
        }
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
            *bktout++ = find_bkt(key);
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
        TreeBuilderLevelOrder<numsplitters>(
            splitter_tree, splitter_lcp, samples, samplesize);
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
ClassifyEqualUnrollTree<3>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<4>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<5>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<6>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<7>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<8>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<9>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<10>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<11>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<12>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<13>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<14>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<15>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<16>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<17>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<18>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<19>::find_bkt(const key_type& key) const
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
ClassifyEqualUnrollTree<20>::find_bkt(const key_type& key) const
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

} // namespace bingmann_sample_sort

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_SAMPLE_SORTBTCE_HEADER

/******************************************************************************/
