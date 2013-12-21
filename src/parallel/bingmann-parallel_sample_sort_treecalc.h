/******************************************************************************
 * src/parallel/bingmann-parallel_sample_sort_treecalc.h
 *
 * Parallel Super Scalar String Sample-Sort, Classifier variant with
 * TreeCalculations instead of separate splitter array.
 *
 * This file is included by bingmann_parallel_sample_sort.cc
 *
 ******************************************************************************
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
 *****************************************************************************/

//! Recursive TreeBuilder for full-descent and unrolled variants, constructs
//! only a binary tree
template <size_t numsplitters>
struct TreeBuilderFullDescentTreeCalc
{
    key_type*       m_tree;
    unsigned char*  m_lcp_iter;
    key_type*       m_samples;

    TreeBuilderFullDescentTreeCalc(key_type* splitter_tree, unsigned char* splitter_lcp,
                                   key_type* samples, size_t samplesize)
        : m_tree( splitter_tree ),
          m_lcp_iter( splitter_lcp ),
          m_samples( samples )
    {

        key_type sentinel = 0;
        recurse(samples, samples + samplesize, 1, sentinel);

        assert(m_lcp_iter == splitter_lcp + numsplitters);
        splitter_lcp[0] &= 0x80; // overwrite sentinel lcp for first < everything bucket
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
            << toHex(*mid));

        key_type mykey = m_tree[treeidx] = *mid;
#if 1
        key_type* midlo = mid;
        while (lo < midlo && *(midlo-1) == mykey) midlo--;

        key_type* midhi = mid;
        while (midhi+1 < hi && *midhi == mykey) midhi++;

        if (midhi - midlo > 1)
            DBG(0, "key range = [" << snum(midlo) << "," << snum(midhi) << ")");
#else
        key_type *midlo = mid, *midhi = mid+1;
#endif
        if (2 * treeidx < numsplitters)
        {
            key_type prevkey = recurse(lo, midlo, 2 * treeidx + 0, rec_prevkey);

            key_type xorSplit = prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " << toHex(prevkey) << " XOR " << toHex(mykey) << " = "
                << toHex(xorSplit) << " - " << count_high_zero_bits(xorSplit) << " bits = "
                << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                | ((mykey & 0xFF) ? 0 : 0x80); // marker for done splitters

            return recurse(midhi, hi, 2 * treeidx + 1, mykey);
        }
        else
        {
            key_type xorSplit = rec_prevkey ^ mykey;

            DBG(debug_splitter, "    lcp: " << toHex(rec_prevkey) << " XOR " << toHex(mykey) << " = "
                << toHex(xorSplit) << " - " << count_high_zero_bits(xorSplit) << " bits = "
                << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

            *m_lcp_iter++ = (count_high_zero_bits(xorSplit) / 8)
                | ((mykey & 0xFF) ? 0 : 0x80); // marker for done splitters

            return mykey;
        }
    }
};

template <size_t treebits>
struct ClassifySimpleTreeCalc
{
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters+1];

    /// binary search on splitter array for bucket number
    inline unsigned int
    find_bkt_tree(const key_type& key) const
    {
        // binary tree traversal without left branch

        unsigned int i = 1;

        while ( i <= numsplitters )
        {
#if 0
            // in gcc-4.6.3 this produces a SETA, LEA sequence
            i = 2 * i + (key <= splitter_tree[i] ? 0 : 1);
#else
            // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which is slightly faster
            if (key <= splitter_tree[i])
                i = 2*i + 0;
            else // (key > splitter_tree[i])
                i = 2*i + 1;
#endif
        }

        i -= numsplitters+1;

        size_t b = i * 2;                                       // < bucket
        if (i < numsplitters && get_splitter(i) == key) b += 1; // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);
            *bktout++ = find_bkt_tree(key);
        }
    }

    //! return a splitter
    inline key_type get_splitter(unsigned int i) const
    {
        return splitter_tree[ TreeCalculations<treebits>::in_to_levelorder(i+1) ];
    }

    /// build tree and splitter array from sample
    inline void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderFullDescentTreeCalc<numsplitters>
            (splitter_tree, splitter_lcp, samples, samplesize);
    }
};

template <size_t treebits>
struct ClassifyUnrollTreeCalc
{
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter_tree[numsplitters+1];

    /// binary search on splitter array for bucket number
    __attribute__((optimize("unroll-all-loops")))
    inline unsigned int
    find_bkt_tree(const key_type& key) const
    {
        // binary tree traversal without left branch

        unsigned int i = 1;

        for (size_t l = 0; l < treebits; ++l)
        {
#if 0
            // in gcc-4.6.3 this produces a SETA, LEA sequence
            i = 2 * i + (key <= splitter_tree[i] ? 0 : 1);
#else
            // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which is slightly faster
            if (key <= splitter_tree[i])
                i = 2*i + 0;
            else // (key > splitter_tree[i])
                i = 2*i + 1;
#endif
        }

        i -= numsplitters+1;

        size_t b = i * 2;                                       // < bucket
        if (i < numsplitters && get_splitter(i) == key) b += 1; // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);

            *bktout++ = find_bkt_tree(key);
        }
    }

    //! return a splitter
    inline key_type get_splitter(unsigned int i) const
    {
        return splitter_tree[ TreeCalculations<treebits>::in_to_levelorder(i+1) ];
    }

    /// build tree and splitter array from sample
    inline void
    build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        TreeBuilderFullDescentTreeCalc<numsplitters>
            (splitter_tree, splitter_lcp, samples, samplesize);
    }
};

template <size_t treebits>
struct ClassifyUnrollBothTreeCalc : public ClassifyUnrollTreeCalc<treebits>
{
   /// search in splitter tree for bucket number, unrolled for U keys at once.
    template <int U>
    __attribute__((optimize("unroll-all-loops")))
    inline void
    find_bkt_tree_unroll(const key_type key[U], uint16_t obkt[U]) const
    {
        // binary tree traversal without left branch

        static const size_t numsplitters = this->numsplitters;
        const key_type* splitter_tree = this->splitter_tree;

        unsigned int i[U];
        std::fill(i+0, i+U, 1);

        for (size_t l = 0; l < treebits; ++l)
        {
            for (int u = 0; u < U; ++u)
            {
#if 0
                // in gcc-4.6.3 this produces a SETA, LEA sequence
                i[u] = 2 * i[u] + (key[u] <= splitter_tree[i[u]] ? 0 : 1);
#else
                // in gcc-4.6.3 this produces two LEA and a CMOV sequence, which is slightly faster
                if (key[u] <= splitter_tree[i[u]])
                    i[u] = 2*i[u] + 0;
                else // (key > splitter_tree[i[u]])
                    i[u] = 2*i[u] + 1;
#endif
            }
        }

        for (int u = 0; u < U; ++u)
            i[u] -= numsplitters+1;

        for (int u = 0; u < U; ++u)
            obkt[u] = i[u] * 2; // < bucket

        for (int u = 0; u < U; ++u)
        {
            if (i[u] < numsplitters && this->get_splitter(i[u]) == key[u]) obkt[u] += 1; // equal bucket
        }
    }

    /// classify all strings in area by walking tree and saving bucket id, unrolled loops
    inline void
    classify(string* strB, string* strE, uint16_t* bktout, size_t depth) const
    {
        for (string* str = strB; str != strE; )
        {
            static const int rollout = 4;
            if (str + rollout < strE)
            {
                key_type key[rollout];
                key[0] = get_char<key_type>(str[0], depth);
                key[1] = get_char<key_type>(str[1], depth);
                key[2] = get_char<key_type>(str[2], depth);
                key[3] = get_char<key_type>(str[3], depth);

                find_bkt_tree_unroll<rollout>(key, bktout);

                str += rollout;
                bktout += rollout;
            }
            else
            {
                // binary search in splitter with equal check
                key_type key = get_char<key_type>(*str++, depth);

                *bktout++ = this->find_bkt_tree(key);
            }
        }
    }
};

void parallel_sample_sortBTCT(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifySimpleTreeCalc>(strings, n, 0);
}

CONTESTANT_REGISTER_PARALLEL_LCP(
    parallel_sample_sortBTCT,
    "bingmann/parallel_sample_sortBTCT",
    "bingmann/parallel_sample_sortBTCT: binary tree, bktcache, tree calc")

void parallel_sample_sortBTCTU1(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyUnrollTreeCalc>(strings, n, 0);
}

CONTESTANT_REGISTER_PARALLEL_LCP(
    parallel_sample_sortBTCTU1,
    "bingmann/parallel_sample_sortBTCTU1",
    "bingmann/parallel_sample_sortBTCTU1: binary tree, bktcache, unroll tree, tree calc")

void parallel_sample_sortBTCTU2(string* strings, size_t n)
{
    parallel_sample_sort_base<ClassifyUnrollBothTreeCalc>(strings, n, 0);
}

CONTESTANT_REGISTER_PARALLEL_LCP(
    parallel_sample_sortBTCTU2,
    "bingmann/parallel_sample_sortBTCTU2",
    "bingmann/parallel_sample_sortBTCTU2: binary tree, bktcache, unroll tree and strings, tree calc")
