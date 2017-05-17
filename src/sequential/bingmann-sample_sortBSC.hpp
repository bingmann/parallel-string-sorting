/*******************************************************************************
 * src/sequential/bingmann-sample_sortBSC.hpp
 *
 * Experiments with sequential Super Scalar String Sample-Sort (S^5).
 *
 * Binary tree search with bucket cache.
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

#ifndef PSS_SRC_SEQUENTIAL_BINGMANN_SAMPLE_SORTBSC_HEADER
#define PSS_SRC_SEQUENTIAL_BINGMANN_SAMPLE_SORTBSC_HEADER

#include "bingmann-sample_sort.hpp"
#include "bingmann-sample_sort_tree_builder.hpp"

namespace bingmann_sample_sort {

// ****************************************************************************
// *** Classification with Binary Search in Splitter Array (old BSC variant)

template <size_t treebits = DefaultTreebits>
class ClassifyBinarySearch
{
public:
    // NOTE: for binary search numsplitters need not be 2^k-1, any size will
    // do, but the tree implementations are always faster, so we keep this only
    // for historical reasons.
    static const size_t numsplitters = (1 << treebits) - 1;

    key_type splitter[numsplitters];

    /// binary search on splitter array for bucket number
    unsigned int find_bkt_tree(const key_type& key) const
    {
        unsigned int lo = 0, hi = numsplitters;

        while (lo < hi)
        {
            unsigned int mid = (lo + hi) >> 1;
            assert(mid < numsplitters);

            if (key <= splitter[mid])
                hi = mid;
            else                                               // (key > splitter[mid])
                lo = mid + 1;
        }

        size_t b = lo * 2;                                     // < bucket
        if (lo < numsplitters && splitter[lo] == key) b += 1;  // equal bucket

        return b;
    }

    /// classify all strings in area by walking tree and saving bucket id
    void classify(string* strB, string* strE, uint16_t* bktout,
                  size_t depth)
    {
        for (string* str = strB; str != strE; )
        {
            key_type key = get_char<key_type>(*str++, depth);
            *bktout++ = find_bkt_tree(key);
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
            *bktout++ = find_bkt_tree(key);
        }
    }

    //! return a splitter
    key_type get_splitter(unsigned int i) const
    { return splitter[i]; }

    /// build tree and splitter array from sample
    void build(key_type* samples, size_t samplesize, unsigned char* splitter_lcp)
    {
        bingmann_sample_sort::TreeBuilderPreorder<numsplitters>(
            splitter, splitter_lcp, samples, samplesize);
    }
};

} // namespace bingmann_sample_sort

#endif // !PSS_SRC_SEQUENTIAL_BINGMANN_SAMPLE_SORTBSC_HEADER

/******************************************************************************/
