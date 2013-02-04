/******************************************************************************
 * src/sequential/bingmann-sample_sort.h
 *
 * Experiments with sequential Super Scalar String Sample-Sort (S^5).
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

namespace bingmann_sample_sort {

static const bool debug_splitter = false;
static const bool debug_bucketsize = false;
static const bool debug_recursion = false;
static const bool debug_splitter_tree = false;

using namespace stringtools;

typedef uint64_t key_type;

static const size_t l2cache = 256*1024;

/// binary search on splitter array for bucket number
inline unsigned int find_bkt(const key_type& key, const key_type* splitter, size_t leaves)
{
    unsigned int lo = 0, hi = leaves;

    while ( lo < hi )
    {
        unsigned int mid = (lo + hi) >> 1;
        assert(mid < leaves);

        if (key <= splitter[mid])
            hi = mid;
        else // (key > splitter[mid])
            lo = mid + 1;
    }

#if 0
    // Verify result of binary search:
    int pos = leaves-1;
    while ( pos >= 0 && key <= splitter[pos] ) --pos;
    pos++;

    //std::cout << "lo " << lo << " hi " << hi << " pos " << pos << "\n";
    assert(lo == pos);
#endif

    size_t b = lo * 2;                              // < bucket
    if (lo < leaves && splitter[lo] == key) b += 1; // equal bucket

    return b;
}

/// Variant 1 of string sample-sort: use binary search on splitters, no caching.
void sample_sortBS(string* strings, size_t n, size_t depth)
{
#if 0
    static const size_t leaves = 32;
#else
    //static const size_t l2cache = 256*1024;

    // bounding equations:
    // splitters            + bktsize
    // n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    static const size_t leaves = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));

#endif

    if (n < 1024*1024)
    {
        return bingmann_radix_sort::msd_CI5(strings, n, depth);
    }

    //std::cout << "leaves: " << leaves << "\n";

    // step 1: select splitters with oversampling

    const size_t oversample_factor = 4;
    size_t samplesize = oversample_factor * leaves;

    key_type* samples = new key_type[ samplesize ];

    LCGRandom rng(9384234);

    for (unsigned int i = 0; i < samplesize; ++i)
    {
        samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
    }

    std::sort(samples, samples + samplesize);

    key_type splitter[leaves];
    unsigned char splitter_lcp[leaves];

    DBG(debug_splitter, "splitter:");
    splitter_lcp[0] = 0; // sentinel for first < everything bucket
    for (size_t i = 0, j = oversample_factor/2; i < leaves; ++i)
    {
        splitter[i] = samples[j];
        DBG(debug_splitter, "key " << toHex(splitter[i]));

        if (i != 0) {
            key_type andSplit = splitter[i-1] ^ splitter[i];

            DBG1(debug_splitter, "    XOR -> " << toHex(andSplit) << " - ");

            DBG3(debug_splitter, count_high_zero_bits(andSplit) << " bits = "
                << count_high_zero_bits(andSplit) / 8 << " chars lcp");

            splitter_lcp[i] = count_high_zero_bits(andSplit) / 8;
        }

        j += oversample_factor;
    }

    delete [] samples;

    // step 2: classify all strings and count bucket sizes

    static const size_t bktnum = 2*leaves+1;

    size_t bktsize[2*leaves+1] = { 0 };

    for (size_t si = 0; si < n; ++si)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(strings[si], depth);

        unsigned int b = find_bkt(key, splitter, leaves);

        assert(b < bktnum);

        ++bktsize[ b ];
    }

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

    size_t bktindex[bktnum];
    bktindex[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned int i=1; i < bktnum; ++i) {
        bktindex[i] = bktindex[i-1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }
    assert(bktindex[bktnum-1] == n);

    // step 4: premute in-place

    for (size_t i=0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        key_type key;
        unsigned int b;

        while (1)
        {
            key = get_char<key_type>(perm, depth);
            b = find_bkt(key, splitter, leaves);

            j = --bktindex[ b ];

            if ( j <= i )
                break;

            std::swap(perm, strings[j]);
        }

        strings[i] = perm;
        i += bktsize[ b ];
    }

    // step 5: recursion

    size_t bsum = 0;
    for (size_t i=0; i < bktnum-1; ++i) {
        // i is even -> bkt[i] is less-than bucket
        if (bktsize[i] > 1)
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bsum << " size " << bktsize[i] << " lcp " << int(splitter_lcp[i/2]));
            sample_sortBS(strings+bsum, bktsize[i], depth + splitter_lcp[i/2]);
        }
        bsum += bktsize[i];
        ++i;
        // i is odd -> bkt[i] is equal bucket
        if (bktsize[i] > 1)
        {
            if ( (splitter[i/2] & 0xFF) == 0 ) { // equal-bucket has NULL-terminated key
                // done.
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " is done!");
            }
            else {
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " lcp keydepth!");
                sample_sortBS(strings+bsum, bktsize[i], depth + sizeof(key_type));
            }
        }
        bsum += bktsize[i];
    }
    if (bktsize[bktnum-1] > 0)
    {
        DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bsum << " size " << bktsize[bktnum-1] << " no lcp");
        sample_sortBS(strings+bsum, bktsize[bktnum-1], depth);
    }
    bsum += bktsize[bktnum-1];
    assert(bsum == n);
}

void bingmann_sample_sortBS(string* strings, size_t n) { return sample_sortBS(strings,n,0); }
CONTESTANT_REGISTER_UCARRAY(bingmann_sample_sortBS, "bingmann/sample_sortBS (binary search, no cache)")

// ------------------------------------------------------------------------------------------------------------------------

/// Variant 2 of string sample-sort: use binary search on splitters, with index caching.
void sample_sortBSC(string* strings, size_t n, size_t depth)
{
#if 0
    static const size_t leaves = 32;
#else
    //static const size_t l2cache = 256*1024;

    // bounding equations:
    // splitters            + bktsize
    // n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    static const size_t leaves = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));

#endif

    if (n < 1024*1024)
    {
        return bingmann_radix_sort::msd_CI5(strings, n, depth);
    }

    // step 1: select splitters with oversampling

    const size_t oversample_factor = 4;
    size_t samplesize = oversample_factor * leaves;

    key_type* samples = new key_type[ samplesize ];

    LCGRandom rng(9384234);

    for (unsigned int i = 0; i < samplesize; ++i)
    {
        samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
    }

    std::sort(samples, samples + samplesize);

    key_type splitter[leaves];
    unsigned char splitter_lcp[leaves];

    DBG(debug_splitter, "splitter:");
    splitter_lcp[0] = 0; // sentinel for first < everything bucket
    for (size_t i = 0, j = oversample_factor/2; i < leaves; ++i)
    {
        splitter[i] = samples[j];
        DBG(debug_splitter, "key " << toHex(splitter[i]));

        if (i != 0) {
            key_type andSplit = splitter[i-1] ^ splitter[i];

            DBG1(debug_splitter, "    XOR -> " << toHex(andSplit) << " - ");

            DBG3(debug_splitter, count_high_zero_bits(andSplit) << " bits = "
                << count_high_zero_bits(andSplit) / 8 << " chars lcp");

            splitter_lcp[i] = count_high_zero_bits(andSplit) / 8;
        }

        j += oversample_factor;
    }

    delete [] samples;

    // step 2: classify all strings and count bucket sizes

    static const size_t bktnum = 2*leaves+1;

    size_t bktsize[2*leaves+1] = { 0 };

    uint16_t* bktcache = new uint16_t[n];

    for (size_t si = 0; si < n; ++si)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(strings[si], depth);

        unsigned int b = find_bkt(key, splitter, leaves);

        assert(b < bktnum);

        bktcache[si] = b;
        ++bktsize[ b ];
    }

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

    size_t bktindex[bktnum];
    bktindex[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned int i=1; i < bktnum; ++i) {
        bktindex[i] = bktindex[i-1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }
    assert(bktindex[bktnum-1] == n);

    // step 4: premute in-place

    for (size_t i=0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        uint16_t permbkt = bktcache[i];

        while ( (j = --bktindex[ permbkt ]) > i )
        {
            std::swap(perm, strings[j]);
            std::swap(permbkt, bktcache[j]);
        }

        strings[i] = perm;
        i += bktsize[ permbkt ];
    }

    delete [] bktcache;

    // step 5: recursion

    size_t bsum = 0;
    for (size_t i=0; i < bktnum-1; ++i) {
        // i is even -> bkt[i] is less-than bucket
        if (bktsize[i] > 1)
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bsum << " size " << bktsize[i] << " lcp " << int(splitter_lcp[i/2]));
            sample_sortBSC(strings+bsum, bktsize[i], depth + splitter_lcp[i/2]);
        }
        bsum += bktsize[i];
        ++i;
        // i is odd -> bkt[i] is equal bucket
        if (bktsize[i] > 1)
        {
            if ( (splitter[i/2] & 0xFF) == 0 ) { // equal-bucket has NULL-terminated key
                // done.
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " is done!");
            }
            else {
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " lcp keydepth!");
                sample_sortBSC(strings+bsum, bktsize[i], depth + sizeof(key_type));
            }
        }
        bsum += bktsize[i];
    }
    if (bktsize[bktnum-1] > 0)
    {
        DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bsum << " size " << bktsize[bktnum-1] << " no lcp");
        sample_sortBSC(strings+bsum, bktsize[bktnum-1], depth);
    }
    bsum += bktsize[bktnum-1];
    assert(bsum == n);
}

void bingmann_sample_sortBSC(string* strings, size_t n) { return sample_sortBSC(strings,n,0); }
CONTESTANT_REGISTER_UCARRAY(bingmann_sample_sortBSC, "bingmann/sample_sortBSC (binary search, bkt cache)")

// ------------------------------------------------------------------------------------------------------------------------

/// binary search on splitter array for bucket number
inline unsigned int find_bkt_splittertree(const key_type& key, const key_type* splitter, const key_type* splitter_tree0, size_t numsplitters)
{
#if 1
    // binary tree traversal without left branch

    const key_type* splitter_tree = splitter_tree0 - 1;
    unsigned int i = 1;

    while ( i <= numsplitters )
    {
        if (key <= splitter_tree[i])
            i = 2*i + 0;
        else // (key > splitter_tree[i])
            i = 2*i + 1;
    }

    i -= numsplitters+1;

    size_t b = i * 2;                                   // < bucket
    if (i < numsplitters && splitter[i] == key) b += 1; // equal bucket

#else
    // binary search variant with keeping the last left branch
    const key_type* splitter_tree = splitter_tree0 - 1;
    unsigned int i = 1;
    unsigned int ll = 1;        // last left branch taken

    while ( i <= numsplitters )
    {
        if (key <= splitter_tree[i]) {
            ll = i;
            i = 2*i + 0;
        }
        else // (key > splitter_tree[i])
            i = 2*i + 1;
    }

    i -= numsplitters+1;

#if 0
    // Verify result of binary search:
    int pos = numsplitters-1;
    while ( pos >= 0 && key <= splitter[pos] ) --pos;
    pos++;

    std::cout << "i " << i << " pos " << pos << "\n";
    assert(i == pos);
#endif

    assert(i >= numsplitters || splitter_tree[ll] == splitter[i]);

    size_t b = i * 2;                                   // < bucket
    if (i < numsplitters && splitter_tree[ll] == key) b += 1; // equal bucket

#endif

    return b;
}

/// Variant 3 of string sample-sort: use super-scalar binary search on splitters, without index caching.
void sample_sortBT(string* strings, size_t n, size_t depth)
{
#if 0
    static const size_t numsplitters = 31;
#else
    //static const size_t l2cache = 256*1024;

    // bounding equations:
    // splitters            + bktsize
    // n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));
    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / ( sizeof(key_type) );

    static const size_t numsplitters = (1 << logfloor_<numsplitters2>::value) - 1;
#endif

    if (depth != 0 && n < 1024*1024)
    {
        return bingmann_radix_sort::msd_CI5(strings, n, depth);
    }

    //std::cout << "numsplitters: " << numsplitters << "\n";

    // step 1: select splitters with oversampling

    const size_t oversample_factor = 4;
    size_t samplesize = oversample_factor * numsplitters;

    key_type* samples = new key_type[ samplesize ];

    LCGRandom rng(9384234);

    for (unsigned int i = 0; i < samplesize; ++i)
    {
        samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
    }

    std::sort(samples, samples + samplesize);

    key_type splitter[numsplitters];
    unsigned char splitter_lcp[numsplitters];

    DBG(debug_splitter, "splitter:");
    splitter_lcp[0] = 0; // sentinel for first < everything bucket
    for (size_t i = 0, j = oversample_factor/2; i < numsplitters; ++i)
    {
        splitter[i] = samples[j];
        DBG(debug_splitter, "key " << toHex(splitter[i]));

        if (i != 0) {
            key_type andSplit = splitter[i-1] ^ splitter[i];

            DBG1(debug_splitter, "    XOR -> " << toHex(andSplit) << " - ");

            DBG3(debug_splitter, count_high_zero_bits(andSplit) << " bits = "
                << count_high_zero_bits(andSplit) / 8 << " chars lcp");

            splitter_lcp[i] = count_high_zero_bits(andSplit) / 8;
        }

        j += oversample_factor;
    }

    delete [] samples;

    // step 2.1: construct splitter tree to perform binary search

    key_type splitter_tree[numsplitters];

    {
        size_t t = 0;
        size_t highbit = (numsplitters+1) / 2;

        while ( highbit > 0 )
        {
            DBG(debug_splitter_tree, "highbit = " << highbit);

            size_t p = highbit-1;
            size_t inc = highbit << 1;

            while ( p <= numsplitters )
            {
                DBG(debug_splitter_tree, "p = " << p);

                splitter_tree[t++] = splitter[p];

                p += inc;
            }

            highbit >>= 1;
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

    static const size_t bktnum = 2*numsplitters+1;

    size_t bktsize[bktnum] = { 0 };

    for (size_t si = 0; si < n; ++si)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(strings[si], depth);

        unsigned int b = find_bkt_splittertree(key, splitter, splitter_tree, numsplitters);

        assert(b < bktnum);

        ++bktsize[ b ];
    }

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

    size_t bktindex[bktnum];
    bktindex[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned int i=1; i < bktnum; ++i) {
        bktindex[i] = bktindex[i-1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }
    assert(bktindex[bktnum-1] == n);

    // step 4: premute in-place

    for (size_t i=0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        key_type key;
        unsigned int b;

        while (1)
        {
            key = get_char<key_type>(perm, depth);
            b = find_bkt_splittertree(key, splitter, splitter_tree, numsplitters);

            j = --bktindex[ b ];

            if ( j <= i )
                break;

            std::swap(perm, strings[j]);
        }

        strings[i] = perm;
        i += bktsize[ b ];
    }

    // step 5: recursion

    size_t bsum = 0;
    for (size_t i=0; i < bktnum-1; ++i) {
        // i is even -> bkt[i] is less-than bucket
        if (bktsize[i] > 1)
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bsum << " size " << bktsize[i] << " lcp " << int(splitter_lcp[i/2]));
            sample_sortBT(strings+bsum, bktsize[i], depth + splitter_lcp[i/2]);
        }
        bsum += bktsize[i];
        ++i;
        // i is odd -> bkt[i] is equal bucket
        if (bktsize[i] > 1)
        {
            if ( (splitter[i/2] & 0xFF) == 0 ) { // equal-bucket has NULL-terminated key
                // done.
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " is done!");
            }
            else {
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " lcp keydepth!");
                sample_sortBT(strings+bsum, bktsize[i], depth + sizeof(key_type));
            }
        }
        bsum += bktsize[i];
    }
    if (bktsize[bktnum-1] > 0)
    {
        DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bsum << " size " << bktsize[bktnum-1] << " no lcp");
        sample_sortBT(strings+bsum, bktsize[bktnum-1], depth);
    }
    bsum += bktsize[bktnum-1];
    assert(bsum == n);
}

void bingmann_sample_sortBT(string* strings, size_t n) { return sample_sortBT(strings,n,0); }
CONTESTANT_REGISTER_UCARRAY(bingmann_sample_sortBT, "bingmann/sample_sortBT (binary tree, no cache)")

// ------------------------------------------------------------------------------------------------------------------------

/// Variant 4 of string sample-sort: use super-scalar binary search on splitters, without index caching.
void sample_sortBTC(string* strings, size_t n, size_t depth)
{
#if 0
    static const size_t numsplitters = 31;
#else
    //static const size_t l2cache = 256*1024;

    // bounding equations:
    // splitters            + bktsize
    // n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));
    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / ( sizeof(key_type) );

    static const size_t numsplitters = (1 << logfloor_<numsplitters2>::value) - 1;
#endif

    if (depth != 0 && n < 1024*1024)
    {
        return bingmann_radix_sort::msd_CI5(strings, n, depth);
    }

    //std::cout << "numsplitters: " << numsplitters << "\n";

    // step 1: select splitters with oversampling

    const size_t oversample_factor = 4;
    size_t samplesize = oversample_factor * numsplitters;

    key_type* samples = new key_type[ samplesize ];

    LCGRandom rng(9384234);

    for (unsigned int i = 0; i < samplesize; ++i)
    {
        samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
    }

    std::sort(samples, samples + samplesize);

    key_type splitter[numsplitters];
    unsigned char splitter_lcp[numsplitters];

    DBG(debug_splitter, "splitter:");
    splitter_lcp[0] = 0; // sentinel for first < everything bucket
    for (size_t i = 0, j = oversample_factor/2; i < numsplitters; ++i)
    {
        splitter[i] = samples[j];
        DBG(debug_splitter, "key " << toHex(splitter[i]));

        if (i != 0) {
            key_type andSplit = splitter[i-1] ^ splitter[i];

            DBG1(debug_splitter, "    XOR -> " << toHex(andSplit) << " - ");

            DBG3(debug_splitter, count_high_zero_bits(andSplit) << " bits = "
                << count_high_zero_bits(andSplit) / 8 << " chars lcp");

            splitter_lcp[i] = count_high_zero_bits(andSplit) / 8;
        }

        j += oversample_factor;
    }

    delete [] samples;

    // step 2.1: construct splitter tree to perform binary search

    key_type splitter_tree[numsplitters];

    {
        size_t t = 0;
        size_t highbit = (numsplitters+1) / 2;

        while ( highbit > 0 )
        {
            DBG(debug_splitter_tree, "highbit = " << highbit);

            size_t p = highbit-1;
            size_t inc = highbit << 1;

            while ( p <= numsplitters )
            {
                DBG(debug_splitter_tree, "p = " << p);

                splitter_tree[t++] = splitter[p];

                p += inc;
            }

            highbit >>= 1;
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

    uint16_t* bktcache = new uint16_t[n];

    static const size_t bktnum = 2*numsplitters+1;

    size_t bktsize[bktnum] = { 0 };

    for (size_t si = 0; si < n; ++si)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(strings[si], depth);

        unsigned int b = find_bkt_splittertree(key, splitter, splitter_tree, numsplitters);

        assert(b < bktnum);

        bktcache[si] = b;
        ++bktsize[ b ];
    }

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

    size_t bktindex[bktnum];
    bktindex[0] = bktsize[0];
    size_t last_bkt_size = bktsize[0];
    for (unsigned int i=1; i < bktnum; ++i) {
        bktindex[i] = bktindex[i-1] + bktsize[i];
        if (bktsize[i]) last_bkt_size = bktsize[i];
    }
    assert(bktindex[bktnum-1] == n);

    // step 4: premute in-place

    for (size_t i=0, j; i < n - last_bkt_size; )
    {
        string perm = strings[i];
        uint16_t permbkt = bktcache[i];

        while ( (j = --bktindex[ permbkt ]) > i )
        {
            std::swap(perm, strings[j]);
            std::swap(permbkt, bktcache[j]);
        }

        strings[i] = perm;
        i += bktsize[ permbkt ];
    }

    delete [] bktcache;

    // step 5: recursion

    size_t bsum = 0;
    for (size_t i=0; i < bktnum-1; ++i) {
        // i is even -> bkt[i] is less-than bucket
        if (bktsize[i] > 1)
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bsum << " size " << bktsize[i] << " lcp " << int(splitter_lcp[i/2]));
            sample_sortBTC(strings+bsum, bktsize[i], depth + splitter_lcp[i/2]);
        }
        bsum += bktsize[i];
        ++i;
        // i is odd -> bkt[i] is equal bucket
        if (bktsize[i] > 1)
        {
            if ( (splitter[i/2] & 0xFF) == 0 ) { // equal-bucket has NULL-terminated key
                // done.
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " is done!");
            }
            else {
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " lcp keydepth!");
                sample_sortBTC(strings+bsum, bktsize[i], depth + sizeof(key_type));
            }
        }
        bsum += bktsize[i];
    }
    if (bktsize[bktnum-1] > 0)
    {
        DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bsum << " size " << bktsize[bktnum-1] << " no lcp");
        sample_sortBTC(strings+bsum, bktsize[bktnum-1], depth);
    }
    bsum += bktsize[bktnum-1];
    assert(bsum == n);
}

void bingmann_sample_sortBTC(string* strings, size_t n) { return sample_sortBTC(strings,n,0); }
CONTESTANT_REGISTER_UCARRAY(bingmann_sample_sortBTC, "bingmann/sample_sortBTC (binary tree, bkt cache)")

} // namespace bingmann_sample_sort
