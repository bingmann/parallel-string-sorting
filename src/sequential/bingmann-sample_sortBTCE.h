/******************************************************************************
 * src/sequential/bingmann-sample_sortBTCE.h
 *
 * Experiments with sequential Super Scalar String Sample-Sort (S^5).
 *
 * Binary tree search with equality branch and bucket cache.
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

namespace bingmann_sample_sortBTCE {

using namespace bingmann_sample_sort;

// ----------------------------------------------------------------------------

/// Variant of string sample-sort: use super-scalar binary search on splitters with equality check and index caching.

static inline std::string binary(uint16_t v) {
    char binstr[17];
    binstr[16] = 0;
    for (int i = 0; i < 16; i++) {
        binstr[15-i] = (v & 1) ? '1' : '0';
        v /= 2;
    }
    return binstr;
}

static inline unsigned int
treeid_to_bkt(unsigned int id, size_t treebits, size_t numsplitters)
{
    assert(id > 0);
    //std::cout << "index: " << id << " = " << binary(id) << "\n";

    //int treebits = 4;
    //int bitmask = ((1 << treebits)-1);
    static const int bitmask = numsplitters;

    int hi = treebits-32 + count_high_zero_bits<uint32_t>(id); // sdfsdf
    //std::cout << "high zero: " << hi << "\n";

    unsigned int bkt = ((id << (hi+1)) & bitmask) | (1 << hi);

    //std::cout << "bkt: " << bkt << " = " << binary(bkt) << "\n";

    return bkt;
}

/// search in splitter tree for bucket number
inline unsigned int
find_bkt_tree_equal(const key_type& key, const key_type* /* splitter */, const key_type* splitter_tree0, size_t treebits, size_t numsplitters)
{
    // binary tree traversal without left branch

    const key_type* splitter_tree = splitter_tree0 - 1;
    unsigned int i = 1;

    while ( i <= numsplitters )
    {
        if (key == splitter_tree[i])
            return 2 * treeid_to_bkt(i,treebits,numsplitters) - 1;
        else if (key < splitter_tree[i])
            i = 2*i + 0;
        else // (key > splitter_tree[i])
            i = 2*i + 1;
    }

    i -= numsplitters+1;

    return 2 * i; // < or > bucket
}

/// binary search on splitter array for bucket number
inline unsigned int
find_bkt_tree_asmequal(const key_type& key, const key_type* /* splitter */, const key_type* splitter_tree0, size_t treebits, size_t numsplitters)
{
    const key_type* splitter_tree = splitter_tree0 - 1;
    unsigned int i;

    // hand-coded assembler binary tree traversal with equality
    asm("mov    $1, %%rax \n"             // rax = i
        // body of while loop
        "1: \n"
        "cmpq	(%[splitter_tree],%%rax,8), %[key] \n"
        "je     2f \n"
        "lea    (%%rax,%%rax), %%rax \n"
        "lea    1(%%rax), %%rcx \n"
        "cmova  %%rcx, %%rax \n"             // CMOV rax = 2 * i + 1
        "cmp    %[numsplitters1], %%rax \n"  // i < numsplitters+1
        "jb     1b \n"
        "sub    %[numsplitters1], %%rax \n"  // i -= numsplitters+1;
        "lea    (%%rax,%%rax), %%rax \n"     // i = i*2
        "jmp    3f \n"
        "2: \n"
        "bsr    %%rax, %%rdx \n"             // dx = bit number of highest one
        "mov    %[treebits], %%rcx \n"
        "sub    %%rdx, %%rcx \n"             // cx = treebits - highest
        "shl    %%cl, %%rax \n"              // shift ax to left
        "and    %[numsplitters], %%rax \n"   // mask off other bits
        "lea    -1(%%rcx), %%rcx \n"
        "mov    $1, %%rdx \n"                // dx = (1 << (hi-1))
        "shl    %%cl, %%rdx \n"              //
        "or     %%rdx, %%rax \n"             // ax = OR of both
        "lea    -1(%%rax,%%rax), %%rax \n"    // i = i * 2 - 1
        "3: \n"
        : "=&a" (i)
        : [key] "r" (key), [splitter_tree] "r" (splitter_tree),
          [numsplitters1] "g" (numsplitters+1),
          [treebits] "g" (treebits),
          [numsplitters] "g" (numsplitters)
        : "rcx", "rdx");

    return i;
}

/// Variant of string sample-sort: use super-scalar binary search on splitters, with index caching.
template <unsigned int (*find_bkt)(const key_type&, const key_type*, const key_type*, size_t, size_t)>
void sample_sortBTCE2(string* strings, size_t n, size_t depth)
{
#if 0
    static const size_t numsplitters2 = 16;
    static const size_t treebits = logfloor_<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;
#else
    //static const size_t l2cache = 256*1024;

    // bounding equations:
    // splitters            + bktsize
    // n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));
    //static const size_t numsplitters2 = l2cache / sizeof(key_type);
    static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (2 * sizeof(size_t));

    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / ( sizeof(key_type) );

    static const size_t treebits = logfloor_<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;
#endif

    //std::cout << "l2cache : " << l2cache << " - numsplitter " << numsplitters << "\n";

    if (n < g_samplesort_smallsort)
    {
        g_rs_steps++;
        //return inssort::inssort(strings, n, depth);
        return bingmann_radix_sort::msd_CI5(strings, n, depth);
        //return bs_mkqs::ssort2(strings, n, depth);
    }
    g_ss_steps++;

    //std::cout << "numsplitters: " << numsplitters << "\n";

    // step 1: select splitters with oversampling

    const size_t samplesize = oversample_factor * numsplitters;

    static key_type samples[ samplesize ];

    LCGRandom rng(&strings);

    for (unsigned int i = 0; i < samplesize; ++i)
    {
        samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
    }

    std::sort(samples, samples + samplesize);

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

        key_type* splitter_tree1 = splitter_tree-1;
        key_type prevsplitter = 0;

        DBG(debug_splitter, "splitter:");
        splitter_lcp[0] = 0; // sentinel for first < everything bucket
        for (size_t i = 0, j = oversample_factor/2; i < numsplitters; ++i)
        {
            const key_type& splitter = samples[j];

            int l = __builtin_ctz(i+1);
            DBG(debug_splitter, "splitter[" << i << "] on level " << l
                << " = tree[" << treelvl[l] << "] = key " << toHex(splitter));
            splitter_tree1[ treelvl[l]++ ] = samples[j];

            if (i != 0) {
                key_type xorSplit = prevsplitter ^ splitter;

                DBG1(debug_splitter, "    XOR -> " << toHex(xorSplit) << " - ");

                DBG3(debug_splitter, count_high_zero_bits(xorSplit) << " bits = "
                     << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

                splitter_lcp[i] = (count_high_zero_bits(xorSplit) / 8)
                    | ((splitter & 0xFF) ? 0 : 0x80); // marker for done splitters
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

    uint16_t* bktcache = new uint16_t[n];

    static const size_t bktnum = 2*numsplitters+1;

    for (size_t si = 0; si < n; ++si)
    {
        // binary search in splitter with equal check
        key_type key = get_char<key_type>(strings[si], depth);

        unsigned int b = find_bkt(key, NULL, splitter_tree, treebits, numsplitters);

        assert(b < bktnum);

        bktcache[si] = b;
    }

    delete [] splitter_tree;

    size_t* bktsize = new size_t[bktnum];
    memset(bktsize, 0, bktnum * sizeof(size_t));

    for (size_t si = 0; si < n; ++si)
        ++bktsize[ bktcache[si] ];

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

    size_t i = 0, bsum = 0;
    while (i < bktnum-1)
    {
        // i is even -> bkt[i] is less-than bucket
        if (bktsize[i] > 1)
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bsum << " size " << bktsize[i] << " lcp " << int(splitter_lcp[i/2] & 0x7F));
            sample_sortBTCE2<find_bkt>(strings+bsum, bktsize[i], depth + (splitter_lcp[i/2] & 0x7F));
        }
        bsum += bktsize[i++];

        // i is odd -> bkt[i] is equal bucket
        if (bktsize[i] > 1)
        {
            if ( splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key, done.
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " is done!");
            }
            else {
                DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " lcp keydepth!");
                sample_sortBTCE2<find_bkt>(strings+bsum, bktsize[i], depth + sizeof(key_type));
            }
        }
        bsum += bktsize[i++];
    }
    if (bktsize[i] > 0)
    {
        DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bsum << " size " << bktsize[i] << " no lcp");
        sample_sortBTCE2<find_bkt>(strings+bsum, bktsize[i], depth);
    }
    bsum += bktsize[i++];
    assert(i == bktnum && bsum == n);

    delete [] splitter_lcp;
    delete [] bktsize;
}

void bingmann_sample_sortBTCE2(string* strings, size_t n)
{
    sample_sort_pre();
    sample_sortBTCE2<find_bkt_tree_equal>(strings,n,0);
    sample_sort_post();
}

CONTESTANT_REGISTER(bingmann_sample_sortBTCE2, "bingmann/sample_sortBTCE2",
                    "bingmann/sample_sortBTCE2 (binary tree equal, bkt cache)")

void bingmann_sample_sortBTCE2A(string* strings, size_t n)
{
    sample_sort_pre();
    sample_sortBTCE2<find_bkt_tree_asmequal>(strings,n,0);
    sample_sort_post();
}

CONTESTANT_REGISTER(bingmann_sample_sortBTCE2A, "bingmann/sample_sortBTCE2A",
                    "bingmann/sample_sortBTCE2A (binary tree equal, asm CMOV, bkt cache)")

// ----------------------------------------------------------------------------

// Variant of string sample-sort: use recursive binary function on sample array
// to generate splitters.

class SampleSortBTCE3
{
public:

#if 0
    static const size_t numsplitters2 = 16;
    static const size_t treebits = logfloor_<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;
#else
    //static const size_t l2cache = 256*1024;

    // bounding equations:
    // splitters            + bktsize
    // n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));
    //static const size_t numsplitters2 = l2cache / sizeof(key_type);
    static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (2 * sizeof(size_t));

    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / ( sizeof(key_type) );

    static const size_t treebits = logfloor_<numsplitters2>::value;
    static const size_t numsplitters = (1 << treebits) - 1;
#endif

    struct SplitterTree
    {
        key_type splitter_tree[numsplitters+1];
        unsigned char splitter_lcp[numsplitters];

        key_type rec_buildtree(key_type* samples, size_t lo, size_t hi, unsigned int treeidx,
                               key_type& rec_prevkey, size_t& iter)
        {
            DBG(debug_splitter, "rec_buildtree(" << lo << "," << hi << ", treeidx=" << treeidx << ")");

            // pick middle element as splitter
            size_t mid = (lo + hi) >> 1;

            DBG(debug_splitter, "tree[" << treeidx << "] = samples[" << mid << "] = "
                << toHex(samples[mid]));

            key_type mykey = splitter_tree[treeidx] = samples[mid];
#if 0
            size_t midhi = mid;
            while (lo < midhi && samples[midhi-1] == mykey) midhi--;

            size_t midlo = mid+1;
            while (midlo < hi && samples[midlo] == mykey) midlo++;

            if (midlo - midhi > 1)
                DBG(0, "key range = [" << midhi << "," << midlo << ")");
#else
            const size_t midhi = mid, midlo = mid+1;
#endif
            if (2 * treeidx < numsplitters)
            {
                key_type prevkey = rec_buildtree(samples, lo, midhi, 2 * treeidx + 0, rec_prevkey, iter);

                key_type xorSplit = prevkey ^ mykey;

                DBG(debug_splitter, "    lcp: " << toHex(prevkey) << " XOR " << toHex(mykey) << " = "
                    << toHex(xorSplit) << " - " << count_high_zero_bits(xorSplit) << " bits = "
                    << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

                splitter_lcp[iter++] = (count_high_zero_bits(xorSplit) / 8)
                    | ((mykey & 0xFF) ? 0 : 0x80); // marker for done splitters

                return rec_buildtree(samples, midlo, hi, 2 * treeidx + 1, mykey, iter);
            }
            else
            {
                key_type xorSplit = rec_prevkey ^ mykey;

                DBG(debug_splitter, "    lcp: " << toHex(rec_prevkey) << " XOR " << toHex(mykey) << " = "
                    << toHex(xorSplit) << " - " << count_high_zero_bits(xorSplit) << " bits = "
                    << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

                splitter_lcp[iter++] = (count_high_zero_bits(xorSplit) / 8)
                    | ((mykey & 0xFF) ? 0 : 0x80); // marker for done splitters

                return mykey;
            }
        }
    };

    /// Variant of string sample-sort: use super-scalar binary search on splitters, with index caching.
    template <unsigned int (*find_bkt)(const key_type&, const key_type*, const key_type*, size_t, size_t)>
    static void sort(string* strings, size_t n, size_t depth)
    {
        if (n < g_samplesort_smallsort)
        {
            g_rs_steps++;
            //return inssort::inssort(strings, n, depth);
            return bingmann_radix_sort::msd_CI5(strings, n, depth);
            //return bs_mkqs::ssort2(strings, n, depth);
        }
        g_ss_steps++;

        // step 1: select splitters with oversampling

        //const size_t oversample_factor = 1;
        const size_t samplesize = oversample_factor * numsplitters;

        static key_type samples[ samplesize ];

        LCGRandom rng(&strings);

        for (unsigned int i = 0; i < samplesize; ++i)
        {
            samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
        }

        std::sort(samples, samples + samplesize);

        SplitterTree tree;
        {
            key_type sentinel = 0;
            size_t iter = 0;
            tree.rec_buildtree(samples, 0, samplesize, 1, sentinel, iter);
            tree.splitter_lcp[0] = 0; // overwrite sentinel for first < everything bucket
            assert(iter == numsplitters);
        }

#if 0
        key_type splitter_tree[numsplitters+1];
        unsigned char splitter_lcp[numsplitters];

        {
            unsigned int treelvl[treebits]; // treebits == treelevels

            // fill in initial level indexes
            int idx = 1;
            for (unsigned int i = treebits; i > 0; ) {
                treelvl[--i] = idx;
                idx <<= 1;
            }

            key_type prevsplitter;

            DBG(debug_splitter, "splitter:");
            splitter_lcp[0] = 0; // sentinel for first < everything bucket
            for (size_t i = 0, j = oversample_factor/2; i < numsplitters; ++i)
            {
                const key_type& splitter = samples[j];

                int l = __builtin_ctz(i+1);
                DBG(debug_splitter, "splitter[" << i << "] on level " << l
                    << " = tree[" << treelvl[l] << "] = key " << toHex(splitter));
                splitter_tree[ treelvl[l]++ ] = samples[j];

                if (i != 0) {
                    key_type xorSplit = prevsplitter ^ splitter;

                    DBG1(debug_splitter, "    XOR -> " << toHex(xorSplit) << " - ");

                    DBG3(debug_splitter, count_high_zero_bits(xorSplit) << " bits = "
                         << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

                    splitter_lcp[i] = (count_high_zero_bits(xorSplit) / 8)
                        | ((splitter & 0xFF) ? 0 : 0x80); // marker for done splitters
                }

                prevsplitter = splitter;
                j += oversample_factor;
            }
        }

        for(size_t i = 0; i < numsplitters; ++i)
        {
            DBG(1, "splitter_tree[" << i+1 << "] = " << tree.splitter_tree[i+1] << " -? " << splitter_tree[i+1]);
            assert(tree.splitter_tree[i+1] == splitter_tree[i+1]);

            DBG(1, "splitter_lcp[" << i << "] = " << int(tree.splitter_lcp[i]) << " -? " << int(splitter_lcp[i]));
            assert(tree.splitter_lcp[i] == splitter_lcp[i]);
        }
#endif

        // step 2.2: classify all strings and count bucket sizes

        uint16_t* bktcache = new uint16_t[n];

        static const size_t bktnum = 2*numsplitters+1;

        for (size_t si = 0; si < n; ++si)
        {
            // binary search in splitter with equal check
            key_type key = get_char<key_type>(strings[si], depth);

            unsigned int b = find_bkt(key, NULL, tree.splitter_tree, treebits, numsplitters);

            assert(b < bktnum);

            bktcache[si] = b;
        }

        size_t* bktsize = new size_t[bktnum];
        memset(bktsize, 0, bktnum * sizeof(size_t));

        for (size_t si = 0; si < n; ++si)
            ++bktsize[ bktcache[si] ];

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

        size_t i = 0, bsum = 0;
        while (i < bktnum-1)
        {
            // i is even -> bkt[i] is less-than bucket
            if (bktsize[i] > 1)
            {
                DBG(debug_recursion, "Recurse[" << depth << "]: < bkt " << bsum << " size " << bktsize[i] << " lcp " << int(tree.splitter_lcp[i/2] & 0x7F));
                sort<find_bkt>(strings+bsum, bktsize[i], depth + (tree.splitter_lcp[i/2] & 0x7F));
            }
            bsum += bktsize[i++];

            // i is odd -> bkt[i] is equal bucket
            if (bktsize[i] > 1)
            {
                if ( tree.splitter_lcp[i/2] & 0x80 ) { // equal-bucket has NULL-terminated key, done.
                    DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " is done!");
                }
                else {
                    DBG(debug_recursion, "Recurse[" << depth << "]: = bkt " << bsum << " size " << bktsize[i] << " lcp keydepth!");
                    sort<find_bkt>(strings+bsum, bktsize[i], depth + sizeof(key_type));
                }
            }
            bsum += bktsize[i++];
        }
        if (bktsize[i] > 0)
        {
            DBG(debug_recursion, "Recurse[" << depth << "]: > bkt " << bsum << " size " << bktsize[i] << " no lcp");
            sort<find_bkt>(strings+bsum, bktsize[i], depth);
        }
        bsum += bktsize[i++];
        assert(i == bktnum && bsum == n);

        delete [] bktsize;
    }
};

void bingmann_sample_sortBTCE3(string* strings, size_t n)
{
    sample_sort_pre();
    SampleSortBTCE3::sort<find_bkt_tree_equal>(strings,n,0);
    sample_sort_post();
}

CONTESTANT_REGISTER(bingmann_sample_sortBTCE3, "bingmann/sample_sortBTCE3",
                    "bingmann/sample_sortBTCE3 (adapt binary tree equal, bkt cache)")

void bingmann_sample_sortBTCE3A(string* strings, size_t n)
{
    sample_sort_pre();
    SampleSortBTCE3::sort<find_bkt_tree_asmequal>(strings,n,0);
    sample_sort_post();
}

CONTESTANT_REGISTER(bingmann_sample_sortBTCE3A, "bingmann/sample_sortBTCE3A",
                    "bingmann/sample_sortBTCE3A (adapt binary tree equal, asm CMOV, bkt cache)")

// ----------------------------------------------------------------------------

} // namespace bingmann_sample_sortBTCE
