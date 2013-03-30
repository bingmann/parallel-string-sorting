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

static const size_t g_samplesort_smallsort = 128;

static const size_t oversample_factor = 1;

static size_t g_ss_steps, g_rs_steps;

static inline void sample_sort_pre()
{
    g_ss_steps = g_rs_steps = 0;
}

static inline void sample_sort_post()
{
    g_statscache >> "l2cache" << l2cache
                 >> "steps_sample_sort" << g_ss_steps
                 >> "steps_base_sort" << g_rs_steps;
}

// ----------------------------------------------------------------------------

#if 0
/// Variant 5 of string sample-sort: use super-scalar binary search on splitters, with index caching non-recursive.
struct SampleSortBTCnr
{
#if 0
    static const size_t numsplitters = 31;
#else
    //static const size_t l2cache = 256*1024;

    // bounding equations:
    // splitters            + bktsize
    // n * sizeof(key_type) + (2*n+1) * sizeof(size_t) <= l2cache

    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (sizeof(key_type) + 2 * sizeof(size_t));
    //static const size_t numsplitters2 = l2cache / sizeof(key_type);
    static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / (2 * sizeof(size_t));

    //static const size_t numsplitters2 = ( l2cache - sizeof(size_t) ) / ( sizeof(key_type) );

    static const size_t numsplitters = (1 << logfloor_<numsplitters2>::value) - 1;
#endif

    static const size_t bktnum = 2 * numsplitters + 1;

    string* str;
    size_t idx;
    size_t depth;
    size_t bktsize[bktnum];

    key_type splitter[numsplitters];
    unsigned char splitter_lcp[numsplitters+1];

    SampleSortBTCnr(string* strings, size_t n, size_t depth, uint16_t* bktcache)
    {
        // step 1: select splitters with oversampling

        const size_t samplesize = oversample_factor * numsplitters;

        static key_type samples[ samplesize ];

        LCGRandom rng(&strings);

        for (unsigned int i = 0; i < samplesize; ++i)
        {
            samples[i] = get_char<key_type>(strings[ rng() % n ], depth);
        }

        std::sort(samples, samples + samplesize);

        DBG(debug_splitter, "splitter:");
        splitter_lcp[0] = 0; // sentinel for first <-everything bucket
        for (size_t i = 0, j = oversample_factor/2; i < numsplitters; ++i)
        {
            splitter[i] = samples[j];
            DBG(debug_splitter, "key " << toHex(splitter[i]));

            if (i != 0) {
                key_type xorSplit = splitter[i-1] ^ splitter[i];

                DBG1(debug_splitter, "    XOR -> " << toHex(xorSplit) << " - ");

                DBG3(debug_splitter, count_high_zero_bits(xorSplit) << " bits = "
                     << count_high_zero_bits(xorSplit) / 8 << " chars lcp");

                splitter_lcp[i] = count_high_zero_bits(xorSplit) / 8;
            }

            j += oversample_factor;
        }
        splitter_lcp[numsplitters] = 0; // sentinel for last >-everything bucket

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
#if 0
        //uint16_t* bktcache = new uint16_t[n];

        memset(bktsize, 0, bktnum * sizeof(size_t));

        for (size_t si = 0; si < n; ++si)
        {
            // binary search in splitter with equal check
            key_type key = get_char<key_type>(strings[si], depth);

            unsigned int b = find_bkt_tree(key, splitter, splitter_tree, numsplitters);

            assert(b < bktnum);

            bktcache[si] = b;
            ++bktsize[ b ];
        }
#else
        //uint16_t* bktcache = new uint16_t[n];

        for (size_t si = 0; si < n; ++si)
        {
            // binary search in splitter with equal check
            key_type key = get_char<key_type>(strings[si], depth);

            unsigned int b = find_bkt_tree(key, splitter, splitter_tree, numsplitters);

            assert(b < bktnum);

            bktcache[si] = b;
        }

        memset(bktsize, 0, bktnum * sizeof(size_t));

        for (size_t si = 0; si < n; ++si)
            ++bktsize[ bktcache[si] ];
#endif

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

        //delete [] bktcache;

        str = strings;
        idx = 0; // will increment to 1 on first process, bkt 0 is not sorted further
        this->depth = depth;
    }

};

void bingmann_sample_sortBTCnr(string* strings, size_t n)
{
    g_statscache >> "l2cache" << l2cache;

    if (n < g_samplesort_smallsort && 0)
        return bingmann_radix_sort::msd_CI5(strings,n,0);

    typedef SampleSortBTCnr Step;

    uint16_t* bktcache = new uint16_t[n];

    std::stack< Step, std::vector<Step> > stack;
    stack.push( Step(strings,n,0,bktcache) );

    size_t ss_steps = 0, rs_steps = 0;

    // step 5: "recursion"

    while ( stack.size() )
    {
        while ( stack.top().idx < Step::bktnum )
        {
            Step& s = stack.top();
            size_t i = s.idx++; // process the bucket s.idx

            // i is even -> bkt[i] is less-than bucket
            if ((i & 1) == 0)
            {
                if (s.bktsize[i] == 0)
                    ;
                else if (s.bktsize[i] < g_samplesort_smallsort)
                {
                    if (i == Step::bktnum-1)
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << s.bktsize[i] << " no lcp");
                    else
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << s.bktsize[i] << " lcp " << int(s.splitter_lcp[i/2]));

                    ++rs_steps;
                    bingmann_radix_sort::msd_CI5(s.str, s.bktsize[i], s.depth + s.splitter_lcp[i/2]);
                    s.str += s.bktsize[i];
                }
                else
                {
                    if (i == Step::bktnum-1)
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: > bkt " << i << " size " << s.bktsize[i] << " no lcp");
                    else
                        DBG(debug_recursion, "Recurse[" << s.depth << "]: < bkt " << i << " size " << s.bktsize[i] << " lcp " << int(s.splitter_lcp[i/2]));

                    ++ss_steps;
                    s.str += s.bktsize[i];
                    stack.push( Step(s.str - s.bktsize[i], s.bktsize[i], s.depth + s.splitter_lcp[i/2], bktcache) );
                    // cannot add here, because s may have invalidated
                }
            }
            // i is odd -> bkt[i] is equal bucket
            else
            {
                if (s.bktsize[i] == 0)
                    ;
                else if ( (s.splitter[i/2] & 0xFF) == 0 ) { // equal-bucket has NULL-terminated key, done.
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << s.bktsize[i] << " is done!");
                    s.str += s.bktsize[i];
                }
                else if (s.bktsize[i] < g_samplesort_smallsort)
                {
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << s.bktsize[i] << " lcp keydepth!");

                    ++rs_steps;
                    bingmann_radix_sort::msd_CI5(s.str, s.bktsize[i], s.depth + sizeof(key_type));
                    s.str += s.bktsize[i];
                }
                else
                {
                    DBG(debug_recursion, "Recurse[" << s.depth << "]: = bkt " << i << " size " << s.bktsize[i] << " lcp keydepth!");

                    ++ss_steps;
                    s.str += s.bktsize[i];
                    stack.push( Step(s.str - s.bktsize[i], s.bktsize[i], s.depth + sizeof(key_type), bktcache) );
                    // cannot add here, because s may have invalidated
                }
            }
        }

        stack.pop();
    }

    g_statscache >> "ss_steps" << ss_steps
                 >> "rs_steps" << rs_steps;
}

CONTESTANT_REGISTER(bingmann_sample_sortBTCnr, "bingmann/sample_sortBTCnr",
                    "bingmann/sample_sortBTCnr (binary tree, bkt cache, non-recursive)")
#endif

// ------------------------------------------------------------------------------------------------------------------------

} // namespace bingmann_sample_sort
