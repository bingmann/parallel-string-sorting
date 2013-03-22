/******************************************************************************
 * src/tools/checker.h
 *
 * Tools to check correctness of sort.
 *
 ******************************************************************************
 * Copyright (C) 2012 Timo Bingmann <tb@panthema.net>
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

class PermutationCheck
{
private:
    // polynomial evaluation of (z-0)*(z-1)*(z-2)...(z-(n-1)) mod (p)

    mpz_class v;        // evaluation result
    size_t iter;        // offset / iteration

public:
    PermutationCheck(const std::vector<unsigned char*>& stringptr)
    {
        mpz_class p("18446744073709551557", 10); // largest prime < 2^64

        for (iter = 1; iter < 100000; ++iter)
        {
            // use address after end of string area as evalution point
            uint64_t z = (uint64_t)g_string_data + g_string_datasize + iter;
            v = 1;

            for (size_t i = 0; i < stringptr.size(); ++i)
            {
                assert(z >= (uint64_t)stringptr[i]);
                v *= (z - (uint64_t)stringptr[i]);
                v %= p;
            }

            //std::cout << "eval = " << v << "\n";

            if (v != 0) break; // v != 0, thus a true check value

            std::cerr << "PermutationCheck: evaluation returned zero, repeating.\n";
        }
    }

    bool check(const std::vector<unsigned char*>& stringptr) const
    {
        mpz_class p("18446744073709551557", 10); // largest prime < 2^64

        // use address after end of string area as evalution point
        uint64_t z = (uint64_t)g_string_data + g_string_datasize + iter;
        mpz_class w = 1;

        for (size_t i = 0; i < stringptr.size(); ++i)
        {
            assert(z >= (uint64_t)stringptr[i]);
            w *= (z - (uint64_t)stringptr[i]);
            w %= p;
        }

        return (v == w);
    }
};

class PermutationCheck64
{
private:
    // polynomial evaluation of (z-0)*(z-1)*(z-2)...(z-(n-1)) mod (p)
    static const uint64_t p = 4294967291;       // largest prime < 2^32

    uint64_t v;         // evaluation result
    size_t iter;        // offset / iteration

public:
    PermutationCheck64(const std::vector<unsigned char*>& stringptr)
    {
        for (iter = 1; iter < 100000; ++iter)
        {
            // use address after end of string area as evalution point
            uint64_t z = (uint64_t)g_string_data + g_string_datasize + iter;
            v = 1;

            for (size_t i = 0; i < stringptr.size(); ++i)
            {
                assert(z >= (uint64_t)stringptr[i]);
                v *= (z - (uint64_t)stringptr[i]);
                v %= p;
            }

            //std::cout << "eval = " << v << "\n";

            if (v != 0) return;

            std::cerr << "PermutationCheck: evaluation returned zero, repeating.\n";
        }
    }

    bool check(const std::vector<unsigned char*>& stringptr) const
    {
        // use address after end of string area as evalution point
        uint64_t z = (uint64_t)g_string_data + g_string_datasize + iter;
        uint64_t w = 1;

        for (size_t i = 0; i < stringptr.size(); ++i)
        {
            assert(z >= (uint64_t)stringptr[i]);
            w *= (z - (uint64_t)stringptr[i]);
            w %= p;
        }

        return (v == w);
    }
};

bool check_sorted_order(const std::vector<unsigned char*>& stringptr, const PermutationCheck& pc)
{
    // check permutation via polynomial evaluation

    if (!pc.check(stringptr))
    {
        std::cerr << "error: not a permutation\n";
        //return false;
    }

    // check order naively

    for (size_t i = 1; i < stringptr.size(); ++i)
    {
        if (strcmp((const char*)stringptr[i-1], (const char*)stringptr[i]) > 0)
        {
            std::cerr << "error: invalid order at pair " << i-1 << " and " << i << "\n";
            return false;
        }
    }

    return true;
}

size_t calc_distinguishing_prefix(const std::vector<unsigned char*>& stringptr, size_t n)
{
    size_t D = 0;
    size_t pdepth = 0;

    for (size_t i = 1; i < stringptr.size(); ++i)
    {
        size_t depth = 0;
        while ( stringptr[i-1][depth] == stringptr[i][depth] &&
                stringptr[i-1][depth] != 0 ) ++depth;
        // depth == LCP of prev and this
        depth++;

        if (pdepth < depth) {
            D += depth - pdepth; // add extra distinguishing characters
        }
        D += depth;
        pdepth = depth;
    }

    std::cout << "percentage dprefix = " << (D * 100.0 / n) << "\n";

    return D;
}
