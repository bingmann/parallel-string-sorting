/******************************************************************************
 * src/eberle/utils/utility-functions.h
 *
 * Header for utility functions used by src/eberle/ algorithms.
 *
 ******************************************************************************
 * Copyright (C) 2013-2014 Andreas Eberle <email@andreas-eberle.com>
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

#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <stdlib.h>
#include "types.h"

namespace eberle_utils
{

using namespace types;

static inline
void
calculateRanges(std::pair<size_t, size_t>* ranges, unsigned numberOfSplits, size_t lengthToSplit)
{
    const size_t split = lengthToSplit / numberOfSplits;
    for (unsigned i = 0; i < numberOfSplits - 1; ++i)
    {
        ranges[i] = std::make_pair(i * split, split);
    }
    ranges[numberOfSplits - 1] = std::make_pair((numberOfSplits - 1) * split, lengthToSplit - (numberOfSplits - 1) * split);
}

static inline
unsigned
getNextHigherPowerOfTwo(unsigned v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

}
// namespace eberle_utils

#endif // UTILITY_FUNCTIONS_H_
