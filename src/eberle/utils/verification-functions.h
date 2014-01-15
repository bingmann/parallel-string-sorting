/******************************************************************************
 * src/eberle/utils/verification-functions.h
 *
 * Functions to verify sort order and correctnes of lcps.
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

#ifndef VERIFICATION_FUNCTIONS_H_
#define VERIFICATION_FUNCTIONS_H_

#include "types.h"

namespace eberle_utils
{

using namespace types;

static inline
bool
checkLcps(LcpStringPtr output, size_t n, lcp_t expectedFirstLcp)
{
    bool allValid = true;

    if (output.lcp() != expectedFirstLcp)
    {
        std::cout << "output.lcp() " << output.lcp() << " but should be " << expectedFirstLcp << std::endl;
        allValid = false;
    }

    string lastText = output.str();
    ++output;

    for (size_t i = 1; i < n; ++i, ++output)
    {
        string s1 = lastText, s2 = output.str();
        lcp_t lcp = 0;
        while (*s1 != 0 && *s1 == *s2)
            ++lcp, ++s1, ++s2;

        if (lcp != output.lcp())
        {
            std::cout << "output[" << i << "].lcp mismatch " << lcp << " != " << output.lcp() << std::endl;
            allValid = false;
        }

        lastText = output.str();
    }

    if (allValid)
    {
        std::cout << "All LCPs valid!" << std::endl;
    }
    else
    {
        std::cout << "Found invalid LCPS!" << std::endl;
    }

    return allValid;
}

} // namespace eberle_utils

#endif // VERIFICATION_FUNCTIONS_H_
