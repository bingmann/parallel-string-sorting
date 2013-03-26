/******************************************************************************
 * src/tools/contest.h
 *
 * Simple linear congruential random generator for 64-bit pseudo-random numbers.
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

/// Simple linear congruential random generator
class LCGRandom
{
private:
    size_t      xn;
public:
    inline LCGRandom(size_t seed) : xn(seed) { }
    inline LCGRandom(void* ptrseed) : xn((size_t)ptrseed) { }
    inline size_t operator()() {
        xn = 0x27BB2EE687B0B0FDLLU * xn + 0xB504F32DLU;
        return xn;
    }
};
