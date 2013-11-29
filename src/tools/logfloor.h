/******************************************************************************
 * src/tools/logfloor.h
 *
 * Recursive template magic to calculate the floor(log(x)/log(2)) at run-time.
 * 
 * Usage: logfloor_<N>::value
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

#ifndef LOGFLOOT_H_
#define LOGFLOOT_H_

template <size_t N, size_t base=2>
struct logfloor_ {
    enum { value = 1 + logfloor_<N/base, base>::value };
};

template <size_t base>
struct logfloor_<1, base> { 
    enum { value = 0 };
};

template <size_t base>
struct logfloor_<0, base> {
    enum { value = 0 };
};

#endif // LOGFLOOT_H_
