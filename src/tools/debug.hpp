/*******************************************************************************
 * src/tools/debug.hpp
 *
 * Debugging macros and memory debugging utility
 *
 *******************************************************************************
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
 ******************************************************************************/

#ifndef PSS_SRC_TOOLS_DEBUG_HEADER
#define PSS_SRC_TOOLS_DEBUG_HEADER

#include <sstream>
#include <iostream>

#define DBGX(dbg, X)   do { if (dbg) { std::cout << X; } \
} while (0)

#define DBG(dbg, X)    DBGX(dbg, __FUNCTION__ << "() " << X << std::endl)

#define DBG1(dbg, X)   DBGX(dbg, __FUNCTION__ << "() " << X)
#define DBG2(dbg, X)   DBGX(dbg, X)
#define DBG3(dbg, X)   DBGX(dbg, X << std::endl)

#define STRINGIFY_X(s) #s
#define STRINGIFY(s)  STRINGIFY_X(s)

#define CONCAT(a, b)     a ## b
#define CONCAT_EXPANDED(a, b) CONCAT(a, b)

#define UNUSED(x) do { (void)x; } while (0)

// *** Support for multi-threaded programs, to activate:
// #undef DBGX
// #define DBGX DBGX_OMP

#define DBGX_OMP(dbg, X) do {                  \
        if (dbg) {                             \
            std::ostringstream os;             \
            os << X; g_debug_output(os.str()); \
        }                                      \
} while (0)

static inline void g_debug_output(const std::string& s)
{
#pragma omp critical
    std::cout << s;
}

// *** an always-on ASSERT

#define ASSERT(expr)  do {                                           \
        if (!(expr)) {                                               \
            fprintf(stderr, "%s:%u %s: Assertion '%s' failed!\n",    \
                    __FILE__, __LINE__, __PRETTY_FUNCTION__, #expr); \
            abort();                                                 \
        }                                                            \
} while (0)

#endif // !PSS_SRC_TOOLS_DEBUG_HEADER

/******************************************************************************/
