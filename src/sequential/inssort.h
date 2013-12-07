/******************************************************************************
 * src/sequential/inssort.h
 *
 * Base insertion sort from Rantala string-sorting/external/
 *
 ******************************************************************************
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

#ifndef INSSORT_H_
#define INSSORT_H_

namespace inssort {

typedef unsigned char* string;

static inline void
inssort(string* str, int n, int d)
{
    string *pj, s, t;

    for (string* pi = str + 1; --n > 0; pi++) {
        string tmp = *pi;

        for (pj = pi; pj > str; pj--) {
            for (s = *(pj-1)+d, t = tmp+d; *s == *t && *s != 0; ++s, ++t)
                ;
            if (*s <= *t)
                break;
            *pj = *(pj-1);
        }
        *pj = tmp;
    }
}

static inline
void insertion_sort(string* a, size_t n)
{ inssort(a, n, 0); }

CONTESTANT_REGISTER(insertion_sort, "insertion_sort", "String Insertion-Sort")

} // namespace  inssort

#endif // INSSORT_H_
