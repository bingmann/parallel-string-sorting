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

namespace inssort {

int scmp(unsigned char *s1, unsigned char *s2)
{
    while( *s1 != '\0' && *s1 == *s2 )
        s1++, s2++;
    return( *s1-*s2 );
}

void
inssort(unsigned char** a, int n, int d)
{
    unsigned char **pi, **pj, *s, *t;

    for (pi = a + 1; --n > 0; pi++) {
        unsigned char* tmp = *pi;

        for (pj = pi; pj > a; pj--) {
            for (s=*(pj-1)+d, t=tmp+d; *s==*t && *s!=0; ++s, ++t)
                ;
            if (*s <= *t)
                break;
            *pj = *(pj-1);
        }
        *pj = tmp;
    }
}

} // namespace  inssort
