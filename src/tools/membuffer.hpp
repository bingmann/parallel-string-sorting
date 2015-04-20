/*******************************************************************************
 * src/tools/membuffer.h
 *
 * Dumb memory buffer using new[]/delete[] which does not initialize the
 * values. This is used instead of std::vector when the values need not be
 * initialized.
 *
 *******************************************************************************
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
 ******************************************************************************/

#ifndef PSS_SRC_TOOLS_MEMBUFFER_HEADER
#define PSS_SRC_TOOLS_MEMBUFFER_HEADER

template <typename Type>
class membuffer
{
private:
    /// pointer to allocated memory area
    Type* m_ptr;

    /// size of allocated memory
    size_t m_size;

    /// protect copy-constructor
    membuffer(const membuffer& b);

    /// protect assignment operator
    membuffer& operator = (const membuffer& b);

public:
    /// Invalid memory buffer
    inline membuffer() : m_ptr(NULL), m_size(0) { }

    /// Allocate memory buffer
    inline membuffer(size_t size)
        : m_ptr(new Type[size]),
          m_size(size)
    { }

    /// Deallocate memory buffer
    inline ~membuffer()
    {
        delete[] m_ptr;
    }

    /// Accessor to elements
    inline Type& operator [] (size_t n) const
    {
        return m_ptr[n];
    }

    inline Type * data() { return m_ptr; }
    inline const Type * data() const { return m_ptr; }

    inline size_t size() const { return m_size; }

    inline Type * begin() { return m_ptr; }
    inline Type * end() { return m_ptr + m_size; }

    inline const Type * begin() const { return m_ptr; }
    inline const Type * end() const { return m_ptr + m_size; }

    inline void copy(const membuffer& b)
    {
        if (m_ptr) delete[] m_ptr;
        m_ptr = new Type[b.m_size];
        m_size = b.m_size;
        memcpy(m_ptr, b.m_ptr, b.m_size * sizeof(Type));
    }
};

#endif // !PSS_SRC_TOOLS_MEMBUFFER_HEADER

/******************************************************************************/
