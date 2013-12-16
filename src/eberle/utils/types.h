#ifndef TYPES_H_
#define TYPES_H_

namespace types
{

typedef unsigned char* string;
typedef unsigned int UINT;

struct AS
{
    UINT lcp;
    string text;

    AS()
    {
        lcp = 0;
        text = 0;
    }

    AS(UINT lcp, string text)
    {
        this->lcp = lcp;
        this->text = text;
    }
};

struct LcpStringPtr
{
    string * strings;
    unsigned* lcps;

public:
    LcpStringPtr(string* strings, unsigned* lcps) :
            strings(strings), lcps(lcps)
    {
    }

    inline void
    set(string s, unsigned lcp) const
    {
        *strings = s;
        *lcps = lcp;
    }

    inline string&
    str() const
    {
        return *strings;
    }

    inline unsigned&
    lcp() const
    {
        return *lcps;
    }

    inline void
    copyFrom(LcpStringPtr& other, size_t length) const
    {
        memcpy(strings, other.strings, length * sizeof(string));
        memcpy(lcps, other.lcps, length * sizeof(unsigned));
    }

    // preincrement
    inline LcpStringPtr&
    operator++()
    {
        ++strings;
        ++lcps;
        return *this;
    }

    friend inline LcpStringPtr
    operator+(const LcpStringPtr& ptr, size_t delta);
};

inline LcpStringPtr
operator+(const LcpStringPtr& ptr, size_t delta)
{
    return LcpStringPtr(ptr.strings + delta, ptr.lcps + delta);
}

} // namespace types

#endif // TYPES_H_
