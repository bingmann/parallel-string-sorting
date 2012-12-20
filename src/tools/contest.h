/******************************************************************************
 * src/tools/contest.h
 *
 * Class to register contestants for speed test.
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

class Contestant;

class Contest
{
public:
    typedef std::vector<Contestant*>    list_type;

    list_type           m_list;

    void register_contestant(Contestant* c)
    {
        m_list.push_back(c);
    }

    void run_contest(const char* path); // implemented in main.cc
};

Contest* getContestSingleton()
{
    static Contest* c = NULL;
    if (!c) c = new Contest;
    return c;
}

class Contestant
{
public:
    std::string         m_funcname, m_description;

    Contestant(const std::string& funcname,
               const std::string& description)
        : m_funcname(funcname),
          m_description(description)
    {
        getContestSingleton()->register_contestant(this);
    }

    virtual void run() = 0; // depends on individual sorter's interface
};

class Contestant_UCArray : public Contestant
{
public:
    typedef void (*func_type)(unsigned char** strings, size_t n);

    func_type           m_func;

    Contestant_UCArray(func_type func,
                       const std::string& funcname,
                       const std::string& description)
        : Contestant(funcname,description),
          m_func(func)
    {
    }

    virtual void run(); // implemented in main.cc
};

#define CONTESTANT_REGISTER_UCARRAY(func, desc)                         \
    static const class Contestant* _Contestant_##func##_register =      \
        new Contestant_UCArray(func,#func,desc);
