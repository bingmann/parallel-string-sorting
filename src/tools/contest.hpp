/*******************************************************************************
 * src/tools/contest.hpp
 *
 * Class to register contestants for speed test.
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

#ifndef PSS_SRC_TOOLS_CONTEST_HEADER
#define PSS_SRC_TOOLS_CONTEST_HEADER

#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>

#include "membuffer.hpp"

#if PSS_CONTEST

//! forward declaration
class Contestant;

class Contest
{
public:
    typedef std::vector<Contestant*> list_type;

    list_type m_list;

    void register_contestant(Contestant* c);     // implemented in main.cc

    void run_contest(const char* path);          // implemented in main.cc

    bool exist_contestant(const char* algoname); // implemented in main.cc

    void list_contentants();                     // implemented in main.cc
};

extern Contest * getContestSingleton();

class Contestant
{
public:
    const char* m_algoname;
    const char* m_description;

    Contestant(const char* algoname,
               const char* description)
        : m_algoname(algoname),
          m_description(description)
    {
        getContestSingleton()->register_contestant(this);
    }

    virtual void run() = 0; // depends on individual sorter's interface

    virtual bool is_parallel() const = 0;

    inline bool operator < (const Contestant& b) const
    {
        return (strcmp(m_algoname, b.m_algoname) < 0);
    }
};

static inline bool sort_contestants(const Contestant* a, const Contestant* b)
{
    return (*a < *b);
}

class Contestant_UCArray : public Contestant
{
public:
    typedef void (* func_type)(unsigned char** strings, size_t n);
    typedef void (* func_lcp_type)(unsigned char** strings, uintptr_t* lcp, size_t n);
    typedef void (* func_lcp_cache_type)(
        uint8_t** strings, uintptr_t* lcp, uint8_t* cache, size_t n);

    func_type m_prepare_func;
    func_type m_run_func;
    func_lcp_type m_run_lcp_func;
    func_lcp_cache_type m_run_lcp_cache_func;

    Contestant_UCArray(func_type prepare_func,
                       func_type run_func,
                       const char* algoname,
                       const char* description)
        : Contestant(algoname, description),
          m_prepare_func(prepare_func),
          m_run_func(run_func),
          m_run_lcp_func(nullptr),
          m_run_lcp_cache_func(nullptr)
    { }

    Contestant_UCArray(func_type prepare_func,
                       func_lcp_type run_lcp_func,
                       const char* algoname,
                       const char* description)
        : Contestant(algoname, description),
          m_prepare_func(prepare_func),
          m_run_func(nullptr),
          m_run_lcp_func(run_lcp_func),
          m_run_lcp_cache_func(nullptr)
    { }

    Contestant_UCArray(func_type prepare_func,
                       func_lcp_cache_type run_lcp_cache_func,
                       const char* algoname,
                       const char* description)
        : Contestant(algoname, description),
          m_prepare_func(prepare_func),
          m_run_func(nullptr),
          m_run_lcp_func(nullptr),
          m_run_lcp_cache_func(run_lcp_cache_func)
    { }

    virtual void run();         // implemented in main.cc
    void run_forked();          // implemented in main.cc
    void prepare_run();         // implemented in main.cc

    // implemented in main.cc
    void real_run(
        membuffer<uint8_t*>& stringptr,
        std::vector<uintptr_t>& lcp, std::vector<uint8_t>& charcache);

    virtual bool is_parallel() const { return false; }

    bool is_lcp_func() const { return m_run_lcp_func; }
    bool is_lcp_cache_func() const { return m_run_lcp_cache_func; }
};

#define PSS_CONTESTANT(func, algoname, desc)                         \
    static const class Contestant* _Contestant_ ## func ## _register \
        __attribute__ ((unused)) =                                   \
            new Contestant_UCArray(NULL, func, algoname, desc);

#define PSS_CONTESTANT_PREPARE(pfunc, func, algoname, desc)          \
    static const class Contestant* _Contestant_ ## func ## _register \
        __attribute__ ((unused)) =                                   \
            new Contestant_UCArray(pfunc, func, algoname, desc);

class Contestant_UCArray_Parallel : public Contestant_UCArray
{
public:
    Contestant_UCArray_Parallel(func_type prepare_func,
                                func_type run_func,
                                const char* algoname,
                                const char* description)
        : Contestant_UCArray(prepare_func, run_func, algoname, description)
    { }

    virtual void run(); // implemented in main.cc

    virtual bool is_parallel() const { return true; }
};

#define PSS_CONTESTANT_PARALLEL(func, algoname, desc)                \
    static const class Contestant* _Contestant_ ## func ## _register \
        __attribute__ ((unused)) =                                   \
            new Contestant_UCArray_Parallel(NULL, func, algoname, desc);

#define PSS_CONTESTANT_PARALLEL_PREPARE(pfunc, func, algoname, desc) \
    static const class Contestant* _Contestant_ ## func ## _register \
        __attribute__ ((unused)) =                                   \
            new Contestant_UCArray_Parallel(pfunc, func, algoname, desc);

#else // !PSS_CONTEST

#define PSS_CONTESTANT(func, algoname, desc)
#define PSS_CONTESTANT_PREPARE(pfunc, func, algoname, desc)
#define PSS_CONTESTANT_PARALLEL(func, algoname, desc)
#define PSS_CONTESTANT_PARALLEL_PREPARE(pfunc, func, algoname, desc)

#endif // !PSS_CONTEST

#endif // !PSS_SRC_TOOLS_CONTEST_HEADER

/******************************************************************************/
