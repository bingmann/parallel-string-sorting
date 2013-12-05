/******************************************************************************
 * src/tools/input.h
 *
 * Tools to read input files as null-terminated strings.
 *
 ******************************************************************************
 * Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
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

namespace input {

/// Returns true for a valid memory type
bool check_memory_type(const std::string& memtype)
{
    if (memtype == "malloc") return true;
    if (memtype == "mmap") return true;
    if (memtype == "interleaved") return true;
    if (memtype == "node0") return true;
    if (memtype == "segment") return true;
    if (memtype == "mmap_interleaved") return true;
    if (memtype == "mmap_node0") return true;
    if (memtype == "mmap_segment") return true;

    std::cout << "Following --memory types are available:" << std::endl
              << "  malloc        use plain malloc() call (default)" << std::endl
              << "  mmap          use mmap() to allocate private memory" << std::endl
              << "  interleaved   use libnuma to interleave onto nodes" << std::endl
              << "  node0         pin memory to numa node 0" << std::endl
              << "  segment       segment characters equally onto all numa nodes" << std::endl
              << "  mmap_{interleaved,node0,segment} combination of mmap and numa allocation" << std::endl
        ;

    return false;
}

void do_numa_segment(char* buff, size_t buffsize)
{
    int numnodes = numa_num_configured_nodes();
    size_t segsize = (buffsize + numnodes-1) / numnodes;

    std::cout << "Segmenting string characters onto " << numnodes << " NUMA nodes, "
              << segsize << " characters each." << std::endl;

    int pagesize = sysconf(_SC_PAGE_SIZE);

    segsize += pagesize - (segsize % pagesize);
    assert(segsize % pagesize == 0);

    for (int n = 0; n < numnodes; ++n)
    {
        char* p = buff + n * segsize;

        // round p down to pagesize
        p -= (uintptr_t)p % pagesize;

        // segsize need not be page aligned
        numa_tonode_memory(p, segsize, n);
    }
}

/// Free previous data file
void free_stringdata()
{
    if (!g_string_databuff) return;

    if (gopt_memory_type == "mmap" ||
        gopt_memory_type == "mmap_interleaved" ||
        gopt_memory_type == "mmap_node0" ||
        gopt_memory_type == "mmap_segment")
    {
        if (munmap(g_string_databuff, g_string_buffsize)) {
            std::cout << "Error unmapping string data memory: " << strerror(errno) << std::endl;
        }
    }
    else if (gopt_memory_type == "interleaved" ||
             gopt_memory_type == "node0" ||
             gopt_memory_type == "segment")
    {
        numa_free(g_string_databuff, g_string_buffsize);
    }
    else // use plain malloc()/free()
    {
        free(g_string_databuff);
    }

    g_string_databuff = NULL;
}

/// Allocate space in g_string_data
char* allocate_stringdata(size_t size, const std::string& path)
{
    // free previous data file
    free_stringdata();

    // allocate one continuous area of memory
    g_string_buffsize = size + 2 + 8;
    std::cout << "Allocating " << size << " bytes in RAM, reading " << path << std::endl;

    char* stringdata;
    if (gopt_memory_type == "mmap" ||
        gopt_memory_type == "mmap_interleaved" ||
        gopt_memory_type == "mmap_node0" ||
        gopt_memory_type == "mmap_segment")
    {
        stringdata = (char*)mmap(NULL, g_string_buffsize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
        if (stringdata == MAP_FAILED) {
            std::cout << "Error allocating memory: " << strerror(errno) << std::endl;
            return NULL;
        }

        if (gopt_memory_type == "mmap_interleaved")
            numa_interleave_memory(stringdata, g_string_buffsize, numa_all_cpus_ptr);
        if (gopt_memory_type == "mmap_node0")
            numa_tonode_memory(stringdata, g_string_buffsize, 0);
        if (gopt_memory_type == "mmap_segment")
            do_numa_segment(stringdata, g_string_buffsize);
    }
    else if (gopt_memory_type == "interleaved")
    {
        stringdata = (char*)numa_alloc_interleaved(g_string_buffsize);
    }
    else if (gopt_memory_type == "node0")
    {
        stringdata = (char*)numa_alloc_onnode(g_string_buffsize, 0);
    }
    else if (gopt_memory_type == "segment")
    {
        stringdata = (char*)numa_alloc(g_string_buffsize);
        do_numa_segment(stringdata, g_string_buffsize);
    }
    else // use plain malloc()/free()
    {
        stringdata = (char*)malloc(g_string_buffsize);
    }

    // CPL-burstsort needs terminator immediately before and after stringdata
    g_string_databuff = stringdata;
    stringdata[0] = 0;
    stringdata[size+1] = 0;
    ++stringdata;

    g_string_data = stringdata;
    g_string_datasize = size;

    return stringdata;
}

void protect_stringdata()
{
    if (!g_string_databuff) return;

    if (gopt_memory_type == "mmap")
    {
        // protect inherited mmap() area from writes by the algorithms.
        if (mprotect(g_string_databuff, g_string_buffsize, PROT_READ)) {
            std::cout << "Error protecting string data memory: " << strerror(errno) << std::endl;
        }
    }
}

/// Strip data path to just a filename
std::string strip_datapath(const std::string& path)
{
    std::string::size_type slashpos = path.rfind('/');
    std::string name = (slashpos == std::string::npos ? path : path.substr(slashpos+1));

    if ( name.substr(name.size()-3,3) == ".gz" ||
         name.substr(name.size()-4,4) == ".bz2" ||
         name.substr(name.size()-3,3) == ".xz" ||
         name.substr(name.size()-4,4) == ".lzo" )
    {
        // remove compression suffix and size, both separated by dots
        std::string::size_type dotpos = name.rfind('.');
        name.erase(dotpos);
        std::string::size_type dot2pos = name.rfind('.');
        name.erase(dot2pos);
    }

    // check for problems
    if (!name.size()) name = path;

    return name;
}

/// Read a plain file containing newline terminated strings
bool load_plain(const std::string& path)
{
    FILE* file;
    size_t size = 0;

    if (!(file = fopen(path.c_str(), "r"))) {
        std::cout << "Cannot open " << path << ": " << strerror(errno) << std::endl;
        return false;
    }

    if (fseek(file,0,SEEK_END)) {
        std::cout << "Cannot seek in " << path << ": " << strerror(errno) << std::endl;
        fclose(file);
        return false;
    }

    size = ftell(file);
    rewind(file);

    // apply size limit
    if (gopt_inputsize && size > gopt_inputsize)
        size = gopt_inputsize;

    // create memory area
    char* stringdata = allocate_stringdata(size, path);
    if (!stringdata) {
        fclose(file);
        return false;
    }

    // mark offsets of '\n' -> '\0' terminated line starts
    g_string_count = 1;

    // read complete file
    size_t rpos = 0;
    while ( rpos < size )
    {
        size_t batch = std::min<size_t>(8*1024*1024, size - rpos);
        if (batch + rpos > size) batch = size - rpos;

        ssize_t rb = fread(stringdata+rpos, sizeof(char), batch, file);

        if (rb < 0) {
            std::cout << "Cannot read from " << path << ": " << strerror(errno) << std::endl;
            fclose(file);
            return false;
        }

        // iterate over read buffer, identify lines and replace \n -> \0
        if (!gopt_suffixsort)
        {
            for (size_t i = rpos; i < rpos + rb; ++i)
            {
                if (stringdata[i] == '\n' || stringdata[i] == 0) {
                    stringdata[i] = 0;
                    if (i+1 < size) g_string_count++;
                }
            }
        }

        rpos += rb;
    }

    if (gopt_suffixsort) g_string_count = size;

    // force terminatation of last string
    stringdata[ size-1 ] = 0;

    // add more termination
    for (size_t i = size; i < size+9; ++i)
        stringdata[i] = 0;

    fclose(file);

    g_dataname = strip_datapath(path);

    return true;
}

/// Read a compressed file containing newline terminated strings
bool load_compressed(const std::string& path)
{
    if (path.size() < 4) return false;

    const char* decompressor = NULL;

    if ( path.substr(path.size()-3,3) == ".gz" )
        decompressor = "gzip";
    else if ( path.substr(path.size()-4,4) == ".bz2" )
        decompressor = "bzip2";
    else if ( path.substr(path.size()-3,3) == ".xz" )
        decompressor = "xz";
    else if ( path.substr(path.size()-4,4) == ".lzo" )
        decompressor = "lzop";

    if (!decompressor) return false;

    size_t size = 0;

    // extract filesize from filename
    std::string::size_type i = path.rfind('.')-1;
    size_t v = 1;

    while ( isdigit(path[i]) ) {
        size += (path[i] - '0') * v;
        v *= 10; --i;
    }
    if (size == 0 || path[i] != '.') {
        std::cout << "Could not find decompressed size in filename " << path << std::endl;
        return false;
    }

    // apply size limit
    if (gopt_inputsize && size > gopt_inputsize)
        size = gopt_inputsize;

    // create pipe, fork and call decompressor as child
    int pipefd[2]; // pipe[0] = read, pipe[1] = write
    if (pipe(pipefd) != 0) {
        std::cout << "Error creating pipe: " << strerror(errno) << std::endl;
        exit(-1);
    }

    pid_t pid = fork();
    if (pid == 0)
    {
        close(pipefd[0]); // close read end
        dup2(pipefd[1], STDOUT_FILENO); // replace stdout with pipe

        execlp(decompressor, decompressor, "-dc", path.c_str(), NULL);

        std::cout << "Pipe execution failed: " << strerror(errno) << std::endl;
        close(pipefd[1]); // close write end
        exit(-1);
    }

    close(pipefd[1]); // close write end

    // create memory area
    char* stringdata = allocate_stringdata(size, path);
    if (!stringdata) {
        exit(-1);
    }
    g_string_count = 1;

    // read complete file from decompressor child's pipe
    size_t rpos = 0;
    while ( rpos < size )
    {
        size_t batch = std::min<size_t>(8*1024*1024, size - rpos);
        if (batch + rpos > size) batch = size - rpos;

        ssize_t rb = read(pipefd[0], stringdata+rpos, batch);

        if (rb <= 0) {
            std::cout << "Error reading pipe: " << strerror(errno) << std::endl;
            close(pipefd[1]);
            exit(-1);
        }

        // iterate over read buffer, identify lines and replace \n -> \0
        if (!gopt_suffixsort)
        {
            for (size_t i = rpos; i < rpos + rb; ++i)
            {
                if (stringdata[i] == '\n' || stringdata[i] == 0) {
                    stringdata[i] = 0;
                    if (i+1 < size) g_string_count++;
                }
            }
        }

        rpos += rb;
    }

    if (gopt_suffixsort) g_string_count = size;

    // force terminatation of last string
    stringdata[ size-1 ] = 0;

    // add more termination
    for (size_t i = size; i < size+9; ++i)
        stringdata[i] = 0;

    // kill and reap child program
    close(pipefd[1]);

    kill(pid, SIGTERM);
    int status;
    wait(&status);

    g_dataname = strip_datapath(path);
    return true;
}

/// Generate artificial random input with given base letter set.
bool generate_random(const std::string& path, const std::string& letters)
{
    if (!gopt_inputsize) {
        std::cout << "Random input size must be specified via '-s <size>'" << std::endl;
        return false;
    }

    size_t size = gopt_inputsize;

    // create memory area
    char* stringdata = allocate_stringdata(size, path);

    if (!stringdata) {
        return false;
    }

    g_string_count = 0;
    LCGRandom rng(1234567);
    size_t slen = 0;

    for (size_t i = 0; i < size; ++i)
    {
        if (i == slen) { // start of string
            g_string_count++;
            slen += (rng() % 3) + 16;
        }

        if (i+1 == slen) // end of string
            stringdata[i] = 0;
        else
            stringdata[i] = letters[ (rng() / 100) % letters.size() ];
    }

    if (gopt_suffixsort) g_string_count = size;

    // force terminatation of last string
    stringdata[ size-1 ] = 0;

    // add more termination
    for (size_t i = size; i < size+9; ++i)
        stringdata[i] = 0;

    return true;
}

/// Generate artificial random input in exactly the same way as sinha's
/// randomstrings.c does.
bool generate_sinha_randomASCII()
{
    if (!gopt_inputsize) {
        std::cout << "Random input size must be specified via '-s <size>'" << std::endl;
        return false;
    }

    size_t size = gopt_inputsize;

    // create memory area
    char* stringdata = allocate_stringdata(size, "randomASCII");
    if (!stringdata) return false;

    g_string_count = 0;
    srandom(73802);

    size_t slen = (rand() % 20); // excludes zero terminator
    g_string_count++;

    for (size_t i = 0; i < size; ++i)
    {
        if (i == slen) // end of string
        {
            stringdata[i] = 0;
            slen += 1 + (rand() % 20); // includes zero terminator of this
            if (i+1 < size)
                g_string_count++;
        }
        else
        {
            int value = rand() % 127;
            if (value > 32 && value < 127)
                stringdata[i] = value;
            else
                i--;
        }
    }

    if (gopt_suffixsort) g_string_count = size;

    // force terminatation of last string
    stringdata[ size-1 ] = 0;

    // add more termination
    for (size_t i = size; i < size+9; ++i)
        stringdata[i] = 0;

    return true;
}

/// Run through a list of artificial inputs and maybe generate one.
bool load_artifical(const std::string& path)
{
    if (path == "random2") {
        return generate_random("random2", "01");
    }
    else if (path == "random4") {
        return generate_random("random4", "ACGT");
    }
    else if (path == "random10") {
        return generate_random("random10", "0123456789");
    }
    else if (path == "random62") {
        return generate_random("random62", "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
    }
    else if (path == "random255")
    {
        std::string letters(255,0);
        for (int i = 0; i < 255; ++i) letters[i] = (char)(i+1);
        return generate_random("random255", letters);
    }
    else if (path == "randomASCII") {
        return generate_sinha_randomASCII();
    }
    else
        return false;
}

/// Parse a string like "343KB" or "44g" into the corresponding size in bytes
bool parse_filesize(const char* str, size_t& outsize)
{
    char* endptr;
    outsize = strtoul(str,&endptr,10);
    if (!endptr) return false;

    if ( *endptr == 0 || ( (*endptr == 'b' || *endptr == 'B') && *(endptr+1) == 0) )
        outsize *= 1;
    else if ( (*endptr == 'k' || *endptr == 'K') &&
              (*(endptr+1) == 0 || ( (*(endptr+1) == 'b' || *(endptr+1) == 'B') && *(endptr+2) == 0) ) )
        outsize *= 1024;
    else if ( (*endptr == 'm' || *endptr == 'M') &&
              (*(endptr+1) == 0 || ( (*(endptr+1) == 'b' || *(endptr+1) == 'B') && *(endptr+2) == 0) ) )
        outsize *= 1024*1024;
    else if ( (*endptr == 'g' || *endptr == 'G') &&
              (*(endptr+1) == 0 || ( (*(endptr+1) == 'b' || *(endptr+1) == 'B') && *(endptr+2) == 0) ) )
        outsize *= 1024*1024*1024;
    else
        return false;

    return true;
}

/// Load an input set into memory
bool load(const std::string& path)
{
    double ts1 = omp_get_wtime();

    if (load_artifical(path)) {
        g_dataname = path;
    }
    else if (load_compressed(path)) {
    }
    else if (load_plain(path)) {
    }
    else {
        return false;
    }

    double ts2 = omp_get_wtime();

    std::cout << "Loaded input in " << ts2-ts1 << " sec with "
              << (g_string_datasize / (ts2-ts1) / 1024.0 / 1024.0) << " MiB/s" << std::endl;

    protect_stringdata();

    return true;
}

} // namespace input
