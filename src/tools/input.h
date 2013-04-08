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

/// Allocate space in g_string_data
char* allocate_stringdata(size_t size, const std::string& path)
{
    // free previous data file
    if (g_string_databuff)
        free(g_string_databuff);

    // allocate one continuous area of memory
    std::cout << "Allocating " << size << " bytes in RAM, reading " << path << "\n";
    char* stringdata = (char*)malloc(size + 2 + 8);

    // CPL-burstsort needs terminator immediately before and after stringdata
    g_string_databuff = stringdata;
    stringdata[0] = 0;
    stringdata[size+1] = 0;
    ++stringdata;

    g_string_data = stringdata;
    g_string_datasize = size;

    if (!stringdata) {
        std::cout << "Error allocating memory: " << strerror(errno) << "\n";
        return NULL;
    }

    return stringdata;
}

/// Read a plain file containing newline terminated strings
bool load_plain(const std::string& path)
{
    FILE* file;
    size_t size = 0;

    if (!(file = fopen(path.c_str(), "r"))) {
        std::cout << "Cannot open file: " << strerror(errno) << "\n";
        return false;
    }

    if (fseek(file,0,SEEK_END)) {
        std::cout << "Cannot seek in file: " << strerror(errno) << "\n";
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
    g_string_count = 0;

    // read complete file
    size_t rpos = 0;
    while ( rpos < size )
    {
        size_t batch = std::min<size_t>(8*1024*1024, size - rpos);
        if (batch + rpos > size) batch = size - rpos;

        ssize_t rb = fread(stringdata+rpos, sizeof(char), batch, file);

        if (rb < 0) {
            std::cout << "Cannot read from file: " << strerror(errno) << "\n";
            fclose(file);
            return false;
        }

        // iterate over read buffer, identify lines and replace \n -> \0
        for (size_t i = rpos; i < rpos + rb; ++i)
        {
            if (gopt_suffixsort) {
                g_string_count++;
            }
            else
            {
                if (stringdata[i] == '\n') {
                    stringdata[i] = 0;
                    if (i+1 < size)
                        g_string_count++;
                }
            }
        }

        rpos += rb;
    }

    // force terminatation of last string
    stringdata[ size-1 ] = 0;

    // add more termination
    for (size_t i = size; i < size+9; ++i)
        stringdata[i] = 0;
  
    fclose(file);

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
        std::cout << "Could not find decompressed size in filename " << path << "\n";
        return false;
    }

    // apply size limit
    if (gopt_inputsize && size > gopt_inputsize)
        size = gopt_inputsize;

    // create pipe, fork and call decompressor as child
    int pipefd[2]; // pipe[0] = read, pipe[1] = write
    if (pipe(pipefd) != 0) {
        std::cout << "Error creating pipe: " << strerror(errno) << "\n";
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
    g_string_count = 0;

    // read complete file from decompressor child's pipe
    size_t rpos = 0;
    while ( rpos < size )
    {
        size_t batch = std::min<size_t>(8*1024*1024, size - rpos);
        if (batch + rpos > size) batch = size - rpos;

        ssize_t rb = read(pipefd[0], stringdata+rpos, batch);

        if (rb <= 0) {
            std::cout << "Error reading pipe: " << strerror(errno) << "\n";
            close(pipefd[1]);
            exit(-1);
        }

        // iterate over read buffer, identify lines and replace \n -> \0
        for (size_t i = rpos; i < rpos + rb; ++i)
        {
            if (gopt_suffixsort) {
                g_string_count++;
            }
            else
            {
                if (stringdata[i] == '\n') {
                    stringdata[i] = 0;
                    if (i+1 < size)
                        g_string_count++;
                }
            }
        }

        rpos += rb;
    }

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

    return true;
}

/// Generate artificial random input with given base letter set.
bool generate_random(const std::string& path, const std::string& letters)
{
    if (!gopt_inputsize) {
        std::cout << "Random input size must be specified via '-s <size>'\n";
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
            stringdata[i] = letters[ rng() % letters.size() ];
    }

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
        std::cout << "Random input size must be specified via '-s <size>'\n";
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
    if (path == "random4") {
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

    return true;
}

} // namespace input
