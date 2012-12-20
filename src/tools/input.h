/******************************************************************************
 * src/tools/input.h
 *
 * Tools to read input files as null-terminated strings.
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

/// Read a compressed file containing newline terminated strings
bool readfile_lines(const std::string& path)
{
    FILE* file;
    size_t filesize = 0;

    if (path.size() >= 4 &&
        ( path.substr(path.size()-3,3) == ".gz" ||
          path.substr(path.size()-4,4) == ".bz2" ||
          path.substr(path.size()-3,3) == ".xz" )
        )
    {
        // extract filesize from filename
        std::string::size_type i = path.rfind('.')-1;
        size_t v = 1;

        while ( isdigit(path[i]) ) {
            filesize += (path[i] - '0') * v;
            v *= 10; --i;
        }
        if (filesize == 0 || path[i] != '.') {
            std::cerr << "\nCould not find decompressed size in filename " << path << "\n";
            return false;
        }

        if (!(file = fzopen(path.c_str(), "r"))) {
            std::cerr << "\n" << strerror(errno) << "\n";
            return false;
        }
    }
    else
    {
        if (!(file = fopen(path.c_str(), "r"))) {
            std::cerr << "\n" << strerror(errno) << "\n";
            return false;
        }

        if (fseek(file,0,SEEK_END)) {
            std::cerr << "\n" << strerror(errno) << "\n";
            fclose(file);
            return false;
        }

        filesize = ftell(file);

        rewind(file);
    }

    // apply size limit
    if (gopt_inputsizelimit && filesize > gopt_inputsizelimit)
        filesize = gopt_inputsizelimit;

    // free previous data file
    if (g_stringdata) free( (char*)g_stringdata );

    // allocate one continuous area of memory
    std::cerr << "Allocating " << filesize << " bytes in RAM, reading " << path << "\n";
    char* stringdata = (char*)malloc(filesize);

    g_stringdata = stringdata;
    g_stringdatasize = filesize;

    if (!stringdata) {
        std::cerr << "\n" << strerror(errno) << "\n";
        fclose(file);
        return false;
    }

    // mark offsets of '\n' -> '\0' terminated line starts
    g_stringoffsets.clear();
    g_stringoffsets.push_back(0);

    // read complete file
    size_t rpos = 0;
    while ( rpos < filesize )
    {
        size_t batch = std::min<size_t>(8*1024*1024, filesize - rpos);
        if (batch + rpos > filesize) batch = filesize - rpos;

        ssize_t rb = fread(stringdata+rpos, sizeof(char), batch, file);

        if (rb < 0) {
            std::cerr << "\n" << strerror(errno) << "\n";
            fclose(file);
            return false;
        }

        // iterate over read buffer, identify lines and replace \n -> \0
        for (size_t i = rpos; i < rpos + rb; ++i)
        {
            if (stringdata[i] == '\n') {
                stringdata[i] = 0;
                if (i+1 < filesize)
                    g_stringoffsets.push_back(i+1);
            }
        }

        rpos += rb;
    }

    // force terminatation of last string
    stringdata[ filesize-1 ] = 0;
    
    fclose(file);

    return true;
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
