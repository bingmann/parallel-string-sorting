# CMake targets for libzio

# Use zlib
include(FindZLIB)
if(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIR})
  add_definitions(-DHAS_ZLIB_H)
  list(APPEND COMPRESSION_LIBRARIES ${ZLIB_LIBRARIES})
endif(ZLIB_FOUND)

# Use bzip2
include(FindBZip2)
if(BZIP2_FOUND)
  include_directories(${BZIP2_INCLUDE_DIR})
  add_definitions(-DHAS_BZLIB_H)
  list(APPEND COMPRESSION_LIBRARIES ${BZIP2_LIBRARIES})
endif(BZIP2_FOUND)

# Use liblzma
include(FindLibLZMA)
if(LIBLZMA_FOUND)
  include_directories(${LIBLZMA_INCLUDE_DIR})
  add_definitions(-DHAS_LZMA_H)
  list(APPEND COMPRESSION_LIBRARIES ${LIBLZMA_LIBRARIES})
endif(LIBLZMA_FOUND)

# Define sources
set(LIBZIO_SOURCES

  zio.c
  zioP.h
  zio.h
  lzw.h
  unlzw.c
)

add_definitions(-DHAVE_FOPENCOOKIE)

add_library(zio ${LIBZIO_SOURCES})
set_target_properties(zio PROPERTIES VERSION 1.00)

target_link_libraries(zio ${COMPRESSION_LIBRARIES}) # append dependencies

# test programs

add_executable(testt testt.c)
target_link_libraries(testt zio)

add_executable(tests tests.c)
target_link_libraries(tests zio)
