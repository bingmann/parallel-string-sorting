# CMake file for parallel-string-sorting project

project(parallel-string-sorting)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(WITH_MALLOC_COUNT "Use malloc_count for memory profiling" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# disallow in-source builds

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# read .git directory (if it exists) and find git sha

include(GetGitRevisionDescription)
get_git_head_revision(GIT_VERSION_REFSPEC GIT_VERSION_SHA1)
if(GIT_VERSION_REFSPEC)
  message(STATUS "Detected git refspec ${GIT_VERSION_REFSPEC} sha ${GIT_VERSION_SHA1}")
  add_definitions("-DGIT_VERSION_SHA1=\"${GIT_VERSION_SHA1}\"")
endif()

# compilation options 

option(BUILD_ALL "Build different test programs" OFF)

# Enable warnings

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Drestrict=__restrict__ -std=c++0x")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -g")

# Use OpenMP

include(FindOpenMP)

if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(FATAL_ERROR "OpenMP not found.")
endif()

# Find <atomic> header

include(CheckIncludeFileCXX)
check_include_file_cxx(atomic HAVE_ATOMIC_H)
check_include_file_cxx(cstdatomic HAVE_CSTDATOMIC_H)

# Find GMP libraries and header

find_path(GMP_INCLUDE_DIR NAMES gmpxx.h)
find_library(GMP_LIBRARIES NAMES gmp gmp)
find_library(GMPXX_LIBRARIES NAMES gmpxx gmpxx)
if(NOT GMP_LIBRARIES)
  message(FATAL_ERROR "gmpxx library not found!")
endif()

message(STATUS "GMP library: " ${GMP_LIBRARIES} " " ${GMPXX_LIBRARIES})
set(GMP_LIBRARIES ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARIES)

# check for libnuma

find_path(NUMA_INCLUDE_DIR NAMES numa.h)
find_library(NUMA_LIBRARIES NAMES numa_available numa)
if(NOT NUMA_LIBRARIES)
  message(FATAL_ERROR "libnuma library not found!")
endif()
message(STATUS "libnuma library: " ${NUMA_LIBRARIES})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NUMA DEFAULT_MSG NUMA_INCLUDE_DIR NUMA_LIBRARIES)

# build subset of Intel TBB
add_subdirectory(minitbb)

# descend into source
add_subdirectory(src)

# write platform specific configuration
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)
