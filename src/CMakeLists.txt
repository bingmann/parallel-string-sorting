################################################################################
# src/CMakeLists.txt
#
# Simple CMake file
################################################################################

include_directories("${PROJECT_SOURCE_DIR}/minitbb/")
include_directories("${PROJECT_BINARY_DIR}/")
include_directories("${GMP_INCLUDE_DIR}")
include_directories("${NUMA_INCLUDE_DIR}")

set(PSSLIB_SOURCES
  tools/globals.cpp
  sequential/nilsson-adaptive.cpp
  sequential/nilsson-forward8.cpp
  sequential/nilsson-forward16.cpp
  sequential/bingmann-qsort.cpp
  sequential/bingmann-lcp-mergesort.cpp
  rantala/burstsort.cpp
  rantala/burstsort2.cpp
  rantala/burstsort_mkq.cpp
  rantala/funnelsort.cpp
  rantala/mergesort.cpp
  rantala/mergesort_lcp.cpp
  rantala/mergesort_losertree.cpp
  rantala/mergesort_unstable.cpp
  rantala/multikey_block.cpp
  rantala/multikey_cache.cpp
  rantala/multikey_dynamic.cpp
  rantala/multikey_multipivot.cpp
  rantala/multikey_simd.cpp
  parallel/akiba-parallel_string_radix_sort.cpp
  parallel/bingmann-parallel_radix_sort.cpp
  parallel/bingmann-parallel_mkqs.cpp
  parallel/shamsundar-lcp-merge-string-sort.cpp
  )

set(PSSBIN_SOURCES
  psstest.cpp
  ${PSSLIB_SOURCES}
  sinha-copy-burstsort/glue.cpp
  )

set(PSSBIN_LIBRARIES ${TBB_LIBRARIES} ${GMP_LIBRARIES} ${NUMA_LIBRARIES} rt dl)

set(PSSBIN_LIBRARIES ${PSSBIN_LIBRARIES} CACHE STRING "Pssbin Libraries" FORCE)

if(WITH_MALLOC_COUNT)

  add_definitions(-DMALLOC_COUNT)
  list(APPEND PSSBIN_SOURCES
    tools/malloc_count.c
    tools/stack_count.c
  )

endif(WITH_MALLOC_COUNT)

# build benchmark binary from sources (without library)
if(PSS_BINARY)

  add_executable(psstest ${PSSBIN_SOURCES})
  target_link_libraries(psstest ${PSSBIN_LIBRARIES})

  # enable compilation of hooks for pss contest list
  set_property(TARGET psstest
    PROPERTY COMPILE_DEFINITIONS "PSS_CONTEST=1")

  if(BUILD_ALL)
    # build parallel string sorting with additional parameters

    foreach(BPARAM 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576)

      add_executable(psstest-${BPARAM} ${PSSBIN_SOURCES})
      target_link_libraries(psstest-${BPARAM} ${PSSBIN_LIBRARIES})
      set_property(TARGET psstest-${BPARAM} PROPERTY COMPILE_DEFINITIONS "BPARAM=${BPARAM}")

    endforeach(BPARAM)

  endif(BUILD_ALL)

endif(PSS_BINARY)

# build parallel string sorting static library from sources (possibly again)
if(PSS_STATIC)

  add_library(pss_static STATIC ${PSSLIB_SOURCES})

  set_target_properties(pss_static PROPERTIES
    VERSION   "${PSSLIB_VERSION}"
    SOVERSION "${PSSLIB_SOVERSION}")

endif()

# build parallel string sorting shared library from sources (possibly again)
if(PSS_SHARED)

  add_library(pss_shared SHARED ${PSSLIB_SOURCES})

  set_target_properties(pss_shared PROPERTIES
    VERSION   "${PSSLIB_VERSION}"
    SOVERSION "${PSSLIB_SOVERSION}")

endif()

# silence some warnings from legacy code
set_property(SOURCE sinha-copy-burstsort/glue.cpp
  PROPERTY COMPILE_FLAGS "-Wno-write-strings -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-sign-compare -Wno-sequence-point -Wno-uninitialized")
set_property(SOURCE parallel/shamsundar-lcp-merge-string-sort.cpp
  PROPERTY COMPILE_FLAGS "-Wno-sign-compare")

################################################################################
