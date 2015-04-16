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
  sequential/nilsson-adaptive.cc
  sequential/nilsson-forward8.cc
  sequential/nilsson-forward16.cc
  sequential/bingmann-qsort.cc
  sequential/bingmann-lcp-mergesort.cc
  rantala/burstsort.cc
  rantala/burstsort2.cc
  rantala/burstsort_mkq.cc
  rantala/funnelsort.cc
  rantala/mergesort.cc
  rantala/mergesort_lcp.cc
  rantala/mergesort_losertree.cc
  rantala/mergesort_unstable.cc
  rantala/multikey_block.cc
  rantala/multikey_cache.cc
  rantala/multikey_dynamic.cc
  rantala/multikey_multipivot.cc
  rantala/multikey_simd.cc
  sinha-copy-burstsort/glue.cc
  parallel/akiba-parallel_string_radix_sort.cc
  parallel/bingmann-parallel_radix_sort.cc
  parallel/bingmann-parallel_mkqs.cc
  parallel/bingmann-parallel_sample_sort.cc
  parallel/bingmann-parallel_sample_sort_lcp.cc
  parallel/shamsundar-lcp-merge-string-sort.cc
  )

set(PSSBIN_SOURCES
  psstest.cc
  ${PSSLIB_SOURCES}
  )

set(LIBS minitbb ${GMP_LIBRARIES} ${NUMA_LIBRARIES} rt)

if(WITH_MALLOC_COUNT)

  add_definitions(-DMALLOC_COUNT)
  list(APPEND PSSBIN_SOURCES
    tools/malloc_count.c
    tools/stack_count.c
  )
  list(APPEND LIBS "dl")

endif(WITH_MALLOC_COUNT)

# build benchmark binary from sources (without library)
if(PSS_BINARY)

  add_executable(psstest ${PSSBIN_SOURCES})
  target_link_libraries(psstest ${LIBS})

  if(BUILD_ALL)
    # build parallel string sorting with additional parameters

    foreach(BPARAM 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576)

      add_executable(psstest-${BPARAM} ${PSSBIN_SOURCES})
      target_link_libraries(psstest-${BPARAM} ${LIBS})
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
set_property(SOURCE sinha-copy-burstsort/glue.cc
  PROPERTY COMPILE_FLAGS "-Wno-write-strings -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-sign-compare -Wno-sequence-point -Wno-uninitialized")
set_property(SOURCE parallel/shamsundar-lcp-merge-string-sort.cc
  PROPERTY COMPILE_FLAGS "-Wno-sign-compare")

################################################################################
