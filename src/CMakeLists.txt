# Simple CMake file

include_directories("${PROJECT_SOURCE_DIR}/minitbb/")
include_directories("${PROJECT_BINARY_DIR}/")
include_directories("${GMP_INCLUDE_DIR}")
include_directories("${NUMA_INCLUDE_DIR}")

set(SOURCES
  psstest.cc
  sequential/nilsson-adaptive.cc
  sequential/nilsson-forward8.cc
  sequential/nilsson-forward16.cc
  sequential/bingmann-qsort.cc
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
  parallel/shamsundar-lcp-merge-string-sort.cc
)

set(LIBS minitbb ${GMP_LIBRARIES} ${NUMA_LIBRARIES})

if(WITH_MALLOC_COUNT)

  add_definitions(-DMALLOC_COUNT)
  list(APPEND SOURCES
    tools/malloc_count.c
    tools/stack_count.c
  )
  list(APPEND LIBS "dl")

endif(WITH_MALLOC_COUNT)

add_executable(psstest ${SOURCES})
target_link_libraries(psstest ${LIBS})

# silence some warnings from legacy code
set_property(SOURCE sinha-copy-burstsort/glue.cc PROPERTY COMPILE_FLAGS "-Wno-write-strings -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-sign-compare -Wno-sequence-point -Wno-uninitialized")
set_property(SOURCE parallel/shamsundar-lcp-merge-string-sort.cc PROPERTY COMPILE_FLAGS "-Wno-sign-compare")

if(BUILD_ALL)

  foreach(BPARAM 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576)

    add_executable(psstest-${BPARAM} ${SOURCES})
    target_link_libraries(psstest-${BPARAM} ${LIBS})
    set_property(TARGET psstest-${BPARAM} PROPERTY COMPILE_DEFINITIONS "BPARAM=${BPARAM}")

  endforeach(BPARAM)

endif(BUILD_ALL)
