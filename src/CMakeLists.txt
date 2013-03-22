# Simple CMake file

include_directories("${PROJECT_SOURCE_DIR}/libzio/")
include_directories("${PROJECT_SOURCE_DIR}/liblfds/inc/")
include_directories("${PROJECT_BINARY_DIR}/")
include_directories("${GMP_INCLUDE_DIR}")

set(SOURCES
  psstest.cc
  sequential/nilsson-adaptive.cc
  sequential/nilsson-forward8.cc
  sequential/nilsson-forward16.cc
  rantala/funnelsort.cc
  sinha-copy-burstsort/glue.cc
  parallel/akiba-parallel_string_radix_sort.cc
  parallel/bingmann-parallel_radix_sort1.cc
  parallel/bingmann-parallel_radix_sort2.cc
  parallel/bingmann-parallel_radix_sort3.cc
  parallel/bingmann-parallel_sample_sortBS.cc
  parallel/bingmann-parallel_sample_sortBSC.cc
  parallel/bingmann-parallel_sample_sortBT.cc
  parallel/bingmann-parallel_sample_sortBTC.cc
)

set(LIBS zio lfds ${GMP_LIBRARIES})

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
