# Simple CMake file

include_directories("${PROJECT_SOURCE_DIR}/libzio/")
include_directories("${PROJECT_SOURCE_DIR}/liblfds/inc/")
include_directories("${PROJECT_BINARY_DIR}/")

set(SOURCES
  psstest.cc
  rantala/funnelsort.cc
  sequential/nilsson-adaptive.cc
  sequential/nilsson-forward8.cc
  sequential/nilsson-forward16.cc
  parallel/akiba-parallel_string_radix_sort.cc
  parallel/bingmann-parallel_radix_sort1.cc
  parallel/bingmann-parallel_radix_sort2.cc
  parallel/bingmann-parallel_radix_sort3.cc
  parallel/bingmann-parallel_sample_sortBS.cc
  parallel/bingmann-parallel_sample_sortBSC.cc
  parallel/bingmann-parallel_sample_sortBT.cc
  parallel/bingmann-parallel_sample_sortBTC.cc
)

set(LIBS zio lfds)

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
