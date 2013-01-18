# Simple CMake file

include_directories("${PROJECT_SOURCE_DIR}/libzio/")
include_directories("${PROJECT_SOURCE_DIR}/liblfds/inc")

add_executable(psstest
  psstest.cc
  rantala/funnelsort.cc
  sequential/nilsson-adaptive.cc
  sequential/nilsson-forward8.cc
  sequential/nilsson-forward16.cc
  parallel/akiba-parallel_string_radix_sort.cc
)
target_link_libraries(psstest zio)

if(WITH_DMALLOC)

  add_definitions(-DDMALLOC)
  target_link_libraries(psstest "dmallocthcxx")

endif(WITH_DMALLOC)
