# Simple CMake file

include_directories("${PROJECT_SOURCE_DIR}/libzio/")

add_executable(psstest
  psstest.cc
  sequential/nilsson-adaptive.cc
  sequential/nilsson-forward8.cc
  sequential/nilsson-forward16.cc
)
target_link_libraries(psstest zio)

if(WITH_DMALLOC)

  add_definitions(-DDMALLOC)
  target_link_libraries(psstest "dmallocthcxx")

endif(WITH_DMALLOC)
