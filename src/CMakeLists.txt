# Simple CMake file

include_directories("${PROJECT_SOURCE_DIR}/libzio/")

add_executable(psstest psstest.cc)
target_link_libraries(psstest zio)
