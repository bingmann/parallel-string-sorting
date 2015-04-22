################################################################################
# tests/CMakeLists.txt
#
# CMake file for parallel-string-sorting testsuite
################################################################################

if(PSS_SHARED)
  set(PSS_TESTLIB pss_shared)
  set(ALL_LIBRARIES pss_shared ${PSSBIN_LIBRARIES})
elseif(PSS_STATIC)
  set(PSS_TESTLIB pss_static)
  set(ALL_LIBRARIES pss_static ${PSSBIN_LIBRARIES})
else()
  message(STATUS "No parallel string sorting library built, "
    "hence, no tests can be build.")
endif()

# macro for building test programs
macro(build_only PROGNAME)

  string(REPLACE "/" "_" TESTNAME "${PROGNAME}") # replace slashes

  add_executable(${TESTNAME} ${PROGNAME}.cpp ${ARGN})
  target_link_libraries(${TESTNAME} ${ALL_LIBRARIES})

endmacro(build_only)

# macro for building and running test programs
macro(build_test PROGNAME)

  build_only(${PROGNAME})

  string(REPLACE "/" "_" TESTNAME "${PROGNAME}") # replace slashes
  add_test(${TESTNAME} ${TESTNAME} ${ARGN})

endmacro(build_test)

################################################################################

if(PSS_TESTLIB)
  
  include_directories(${PROJECT_SOURCE_DIR}/src)

  build_test(test-bingmann-simple)
  build_test(test-bingmann-parallel_sample_sort)

endif()

################################################################################