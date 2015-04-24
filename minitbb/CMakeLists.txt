# CMake targets for subset of Intel TBB

include_directories("${PROJECT_SOURCE_DIR}/minitbb/")
include_directories("${PROJECT_SOURCE_DIR}/minitbb/src/")

# Define sources
set(LIBMINITBB_SOURCES

  src/tbb/cache_aligned_allocator.cpp
  src/tbb/concurrent_monitor.cpp
  src/tbb/concurrent_monitor.h
  src/tbb/concurrent_queue.cpp
  src/tbb/dynamic_link.cpp
  src/tbb/dynamic_link.h
  src/tbb/itt_notify.cpp
  src/tbb/itt_notify.h
  src/tbb/semaphore.cpp
  src/tbb/semaphore.h
  src/tbb/tbb_assert_impl.h
  src/tbb/tbb_misc.cpp
  src/tbb/tbb_misc.h
  src/tbb/tbb_version.h
  src/tbb/version_string.ver

  tbb/aligned_space.h
  tbb/atomic.h
  tbb/cache_aligned_allocator.h
  tbb/concurrent_queue.h
  tbb/internal/_concurrent_queue_impl.h
  tbb/internal/_mutex_padding.h
  tbb/internal/_tbb_strings.h
  tbb/internal/_x86_eliding_mutex_impl.h
  tbb/machine/gcc_armv7.h
  tbb/machine/gcc_generic.h
  tbb/machine/gcc_ia32_common.h
  tbb/machine/gcc_itsx.h
  tbb/machine/ibm_aix51.h
  tbb/machine/icc_generic.h
  tbb/machine/linux_common.h
  tbb/machine/linux_ia32.h
  tbb/machine/linux_ia64.h
  tbb/machine/linux_intel64.h
  tbb/machine/mac_ppc.h
  tbb/machine/macos_common.h
  tbb/machine/mic_common.h
  tbb/machine/msvc_armv7.h
  tbb/machine/msvc_ia32_common.h
  tbb/machine/sunos_sparc.h
  tbb/machine/windows_api.h
  tbb/machine/windows_ia32.h
  tbb/machine/windows_intel64.h
  tbb/machine/xbox360_ppc.h
  tbb/spin_mutex.h
  tbb/tbb_allocator.h
  tbb/tbb_config.h
  tbb/tbb_exception.h
  tbb/tbb_machine.h
  tbb/tbb_profiling.h
  tbb/tbb_stddef.h
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__TBB_BUILD=1 -DUSE_PTHREAD")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DTBB_USE_DEBUG")

add_library(minitbb ${LIBMINITBB_SOURCES})
set_target_properties(minitbb PROPERTIES VERSION 20121003)
