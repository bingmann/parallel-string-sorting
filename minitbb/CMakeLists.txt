# CMake targets for subset of Intel TBB

include_directories("${PROJECT_SOURCE_DIR}/minitbb/")

# Define sources
set(LIBMINITBB_SOURCES

  src/concurrent_monitor.cpp
  src/concurrent_monitor.h
  src/concurrent_queue.cpp
  src/itt_notify.cpp
  src/itt_notify.h
  src/semaphore.cpp
  src/semaphore.h
  src/version_string.ver

  tbb/aligned_space.h
  tbb/atomic.h
  tbb/cache_aligned_allocator.h
  tbb/concurrent_queue.h
  tbb/internal/_concurrent_queue_impl.h
  tbb/machine/gcc_generic.h
  tbb/machine/gcc_ia32_common.h
  tbb/machine/ibm_aix51.h
  tbb/machine/icc_generic.h
  tbb/machine/linux_common.h
  tbb/machine/linux_ia32.h
  tbb/machine/linux_ia64.h
  tbb/machine/linux_intel64.h
  tbb/machine/mac_ppc.h
  tbb/machine/macos_common.h
  tbb/machine/mic_common.h
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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__TBB_BUILD=1 -DUSE_PTHREAD -D__TBB_DYNAMIC_LOAD_ENABLED=0")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DTBB_USE_DEBUG")

add_library(minitbb ${LIBMINITBB_SOURCES})
set_target_properties(minitbb PROPERTIES VERSION 20121003)
