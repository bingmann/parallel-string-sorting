# CMake targets for liblfds

include_directories("${PROJECT_SOURCE_DIR}/liblfds/inc" "${PROJECT_SOURCE_DIR}/liblfds/src")

# Define sources
set(LIBLFDS_SOURCES

  inc/liblfds611.h
  src/liblfds611_internal.h
  #src/lfds611_abstraction/lfds611_abstraction_cas.c
  #src/lfds611_abstraction/lfds611_abstraction_dcas.c
  src/lfds611_abstraction/lfds611_abstraction_free.c
  #src/lfds611_abstraction/lfds611_abstraction_increment.c
  src/lfds611_abstraction/lfds611_abstraction_internal_body.h
  src/lfds611_abstraction/lfds611_abstraction_internal_wrapper.h
  src/lfds611_abstraction/lfds611_abstraction_malloc.c
  src/lfds611_freelist/lfds611_freelist_delete.c
  src/lfds611_freelist/lfds611_freelist_get_and_set.c
  src/lfds611_freelist/lfds611_freelist_internal.h
  src/lfds611_freelist/lfds611_freelist_new.c
  src/lfds611_freelist/lfds611_freelist_pop_push.c
  src/lfds611_freelist/lfds611_freelist_query.c
  src/lfds611_liblfds/lfds611_liblfds_abstraction_test_helpers.c
  src/lfds611_liblfds/lfds611_liblfds_aligned_free.c
  src/lfds611_liblfds/lfds611_liblfds_aligned_malloc.c
  src/lfds611_liblfds/lfds611_liblfds_internal.h
  src/lfds611_queue/lfds611_queue_delete.c
  src/lfds611_queue/lfds611_queue_internal.h
  src/lfds611_queue/lfds611_queue_new.c
  src/lfds611_queue/lfds611_queue_query.c
  src/lfds611_queue/lfds611_queue_queue.c
  src/lfds611_ringbuffer/lfds611_ringbuffer_delete.c
  src/lfds611_ringbuffer/lfds611_ringbuffer_get_and_put.c
  src/lfds611_ringbuffer/lfds611_ringbuffer_internal.h
  src/lfds611_ringbuffer/lfds611_ringbuffer_new.c
  src/lfds611_ringbuffer/lfds611_ringbuffer_query.c
  src/lfds611_slist/lfds611_slist_delete.c
  src/lfds611_slist/lfds611_slist_get_and_set.c
  src/lfds611_slist/lfds611_slist_internal.h
  src/lfds611_slist/lfds611_slist_link.c
  src/lfds611_slist/lfds611_slist_new.c
  src/lfds611_stack/lfds611_stack_delete.c
  src/lfds611_stack/lfds611_stack_internal.h
  src/lfds611_stack/lfds611_stack_new.c
  src/lfds611_stack/lfds611_stack_push_pop.c
  src/lfds611_stack/lfds611_stack_query.c
)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -march=core2")

add_library(lfds ${LIBLFDS_SOURCES})
set_target_properties(lfds PROPERTIES VERSION 6.11)
