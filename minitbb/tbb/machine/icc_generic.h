/*
    Copyright 2005-2012 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#if !defined(__TBB_machine_H) || defined(__TBB_machine_icc_generic_H)
#error Do not #include this internal file directly; use public TBB headers instead.
#endif

#if ! __TBB_ICC_BUILTIN_ATOMICS_PRESENT
    #error "Intel C++ Compiler of at least 12.1 version is needed to use ICC intrinsics port"
#endif

#define __TBB_machine_icc_generic_H

//ICC mimics the "native" target compiler
#if _MSC_VER
    #include "msvc_ia32_common.h"
#else
    #include "gcc_ia32_common.h"
#endif

//TODO: Make __TBB_WORDSIZE macro optional for ICC intrinsics port.
//As compiler intrinsics are used for all the operations it is possible to do.

#if __TBB_x86_32
    #define __TBB_WORDSIZE 4
#else
    #define __TBB_WORDSIZE 8
#endif
#define __TBB_BIG_ENDIAN 0

//__TBB_compiler_fence() defined just in case, as it seems not to be used on its own anywhere else
#if _MSC_VER
    //TODO: any way to use same intrinsics on windows and linux?
    #pragma intrinsic(_ReadWriteBarrier)
    #pragma intrinsic(_mm_mfence)
    #define __TBB_compiler_fence()    _ReadWriteBarrier()
    #define __TBB_full_memory_fence() _mm_mfence()
#else
    #define __TBB_compiler_fence()    __asm__ __volatile__("": : :"memory")
    #define __TBB_full_memory_fence() __asm__ __volatile__("mfence": : :"memory")
#endif

#define __TBB_control_consistency_helper() __TBB_compiler_fence()

namespace tbb { namespace internal {
//TODO: is there any way to reuse definition of memory_order enum from ICC instead of copy paste.
//however it seems unlikely that ICC will silently change exact enum values, as they are defined
//in the ISO exactly like this.
//TODO: add test that exact values of the enum are same as in the ISO C++11
typedef enum memory_order {
    memory_order_relaxed, memory_order_consume, memory_order_acquire,
    memory_order_release, memory_order_acq_rel, memory_order_seq_cst
} memory_order;

template <typename T, size_t S>
struct machine_load_store {
    static T load_with_acquire ( const volatile T& location ) {
        return __atomic_load_explicit(&location, memory_order_acquire);
    }
    static void store_with_release ( volatile T &location, T value ) {
        __atomic_store_explicit(&location, value, memory_order_release);
    }
};
//The specializations below are needed to have explicit conversion of pointer to void* in argument list.
//compiler bug?
//TODO: move to workaround to separate layer by wrapping the call to __atomic_store_explicit in overloaded template function
template <typename T, size_t S>
struct machine_load_store<T*, S> {
    static T* load_with_acquire ( T* const volatile & location ) {
        return __atomic_load_explicit(&location, memory_order_acquire);
    }
    static void store_with_release ( T* volatile &location, T* value ) {
        __atomic_store_explicit(&location, (void*) value, memory_order_release);
    }
};
template <typename T, size_t S>
struct machine_load_store_relaxed {
    static inline T load ( const T& location ) {
        return __atomic_load_explicit(&location, memory_order_relaxed);
    }
    static inline void store (  T& location, T value ) {
        __atomic_store_explicit(&location, value, memory_order_relaxed);
    }
};

template <typename T, size_t S>
struct machine_load_store_relaxed<T*, S> {
    static inline T* load ( T* const& location ) {
        return (T*)__atomic_load_explicit(&location, memory_order_relaxed);
    }
    static inline void store (  T*& location, T* value ) {
        __atomic_store_explicit(&location, (void*)value, memory_order_relaxed);
    }
};
template <typename T, size_t S>
struct machine_load_store_seq_cst {
    static T load ( const volatile T& location ) {
        return __atomic_load_explicit(&location, memory_order_seq_cst);
    }

    static void store ( volatile T &location, T value ) {
        __atomic_store_explicit(&location, value, memory_order_seq_cst);
    }
};
}} // namespace tbb::internal

namespace tbb{ namespace internal { namespace icc_intrinsics_port{
    typedef enum memory_order_map {
        relaxed = memory_order_relaxed,
        acquire = memory_order_acquire,
        release = memory_order_release,
        full_fence=  memory_order_seq_cst
    } memory_order_map;
}}}// namespace tbb::internal

#define __TBB_MACHINE_DEFINE_ATOMICS(S,T,M)                                                     \
inline T __TBB_machine_cmpswp##S##M( volatile void *ptr, T value, T comparand ) {               \
    __atomic_compare_exchange_strong_explicit(                                                  \
            (T*)ptr                                                                             \
            ,&comparand                                                                         \
            ,value                                                                              \
            , tbb::internal::icc_intrinsics_port::M                                             \
            , tbb::internal::icc_intrinsics_port::M);                                           \
    return comparand;                                                                           \
}                                                                                               \
                                                                                                \
inline T __TBB_machine_fetchstore##S##M(volatile void *ptr, T value) {                          \
    return __atomic_exchange_explicit((T*)ptr, value, tbb::internal::icc_intrinsics_port::M);   \
}                                                                                               \
                                                                                                \
inline T __TBB_machine_fetchadd##S##M(volatile void *ptr, T value) {                            \
    return __atomic_fetch_add_explicit((T*)ptr, value, tbb::internal::icc_intrinsics_port::M);  \
}                                                                                               \

__TBB_MACHINE_DEFINE_ATOMICS(1,tbb::internal::int8_t, full_fence)
__TBB_MACHINE_DEFINE_ATOMICS(1,tbb::internal::int8_t, acquire)
__TBB_MACHINE_DEFINE_ATOMICS(1,tbb::internal::int8_t, release)
__TBB_MACHINE_DEFINE_ATOMICS(1,tbb::internal::int8_t, relaxed)

__TBB_MACHINE_DEFINE_ATOMICS(2,tbb::internal::int16_t, full_fence)
__TBB_MACHINE_DEFINE_ATOMICS(2,tbb::internal::int16_t, acquire)
__TBB_MACHINE_DEFINE_ATOMICS(2,tbb::internal::int16_t, release)
__TBB_MACHINE_DEFINE_ATOMICS(2,tbb::internal::int16_t, relaxed)

__TBB_MACHINE_DEFINE_ATOMICS(4,tbb::internal::int32_t, full_fence)
__TBB_MACHINE_DEFINE_ATOMICS(4,tbb::internal::int32_t, acquire)
__TBB_MACHINE_DEFINE_ATOMICS(4,tbb::internal::int32_t, release)
__TBB_MACHINE_DEFINE_ATOMICS(4,tbb::internal::int32_t, relaxed)

__TBB_MACHINE_DEFINE_ATOMICS(8,tbb::internal::int64_t, full_fence)
__TBB_MACHINE_DEFINE_ATOMICS(8,tbb::internal::int64_t, acquire)
__TBB_MACHINE_DEFINE_ATOMICS(8,tbb::internal::int64_t, release)
__TBB_MACHINE_DEFINE_ATOMICS(8,tbb::internal::int64_t, relaxed)


#undef __TBB_MACHINE_DEFINE_ATOMICS

#define __TBB_USE_FENCED_ATOMICS                            1

template <typename T>
inline void __TBB_machine_OR( T *operand, T addend ) {
    __atomic_fetch_or_explicit(operand, addend, tbb::internal::memory_order_seq_cst);
}

template <typename T>
inline void __TBB_machine_AND( T *operand, T addend ) {
    __atomic_fetch_and_explicit(operand, addend, tbb::internal::memory_order_seq_cst);
}
