#pragma once

#include <cstdlib>
#include <climits>
#include <cstdint>
#include <thread>
#include <cassert>
#include <unistd.h>
#include <random>

#define __attribute__(x)

#define NAMESPACE_PMT namespace pmt {
#define NAMESPACE_PMT_END }

#define ALWAYS_INLINE __attribute__((always_inline))
#if defined(__clang__)
#  define ALWAYS_INL_L(...) ALWAYS_INLINE -> __VA_ARGS__
#else
#  define ALWAYS_INL_L(...) -> __VA_ARGS__ ALWAYS_INLINE
#endif

#define INLINE ALWAYS_INLINE inline
#define NO_INLINE __attribute__ ((noinline))

#define PACKED __attribute__((packed))

#define RESTRICT __restrict__

#define PMT_ASSERT(x__) assert(x__)
//#define PMT_ASSERT(x__) x__;

NAMESPACE_PMT

constexpr unsigned parallel_block_size = 1024 * 16;

constexpr unsigned max_threads = 12U;
constexpr unsigned n_threads = 8;
  // std::min(std::thread::hardware_concurrency(), max_threads);

constexpr unsigned const page_size = 4096U; // sysconf(_SC_PAGESIZE);
constexpr unsigned const mem_alignment = 4096U; // std::min(page_size, 4096U);

NAMESPACE_PMT_END

#include "sort_settings.h"
