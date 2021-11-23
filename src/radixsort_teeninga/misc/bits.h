#pragma once

#include <climits>
#include "../pmt_common.h"

NAMESPACE_PMT

constexpr unsigned log2(unsigned int i)
{
  return CHAR_BIT * sizeof(i) - __builtin_clz(i) - 1;
}

constexpr unsigned log2(unsigned long i)
{
  return CHAR_BIT * sizeof(i) - __builtin_clzl(i) - 1;
}

constexpr unsigned log2(unsigned long long i)
{
  return CHAR_BIT * sizeof(i) - __builtin_clzll(i) - 1;
}

template <typename index_t>
constexpr index_t ceil_pow2(index_t i)
{
  return (index_t(1) << (pmt::log2(i) + 1));
}


constexpr unsigned bits_per_word_log2 =
  pmt::log2(sizeof(size_t) * CHAR_BIT);
constexpr unsigned bits_per_word = sizeof(size_t) * CHAR_BIT;
constexpr unsigned bits_word_mask = sizeof(size_t) * CHAR_BIT - 1;

template <typename U>
INLINE void set_bit(U& value, unsigned bit_nr)
{
  value |= U(1) << bit_nr;
}

template <typename U>
INLINE void clear_bit(U& value, unsigned bit_nr)
{
  value &= ~(U(1) << bit_nr);
}

template <typename U>
INLINE bool is_bit_set(U& value, unsigned bit_nr)
{
  return value & (U(1) << bit_nr);
}

NAMESPACE_PMT_END