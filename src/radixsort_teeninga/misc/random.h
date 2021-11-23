#pragma once

#include "../pmt_common.h"
#include "../pcg/pcg_random.hpp"
#include "../misc/bits.h"

NAMESPACE_PMT

pcg32_unique random_type(uint8_t u);
pcg32_unique random_type(uint16_t u);
pcg32_unique random_type(uint32_t u);
pcg64_unique random_type(uint64_t u);

template <typename U>
struct rng
{
  using type = decltype(random_type(U(0U)));
};

template <typename Random, typename U>
U random_uint(Random& rng, U max_val)
{
  U const mask = ceil_pow2(max_val) - 1;
  U next;

  do
  {
    next = rng() & mask;
  } while (next > max_val);

  return next;
}

NAMESPACE_PMT_END
