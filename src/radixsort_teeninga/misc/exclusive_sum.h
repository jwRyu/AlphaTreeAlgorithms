#pragma once

#include "../pmt_common.h"

NAMESPACE_PMT

template <typename Iter>
void exclusive_sum(Iter iter, Iter iter_end)
{
  using Value = typename std::iterator_traits<Iter>::value_type;

  Value sum = 0;
  Value tmp;

  for (; iter < iter_end; ++iter)
  {
    tmp = *iter;
    *iter = sum;
    sum += tmp;
  }
}

NAMESPACE_PMT_END