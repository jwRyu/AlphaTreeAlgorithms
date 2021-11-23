#pragma once

#include "../pmt_common.h"

NAMESPACE_PMT

template <typename RAI, typename Index>
void copy_rai(RAI in, Index n, RAI out)
{
  RAI const past_last = in + n;

  for (; in < past_last;)
  {
    *out++ = *in++;
  }
}

NAMESPACE_PMT_END