#pragma once

#include "../pmt_common.h"

NAMESPACE_PMT

template <typename index_t_, typename sub_t_ = index_t_, size_t n = 1U> // n >= 1
struct Dim
{
  using index_t = index_t_;
  using sub_t = sub_t_;
  using dim_size_t = uint_fast8_t;
  constexpr static dim_size_t n_dims = n;

  static Dim from_index(index_t index, Dim const& dims)
  {
    Dim result;

    for (dim_size_t d = 0; d < n_dims; ++d)
    {
      index_t dsz = dim_sz(dims, d);

      result.xs_[d] = index % dsz;

      index /= dsz;
    }

    return result;
  }

  index_t length() const
  {
    index_t x = 1U;

    for (dim_size_t d = 0U; d < n_dims; ++d)
    {
      x *= dim_sz(d);
    }

    return x;
  }

  index_t index(Dim const& dims) const
  {
    index_t result = 0;
    index_t factor = 1;

    for (dim_size_t d = 0; d < n_dims; ++d)
    {
      result += factor * xs_[d];
      factor *= dim_sz(dims, d);
    }

    return result;
  }

  constexpr index_t dim_sz(dim_size_t d) const
  {
    return dim_sz(*this, d);
  }

  static index_t dim_sz(Dim const& dims, dim_size_t d)
  {
    if (sizeof(index_t) > sizeof(sub_t) && dims.xs_[d] == 0U)
    {
      return index_t(1U) << (sizeof(sub_t) * CHAR_BIT);
    }

    return dims.xs_[d];
  }

  sub_t xs_[n_dims];
};

template <typename index_t, typename sub_t = index_t>
using Dim2 = Dim<index_t, sub_t, 2>;
template <typename index_t, typename sub_t = index_t>
using Dim3 = Dim<index_t, sub_t, 3>;

NAMESPACE_PMT_END
