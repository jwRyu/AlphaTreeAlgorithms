#pragma once

#include "dim.h"
 
NAMESPACE_PMT

template <typename D>
struct BlockRange
{
  using dim_t = D;
  using index_t = typename D::index_t;
  using sub_t = typename D::sub_t;
  using dim_size_t = typename D::dim_size_t;

  constexpr static dim_size_t n_dims = D::n_dims;

  BlockRange() {}
  BlockRange(D const& dims)
  {
    for (dim_size_t d = 0; d < n_dims; ++d)
    {
      first_.xs_[d] = 0U;
      last_.xs_[d] = dims.xs_[d] - sub_t(1U);
    }
  }

  index_t length(D const& dims) const
  {
    return last_.index(dims) - first_.index(dims) + index_t(1U);
  }

  void split(BlockRange& o, D const& dims)
  {
    index_t index_begin = first_.index(dims);
    index_t index_end = last_.index(dims);
    index_t n_elems_half = (index_end - index_begin) / index_t(2U) + 1U;

    split_n(o, dims, n_elems_half);
  }

  void split_n(BlockRange& o, D const& dims, index_t n_elems)
  {
    index_t index_begin = first_.index(dims);
    index_t mid = index_begin + n_elems - index_t(1U);

    o.first_ = first_;
    o.last_ = D::from_index(mid, dims);
    first_ = o.last_;

      // avoid an overlap of one.
    take_first(dims);
  }  

  D take_first(D const& dims)
  {
    D result = first_;
    bool carry = true;
    
    for (dim_size_t d = 0; carry; ++d)
    {
      ++first_.xs_[d];
      carry = d < n_dims - 1 && first_.xs_[d] == dims.xs_[d];
      if (carry)
      {
        first_.xs_[d] = 0U;
      }
    }

    return result;
  }

private:
  D first_;
  D last_;
};

NAMESPACE_PMT_END