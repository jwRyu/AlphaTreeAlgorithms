#pragma once

#include "../pmt_common.h"

NAMESPACE_PMT

template <typename UValue, typename Data>
struct SortPair
{
  SortPair() {}
  SortPair(UValue val, Data data) : val_(val), data_(data) {}

  INLINE UValue value() const
  {
    return val_;
  }

  INLINE Data data() const
  {
    return data_;
  }

  UValue val_;
  Data data_;
} PACKED;

template <typename UValue>
struct SortValue
{
  INLINE UValue value() const
  {
    return val_;
  }

  //INLINE UValue data() const
  //{
    //return data_;
  //}

  union
  {
    UValue val_;
    //UValue data_;
  };
};

NAMESPACE_PMT_END
