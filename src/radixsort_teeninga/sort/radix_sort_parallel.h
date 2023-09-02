#pragma once

#include "../pmt_common.h"
#include "../misc/exclusive_sum.h"
#include "sort_item.h"
#include "../misc/parallel_for.h"
#include "../misc/range.h"
#include "../misc/timer.h"
#include "../misc/alloc.h"
#include "../misc/copy.h"

NAMESPACE_PMT

template<
  unsigned histo_sz_log2,
  typename DataRAI,
  typename ItemRAI,
  typename Index,
  typename InitialItemF,
  typename LastItemF>
class RadixSortParallel
{
public:
  static constexpr signed histo_sz = 1U << histo_sz_log2;
  static constexpr signed histo_mask = histo_sz - 1U;

  using Histogram = Index[histo_sz];
  using Data = typename std::iterator_traits<DataRAI>::value_type;
  using Item = typename std::iterator_traits<ItemRAI>::value_type;
  using UValue = typename std::remove_const<decltype(Item().value())>::type;

  RadixSortParallel(
    DataRAI const data,
    ItemRAI const items1,
    ItemRAI const items2,
    Index const len,
    uint_fast8_t bit_start,
    uint_fast8_t bit_end,
    InitialItemF const& f_initial_item,
    LastItemF const& f_last_item,
    unsigned const n_threads) :
    data_(data),
    items1_(items1),
    items2_(items2),
    len_(len),
    f_initial_item_(f_initial_item),
    f_last_item_(f_last_item),
    n_threads_(n_threads),
    hs_(new Histogram[2 * n_blocks()]),
    global_offsets_(hs_ + n_blocks()),
    buffers_(new Item[n_threads * block_sz()]),
    bits_{bit_start, bit_end}
  {
    sort_digits();
  }

  ~RadixSortParallel()
  {
    delete[] buffers_;
    delete[] hs_;
  }

private:
  template <typename ItemF>
  void create_histograms(uint8_t shift, ItemF const& f_item)
  {
    parallel_for(n_blocks(), [=](Index b, unsigned thread_nr)
    {
      Range<Index> const range = make_range(b);
      Histogram& h = hs_[b];

      std::fill(h, h + histo_sz, 0);

      for (Index i = range.begin_; i < range.end_; ++i)
      {
        UValue u = f_item(i).val_ >> shift;

        ++h[u & histo_mask];
      }
    }, Index(1), n_threads_);
  }

  void make_offsets()
  {
    Index sums[histo_sz + 1];

    std::fill(sums, sums + histo_sz + 1, 0);

    for (Index b = 0; b < n_blocks(); ++b)
    {
      Histogram& h = hs_[b];
      Histogram& g = global_offsets_[b];

      for (Index i = 0; i < histo_sz; ++i)
      {
        g[i] = sums[i];
        sums[i] += h[i];
      }
    }

    exclusive_sum(sums, sums + histo_sz + 1);

    assert(sums[histo_sz] == len_);

    for (Index b = 0; b < n_blocks(); ++b)
    {
      Histogram& g = global_offsets_[b];

      for (Index i = 0; i < histo_sz; ++i)
      {
        g[i] += sums[i];
      }
    }
  }

  template <typename ItemOutRAI, typename ItemF, typename OutF>
  void scatter_digit(
    ItemOutRAI const items_out,
    uint8_t shift,
    ItemF const& f_item,
    OutF const& f_out)
  {
    using ItemOut = typename std::iterator_traits<ItemOutRAI>::value_type;

    parallel_for(n_blocks(), [=](Index b, unsigned thread_nr)
    {
      Range<Index> const range = make_range(b);
      Histogram& h = hs_[b];
      Histogram buffer_offsets;
      ItemOut* buffer = (ItemOut*)buffers_ + thread_nr * block_sz();

      copy_rai(h, histo_sz, buffer_offsets);
      //std::copy(h, h + histo_sz, buffer_offsets);
      exclusive_sum(buffer_offsets, buffer_offsets + histo_sz);

      for (Index i = range.begin_; i < range.end_; ++i)
      {
        Item const& item = f_item(i);

        UValue u = item.val_ >> shift;

        Index out = buffer_offsets[u & histo_mask]++;

        f_out(buffer[out], item);
      }

      Histogram& g = global_offsets_[b];

      ItemOut* buffer_ptr = buffer;
      for (Index k = 0; k < histo_sz; ++k)
      {
        copy_rai(buffer_ptr, h[k], items_out + g[k]);
        //std::copy(buffer_ptr, buffer_ptr + h[k], items_out + g[k]);
        buffer_ptr += h[k];
      }
    }, Index(1), n_threads_);
  }

  void sort_digits()
  {
    uint_fast8_t n_digits = (bits_[1] - bits_[0] + histo_sz_log2 - 1U) / histo_sz_log2;

    if (n_digits == 1)
    {
      sort_digit(data_, 0, f_initial_item_, f_last_item_);
      return;
    }

    auto const& f_item =
      [=](Index i) ALWAYS_INL_L(Item const&)
      {
        return items1_[i];
      };

    auto const& f_out =
      [=](Item& out, Item const& item) ALWAYS_INLINE
      {
        out = item;
      };

    sort_digit(items2_, bits_[0], f_initial_item_, f_out);

    std::swap(items1_, items2_);

    for (uint_fast8_t d = 1; d < n_digits - 1U; ++d)
    {
      sort_digit(items2_, bits_[0] + d * histo_sz_log2, f_item, f_out);

      std::swap(items1_, items2_);
    }

    sort_digit(data_, bits_[0] + (n_digits - 1U) * histo_sz_log2, f_item, f_last_item_);
  }

  template <typename ItemOutRAI, typename ItemF, typename OutF>
  void sort_digit(
    ItemOutRAI items_out,
    uint8_t shift,
    ItemF const& f_item,
    OutF const& f_out)
  {
    create_histograms(shift, f_item);
    make_offsets();
    scatter_digit(items_out, shift, f_item, f_out);
  }

  constexpr Range<Index> const make_range(Index block_nr) const
  {
    Index const begin = block_nr * block_sz();
    Index const end = begin + std::min(len_ - begin, block_sz());

    return {begin, end};
  }

  constexpr Index block_sz() const
  {
    return block_size_sorting / sizeof(Item);
  }

  constexpr Index n_blocks() const
  {

    return (len_ + block_sz() - 1) / block_sz();
  }

  DataRAI const data_;
  ItemRAI items1_;
  ItemRAI items2_;
  Index const len_;
  InitialItemF const& f_initial_item_;
  LastItemF const& f_last_item_;
  unsigned const n_threads_;
  Histogram* const hs_;
  Histogram* const global_offsets_;
  Item* buffers_;
  uint8_t bits_[2];
};

template<
  unsigned histo_sz_log2 = pmt::histo_sz_log2,
  typename DataRAI,
  typename ItemRAI,
  typename Index,
  typename InitialItemF,
  typename LastItemF>
void radix_sort_parallel(
  DataRAI const data,
  ItemRAI const items1,
  ItemRAI const items2,
  Index const len,
  uint_fast8_t bit_start,
  uint_fast8_t bit_end,
  InitialItemF const& f_initial_item,
  LastItemF const& f_last_item,
  unsigned const n_threads)
{
  if (bit_end <= bit_start) return;

  RadixSortParallel<histo_sz_log2, DataRAI, ItemRAI, Index, InitialItemF, LastItemF>(
    data, items1, items2, len, bit_start, bit_end, f_initial_item, f_last_item, n_threads);
}

template<
  unsigned histo_sz_log2 = pmt::histo_sz_log2,
  typename DataRAI,
  typename ItemRAI,
  typename Index,
  typename InitialItemF,
  typename LastItemF>
void radix_sort_parallel(
  DataRAI const data,
  ItemRAI const items1,
  ItemRAI const items2,
  Index const len,
  InitialItemF const& f_initial_item,
  LastItemF const& f_last_item,
  unsigned const n_threads)
{
  using UValue = typename std::remove_const<decltype(f_initial_item(0U).value())>::type;

  radix_sort_parallel(data, items1, items2, len, 0U, sizeof(UValue) * CHAR_BIT, f_initial_item, f_last_item, n_threads);
}

template<
  unsigned histo_sz_log2 = pmt::histo_sz_log2,
  typename ItemRAI,
  typename Index,
  typename IndexRAI,
  typename SortItemRAI>
void rank_to_index(
  ItemRAI const items,
  Index const len,
  IndexRAI const out,
  uint_fast8_t bit_start,
  uint_fast8_t bit_end,
  SortItemRAI sort_space,
  unsigned n_threads = 0U)
{
  using UValue = typename std::remove_const<decltype(items[0].value())>::type;
  using SortPair = SortPair<UValue, Index>;

  if (n_threads == 0)
  {
    n_threads = pmt::n_threads;
  }

  auto const& f_initial_item = [=](Index i) ALWAYS_INL_L(SortPair)
  {
    return {items[i].value(), i};
  };

  auto const& f_out = [=](Index& out, SortPair const& item) ALWAYS_INLINE
  {
    out = item.data_;
  };

  radix_sort_parallel<histo_sz_log2>(out, sort_space + len, sort_space, len, bit_start, bit_end, f_initial_item, f_out, n_threads);
}

NAMESPACE_PMT_END
