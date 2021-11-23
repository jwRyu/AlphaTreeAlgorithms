#include "../pmt_common.h"
#include <algorithm>

NAMESPACE_PMT

template <typename Value, typename Index>
void validate_sort(std::string const& info, Value* values, Index* rank_to_index, Index n)
{
  Index* counts = new Index[n];
  std::fill(counts, counts + n, Index(0));

  for (Index i = 0; i < n; ++i)
  {
    ++counts[rank_to_index[i]];
  }

  for (Index i = 0; i < n; ++i)
  {
    if (counts[i] != 1U)
    {
      std::cout << "sortmap is not a permutation @ "<< info << '\n';
      assert(false);
    }

    assert(counts[i] == 1U);
  }

  delete[] counts;

  // check order

  bool fail = false;

  for (Index i = 1; i < n; ++i)
  {
    bool chk = 
      values[rank_to_index[i]] >= values[rank_to_index[i - 1]];

    if (!chk)
    {
      fail = true;
      break;
    }
  }

  if (fail)
  {
    std::cout << "invalid order @ "<< info << '\n';
    assert(false);
  }
}

NAMESPACE_PMT_END