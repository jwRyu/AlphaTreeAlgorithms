extern "C" {
#include <sys/mman.h>
}

#include "../pmt_common.h"
#include <iostream>

template <typename Elem>
struct AllocatedMemory
{
  void free()
  {
    if (ptr_ == nullptr) return;

    if (munmap(ptr_, length_) != 0)
    {
      std::cerr << "Error: unable to free memory\n";
      assert(false);
    }
  }

  Elem* ptr_ = nullptr;
  size_t length_ = 0U;
};

size_t aligned_length(size_t length, unsigned n_bits)
{
  size_t mask = (size_t(1) << n_bits) - size_t(1);
  size_t aligned = length & ~size_t(mask);

  if (length & mask)
  {
    aligned += 1U << n_bits;    
  }

  return aligned;
}

template <typename Elem>
AllocatedMemory<Elem> allocate(size_t length)
{
  AllocatedMemory<Elem> am;

  am.length_ = aligned_length(length * sizeof(Elem), 21U);
  //printf("%lu %lu\n", length, am.length_);
  am.ptr_ = (Elem*)mmap(nullptr, am.length_, PROT_READ|PROT_WRITE, MAP_ANON|MAP_PRIVATE|MAP_HUGETLB, 0, 0);

  if (am.ptr_ != MAP_FAILED)
  {
    return am;
  }

  std::cerr << "Warning: unable to allocate 2MB pages\n";

  am.length_ = aligned_length(length * sizeof(Elem), 12U);
  am.ptr_ = (Elem*)mmap(nullptr, am.length_, PROT_READ|PROT_WRITE, MAP_ANON|MAP_PRIVATE, 0, 0);

  if (am.ptr_ == MAP_FAILED)
  {
    std::cerr << "Error: unable to allocate memory\n";
    assert(false);
  }

  return am;
}


