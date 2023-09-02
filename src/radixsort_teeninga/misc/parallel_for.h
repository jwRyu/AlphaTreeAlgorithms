#pragma once

#include <pthread.h>
#include "../misc/random.h"
#include "../misc/block_range.h"
#include <algorithm>

NAMESPACE_PMT

template <typename D>
struct ParallelSlot
{
  ParallelSlot()
  {
    PMT_ASSERT(!pthread_mutex_init(&lock_, nullptr));
  }

  ~ParallelSlot()
  {
    PMT_ASSERT(!pthread_mutex_destroy(&lock_));
  }

  pthread_mutex_t lock_;
  BlockRange<D> range_;
};

template <typename F, typename D>
struct SharedThreadData
{
  SharedThreadData(
    F const& f,
    ParallelSlot<D>* slots,
    uint16_t n_slots,
    D const& dims) :
    f_(f), slots_(slots), n_slots_(n_slots), dims_(dims)
  {
  }

  F const& f_;
  ParallelSlot<D>* const slots_;
  uint16_t const n_slots_;
  D const dims_;
};

template <typename F, typename D>
struct ThreadData
{
  ThreadData()
  {

  }

  ThreadData(SharedThreadData<F, D>* shared, uint16_t thread_nr) :
    shared_(shared), thread_nr_(thread_nr)
  {
  }

  SharedThreadData<F, D>* shared_;
  typename rng<unsigned>::type rng_;
  uint16_t thread_nr_;
};

template <typename F, typename R>
unsigned next_nonempty_slot(ThreadData<F, R>& data, unsigned slot_nr)
{
  auto const shared = data.shared_;
  unsigned const n_slots = shared->n_slots_;
  unsigned x = random_uint(data.rng_, n_slots - 1U);
  //unsigned x = slot_nr + 1U;

  bool available_work = true;
  while (available_work)
  {
    available_work = false;

    for (unsigned ctr = 0U; ctr < n_slots; ++ctr, ++x)
    {
      if (x == n_slots) x = 0U;
      if (x == slot_nr) continue;

      auto& slot = shared->slots_[x];

      if (!pthread_mutex_trylock(&slot.lock_))
      {
        if (slot.range_.length(shared->dims_) != 0U)
        {
          return x;
        }

        pthread_mutex_unlock(&slot.lock_);
        continue;
      }

      if (slot.range_.length(shared->dims_) != 0U)
      {
        available_work = true;
      }
    }
  }

  return slot_nr;
}

template <typename F, typename D>
bool select_block(ThreadData<F, D>& data, D& block, unsigned thread_nr)
{
  auto const& shared = data.shared_;
  //unsigned slot_nr = random_uint(data.rng_, shared->n_slots_ - 1U);
  unsigned slot_nr = thread_nr;
  auto& slot = shared->slots_[slot_nr];

  pthread_mutex_lock(&slot.lock_);

  if (slot.range_.length(shared->dims_) == 0)
  {
    unsigned next_slot_nr = next_nonempty_slot(data, slot_nr);

    if (next_slot_nr == slot_nr)
    {
      pthread_mutex_unlock(&slot.lock_);
      return false;
    }

    shared->slots_[next_slot_nr].range_.split(shared->slots_[slot_nr].range_, shared->dims_);
    pthread_mutex_unlock(&shared->slots_[next_slot_nr].lock_);
  }

  block = slot.range_.take_first(shared->dims_);
  pthread_mutex_unlock(&slot.lock_);
  return true;
}

template <typename F, typename D>
void* init_f(void* data_ptr)
{
  ThreadData<F, D>& data = *reinterpret_cast<ThreadData<F, D>*>(data_ptr);
  D block;

  while (select_block(data, block, data.thread_nr_))
  {
    data.shared_->f_(block, data.thread_nr_);
  }

  return NULL;
}

template <typename F, typename R, typename D>
void iterate_block_range(F const& f, R range, D dims, unsigned thread_nr)
{
  while (range.length(dims))
  {
    f(range.take_first(dims), thread_nr);
  }
}

template <typename D>
void distribute_range_to_slots(BlockRange<D> range, D const&dims, ParallelSlot<D>* slots, signed n_slots)
{
  using index_t = typename D::index_t;

  index_t n_blocks = range.length(dims);
  index_t blocks_per_slot = n_blocks / n_slots;
  index_t blocks_remainder = n_blocks % n_slots;

  for (signed s = 0; s < n_slots; ++s)
  {
    index_t n_split = blocks_per_slot + (blocks_remainder > s);
    range.split_n(slots[s].range_, dims, n_split);
  }
}

template <typename F, typename D>
void parallel_for_blocks(
  D dims,
  F const& f,
  unsigned n_threads = 0)
{
  using index_t = typename D::index_t;

  index_t n = dims.length();

  if (n == 0)
  {
    return;
  }

  if (n_threads == 0)
  {
    n_threads = pmt::n_threads;
  }

  n_threads = std::min(index_t(n_threads), n);

  using Range = BlockRange<D>;
  Range range(dims);

  if (n_threads == 1)
  {
    iterate_block_range(f, range, dims, 0U);
    return;
  }

  unsigned n_slots = n_threads;

  pthread_t threads[n_threads];
  ParallelSlot<D> slots[n_slots];
  SharedThreadData<F, D> shared(f, slots, n_slots, dims);
  ThreadData<F, D> t_data[n_threads];

  distribute_range_to_slots(range, dims, slots, n_slots);

  for (unsigned t = 0; t < n_threads; ++t)
  {
    t_data[t] = {&shared, static_cast<uint16_t>(t)};
    assert(!pthread_create(&threads[t], nullptr, init_f<F, D>, (void*)&t_data[t]));
  }

  for (unsigned t = 0; t < n_threads; ++t)
  {
    assert(!pthread_join(threads[t], NULL));
  }
}

template <typename F, typename index_t>
NO_INLINE void parallel_for_iterate(F const& f, index_t begin, index_t end, unsigned thread_nr)
{
  for (; begin < end; ++begin)
  {
    f(begin, thread_nr);
  }
}

template <typename index_t, typename F>
void parallel_for(
  index_t n,
  F const& f,
  unsigned block_size = parallel_block_size,
  unsigned n_threads = 0)
{
  using D = Dim<index_t>;

  index_t n_blocks = (n + block_size - 1U) / block_size;
  D const dims({n_blocks});

  parallel_for_blocks(dims, [=](D block, unsigned thread_nr)
  {
    index_t begin = block.index(dims) * block_size;
    index_t end = begin + block_size;

    if (block.index(dims) == n_blocks - 1)
    {
      end = n;
    }

    parallel_for_iterate(f, begin, end, thread_nr);
  }, n_threads);
}

NAMESPACE_PMT_END
