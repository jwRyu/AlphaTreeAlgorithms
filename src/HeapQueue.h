#pragma once
#include "defines.h"
#include <cfloat>
#define PROFILE 0
#define HEAPQUEUE_DEBUG 1

#if PROFILE
#include <vector>
#endif

template <class Pixel> class HQentry {
  public:
    HQentry() {}
    ~HQentry() {}
    ImgIdx pidx;
    Pixel alpha;

#if HEAPQUEUE_DEBUG
    ImgIdx edge = 0;
#endif

    inline void operator=(const HQentry &item) {
        this->pidx = item.pidx;
        this->alpha = item.alpha;
    }
};

// Heap-based priority queue
// Use push buffer to minimize push() and pop() call and time
template <class Pixel> //
class HeapQueue {
    ImgIdx maxsize;
    HQentry<Pixel> *arr;
    Pixel pop_level;
    Pixel max_level;
    HQentry<Pixel> pushlist[7]; // Connectivity - 1
    _uint8 flags[7];
    _int8 pushlistidx;

#if PROFILE
    std::vector<_uint64> num_memmove_push;
    std::vector<_uint64> num_memmove_pop;
    std::vector<_uint64> num_items_push;
    std::vector<_uint64> num_items_pop;

    _uint64 num_memmove_push_i;
    _uint64 num_memmove_pop_i;
#endif

  public:
    ImgIdx cursize;
    double qtime, timing;

    HeapQueue(ImgIdx maxsize_in);
    ~HeapQueue();
    int minidx_pushlist();
    void find_minlev();
    ImgIdx pop();
    void push(ImgIdx pidx, Pixel alpha);
    void push_run(ImgIdx pidx, Pixel alpha);
    inline Pixel get_minlev() { return arr[1].alpha; }
    inline ImgIdx top() { return arr[1].pidx; }
};

// Heap-based priority queue
template <class Pixel> class HeapQueue_naive {
    ImgIdx cursize;
    ImgIdx maxsize;
    HQentry<Pixel> *arr;
    Pixel pop_level;

#if PROFILE
    std::vector<_uint64> num_memmove_push;
    std::vector<_uint64> num_memmove_pop;
    std::vector<_uint64> num_items_push;
    std::vector<_uint64> num_items_pop;

    _uint64 num_memmove_push_i;
    _uint64 num_memmove_pop_i;
#endif

  public:
    double qtime; // tmp
    HeapQueue_naive(ImgIdx maxsize_in);
    ~HeapQueue_naive();
    ImgIdx pop();
    void push(ImgIdx pidx, Pixel alpha);
    inline Pixel get_minlev() { return arr[1].alpha; }
    inline ImgIdx top() { return arr[1].pidx; }
    inline Pixel top_alpha() { return arr[1].alpha; }
};

// Heap-based priority queue
template <class Pixel> class HeapQueue_naive_quad {
    ImgIdx cursize;
    ImgIdx maxsize;
    HQentry<Pixel> *arr;
    Pixel pop_level;

  public:
    double qtime; // tmp
    inline ImgIdx get_cursize() { return cursize; }
    HeapQueue_naive_quad(ImgIdx maxsize_in);
    ~HeapQueue_naive_quad();

#if PROFILE
    _uint64 pop();
    _uint64 push(ImgIdx pidx, Pixel alpha);
#else
    ImgIdx pop();
    void push(ImgIdx pidx, Pixel alpha);
#endif

    inline Pixel get_minlev() { return arr[1].alpha; }
    inline ImgIdx top() { return arr[1].pidx; }
    inline Pixel top_alpha() { return arr[1].alpha; }
};
