#pragma once
#include "HeapQueue.h"
#include "HierarQueue.h"
#include "Trie.h"
#include "allocator.h"
#include "defines.h"
#include "walltime.h" //tmp

#include <cfloat>
#include <cmath>
#include <iostream>
#include <vector>

#define LISTSIZE_DEFAULT 12
#define HEAPSIZE_DEFAULT 128

#define PROFILE 0

template <class Pixel> class HierarHeapQueue_HEQ {
    HQentry<Pixel> *list;
    HeapQueue_naive_quad<Pixel> **hqueue;
    HQentry<Pixel> **storage;
    ImgIdx *storage_cursize;
    ImgIdx *qsizes;
    _uint32 *histeqmap;

    ImgIdx thr_hqueue, curthr, numlevels;
    double a;
    ImgIdx queue_minlev;

    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;
    int emptytop;

  public:
    void initHQ(ImgIdx *dhist, _uint32 *histeqmap_in, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize);
    HierarHeapQueue_HEQ(ImgIdx *dhist, _uint32 *histeqmap_in, ImgIdx numlevels_in, double a_in, ImgIdx size);
    HierarHeapQueue_HEQ(ImgIdx *dhist, _uint32 *histeqmap_in, ImgIdx numlevels_in, ImgIdx size, double a_in,
                        int listsize);
    ~HierarHeapQueue_HEQ();

    void push_1stitem(ImgIdx idx, Pixel alpha);
    void end_pushes(_uint8 *isVisited);

    inline void start_pushes() { emptytop = 1; }
    inline Pixel get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }

    void push(ImgIdx idx, Pixel alpha);
    void push_queue(ImgIdx idx, Pixel alpha);
    ImgIdx pop(_uint8 *isVisited);
    int check_queue_level(_uint8 *isVisited);
    void pop_queue(_uint8 *isVisited);
};

template <class Pixel> class HierarHeapQueue {
    HeapQueue_naive_quad<Pixel> **hqueue;
    HQentry<Pixel> **storage;
    ImgIdx *storage_cursize;
    ImgIdx *qsizes;
    ImgIdx thr_hqueue, curthr, numlevels;
    double a;
    ImgIdx queue_minlev;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;

  public:
    HierarHeapQueue(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity,
                    double r);
    ~HierarHeapQueue();
    inline ImgIdx top() { return hqueue[queue_minlev]->top(); }
    inline Pixel top_alpha() { return hqueue[queue_minlev]->top_alpha(); }
    inline void push_1stitem(ImgIdx idx, Pixel alpha) { push_queue(idx, alpha); }
    inline void end_pushes(_uint8 *isVisited) { pop(isVisited); }
    inline void push(ImgIdx idx, Pixel alpha) { push_queue(idx, alpha); }
    void push_queue(ImgIdx idx, Pixel alpha);

    ImgIdx pop(_uint8 *isVisited);
    int check_queue_level(_uint8 *isVisited);
    void pop_queue(_uint8 *isVisited);
};

// hhpq
template <class Pixel> class HierarHeapQueue_cache {
    HQentry<Pixel> *list;
    HeapQueue_naive_quad<Pixel> **hqueue;
    HQentry<Pixel> **storage;
    ImgIdx *storage_cursize;
    ImgIdx *qsizes;

    ImgIdx thr_hqueue, curthr, numlevels;
    double a;
    ImgIdx queue_minlev;

    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;
    int emptytop;

    ImgIdx maxSize;

  public:
    _uint8 *edgeLabels;

#if PROFILE
    double t0 = get_cpu_time();
    double tconv = 0;
    double tcache = 0;
    double tqueue = 0;

    ImgIdx num_cache = 0;
    ImgIdx num_cache_ovfl = 0;
    ImgIdx num_hq = 0;
    ImgIdx num_store = 0;
    ImgIdx num_conv = 0;

    std::vector<_uint64> num_memmove_push;
    std::vector<_uint64> num_memmove_pop;
    std::vector<_uint64> num_items_push;
    std::vector<_uint64> num_items_pop;
    _uint64 curSize = 0;

    _uint64 num_memmove_push_i;
    _uint64 num_memmove_pop_i;

    void decrement_curSize() { curSize--; }

#endif

    void initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity, double r);
    HierarHeapQueue_cache(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize,
                          ImgIdx connectivity = 4, double r = 0.2);
    ~HierarHeapQueue_cache();
    inline void start_pushes() { emptytop = 1; }
    inline Pixel get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }
    inline Pixel top_edge() { return list[0].edge; }
    void push_1stitem(ImgIdx idx, Pixel alpha, ImgIdx edgeIdx);
    void end_pushes(_uint8 *isVisited);
    void push(ImgIdx idx, Pixel alpha, ImgIdx edgeIdx);
    void push_queue(ImgIdx idx, Pixel alpha, ImgIdx edgeIdx);
    ImgIdx pop(_uint8 *isVisited);
    int check_queue_level(_uint8 *isVisited);
    void pop_queue(_uint8 *isVisited);
};

template <class Pixel> class Cache_Heapqueue {
    HQentry<Pixel> *list;
    HeapQueue_naive<Pixel> *hqueue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;
    int emptytop;

    void initHQ(ImgIdx size, size_t listsize);

  public:
    double qtime; // tmp

    Cache_Heapqueue(ImgIdx size);
    Cache_Heapqueue(ImgIdx size, size_t listsize);
    ~Cache_Heapqueue();

    inline void start_pushes() { emptytop = 1; }
    inline ImgIdx get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }
    inline void end_pushes() {
        if (emptytop)
            pop();
    }
    inline void push_queue(ImgIdx idx, Pixel alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() { hqueue->pop(); }
    void push_1stitem(ImgIdx idx, Pixel alpha);
    void push(ImgIdx idx, Pixel alpha);
    ImgIdx pop();
};

template <class Pixel> class Cache_Quad_Heapqueue {
    // MinList1 *list, *list_end, *head, *tail;
    HQentry<Pixel> *list;
    HeapQueue_naive_quad<Pixel> *hqueue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;
    int emptytop;

    void initHQ(ImgIdx size, size_t listsize);

  public:
    double qtime; // tmp

    Cache_Quad_Heapqueue(ImgIdx size);
    Cache_Quad_Heapqueue(ImgIdx size, size_t listsize);
    ~Cache_Quad_Heapqueue();

    inline void start_pushes() { emptytop = 1; }
    inline Pixel get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }
    inline void end_pushes() {
        if (emptytop)
            pop();
    }
    inline void push_queue(ImgIdx idx, Pixel alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() { hqueue->pop(); }

    void push_1stitem(ImgIdx idx, Pixel alpha);
    void push(ImgIdx idx, Pixel alpha);
    ImgIdx pop();
};

template <class Pixel> class CirCache_Hierqueue {
    HQentry<_int32> *list;
    HierarQueue *hqueue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list, liststart, mask;
    ImgIdx maxSize_queue, mask_field;
    int emptytop;
    void initHQ(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);

  public:
    double qtime; // tmp

    CirCache_Hierqueue(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels);
    CirCache_Hierqueue(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);
    ~CirCache_Hierqueue();

    inline void start_pushes() { emptytop = 1; }
    inline _int32 get_minlev() { return list[liststart].alpha; }
    inline ImgIdx top() { return list[liststart].pidx; }
    inline _int32 top_alpha() { return list[liststart].alpha; }
    inline void end_pushes() {
        if (emptytop)
            pop();
    }
    inline void push_queue(ImgIdx idx, _int32 alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() {
        hqueue->pop();
        hqueue->find_minlev();
    }

    void push_1stitem(ImgIdx idx, _int32 alpha);
    void push(ImgIdx idx, _int32 alpha);
    ImgIdx pop();
};

template <class Pixel> class HierarQueueCache {
    HQentry<_int32> *list;
    HierarQueue *hqueue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;
    int emptytop;
    void initHQ(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize = LISTSIZE_DEFAULT);

  public:
    double qtime; // tmp

    HierarQueueCache(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels);
    HierarQueueCache(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);
    ~HierarQueueCache();

    inline void start_pushes() { emptytop = 1; }
    inline _int32 get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline _int32 top_alpha() { return list[0].alpha; }
    inline void end_pushes() {
        if (emptytop)
            pop();
    }
    inline void push_queue(ImgIdx idx, _int32 alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() {
        hqueue->pop();
        hqueue->find_minlev();
    }

    void push_1stitem(ImgIdx idx, _int32 alpha);
    void push(ImgIdx idx, _int32 alpha);
    ImgIdx pop();
};

template <class Pixel> class Cache_Hierqueue_l1 {
    HQentry<_int32> *list;
    HQueue_l1idx *hqueue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;
    int emptytop;

    void initHQ(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);

  public:
    double qtime; // tmp

    Cache_Hierqueue_l1(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels);
    Cache_Hierqueue_l1(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);
    ~Cache_Hierqueue_l1();

    inline void start_pushes() { emptytop = 1; }
    inline _int32 get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline _int32 top_alpha() { return list[0].alpha; }
    inline void end_pushes() {
        if (emptytop)
            pop();
    }
    inline void push_queue(ImgIdx idx, _int32 alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() {
        hqueue->pop();
        hqueue->find_minlev();
    }

    void push_1stitem(ImgIdx idx, _int32 alpha);
    void push(ImgIdx idx, _int32 alpha);
    ImgIdx pop();
};

template <class Pixel> class Cache_Hierqueue_l2 {
    HQentry<_int32> *list;
    HQueue_l2idx *hqueue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;
    int emptytop;

    void initHQ(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);

  public:
    double qtime; // tmp

    Cache_Hierqueue_l2(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels);
    Cache_Hierqueue_l2(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);
    ~Cache_Hierqueue_l2();

    inline void start_pushes() { emptytop = 1; }
    inline _int32 get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline _int32 top_alpha() { return list[0].alpha; }
    inline void end_pushes() {
        if (emptytop)
            pop();
    }
    inline void push_queue(ImgIdx idx, _int32 alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() {
        hqueue->pop();
        hqueue->find_minlev();
    }

    void push_1stitem(ImgIdx idx, _int32 alpha);
    void push(ImgIdx idx, _int32 alpha);
    ImgIdx pop();
};

class Trie_Cache {
    ImgIdx *list;
    Trie<_uint64> *trie;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;

    void initHQ(ImgIdx size, size_t listsize);

  public:
    Trie_Cache(ImgIdx size);
    Trie_Cache(ImgIdx size, size_t listsize);
    ~Trie_Cache();

    inline ImgIdx get_minlev() { return list[0]; }
    inline ImgIdx top() { return list[0]; }
    inline void push_queue(ImgIdx idx) { trie->push(idx); }
    inline void pop_queue() { trie->pop(); }

    void push(ImgIdx idx);
    void pop();
};

class HybridQueue_HQueue_Rank {
    ImgIdx *list;
    HQueue_l1idx_rank *queue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;

    void initHQ(ImgIdx size, size_t listsize);

  public:
    HybridQueue_HQueue_Rank(ImgIdx size);
    HybridQueue_HQueue_Rank(ImgIdx size, size_t listsize);
    ~HybridQueue_HQueue_Rank();

    inline ImgIdx get_minlev() { return list[0]; }
    inline ImgIdx top() { return list[0]; }
    inline void push_queue(ImgIdx idx) { queue->push(idx); }
    void push(ImgIdx idx);
    void pop();
    inline void pop_queue() { queue->pop(); }
};

class HybridQueue_HQueue_Rank1 {
    ImgIdx *list;
    HQueue_l1idx_rank *queue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 shamt, nbit;
    _int16 l0, mask;

    void initHQ(ImgIdx size, size_t listsize);

  public:
    HybridQueue_HQueue_Rank1(ImgIdx size);
    HybridQueue_HQueue_Rank1(ImgIdx size, size_t listsize);
    ~HybridQueue_HQueue_Rank1();

    inline ImgIdx get_minlev() { return list[l0]; }
    inline ImgIdx top() { return list[l0]; }
    inline void push_queue(ImgIdx idx) { queue->push(idx); }
    inline void pop_queue() { queue->pop(); }

    void push(ImgIdx idx);
    void pop();
};

class HybridQueue_HQueue {
    ImgIdx *list;
    _int64 *levels;
    HQueue_l1idx *queue;
    ImgIdx minidx_queue;
    _int16 curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    _int8 minlevnotfixed;

    void initHQ(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);

  public:
    HybridQueue_HQueue(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels);
    HybridQueue_HQueue(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels, size_t listsize);
    ~HybridQueue_HQueue();

    inline ImgIdx top() { return list[0]; }
    inline _int64 get_minlev() { return levels[0]; }
    inline void find_minlev() { queue->find_minlev(); }
    inline void push_queue(ImgIdx idx, _int64 level) {
        if (queue->push(idx, level))
            minlevnotfixed = 0;
    }
    inline void pop_queue() {
        minlevnotfixed = 1;
        queue->pop();
    }

    void push(ImgIdx idx, _int64 level);
    ImgIdx pop();
};
