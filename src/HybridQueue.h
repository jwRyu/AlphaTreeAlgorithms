#pragma once
#include "HeapQueue.h"
#include "HierarQueue.h"
#include "Trie.h"
#include "allocator.h"
#include "defines.h"
#include "walltime.h" //tmp

#define LISTSIZE_DEFAULT 12
#define HEAPSIZE_DEFAULT 128

#define PROFILE 0

template <class Pixel> class HierarHeapQueue_HEQ {
    HQentry<Pixel> *list;
    HeapQueue_naive_quad<Pixel> **hqueue;
    HQentry<Pixel> **storage;
    ImgIdx *storage_cursize;
    ImgIdx *qsizes;
    uint32_t *histeqmap;

    ImgIdx thr_hqueue, curthr, numlevels;
    double a;
    ImgIdx queue_minlev;

    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;
    int _emptyTop;

  public:
    void initHQ(ImgIdx *dhist, uint32_t *histeqmap_in, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize);
    HierarHeapQueue_HEQ(ImgIdx *dhist, uint32_t *histeqmap_in, ImgIdx numlevels_in, double a_in, ImgIdx size);
    HierarHeapQueue_HEQ(ImgIdx *dhist, uint32_t *histeqmap_in, ImgIdx numlevels_in, ImgIdx size, double a_in,
                        int listsize);
    ~HierarHeapQueue_HEQ();

    void push_1stitem(ImgIdx idx, Pixel alpha);
    void endPushes(uint8_t *isVisited);

    inline void startPushes() { _emptyTop = 1; }
    inline Pixel get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }

    void push(ImgIdx idx, Pixel alpha);
    void push_queue(ImgIdx idx, Pixel alpha);
    ImgIdx pop(uint8_t *isVisited);
    int check_queue_level(uint8_t *isVisited);
    void pop_queue(uint8_t *isVisited);
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
    int8_t shamt, nbit;

  public:
    HierarHeapQueue(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity,
                    double r);
    ~HierarHeapQueue();
    inline ImgIdx top() { return hqueue[queue_minlev]->top(); }
    inline Pixel top_alpha() { return hqueue[queue_minlev]->top_alpha(); }
    inline void push_1stitem(ImgIdx idx, Pixel alpha) { push_queue(idx, alpha); }
    inline void endPushes(uint8_t *isVisited) { pop(isVisited); }
    inline void push(ImgIdx idx, Pixel alpha) { push_queue(idx, alpha); }
    void push_queue(ImgIdx idx, Pixel alpha);

    ImgIdx pop(uint8_t *isVisited);
    int check_queue_level(uint8_t *isVisited);
    void pop_queue(uint8_t *isVisited);
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

    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;
    int _emptyTop;

    ImgIdx maxSize;

  public:
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

    std::vector<uint64_t> num_memmove_push;
    std::vector<uint64_t> num_memmove_pop;
    std::vector<uint64_t> num_items_push;
    std::vector<uint64_t> num_items_pop;
    uint64_t curSize = 0;

    uint64_t num_memmove_push_i;
    uint64_t num_memmove_pop_i;

    void decrement_curSize() { curSize--; }

#endif

    void initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity, double r);
    HierarHeapQueue_cache(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize,
                          ImgIdx connectivity = 4, double r = 0.2);
    ~HierarHeapQueue_cache();
    inline void startPushes() { _emptyTop = 1; }
    inline Pixel get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }
    void push_1stitem(ImgIdx idx, Pixel alpha);
    void endPushes(uint8_t *isVisited);
    void push(ImgIdx idx, Pixel alpha);
    void push_queue(ImgIdx idx, Pixel alpha);
    ImgIdx pop(uint8_t *isVisited);
    int check_queue_level(uint8_t *isVisited);
    void pop_queue(uint8_t *isVisited);
};

template <class Pixel> class Cache_Heapqueue {
    HQentry<Pixel> *list;
    HeapQueue_naive<Pixel> *hqueue;
    ImgIdx minidx_queue;
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;
    int _emptyTop;

    void initHQ(ImgIdx size, size_t listsize);

  public:
    double qtime; // tmp

    Cache_Heapqueue(ImgIdx size);
    Cache_Heapqueue(ImgIdx size, size_t listsize);
    ~Cache_Heapqueue();

    inline void startPushes() { _emptyTop = 1; }
    inline ImgIdx get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }
    inline void endPushes() {
        if (_emptyTop)
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
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;
    int _emptyTop;

    void initHQ(ImgIdx size, size_t listsize);

  public:
    double qtime; // tmp

    Cache_Quad_Heapqueue(ImgIdx size);
    Cache_Quad_Heapqueue(ImgIdx size, size_t listsize);
    ~Cache_Quad_Heapqueue();

    inline void startPushes() { _emptyTop = 1; }
    inline Pixel get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }
    inline void endPushes() {
        if (_emptyTop)
            pop();
    }
    inline void push_queue(ImgIdx idx, Pixel alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() { hqueue->pop(); }

    void push_1stitem(ImgIdx idx, Pixel alpha);
    void push(ImgIdx idx, Pixel alpha);
    ImgIdx pop();
};

template <class Pixel> class CirCache_Hierqueue {
    HQentry<int32_t> *list;
    HierarQueue *hqueue;
    ImgIdx minidx_queue;
    int16_t curSize_list, maxSize_list, liststart, mask;
    ImgIdx maxSize_queue, mask_field;
    int _emptyTop;
    void initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);

  public:
    double qtime; // tmp

    CirCache_Hierqueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels);
    CirCache_Hierqueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);
    ~CirCache_Hierqueue();

    inline void startPushes() { _emptyTop = 1; }
    inline int32_t get_minlev() { return list[liststart].alpha; }
    inline ImgIdx top() { return list[liststart].pidx; }
    inline int32_t top_alpha() { return list[liststart].alpha; }
    inline void endPushes() {
        if (_emptyTop)
            pop();
    }
    inline void push_queue(ImgIdx idx, int32_t alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() {
        hqueue->pop();
        hqueue->find_minlev();
    }

    void push_1stitem(ImgIdx idx, int32_t alpha);
    void push(ImgIdx idx, int32_t alpha);
    ImgIdx pop();
};

template <class Pixel> class HierarQueueCache {
    HQentry<int32_t> *list;
    HierarQueue *hqueue;
    ImgIdx minidx_queue;
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;
    int _emptyTop;
    void initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize = LISTSIZE_DEFAULT);

  public:
    double qtime; // tmp

    HierarQueueCache(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels);
    HierarQueueCache(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);
    ~HierarQueueCache();

    inline void startPushes() { _emptyTop = 1; }
    inline int32_t get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline int32_t top_alpha() { return list[0].alpha; }
    inline void endPushes() {
        if (_emptyTop)
            pop();
    }
    inline void push_queue(ImgIdx idx, int32_t alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() {
        hqueue->pop();
        hqueue->find_minlev();
    }

    void push_1stitem(ImgIdx idx, int32_t alpha);
    void push(ImgIdx idx, int32_t alpha);
    ImgIdx pop();
};

template <class Pixel> class Cache_Hierqueue_l1 {
    HQentry<int32_t> *list;
    HQueue_l1idx *hqueue;
    ImgIdx minidx_queue;
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;
    int _emptyTop;

    void initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);

  public:
    double qtime; // tmp

    Cache_Hierqueue_l1(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels);
    Cache_Hierqueue_l1(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);
    ~Cache_Hierqueue_l1();

    inline void startPushes() { _emptyTop = 1; }
    inline int32_t get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline int32_t top_alpha() { return list[0].alpha; }
    inline void endPushes() {
        if (_emptyTop)
            pop();
    }
    inline void push_queue(ImgIdx idx, int32_t alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() {
        hqueue->pop();
        hqueue->find_minlev();
    }

    void push_1stitem(ImgIdx idx, int32_t alpha);
    void push(ImgIdx idx, int32_t alpha);
    ImgIdx pop();
};

template <class Pixel> class Cache_Hierqueue_l2 {
    HQentry<int32_t> *list;
    HQueue_l2idx *hqueue;
    ImgIdx minidx_queue;
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;
    int _emptyTop;

    void initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);

  public:
    double qtime; // tmp

    Cache_Hierqueue_l2(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels);
    Cache_Hierqueue_l2(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);
    ~Cache_Hierqueue_l2();

    inline void startPushes() { _emptyTop = 1; }
    inline int32_t get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline int32_t top_alpha() { return list[0].alpha; }
    inline void endPushes() {
        if (_emptyTop)
            pop();
    }
    inline void push_queue(ImgIdx idx, int32_t alpha) { hqueue->push(idx, alpha); }
    inline void pop_queue() {
        hqueue->pop();
        hqueue->find_minlev();
    }

    void push_1stitem(ImgIdx idx, int32_t alpha);
    void push(ImgIdx idx, int32_t alpha);
    ImgIdx pop();
};

class Trie_Cache {
    ImgIdx *list;
    Trie<uint64_t> *trie;
    ImgIdx minidx_queue;
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;

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
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;

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
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t shamt, nbit;
    int16_t l0, mask;

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
    int64_t *levels;
    HQueue_l1idx *queue;
    ImgIdx minidx_queue;
    int16_t curSize_list, maxSize_list;
    ImgIdx maxSize_queue, mask_field;
    int8_t minlevnotfixed;

    void initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);

  public:
    HybridQueue_HQueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels);
    HybridQueue_HQueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize);
    ~HybridQueue_HQueue();

    inline ImgIdx top() { return list[0]; }
    inline int64_t get_minlev() { return levels[0]; }
    inline void find_minlev() { queue->find_minlev(); }
    inline void push_queue(ImgIdx idx, int64_t level) {
        if (queue->push(idx, level))
            minlevnotfixed = 0;
    }
    inline void pop_queue() {
        minlevnotfixed = 1;
        queue->pop();
    }

    void push(ImgIdx idx, int64_t level);
    ImgIdx pop();
};
