#pragma once

#include "allocator.h"
#include "defines.h"

// Do not use beyond 20-bit image
class HierarQueue {
  public:
    ImgIdx *queue;
    ImgIdx *bottom, *cur;
    int32_t numlevel;
    int64_t qsize;
    int64_t min_level, max_level;

    void print();
    HierarQueue(uint64_t qsize_in, int32_t numlevels);
    HierarQueue(uint64_t qsize_in);
    void reset_queue();
    ImgIdx set_queue(ImgIdx *dhist);
    ImgIdx set_queue(ImgIdx *dhist, int32_t maxpix);

    HierarQueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels);
    HierarQueue(ImgIdx *dhist, int32_t numlevels);
    HierarQueue(int32_t numlevels, ImgIdx binsize);
    ~HierarQueue();

    int8_t push(ImgIdx pidx, int64_t level);
    void find_minlev();

    inline ImgIdx pop() { return queue[bottom[min_level]++]; }
    inline ImgIdx top() { return queue[bottom[min_level]]; }
    inline int64_t get_minlev() { return min_level; }
};

class HQueue_l1idx {
  public:
    ImgIdx *queue;
    ImgIdx *bottom, *cur;
    uint64_t *seeker;

    int64_t qsize, seekersize;
    int32_t min_level;
    HQueue_l1idx(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels);
    ~HQueue_l1idx();

    int push(ImgIdx pidx, int32_t level);
    ImgIdx pop();
    void find_minlev();

    inline ImgIdx top() { return queue[cur[min_level] - 1]; }
    inline int32_t get_minlev() { return min_level; }
};

class HQueue_l2idx {
    ImgIdx *queue;
    ImgIdx *bottom, *cur;
    uint64_t *seeker, *seeker2;

  public:
    int64_t qsize;
    int64_t min_level;
    HQueue_l2idx(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels);
    ~HQueue_l2idx();

    void push(ImgIdx pidx, int64_t level);
    ImgIdx pop();
    void find_minlev();

    inline ImgIdx top() { return queue[cur[min_level] - 1]; }
    inline int64_t get_minlev() { return min_level; }
};

class HQueue_l1idx_rank {
    struct hqueue_word {
        int64_t qword[64];
        int64_t seeker;
    };

    hqueue_word *queue;

  public:
    int64_t qsize, seekersize;
    int64_t min_level;
    HQueue_l1idx_rank(int64_t qsize_in);
    ~HQueue_l1idx_rank();

    void push(ImgIdx pidx);
    void pop();
    void find_minlev();

    inline ImgIdx top() { return min_level; }
    inline int64_t get_minlev() { return min_level; }
};