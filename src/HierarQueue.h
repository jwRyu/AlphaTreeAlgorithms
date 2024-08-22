#pragma once

#include "allocator.h"
#include "defines.h"
#include <cstdio>

// Do not use beyond 20-bit image
class HierarQueue {
  public:
    ImgIdx *queue;
    ImgIdx *bottom, *cur;
    _int32 numlevel;
    _int64 qsize;
    _int64 min_level, max_level;

    void print();
    HierarQueue(_uint64 qsize_in, _int32 numlevels);
    HierarQueue(_uint64 qsize_in);
    void reset_queue();
    ImgIdx set_queue(ImgIdx *dhist);
    ImgIdx set_queue(ImgIdx *dhist, _int32 maxpix);

    HierarQueue(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels);
    HierarQueue(ImgIdx *dhist, _int32 numlevels);
    HierarQueue(_int32 numlevels, ImgIdx binsize);
    ~HierarQueue();

    _int8 push(ImgIdx pidx, _int64 level);
    void find_minlev();

    inline ImgIdx pop() { return queue[bottom[min_level]++]; }
    inline ImgIdx top() { return queue[bottom[min_level]]; }
    inline _int64 get_minlev() { return min_level; }
};

class HQueue_l1idx {
  public:
    ImgIdx *queue;
    ImgIdx *bottom, *cur;
    _uint64 *seeker;

    _int64 qsize, seekersize;
    _int32 min_level;
    HQueue_l1idx(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels);
    ~HQueue_l1idx();

    int push(ImgIdx pidx, _int32 level);
    ImgIdx pop();
    void find_minlev();

    inline ImgIdx top() { return queue[cur[min_level] - 1]; }
    inline _int32 get_minlev() { return min_level; }
};

class HQueue_l2idx {
    ImgIdx *queue;
    ImgIdx *bottom, *cur;
    _uint64 *seeker, *seeker2;

  public:
    _int64 qsize;
    _int64 min_level;
    HQueue_l2idx(_uint64 qsize_in, ImgIdx *dhist, _int32 numlevels);
    ~HQueue_l2idx();

    void push(ImgIdx pidx, _int64 level);
    ImgIdx pop();
    void find_minlev();

    inline ImgIdx top() { return queue[cur[min_level] - 1]; }
    inline _int64 get_minlev() { return min_level; }
};

class HQueue_l1idx_rank {
    struct hqueue_word {
        _int64 qword[64];
        _int64 seeker;
    };

    hqueue_word *queue;

  public:
    _int64 qsize, seekersize;
    _int64 min_level;
    HQueue_l1idx_rank(_int64 qsize_in);
    ~HQueue_l1idx_rank();

    void push(ImgIdx pidx);
    void pop();
    void find_minlev();

    inline ImgIdx top() { return min_level; }
    inline _int64 get_minlev() { return min_level; }
};