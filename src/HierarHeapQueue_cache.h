#pragma once

#include <HeapQueue.h>

template <class Pixel> class HierarHeapQueue_cache {
    HQentry<double> *list;
    HeapQueue_naive_quad<double> **hqueue;
    HQentry<double> **storage;
    ImgIdx *storage_cursize;
    ImgIdx *qsizes;

    ImgIdx thr_hqueue, curthr, numlevels;
    double a;
    ImgIdx queue_minlev;

    _int16 curSize_list, maxSize_list;
    int emptytop;

    ImgIdx maxSize;
    ImgIdx _size = 0;

  public:
    void initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity, double r);
    HierarHeapQueue_cache(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize,
                          ImgIdx connectivity = 4, double r = 0.2);
    ~HierarHeapQueue_cache();
    inline void start_pushes() { emptytop = 1; }
    inline double get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline double top_alpha() { return list[0].alpha; }
    void end_pushes(_uint8 *isVisited);
    void push(ImgIdx idx, double alpha);
    void push_queue(ImgIdx idx, double alpha);
    ImgIdx pop(_uint8 *isVisited);
    int check_queue_level(_uint8 *isVisited);
    void pop_queue(_uint8 *isVisited);
    bool empty() const { return _size == 0; }
};