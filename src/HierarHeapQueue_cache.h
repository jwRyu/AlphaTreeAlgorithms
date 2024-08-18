#pragma once

#include <HeapQueue.h>

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
    void initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity, double r);
    HierarHeapQueue_cache(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize,
                          ImgIdx connectivity = 4, double r = 0.2);
    ~HierarHeapQueue_cache();
    inline void start_pushes() { emptytop = 1; }
    inline Pixel get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].pidx; }
    inline Pixel top_alpha() { return list[0].alpha; }
    void push_1stitem(ImgIdx idx, Pixel alpha);
    void end_pushes(_uint8 *isVisited);
    void push(ImgIdx idx, Pixel alpha);
    void push_queue(ImgIdx idx, Pixel alpha);
    ImgIdx pop(_uint8 *isVisited);
    int check_queue_level(_uint8 *isVisited);
    void pop_queue(_uint8 *isVisited);
};