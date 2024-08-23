#pragma once

#include <HeapQueue.h>
#include <QItem.hpp>
#include <QuadHeapQueue.hpp>
#include <cstdio>
#include <defines.h>

template <class Pixel> class HHPQ {
    QItem<Pixel> *list;
    QuadHeapQueue<Pixel> **hqueue;
    QItem<Pixel> **storage;
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
    void print();
    void initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity, double r);
    HHPQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, ImgIdx connectivity = 4,
         double r = 0.2);
    ~HHPQ();

    ImgIdx alphaToLevel(const double &alpha) const;
    inline void start_pushes() { emptytop = 1; }
    inline Pixel get_minlev() { return list[0].alpha; }
    inline ImgIdx top() { return list[0].index; }
    inline Pixel top_alpha() { return list[0].alpha; }
    void push_1stitem(ImgIdx idx);
    void end_pushes(_uint8 *isVisited);
    void push(const ImgIdx &idx, const Pixel &alpha);
    void push_queue(const QItem<Pixel> &item);
    ImgIdx pop(_uint8 *isVisited);
    int check_queue_level(_uint8 *isVisited);
    void pop_queue(_uint8 *isVisited);
};