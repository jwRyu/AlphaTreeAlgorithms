#pragma once

#include <HeapQueue.h>
#include <QItem.hpp>
#include <QuadHeapQueue.hpp>
#include <cstdio>
#include <defines.h>

template <class Pixel> class HHPQ {
    QItem<Pixel> *_cache;
    QuadHeapQueue<Pixel> **_sortedLevels;
    QItem<Pixel> **_unsortedLevels;
    ImgIdx *_unsortedLevelSizes;
    ImgIdx *_levelMaxSizes;

    ImgIdx thr_hqueue, _lowestUnsortedLevel, _numLevels;
    double a;
    ImgIdx _lowestNonemptyLevel;

    _int16 curSizeCache, maxSizeCache;
    int emptytop;

    ImgIdx maxSize;

  public:
    void print();
    void initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int cacheSize, int connectivity,
                double r);
    HHPQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int cacheSize, ImgIdx connectivity = 4,
         double r = 0.2);
    ~HHPQ();

    ImgIdx alphaToLevel(const double &alpha) const;
    void start_pushes() { emptytop = 1; }
    Pixel get_minlev() { return _cache[0].alpha; }
    ImgIdx top() { return _cache[0].index; }
    const QItem<Pixel> &cacheBack() { return _cache[curSizeCache]; }

    Pixel top_alpha() { return _cache[0].alpha; }
    void push_1stitem(ImgIdx idx);
    void end_pushes(_uint8 *isVisited);
    void push(const ImgIdx &idx, const Pixel &alpha);
    void push_queue(const QItem<Pixel> &item, const ImgIdx &level);
    ImgIdx pop(_uint8 *isVisited);
    int check_queue_level(_uint8 *isVisited);
    void pop_queue(_uint8 *isVisited);
};