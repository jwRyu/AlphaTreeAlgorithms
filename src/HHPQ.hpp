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
    double _a;
    ImgIdx _lowestNonemptyLevel;

    _int16 _curSizeCache, _maxSizeCache;
    int emptytop;

    ImgIdx maxSize;

    _uint8 *_isVisited;

  public:
    void print();
    void initHQ(ImgIdx *dhist, ImgIdx size, double r);
    HHPQ(ImgIdx *dhist, ImgIdx numLevels_, ImgIdx size, _uint8 *isVisited_, double a_ = 15.0, int cacheSize = 15,
         double r = 0.2);
    ~HHPQ();

    static ImgIdx alphaToLevel(const double &alpha, const double &a);

    void startPushes() { emptytop = 1; }
    void endPushes();

    Pixel get_minlev() { return _cache[0].alpha; }
    const QItem<Pixel> &top() { return _cache[0]; }
    const QItem<Pixel> &cacheBack() { return _cache[_curSizeCache]; }

    void push(const ImgIdx &idx, const Pixel &alpha = std::numeric_limits<Pixel>::max());
    void push_queue(const QItem<Pixel> &item, const ImgIdx &level);
    void pop();
    int check_queue_level();
    void pop_queue();
};