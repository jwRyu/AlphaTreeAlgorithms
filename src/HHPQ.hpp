#pragma once

#include <HeapQueue.h>
#include <QItem.hpp>
#include <QuadHeapQueue.hpp>
#include <cstdio>
#include <defines.h>

template <class Pixel> class HHPQ {
    QItem<Pixel> *_cache = nullptr;
    QuadHeapQueue<Pixel> **_sortedLevels = nullptr;
    QItem<Pixel> **_unsortedLevels = nullptr;
    ImgIdx *_unsortedLevelSizes = nullptr;
    ImgIdx *_levelMaxSizes = nullptr;

    ImgIdx _lowestUnsortedLevelAllocated = -1;
    ImgIdx _lowestUnsortedLevel = -1;
    ImgIdx _numLevels = -1;
    double _a = 0.0;
    ImgIdx _lowestNonemptyLevel = -1;

    _int16 _curSizeCache = -1;
    _int16 _maxSizeCache = -1;
    int emptytop = 0;

    _uint8 *_isVisited = nullptr;

  public:
    void print();
    void initHQ(ImgIdx *dhist, ImgIdx size, double r);
    HHPQ(ImgIdx *dhist, ImgIdx numLevels_, ImgIdx size, _uint8 *isVisited_, double a_ = 15.0, int cacheSize = 15,
         double r = 0.2);
    ~HHPQ();
    void push_queue(const QItem<Pixel> &item, const ImgIdx &level);

    static ImgIdx alphaToLevel(const double &alpha, const double &a);

    void startPushes() { emptytop = 1; }
    void endPushes();

    const QItem<Pixel> &front() { return _cache[0]; }
    const QItem<Pixel> &cacheBack() { return _cache[_curSizeCache]; }

    void push(const ImgIdx &idx, const Pixel &alpha = std::numeric_limits<Pixel>::max());
    void pop();
    bool isFrontLevelEmptyAfterSort();
    void pop_queue();
};