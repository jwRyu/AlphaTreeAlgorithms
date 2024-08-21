#pragma once

#include <cassert>
#include <cstdio>
#include <defines.h>
#include <limits>

template <class Pixel> class HHPQ {
  public:
    struct QItem {
        ImgIdx index = -1;
        double alpha = 0.0;

        QItem() = default;
        QItem(ImgIdx index_, double alpha_) : index(index_), alpha(alpha_) {}
        bool operator<(const QItem &other) const { return alpha < other.alpha; }
        void print() { printf("(%d, %f) ", (int)index, (double)alpha); }
    };
    HHPQ(const ImgIdx *levelSizes, ImgIdx numLevels, ImgIdx sizeTotal, const _uint8 *isVisited_, double a_ = 32.0,
         int cacheSize = 12, int connectivity = 4);
    ~HHPQ();

    void start_pushes() { _isEmptyFront = true; }
    void end_pushes();
    // inline Pixel frontAlpha();
    inline QItem front() const { return _queue[0]; }
    void push(ImgIdx idx, double alpha);
    ImgIdx pop();
    inline bool empty() const { return _size == 0; }

    void print();

    static inline ImgIdx alphaToLevel(const double &a, const double &alpha);

  private:
    QItem *_queue = nullptr;
    static constexpr ImgIdx _cacheStart = 0;
    const ImgIdx _cacheMaxSize;
    ImgIdx _cacheCurSize = 0;

    ImgIdx *_levelStart = nullptr;
    ImgIdx *_levelCurrent = nullptr;

    ImgIdx _lowestUnsortedLevel = -1;
    ImgIdx _lowestNonemptyLevel = -1;

    const double _a;

    ImgIdx _size = 0;
    const ImgIdx _sizeMax = 0;
    const ImgIdx _numLevels;

    const _uint8 *_isVisited;

    bool _isEmptyFront = false;

    static inline ImgIdx quadHeapFirstChild(ImgIdx parent);
    static inline ImgIdx quadHeapParent(ImgIdx child);

    QItem front(ImgIdx level) const {
        assert(_size > 0 && _cacheCurSize > 0);
        return _queue[_levelStart[level]];
    }
    QItem cacheBack() const { return _queue[_cacheCurSize - 1]; }
    inline double frontAlpha(ImgIdx level) const;
    inline ImgIdx findNextNonemptyLevel(ImgIdx level) const;
    inline void popCache();
    inline void pushToLevel(QItem item);
    inline void pushToQuadHeapQueue(QItem item, ImgIdx level);
    inline void popFromQuadHeapQueue(ImgIdx level);
    inline void sort(ImgIdx level);
    inline ImgIdx size(ImgIdx level) const;

    inline void clear(ImgIdx level);
    inline bool isLevelFull(ImgIdx level) const;
    inline bool isLevelEmpty(ImgIdx level) const;
    inline bool isLevelOverflowed(ImgIdx level) const;
};