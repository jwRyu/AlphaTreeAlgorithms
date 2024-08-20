#pragma once

#include <defines.h>

template <class Pixel> class HHPQ {
    struct QItem {
        ImgIdx index = -1;
        Pixel alpha = (Pixel)0;

        QItem() = default;
        QItem(ImgIdx index_, Pixel alpha_) : index(index_), alpha(alpha_) {}
    };

    QItem *_queue = nullptr;
    static constexpr ImgIdx _cacheStart = 0;
    const ImgIdx _cacheMaxSize;
    ImgIdx _cacheCurSize = 0;

    ImgIdx *_levelStart = nullptr;
    ImgIdx *_levelCurrent = nullptr;

    ImgIdx _lowestUnsortedLevel = -1;
    ImgIdx _lowestNonemptyLevel = -1;

    const double a;

    ImgIdx _size = 0;
    const ImgIdx _sizeMax = 0;
    const ImgIdx _numLevels;

    const _uint8 *_isVisited;

    bool _isEmptyTop = false;

    QItem top(ImgIdx level) const { return _queue[_levelStart[level]]; }
    inline Pixel topAlpha(ImgIdx level) const;
    inline Pixel bottomAlphaInCache() const;
    inline void popCache();
    inline void pushToLevel(QItem item);
    inline void pushToQuadHeapQueue(QItem item, ImgIdx level);
    inline void popFromQuadHeapQueue(ImgIdx level);
    inline void sort(ImgIdx level);

    inline void clear(ImgIdx level);
    inline bool isLevelFull(ImgIdx level);
    inline bool isLevelEmpty(ImgIdx level);
    inline bool isLevelOverflowed(ImgIdx level);

  public:
    HHPQ(const ImgIdx *levelSizes, ImgIdx numLevels, ImgIdx sizeTotal, const _uint8 *isVisited_, double a_ = 32.0,
         int cacheSize = 12, int connectivity = 4);
    ~HHPQ();

    void start_pushes() { _isEmptyTop = true; }
    inline Pixel topAlpha();
    inline ImgIdx top();
    void push(ImgIdx idx, Pixel alpha);
    ImgIdx pop();
    bool isEmptyAfterSort(ImgIdx level);

    inline ImgIdx alphaToLevel(double alpha);

    // void initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity, double
    // r); HierarHeapQueue_cache(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize,
    //                       ImgIdx connectivity = 4, double r = 0.2);
    // ~HierarHeapQueue_cache();
    // void push_1stitem(ImgIdx idx, Pixel alpha);
    // void end_pushes(_uint8 *isVisited);
    // void push(ImgIdx idx, Pixel alpha);
    // void push_queue(ImgIdx idx, Pixel alpha);
    // ImgIdx pop(_uint8 *isVisited);
    // int check_queue_level(_uint8 *isVisited);
    // void pop_queue(_uint8 *isVisited);
};