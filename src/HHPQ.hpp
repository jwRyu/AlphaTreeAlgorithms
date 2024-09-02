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

    ImgIdx _lowestUnsortedLevelInitial = -1;
    ImgIdx _lowestUnsortedLevel = -1;
    const ImgIdx _numLevels = -1;
    double _a = 0.0;
    ImgIdx _lowestNonemptyLevel = -1;

    _int16 _curSizeCache = -1;
    const _int16 _maxSizeCache = -1;
    int _emptyTop = 0;

    _uint8 *_isVisited = nullptr;
    bool *_isRedundant = nullptr;

    ImgIdx _size = 0;

    bool isFrontLevelEmptyAfterSort();
    void pop_queue();
    const QItem<Pixel> &cacheBack() { return _cache[_curSizeCache - 1]; }

    void push_queue(const QItem<Pixel> &item, const ImgIdx &level);

    void initHQ(ImgIdx *dhist, ImgIdx size, double r);

  public:
    /// @brief  Print the queue in terminal
    void print();

    ImgIdx size() const { return _size; }

    /// @brief Create an HHPQ object
    /// @param dhist Histogram of all level items in the queue. Use alphaToLevel() on the image to build the histogram
    /// @param numLevels_ Number of levels in the HHPQ, which equates to the dynamic range of dhist
    /// @param size Total number of items to be processed by HHPQ, which should be equal to accumulate(dhist(:)).
    /// @param isVisited_ isVisited array used during flooding, which helps reducing the computational cost of HHPQ
    /// @param a_ Scaler to adjust the number of levels in HHPQ. Higher this value is, the more number levels in the
    /// HHPQ
    /// @param cacheSize Size of the cache. Cache size Bigger than 20 may decrease speed
    /// @param r Percentile of the number of sorted levels in the initial HHPQ setup. Should correlate to the ratio of
    /// redundant edges, but this parameter does not terribly affect the performace. Using default value is recommended.
    HHPQ(ImgIdx *dhist, ImgIdx numLevels_, ImgIdx size, _uint8 *isVisited_, double a_ = 15.0, int cacheSize = 15,
         double r = 0.2, bool *isRedundant = nullptr);
    ~HHPQ();

    /// @brief Convert alpha value to HHPQ level. Use this outside the class to build the pixel difference histogram
    /// (dhist).
    /// @param alpha Pixel dissimilarity value
    /// @param a Scaler constant for the HHPQ
    /// @return Level index converted from the alpha value
    static ImgIdx alphaToLevel(const double &alpha, const double &a);
    ImgIdx alphaToLevel(const double &alpha) const;

    /// @brief Make the HHPQ ready for pushing multiple items followed by a single pop()
    /// @details This helps reducing memory swaps by simply replacing the top item instead of running separate pop() and
    /// push() operations, when the conditions are met. _emptytop represented that the top item should have been emptied
    /// but its operation has been amortized. If an item with an alpha value <= than the top it will replace the top
    /// (which should happen around 50% of the time), effectively removing computational costs of pop() and push().
    // void startPushes();

    /// @brief Finish pushing multiple items
    /// @details If _emptytop is still true after pushing multiple items, which means the top item hasn't been replaced,
    /// run pop() operation. Otherwise, do nothing.
    // void endPushes();

    /// @brief Get the top item of the queue
    /// @return Const reference of the top item. The behavior is undefined if the queue was empty
    const QItem<Pixel> &front() const { return _cache[0]; }

    /// @brief Push an item into HHPQ
    /// @param idx Index element of the new item
    /// @param alpha Alpha element of the new item
    void push(const ImgIdx &idx, const Pixel &alpha = std::numeric_limits<Pixel>::max());

    /// @brief Check if the queue is empty
    /// @return True if the queue is empty
    bool empty() const { return _size == 0; }

    /// @brief Pop the top item from the queue
    void pop();
};