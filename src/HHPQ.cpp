#include <HHPQ.hpp>

#include <allocator.h>
#include <cassert>
#include <cmath>
#include <cstring>
#include <defines.h>

#define BOUNDARYCHECK 1

template <class Pixel>
HHPQ<Pixel>::HHPQ(const ImgIdx *levelSizes, ImgIdx numLevels, ImgIdx sizeTotal, const _uint8 *isVisited_, double a_,
                  ImgIdx cacheSize, int connectivity)
    : _cacheMaxSize(cacheSize), a(a_), _sizeMax(sizeTotal + cacheSize), _numLevels(numLevels), _isVisited(isVisited_) {

    // Validate input parameters
    assert(numLevels > 0 && a_ > 0.0 && cacheSize > 0 &&
           (connectivity == 4 || connectivity == 8 || connectivity == 12));
    ImgIdx levelSizeSum = 0;
    for (int level = 0; level < numLevels; level++)
        levelSizeSum += levelSizes[level];
    assert(levelSizeSum == sizeTotal);

    // Allocate memory
    _queue = new QItem[_sizeMax];
    _levelStart = new ImgIdx[_numLevels + 1];
    _levelCurrent = new ImgIdx[_numLevels];

    ImgIdx n = _cacheMaxSize;
    for (ImgIdx level = 0; level < _numLevels; level++) {
        _levelStart[level] = _levelCurrent[level] = n;
        n += levelSizes[level];
    }
    _levelStart[_numLevels] = n;
}

template <class Pixel> HHPQ<Pixel>::~HHPQ() {
    delete[] _queue;
    delete[] _levelStart;
    delete[] _levelCurrent;
}

template <class Pixel> Pixel HHPQ<Pixel>::topAlpha() {
#ifdef BOUNDARYCHECK
    assert(_size > 0 && _cacheMaxSize > 0);
#endif
    return _queue[0].alpha;
}

template <class Pixel> ImgIdx HHPQ<Pixel>::top() {
#ifdef BOUNDARYCHECK
    assert(_size > 0 && _cacheMaxSize > 0);
#endif
    return _queue[0].index;
}

template <class Pixel> ImgIdx HHPQ<Pixel>::alphaToLevel(double alpha) { return (ImgIdx)(a * log2(1.0 + alpha)); }

template <class Pixel> Pixel HHPQ<Pixel>::topAlpha(ImgIdx level) const {
#ifdef BOUNDARYCHECK
    assert(_size > 0);
    assert(level >= 0 && level <= _numLevels);
    assert(!isLevelEmpty(level) && !isLevelOverflowed(level));
#endif
    return _queue[_levelStart[level]].alpha;
}

template <class Pixel> Pixel HHPQ<Pixel>::bottomAlphaInCache() const {
#ifdef BOUNDARYCHECK
    assert(_size > 0 && _cacheCurSize > 0 && _cacheCurSize <= _cacheMaxSize);
#endif
    return _queue[_cacheCurSize - 1].alpha;
}

template <class Pixel> void HHPQ<Pixel>::push(ImgIdx index, Pixel alpha) {
#ifdef BOUNDARYCHECK
    assert(_size > 0 && _cacheMaxSize > 0);
    assert(_cacheCurSize <= _cacheMaxSize);
#endif
    if (_cacheCurSize == 0) {
#ifdef BOUNDARYCHECK
        assert(_isEmptyTop == false && _size == 0);
#endif
        _size = _cacheMaxSize = 1;
        _queue[0] = QItem(index, alpha);
        return;
    }

    if (_isEmptyTop && alpha < topAlpha()) {
        _isEmptyTop = false;
        _queue[0] = QItem(index, alpha);
        return;
    }

    const ImgIdx newItemLevel = alphaToLevel(alpha);

    const bool isTopLevelSorted = _lowestNonemptyLevel < _lowestUnsortedLevel;
    const ImgIdx topLevelStart = _levelStart[_lowestNonemptyLevel];
    const bool cacheFull = _cacheCurSize == _cacheMaxSize;
    const bool pushToCache =
        ((newItemLevel < _lowestNonemptyLevel) || (isTopLevelSorted && (alpha < topAlpha(_lowestNonemptyLevel)))) &&
        (!cacheFull || alpha < bottomAlphaInCache());

    if (pushToCache) {
        if (cacheFull) {
#ifdef BOUNDARYCHECK
            assert(alpha < _queue[_cacheCurSize - 1].alpha);
#endif
            pushToLevel();
            _cacheCurSize--;
        }

        ImgIdx i = _cacheCurSize++;
        for (; i > 0 && alpha < _queue[i - 1].alpha; i--)
            _queue[i] = _queue[i - 1];
        _queue[i] = QItem(index, alpha);
    } else
        pushToLevel(QItem(index, alpha));
    _size++;
}

template <class Pixel> void HHPQ<Pixel>::pushToLevel(QItem item) {
    const ImgIdx level = alphaToLevel(item.alpha);

#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
    assert(!isLevelFull(level) && !isLevelOverflowed(level));
#endif

    if (level < _lowestUnsortedLevel) {
        if (level < _lowestNonemptyLevel)
            _lowestNonemptyLevel = level;
        // HEAP QUEUE
        pushToQuadHeapQueue(item, level);
    } else
        _queue[_levelCurrent[level]++] = item;

#ifdef BOUNDARYCHECK
    assert(!isLevelOverflowed(level));
#endif
}

template <class Pixel> void HHPQ<Pixel>::pushToQuadHeapQueue(QItem item, ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
    assert(!isLevelFull(level) && !isLevelOverflowed(level));
#endif

    ImgIdx begin = _levelStart[level];
    ImgIdx curSize = ++_levelCurrent[level];

#ifdef BOUNDARYCHECK
    ImgIdx end = _levelStart[level + 1];
    assert(begin >= 0 && end <= _sizeMax);
    assert(begin < end);
#endif

    QItem *heapQueue = _queue + begin - 1;
    ImgIdx current = curSize - begin;
    for (ImgIdx parent = current >> 2; parent > 0 && heapQueue[parent].alpha > item.alpha; parent >>= 2) {
        heapQueue[current] = heapQueue[parent];
        current = parent;
    }
    heapQueue[current] = item;

#ifdef BOUNDARYCHECK
    assert(!isLevelOverflowed(level));
#endif
}

template <class Pixel> void HHPQ<Pixel>::clear(ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
    assert(!isLevelEmpty(level));
#endif
    _levelCurrent[level] = _levelStart[level];
}

template <class Pixel> bool HHPQ<Pixel>::isLevelFull(ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
#endif
    return _levelCurrent[level] == _levelStart[level + 1];
}

template <class Pixel> bool HHPQ<Pixel>::isLevelEmpty(ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
#endif
    return _levelCurrent[level] == _levelStart[level];
}

template <class Pixel> bool HHPQ<Pixel>::isLevelOverflowed(ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
#endif
    return _levelCurrent[level] > _levelStart[level + 1];
}

template <class Pixel> void HHPQ<Pixel>::sort(ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
#endif
    if (isLevelEmpty(level))
        return;

    ImgIdx begin = _levelStart[level];
    ImgIdx curSize = _levelCurrent[level];

    const ImgIdx levelSize = curSize - begin;
    QItem *copy = new QItem[levelSize];
    memcpy(copy, _queue + begin, sizeof(QItem) * (size_t)levelSize);

    clear(level);

    for (int i = 0; i < levelSize; i++) {
        QItem &item = copy[i];
        if (_isVisited[item.index]) // Skip redundant entry
            continue;
        pushToQuadHeapQueue(item, level);
    }

    delete[] copy;

#ifdef BOUNDARYCHECK
    assert(!isLevelOverflowed(level));
#endif
}

template <class Pixel> bool HHPQ<Pixel>::isEmptyAfterSort(ImgIdx level) {
    while (_lowestUnsortedLevel <= level) {
        sort(_lowestUnsortedLevel);
        _lowestUnsortedLevel++;
    }
    return isLevelEmpty(level);
}

template <class Pixel> void HHPQ<Pixel>::popFromQuadHeapQueue(ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
    assert(!isLevelEmpty(level));
#endif

    _levelCurrent[level]--;

    ImgIdx begin = _levelStart[level];
    const ImgIdx curSize = _levelCurrent[level] - begin;

    if (curSize == 0) // Empty after pop - no need for replacement
        return;

    QItem *heapQueue = _queue + begin - 1;
    QItem lastItem = heapQueue[curSize + 1];
    ImgIdx newSpotForLastItem = 1;
    while (true) {
        const ImgIdx child0 = newSpotForLastItem << 2;
        const ImgIdx child1 = newSpotForLastItem << 2;
        const ImgIdx child2 = newSpotForLastItem << 2;
        const ImgIdx child3 = newSpotForLastItem << 2;
        ImgIdx bestChild = firstChild;
        if (firstChild + 3 <= curSize) {
            // clang-format off
            if (arr[firstChild + 1].alpha < arr[next].alpha)
                next = next0 + 1;
            if (arr[next0 + 2].alpha < arr[next].alpha)
                next = next0 + 2;
            if (arr[next0 + 3].alpha < arr[next].alpha)
                next = next0 + 3;
            // clang-format on
        } else {
            if (next0 > cursize)
                break;
            if (next0 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 1].alpha < arr[next].alpha)
                next = next0 + 1;
            if (next0 + 1 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 2].alpha < arr[next].alpha)
                next = next0 + 2;
            if (next0 + 2 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 3].alpha < arr[next].alpha)
                next = next0 + 3;
        }

    MIN_NEXT_FOUND:
        if (curalpha < arr[next].alpha)
            break;

        arr[current] = arr[next];
#if PROFILE
        nummove++;
#endif
        current = next;
    }

    template <class Pixel> ImgIdx HHPQ<Pixel>::pop() {
        if (_cacheCurSize)
            popCache();
        else {
            while (isEmptyAfterSort(_lowestNonemptyLevel))
                _lowestNonemptyLevel++;
            _queue[0] = top(_lowestNonemptyLevel);
            popFromQuadHeapQueue(_lowestNonemptyLevel);
        }

        return ret;
    }

    template <class Pixel> void HHPQ<Pixel>::popCache() {
        for (int i = 0; i < _cacheCurSize - 1; i++)
            _queue[i] = _queue[i + 1];
        _cacheCurSize--;
    }

    template class HHPQ<_uint8>;
    template class HHPQ<_uint16>;
    template class HHPQ<_uint32>;
    template class HHPQ<_uint64>;
