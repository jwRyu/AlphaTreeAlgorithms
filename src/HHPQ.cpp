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
    : _cacheMaxSize(cacheSize), _a(a_), _sizeMax(sizeTotal + cacheSize), _numLevels(numLevels), _isVisited(isVisited_) {

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
    _levelCurrent = new ImgIdx[_numLevels + 1];

    ImgIdx n = _cacheMaxSize;
    for (ImgIdx level = 0; level < _numLevels; level++) {
        _levelStart[level] = _levelCurrent[level] = n;
        n += levelSizes[level];
    }
    _levelCurrent[_numLevels] = _levelStart[_numLevels] = n;
    _lowestNonemptyLevel = _numLevels;
    _lowestUnsortedLevel = (ImgIdx)(_numLevels / (connectivity / 2));
}

template <class Pixel> HHPQ<Pixel>::~HHPQ() {
    delete[] _queue;
    delete[] _levelStart;
    delete[] _levelCurrent;
}

// template <class Pixel> Pixel HHPQ<Pixel>::frontAlpha() {
// #ifdef BOUNDARYCHECK
//     assert(_size > 0 && _cacheMaxSize > 0);
// #endif
//     return _queue[0].alpha;
// }

// template <class Pixel> ImgIdx HHPQ<Pixel>::front() {
// #ifdef BOUNDARYCHECK
//     assert(_size > 0 && _cacheMaxSize > 0);
// #endif
//     return _queue[0].index;
// }

template <class Pixel> void HHPQ<Pixel>::print() {
    printf("---------- HHPQ<Pixel>::print START -------------\n");
    printf("size = %d\n", _size);
    printf("Cache[%d / %d]: ", _cacheCurSize, _cacheMaxSize);
    for (int i = 0; i < _cacheCurSize; i++)
        _queue[i].print();
    printf("\n");

    for (int level = 0; level < _numLevels; level++) {
        if (isLevelEmpty(level))
            continue;
        printf("level[%d][%d / %d]: ", level, _levelCurrent[level] - _levelStart[level],
               _levelStart[level + 1] - _levelStart[level]);
        for (ImgIdx qIdx = _levelStart[level]; qIdx < _levelCurrent[level]; qIdx++)
            _queue[qIdx].print();
        printf("\n");
    }

    printf("---------- HHPQ<Pixel>::print END -------------\n");
    // std::getchar();
}

template <class Pixel> ImgIdx HHPQ<Pixel>::alphaToLevel(const double &a, const double &alpha) {
    return (ImgIdx)(a * log2(1.0 + alpha));
}

template <class Pixel> double HHPQ<Pixel>::frontAlpha(ImgIdx level) const {
#ifdef BOUNDARYCHECK
    assert(_size > 0);
    assert(level >= 0 && level <= _numLevels);
    assert(!isLevelEmpty(level) && !isLevelOverflowed(level));
#endif
    return _queue[_levelStart[level]].alpha;
}

template <class Pixel> void HHPQ<Pixel>::end_pushes() {
    if (_isEmptyFront)
        pop();
    _isEmptyFront = false;
}

template <class Pixel> void HHPQ<Pixel>::push(ImgIdx index, double alpha, ImgIdx edgeIdx) {
    printf("HHPQ::Pushing %d at %.2f \n", index, (double)alpha);
#ifdef BOUNDARYCHECK
    assert(_sizeMax > 0 && _cacheMaxSize > 0);
    assert(_cacheCurSize <= _cacheMaxSize);
#endif

    // First ever push
    if (_cacheCurSize == 0 && _size == 0) {
        _size = _cacheCurSize = 1;
        _queue[0] = QItem(index, alpha, edgeIdx);
        return;
    }

    if (_isEmptyFront && alpha <= front().alpha) {
        _isEmptyFront = false;
        _queue[0] = QItem(index, alpha, edgeIdx);
        return;
    }

    const ImgIdx newItemLevel = alphaToLevel(_a, alpha);
    const bool isFrontLevelSorted = _lowestNonemptyLevel < _lowestUnsortedLevel;
    const bool cacheFull = (_cacheCurSize == _cacheMaxSize);
    const bool pushToCache = ((newItemLevel < _lowestNonemptyLevel) ||
                              (isFrontLevelSorted && (alpha <= frontAlpha(_lowestNonemptyLevel)))) &&
                             (!cacheFull || alpha < cacheBack().alpha);

    if (pushToCache) {
        ImgIdx i = _cacheCurSize;
        if (cacheFull) {
#ifdef BOUNDARYCHECK
            assert(alpha < _queue[_cacheCurSize - 1].alpha);
#endif
            pushToLevel(cacheBack());
        } else
            _cacheCurSize++;

        for (; i > 0 && alpha < _queue[i - 1].alpha; i--)
            _queue[i] = _queue[i - 1];
        _queue[i] = QItem(index, alpha, edgeIdx);
    } else
        pushToLevel(QItem(index, alpha, edgeIdx));
    _size++;
}

template <class Pixel> void HHPQ<Pixel>::pushToLevel(QItem item) {
    const ImgIdx level = alphaToLevel(_a, item.alpha);

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

template <class Pixel> ImgIdx HHPQ<Pixel>::quadHeapFirstChild(ImgIdx parent) { return (parent << 2) - 2; }
template <class Pixel> ImgIdx HHPQ<Pixel>::quadHeapParent(ImgIdx child) { return (child + 2) >> 2; }

template <class Pixel> void HHPQ<Pixel>::pushToQuadHeapQueue(QItem item, ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
    assert(!isLevelFull(level) && !isLevelOverflowed(level));
#endif

    _levelCurrent[level]++;

    const ImgIdx begin = _levelStart[level];
#ifdef BOUNDARYCHECK
    const ImgIdx end = _levelStart[level + 1];
    assert(begin >= 0 && end <= _sizeMax);
    assert(begin < end);
#endif
    QItem *heapQueue = _queue + begin - 1;
    ImgIdx current = _levelCurrent[level] - begin;
    for (ImgIdx parent = quadHeapParent(current); parent > 0 && heapQueue[parent].alpha > item.alpha;
         parent = quadHeapParent(current)) {
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

template <class Pixel> bool HHPQ<Pixel>::isLevelFull(ImgIdx level) const {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
#endif
    return _levelCurrent[level] == _levelStart[level + 1];
}

template <class Pixel> bool HHPQ<Pixel>::isLevelEmpty(ImgIdx level) const {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
#endif
    return _levelCurrent[level] == _levelStart[level];
}

template <class Pixel> bool HHPQ<Pixel>::isLevelOverflowed(ImgIdx level) const {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
#endif
    return _levelCurrent[level] > _levelStart[level + 1];
}

template <class Pixel> inline ImgIdx HHPQ<Pixel>::size(ImgIdx level) const {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
    assert(!isLevelOverflowed(level));
#endif
    return _levelCurrent[level] - _levelStart[level];
}

template <class Pixel> void HHPQ<Pixel>::sort(ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
    assert(!isLevelOverflowed(level));
#endif
    if (isLevelEmpty(level))
        return;

    const ImgIdx levelSizeBefore = size(level);
    QItem *copy = new QItem[levelSizeBefore];
    memcpy(copy, _queue + _levelStart[level], sizeof(QItem) * (size_t)levelSizeBefore);

    clear(level);
    for (ImgIdx i = 0; i < levelSizeBefore; i++) {
        QItem &item = copy[i];
        if (_isVisited[item.index]) { // Skip redundant entry
            if (edge != nullptr) {    // temp
                edge[item.edgeIdx] = 4;
            }
            continue;
        }
        pushToQuadHeapQueue(item, level);
    }

    const ImgIdx levelSizeAfter = size(level);
    _size = _size - (levelSizeBefore - levelSizeAfter);

    delete[] copy;

#ifdef BOUNDARYCHECK
    assert(!isLevelOverflowed(level));
    assert(_size >= 0);
#endif
}

template <class Pixel> void HHPQ<Pixel>::popFromQuadHeapQueue(ImgIdx level) {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
    assert(!isLevelEmpty(level));
    assert(!isLevelOverflowed(level));
#endif

    _levelCurrent[level]--;

    const ImgIdx curSize = _levelCurrent[level] - _levelStart[level];

    if (curSize == 0) // Empty after pop - no need for replacement
        return;

    QItem *heapQueue = _queue + _levelStart[level] - 1;
    QItem lastItem = heapQueue[curSize + 1];
    ImgIdx newSpotForLastItem = 1;
    while (quadHeapFirstChild(newSpotForLastItem) <= curSize) {
        const ImgIdx child0 = quadHeapFirstChild(newSpotForLastItem);
        const ImgIdx child1 = child0 + 1;
        const ImgIdx child2 = child0 + 2;
        const ImgIdx child3 = child0 + 3;

        ImgIdx bestChild = child0;
        // clang-format off
        if (child1 <= curSize && heapQueue[child1].alpha < heapQueue[bestChild].alpha)   bestChild = child1;
        if (child2 <= curSize && heapQueue[child2].alpha < heapQueue[bestChild].alpha)   bestChild = child2;
        if (child3 <= curSize && heapQueue[child3].alpha < heapQueue[bestChild].alpha)   bestChild = child3;
        // clang-format on

        if (lastItem.alpha < heapQueue[bestChild].alpha)
            break;

        heapQueue[newSpotForLastItem] = heapQueue[bestChild];
        newSpotForLastItem = bestChild;
    }
    heapQueue[newSpotForLastItem] = lastItem;
}

template <class Pixel> ImgIdx HHPQ<Pixel>::findNextNonemptyLevel(ImgIdx level) const {
#ifdef BOUNDARYCHECK
    assert(level >= 0 && level < _numLevels);
#endif
    while (isLevelEmpty(level)) {
        level++;
#ifdef BOUNDARYCHECK
        assert(level < _numLevels);
#endif
    }
    return level;
}

template <class Pixel> ImgIdx HHPQ<Pixel>::pop() {
    printf("HHPQ::pop %d at %.2f \n", front().index, (double)front().alpha);
    if (_size == 0) {
        return -1;
    }

    _size--;
    ImgIdx ret = _queue[0].index;
    if (_cacheCurSize) {
        popCache();

        // Pop one item from the levels if the cache is empty. Cache must always be populated (for front())
        if (_cacheCurSize == 0) {
            do {
                if (_size <= 0) {
                    print();
                    return -1;
                }
                _lowestNonemptyLevel = findNextNonemptyLevel(_lowestNonemptyLevel);
#ifdef BOUNDARYCHECK
                assert(_lowestNonemptyLevel < _numLevels);
#endif
                sort(_lowestNonemptyLevel);
            } while (isLevelEmpty(_lowestNonemptyLevel));
            _lowestUnsortedLevel = _lowestNonemptyLevel + 1;
            _queue[0] = front(_lowestNonemptyLevel);
            _cacheCurSize++;
            popFromQuadHeapQueue(_lowestNonemptyLevel);
            _lowestNonemptyLevel = findNextNonemptyLevel(_lowestNonemptyLevel);
        }
    }

    return ret;
}

template <class Pixel> void HHPQ<Pixel>::popCache() {
    _cacheCurSize--;
    for (int i = 0; i < _cacheCurSize; i++)
        _queue[i] = _queue[i + 1];
}

template class HHPQ<_uint8>;
template class HHPQ<_uint16>;
template class HHPQ<_uint32>;
template class HHPQ<_uint64>;
