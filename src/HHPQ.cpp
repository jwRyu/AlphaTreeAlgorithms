#include <HHPQ.hpp>
#include <allocator.h>
#include <cmath>
#include <cstring>

template <class Pixel> ImgIdx HHPQ<Pixel>::alphaToLevel(const double &alpha, const double &a) {
    return (ImgIdx)(a * log2(1.0 + alpha));
}

template <class Pixel>
HHPQ<Pixel>::HHPQ(ImgIdx *dhist, ImgIdx numLevels_, ImgIdx size, _uint8 *isVisited_, double a_, int cacheSize, double r)
    : _numLevels(numLevels_), _a(a_), _lowestNonemptyLevel(numLevels_), _maxSizeCache(cacheSize - 1),
      _isVisited(isVisited_) {
    initHQ(dhist, size, r);
}

template <class Pixel> HHPQ<Pixel>::~HHPQ() {
    Free(_cache);
    Free(_levelMaxSizes);

    for (int level = 0; level < _numLevels; level++)
        if (_sortedLevels[level])
            delete _sortedLevels[level];
    Free(_sortedLevels);

    if (_unsortedLevels) {
        Free(_unsortedLevelSizes);
        for (int level = _lowestUnsortedLevelAllocated; level < _numLevels; level++)
            if (_unsortedLevels[level])
                Free(_unsortedLevels[level]);
        Free(_unsortedLevels + _lowestUnsortedLevelAllocated);
    }
}

template <class Pixel> void HHPQ<Pixel>::initHQ(ImgIdx *dhist, ImgIdx size, double r) {
    _cache = (QItem<Pixel> *)Malloc((_maxSizeCache + 1) * sizeof(QItem<Pixel>));
    _curSizeCache = -1;
    ImgIdx cumsum = 0;
    _levelMaxSizes = (ImgIdx *)Calloc((size_t)_numLevels * sizeof(ImgIdx));
    memcpy(_levelMaxSizes, dhist, (size_t)_numLevels * sizeof(ImgIdx));
    if (r >= 1) {
        _lowestUnsortedLevelAllocated = _lowestUnsortedLevel = _numLevels;
        _sortedLevels = (QuadHeapQueue<Pixel> **)Calloc(_numLevels * sizeof(QuadHeapQueue<Pixel> *));
        for (int level = 0; level < _lowestUnsortedLevelAllocated; level++)
            _sortedLevels[level] = new QuadHeapQueue<Pixel>(_levelMaxSizes[level]);
        _unsortedLevels = 0;
        _unsortedLevelSizes = 0;
    } else {
        _unsortedLevelSizes = (ImgIdx *)Calloc(_numLevels * sizeof(ImgIdx));
        ImgIdx thr_nonredundantnodes = (ImgIdx)(size * r);
        for (int level = 0; level < _numLevels; level++) {
            cumsum += _levelMaxSizes[level];
            if (cumsum > thr_nonredundantnodes) {
                _lowestUnsortedLevelAllocated = _lowestUnsortedLevel = level;
                break;
            }
        }

        _sortedLevels = (QuadHeapQueue<Pixel> **)Calloc(_numLevels * sizeof(QuadHeapQueue<Pixel> *));
        for (int level = 0; level < _lowestUnsortedLevelAllocated; level++)
            _sortedLevels[level] = new QuadHeapQueue<Pixel>(_levelMaxSizes[level]);

        _unsortedLevels =
            (QItem<Pixel> **)Calloc((_numLevels - _lowestUnsortedLevelAllocated) * sizeof(QItem<Pixel> *));
        _unsortedLevels -= _lowestUnsortedLevelAllocated;
        for (int level = _lowestUnsortedLevelAllocated; level < _numLevels; level++)
            _unsortedLevels[level] = (QItem<Pixel> *)Malloc(_levelMaxSizes[level] * sizeof(QItem<Pixel>));
    }
}

template <class Pixel> void HHPQ<Pixel>::print() {
    printf("---------- HHPQ<Pixel>::print START -------------\n");
    printf("Cache[%d / %d]: ", _curSizeCache + 1, _maxSizeCache + 1);
    size_t size = _curSizeCache;
    for (int i = -1; i < _curSizeCache; i++)
        _cache[i + 1].print();
    printf("\n");

    for (int level = 0; level < _numLevels; level++) {
        if (level < _lowestUnsortedLevel) {
            size += _sortedLevels[level]->size();
            if (_sortedLevels[level]->empty())
                continue;
            printf("Q: level[%d][%d / %d]: ", level, _sortedLevels[level]->size(), _sortedLevels[level]->sizeMax());
            _sortedLevels[level]->print();
        } else {
            if (_unsortedLevelSizes[level] == 0)
                continue;
            size += _unsortedLevelSizes[level];
            printf("S: level[%d][%d / %d]: ", level, _unsortedLevelSizes[level], _levelMaxSizes[level]);
            for (int i = 0; i < _unsortedLevelSizes[level]; i++)
                _unsortedLevels[level][i].print();
            printf("\n");
        }
    }

    printf("size = %d\n", (int)size);
    printf("---------- HHPQ<Pixel>::print END -------------\n");
    // std::getchar();
}

template <class Pixel> void HHPQ<Pixel>::endPushes() {
    if (_emptyTop)
        pop();
}

template <class Pixel> void HHPQ<Pixel>::push(const ImgIdx &idx, const Pixel &alpha) {
    // printf("Pushing %d at %.2f\n", idx, (double)alpha);

    const QItem<Pixel> newItem(idx, alpha);

    if (_curSizeCache == -1) {
        _curSizeCache++;
        _cache[0] = newItem;
        return;
    }

    if (_emptyTop && newItem < _cache[0]) {
        _emptyTop = 0;
        _cache[0] = newItem;
        return;
    }

    const ImgIdx newItemLevel = alphaToLevel(newItem.alpha, _a);
    const bool isCacheNotFull = _curSizeCache < _maxSizeCache;
    const bool isLowerThanCacheBack = newItem < cacheBack();
    const bool isFrontLevelSorted = _lowestNonemptyLevel < _lowestUnsortedLevel;
    const bool isLowerThanLevelFront =
        isFrontLevelSorted ? newItem < _sortedLevels[_lowestNonemptyLevel]->top() : newItemLevel < _lowestNonemptyLevel;

    bool pushToCache = isLowerThanLevelFront && (isCacheNotFull || isLowerThanCacheBack);

    if (pushToCache) {
        int i = _curSizeCache - 1;
        if (isCacheNotFull)
            i = _curSizeCache++;
        else
            push_queue(_cache[_curSizeCache], alphaToLevel(_cache[_curSizeCache].alpha, _a));

        for (; i >= 0 && newItem < _cache[i]; i--)
            _cache[i + 1] = _cache[i];

        _cache[i + 1] = newItem;
    } else
        push_queue(newItem, newItemLevel); // Push to the queue
}

template <class Pixel> void HHPQ<Pixel>::push_queue(const QItem<Pixel> &item, const ImgIdx &level) {
    if (level < _lowestNonemptyLevel)
        _lowestNonemptyLevel = level;

    if (level < _lowestUnsortedLevel) {
        _sortedLevels[level]->push(item);
    } else {
        ImgIdx cur = _unsortedLevelSizes[level]++;
        _unsortedLevels[level][cur] = item;
    }
}

template <class Pixel> void HHPQ<Pixel>::pop() {
    if (_curSizeCache) {
        for (int i = 0; i < _curSizeCache; i++)
            _cache[i] = _cache[i + 1];
        _curSizeCache--;
    } else {
        while (isFrontLevelEmptyAfterSort())
            _lowestNonemptyLevel++;
        _cache[0] = _sortedLevels[_lowestNonemptyLevel]->top();
        pop_queue();
    }
}

template <class Pixel> bool HHPQ<Pixel>::isFrontLevelEmptyAfterSort() {
    if (_lowestNonemptyLevel < _lowestUnsortedLevel)
        return _sortedLevels[_lowestNonemptyLevel]->empty();
    else {
        while (_lowestUnsortedLevel < _lowestNonemptyLevel) {
            _sortedLevels[_lowestUnsortedLevel] = new QuadHeapQueue<Pixel>(_levelMaxSizes[_lowestUnsortedLevel]);
            Free(_unsortedLevels[_lowestUnsortedLevel]);
            _unsortedLevels[_lowestUnsortedLevel] = 0;
            _lowestUnsortedLevel++;
        }
        _lowestUnsortedLevel++;

        _sortedLevels[_lowestNonemptyLevel] = new QuadHeapQueue<Pixel>(_levelMaxSizes[_lowestNonemptyLevel]);

        QItem<Pixel> *level = _unsortedLevels[_lowestNonemptyLevel];
        const ImgIdx cur = _unsortedLevelSizes[_lowestNonemptyLevel];
        QuadHeapQueue<Pixel> *pQ = _sortedLevels[_lowestNonemptyLevel];
        for (ImgIdx p = 0; p < cur; p++) {
            if (!_isVisited[level[p].index])
                pQ->push(level[p]);
        }
        Free(_unsortedLevels[_lowestNonemptyLevel]);
        _unsortedLevels[_lowestNonemptyLevel] = 0;
        return pQ->empty();
    }
}

template <class Pixel> void HHPQ<Pixel>::pop_queue() {
    _sortedLevels[_lowestNonemptyLevel]->pop();
    if (!_sortedLevels[_lowestNonemptyLevel]->get_cursize()) {
        do {
            _lowestNonemptyLevel++;
        } while (_lowestNonemptyLevel < _numLevels && isFrontLevelEmptyAfterSort());
    }
}

template class HHPQ<_uint8>;
template class HHPQ<_uint16>;
template class HHPQ<_uint32>;
template class HHPQ<_uint64>;