#include <HHPQ.hpp>
#include <allocator.h>
#include <cmath>
#include <cstring>

template <class Pixel> ImgIdx HHPQ<Pixel>::alphaToLevel(const double &alpha, const double &a) {
    return (ImgIdx)(a * log2(1.0 + alpha));
}

template <class Pixel> ImgIdx HHPQ<Pixel>::alphaToLevel(const double &alpha) const {
    return std::min<ImgIdx>((ImgIdx)(_a * log2(1.0 + alpha)), _numLevels - 1);
}

template <class Pixel>
HHPQ<Pixel>::HHPQ(ImgIdx *dhist, ImgIdx numLevels_, ImgIdx size, _uint8 *isVisited_, double a_, int cacheSize, double r)
    : _numLevels(numLevels_), _a(a_), _lowestNonemptyLevel(numLevels_), _curSizeCache(0), _maxSizeCache(cacheSize),
      _isVisited(isVisited_), _size(0) {
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
        for (int level = _lowestUnsortedLevelInitial; level < _numLevels; level++)
            if (_unsortedLevels[level])
                Free(_unsortedLevels[level]);
        Free(_unsortedLevels + _lowestUnsortedLevelInitial);
    }
}

template <class Pixel> void HHPQ<Pixel>::initHQ(ImgIdx *dhist, ImgIdx size, double r) {
    _cache = (QItem<Pixel> *)Malloc(_maxSizeCache * sizeof(QItem<Pixel>));
    ImgIdx cumsum = 0;
    _levelMaxSizes = (ImgIdx *)Calloc((size_t)_numLevels * sizeof(ImgIdx));
    memcpy(_levelMaxSizes, dhist, (size_t)_numLevels * sizeof(ImgIdx));
    if (r >= 1) {
        _lowestUnsortedLevelInitial = _lowestUnsortedLevel = _numLevels;
        _sortedLevels = (QuadHeapQueue<Pixel> **)Calloc(_numLevels * sizeof(QuadHeapQueue<Pixel> *));
        for (int level = 0; level < _lowestUnsortedLevelInitial; level++)
            _sortedLevels[level] = new QuadHeapQueue<Pixel>(_levelMaxSizes[level]);
        _unsortedLevels = 0;
        _unsortedLevelSizes = 0;
    } else {
        _unsortedLevelSizes = (ImgIdx *)Calloc(_numLevels * sizeof(ImgIdx));
        ImgIdx thr_nonredundantnodes = (ImgIdx)(size * r);
        for (int level = 0; level < _numLevels; level++) {
            cumsum += _levelMaxSizes[level];
            if (cumsum > thr_nonredundantnodes) {
                _lowestUnsortedLevelInitial = _lowestUnsortedLevel = level;
                break;
            }
        }

        _sortedLevels = (QuadHeapQueue<Pixel> **)Calloc(_numLevels * sizeof(QuadHeapQueue<Pixel> *));
        for (int level = 0; level < _lowestUnsortedLevelInitial; level++)
            _sortedLevels[level] = new QuadHeapQueue<Pixel>(_levelMaxSizes[level]);

        _unsortedLevels = (QItem<Pixel> **)Calloc((_numLevels - _lowestUnsortedLevelInitial) * sizeof(QItem<Pixel> *));
        _unsortedLevels -= _lowestUnsortedLevelInitial;
        for (int level = _lowestUnsortedLevelInitial; level < _numLevels; level++)
            _unsortedLevels[level] = (QItem<Pixel> *)Malloc(_levelMaxSizes[level] * sizeof(QItem<Pixel>));
    }
}

template <class Pixel> void HHPQ<Pixel>::print() {
    printf("---------- HHPQ<Pixel>::print START -------------\n");
    printf("Cache[%d / %d]: ", _curSizeCache, _maxSizeCache);
    size_t size = _curSizeCache - _emptyTop;
    for (int i = 0; i < _curSizeCache; i++) {
        if (_emptyTop && i == 0)
            continue;
        _cache[i].print();
    }
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

    printf("size = %d _size = %d\n", (int)size, (int)_size - _emptyTop);
    printf("---------- HHPQ<Pixel>::print END -------------\n");
    // std::getchar();
}

// template <class Pixel> void HHPQ<Pixel>::startPushes() { _emptyTop = true; }

// template <class Pixel> void HHPQ<Pixel>::endPushes() {
//     if (_emptyTop)
//         pop();
//     else
//         _size--;
// }

template <class Pixel> void HHPQ<Pixel>::push(const ImgIdx &idx, const Pixel &alpha) {
    _size++;
    // printf("***HHPQ*** Pushing %d at %.2f (_size = %d)\n", idx, (double)alpha, (int)_size);

    const QItem<Pixel> newItem(idx, alpha);

    if (_curSizeCache == 0) {
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
        // indexForNewItem is initalized to index of empty spot at the back of the cache
        int indexForNewItem = _curSizeCache - 1;
        if (isCacheNotFull) {
            indexForNewItem = _curSizeCache++;
        } else {
            // Empty last spot if it's full (isLowerThanCacheBack ensures newItem < _cache[_curSizeCache - 1])
            push_queue(_cache[_curSizeCache - 1], alphaToLevel(_cache[_curSizeCache - 1].alpha, _a));
        }

        // Shift indexForNewItem to front if the item at the front has higher alpha
        for (; indexForNewItem > 0 && newItem < _cache[indexForNewItem - 1]; indexForNewItem--)
            _cache[indexForNewItem] = _cache[indexForNewItem - 1];

        _cache[indexForNewItem] = newItem;
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
    if (_size == 0) {
        // printf("***HHPQ*** POPPING EMPTY QUEUE\n");
        return;
    }
    _size--;
    // printf("***HHPQ*** pop %d at %.2f (_size = %d)\n", _cache[0].index, (double)_cache[0].alpha, (int)_size);
    bool popLevel = true;
    if (_curSizeCache > 0) {
        // Overwrite _cache[0] and shift cache to the front
        for (int i = 0; i < _curSizeCache - 1; i++)
            _cache[i] = _cache[i + 1];
        _curSizeCache--;
        if (_curSizeCache > 0 || _size == 0)
            popLevel = false;
    }
    if (popLevel) {
        while (_lowestNonemptyLevel < _numLevels && isFrontLevelEmptyAfterSort())
            _lowestNonemptyLevel++;
        if (_lowestNonemptyLevel < _numLevels) {
            _cache[0] = _sortedLevels[_lowestNonemptyLevel]->top();
            _curSizeCache++;
            pop_queue();
        }
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
        const ImgIdx levelSize = _unsortedLevelSizes[_lowestNonemptyLevel];
        QuadHeapQueue<Pixel> *pQ = _sortedLevels[_lowestNonemptyLevel];
        for (ImgIdx p = 0; p < levelSize; p++) {
            if (!_isVisited[level[p].index])
                pQ->push(level[p]);
        }
        Free(_unsortedLevels[_lowestNonemptyLevel]);
        _unsortedLevels[_lowestNonemptyLevel] = 0;
        _size -= (levelSize - pQ->size());
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