#include <HHPQ.hpp>
#include <allocator.h>
#include <cmath>

template <class Pixel> ImgIdx HHPQ<Pixel>::alphaToLevel(const double &alpha) const {
    return (ImgIdx)(a * log2(1.0 + alpha));
}

template <class Pixel>
void HHPQ<Pixel>::initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int cacheSize, int connectivity,
                         double r) {
    maxSize = size;

    _cache = (QItem<Pixel> *)Malloc((cacheSize) * sizeof(QItem<Pixel>));
    maxSizeCache = cacheSize - 1;
    curSizeCache = -1;

    this->_numLevels = numlevels_in;
    this->a = a_in;
    this->_lowestNonemptyLevel = _numLevels;

    ImgIdx cumsum = 0;
    _levelMaxSizes = dhist; // do not free dhist outside
    if (r >= 1) {
        thr_hqueue = _lowestUnsortedLevel = _numLevels;
        _sortedLevels = (QuadHeapQueue<Pixel> **)Calloc(_numLevels * sizeof(QuadHeapQueue<Pixel> *));
        for (int level = 0; level < thr_hqueue; level++)
            _sortedLevels[level] = new QuadHeapQueue<Pixel>(_levelMaxSizes[level]);
        _unsortedLevels = 0;
        _unsortedLevelSizes = 0;
    } else {
        _unsortedLevelSizes = (ImgIdx *)Calloc(_numLevels * sizeof(ImgIdx));
        ImgIdx thr_nonredundantnodes = (ImgIdx)(size * r);
        for (int level = 0; level < _numLevels; level++) {
            cumsum += _levelMaxSizes[level];
            if (cumsum > thr_nonredundantnodes) {
                thr_hqueue = _lowestUnsortedLevel = level;
                break;
            }
        }

        _sortedLevels = (QuadHeapQueue<Pixel> **)Calloc(_numLevels * sizeof(QuadHeapQueue<Pixel> *));
        for (int level = 0; level < thr_hqueue; level++)
            _sortedLevels[level] = new QuadHeapQueue<Pixel>(_levelMaxSizes[level]);

        _unsortedLevels = (QItem<Pixel> **)Calloc((_numLevels - thr_hqueue) * sizeof(QItem<Pixel> *));
        _unsortedLevels -= thr_hqueue;
        for (int level = thr_hqueue; level < _numLevels; level++)
            _unsortedLevels[level] = (QItem<Pixel> *)Malloc(_levelMaxSizes[level] * sizeof(QItem<Pixel>));
    }
}

template <class Pixel> void QuadHeapQueue<Pixel>::print() {
    for (int i = 0; i < cursize; i++)
        arr[i + 1].print();
    printf("\n");
}

template <class Pixel> void HHPQ<Pixel>::print() {
    printf("---------- HHPQ<Pixel>::print START -------------\n");
    printf("Cache[%d / %d]: ", curSizeCache + 1, maxSizeCache + 1);
    size_t size = curSizeCache;
    for (int i = -1; i < curSizeCache; i++)
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

template <class Pixel>
HHPQ<Pixel>::HHPQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int cacheSize, ImgIdx connectivity,
                  double r) {
    initHQ(dhist, numlevels_in, size, a_in, cacheSize, (int)connectivity, r);
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
        for (int level = thr_hqueue; level < _numLevels; level++)
            if (_unsortedLevels[level])
                Free(_unsortedLevels[level]);
        Free(_unsortedLevels + thr_hqueue);
    }
}

template <class Pixel> void HHPQ<Pixel>::push_1stitem(ImgIdx idx) {
    _cache[0].index = idx;
    _cache[0].alpha = std::numeric_limits<Pixel>::max();
    curSizeCache++;
}

template <class Pixel> void HHPQ<Pixel>::end_pushes(_uint8 *isVisited) {
    if (emptytop)
        pop(isVisited);
}

template <class Pixel> void HHPQ<Pixel>::push(const ImgIdx &idx, const Pixel &alpha) {
    const QItem<Pixel> newItem(idx, alpha);

    // printf("Pushing %d at %.2f\n", idx, (double)alpha);
    if (emptytop && newItem < _cache[0]) {
        emptytop = 0;
        _cache[0] = newItem;
        return;
    }

    const ImgIdx newItemLevel = alphaToLevel(newItem.alpha);
    const bool isCacheNotFull = curSizeCache < maxSizeCache;
    const bool isLowerThanCacheBack = newItem < cacheBack();
    const bool isFrontLevelSorted = _lowestNonemptyLevel < _lowestUnsortedLevel;
    const bool isLowerThanLevelFront =
        isFrontLevelSorted ? newItem < _sortedLevels[_lowestNonemptyLevel]->top() : newItemLevel < _lowestNonemptyLevel;

    bool pushToCache = isLowerThanLevelFront && (isCacheNotFull || isLowerThanCacheBack);

    if (pushToCache) {
        if (curSizeCache < maxSizeCache) { // spare room in the _cache
            int i = curSizeCache++;
            for (; i >= 0 && newItem < _cache[i]; i--)
                _cache[i + 1] = _cache[i];

            _cache[i + 1] = newItem;
        } else if (alpha < _cache[curSizeCache].alpha) { // push to the full _cache
            push_queue(_cache[curSizeCache], alphaToLevel(_cache[curSizeCache].alpha));
            int i = curSizeCache - 1;
            for (; i >= 0 && alpha < _cache[i].alpha; i--)
                _cache[i + 1] = _cache[i];

            _cache[i + 1] = newItem;
        } else {
            push_queue(newItem, newItemLevel); // push to the queue
        }
    } else
        push_queue(newItem, newItemLevel); // push to the queue
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

template <class Pixel> ImgIdx HHPQ<Pixel>::pop(_uint8 *isVisited) {
    ImgIdx ret = top();
    if (curSizeCache) {
        for (int i = 0; i < curSizeCache; i++)
            _cache[i] = _cache[i + 1];
        curSizeCache--;
    } else {
        while (!check_queue_level(isVisited))
            _lowestNonemptyLevel++;
        _cache[0] = _sortedLevels[_lowestNonemptyLevel]->top();
        pop_queue(isVisited);
    }
    return ret;
}

template <class Pixel> int HHPQ<Pixel>::check_queue_level(_uint8 *isVisited) {
    if (_lowestNonemptyLevel < _lowestUnsortedLevel)
        return _sortedLevels[_lowestNonemptyLevel]->get_cursize();
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
            if (!isVisited[level[p].index])
                pQ->push(level[p]);
        }
        Free(_unsortedLevels[_lowestNonemptyLevel]);
        _unsortedLevels[_lowestNonemptyLevel] = 0;
        return pQ->get_cursize();
    }
}

template <class Pixel> void HHPQ<Pixel>::pop_queue(_uint8 *isVisited) {
    _sortedLevels[_lowestNonemptyLevel]->pop();
    if (!_sortedLevels[_lowestNonemptyLevel]->get_cursize()) {
        do {
            _lowestNonemptyLevel++;
        } while (_lowestNonemptyLevel < _numLevels && !check_queue_level(isVisited));
    }
}

template class HHPQ<_uint8>;
template class HHPQ<_uint16>;
template class HHPQ<_uint32>;
template class HHPQ<_uint64>;