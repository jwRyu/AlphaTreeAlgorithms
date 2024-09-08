#include <HHPQ.hpp>
#include <allocator.h>
#include <cmath>
#include <cstring>

ImgIdx HHPQ::alphaToLevel(const double &alpha, const double &a) { return (ImgIdx)(a * log2(1.0 + alpha)); }

ImgIdx HHPQ::alphaToLevel(const double &alpha) const {
    return std::min<ImgIdx>((ImgIdx)(_a * log2(1.0 + alpha)), _numLevels - 1);
}

HHPQ::HHPQ(const ImgIdx *dhist, ImgIdx numLevels_, ImgIdx size, const uint8_t *isVisited_, double a_, int cacheSize,
           double r)
    : _numLevels(numLevels_), _a(a_), _lowestNonemptyLevel(numLevels_), _curSizeCache(0), _maxSizeCache(cacheSize),
      _isVisited(isVisited_), _size(0) {
    initHQ(dhist, size, r);
}

HHPQ::~HHPQ() {
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

void HHPQ::initHQ(const ImgIdx *dhist, ImgIdx size, double r) {
    _cache = (QItem *)Malloc(_maxSizeCache * sizeof(QItem));
    ImgIdx cumsum = 0;
    _levelMaxSizes = (ImgIdx *)Calloc((size_t)_numLevels * sizeof(ImgIdx));
    memcpy(_levelMaxSizes, dhist, (size_t)_numLevels * sizeof(ImgIdx));
    if (r >= 1) {
        _lowestUnsortedLevelInitial = _lowestUnsortedLevel = _numLevels;
        _sortedLevels = (QuadHeapQueue **)Calloc(_numLevels * sizeof(QuadHeapQueue *));
        for (int level = 0; level < _lowestUnsortedLevelInitial; level++)
            _sortedLevels[level] = new QuadHeapQueue(_levelMaxSizes[level]);
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

        _sortedLevels = (QuadHeapQueue **)Calloc(_numLevels * sizeof(QuadHeapQueue *));
        for (int level = 0; level < _lowestUnsortedLevelInitial; level++)
            _sortedLevels[level] = new QuadHeapQueue(_levelMaxSizes[level]);

        _unsortedLevels = (QItem **)Calloc((_numLevels - _lowestUnsortedLevelInitial) * sizeof(QItem *));
        _unsortedLevels -= _lowestUnsortedLevelInitial;
        for (int level = _lowestUnsortedLevelInitial; level < _numLevels; level++)
            _unsortedLevels[level] = (QItem *)Malloc(_levelMaxSizes[level] * sizeof(QItem));
    }
}

void HHPQ::print() {
    printf("---------- HHPQ::print START -------------\n");
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
    printf("---------- HHPQ::print END -------------\n");
}

void HHPQ::push(const ImgIdx &idx, const double &alpha) {
    _size++;
    const QItem newItem(idx, alpha);

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

void HHPQ::push_queue(const QItem &item, const ImgIdx &level) {
    if (level < _lowestNonemptyLevel)
        _lowestNonemptyLevel = level;

    if (level < _lowestUnsortedLevel) {
        _sortedLevels[level]->push(item);
    } else {
        ImgIdx cur = _unsortedLevelSizes[level]++;
        _unsortedLevels[level][cur] = item;
    }
}

void HHPQ::pop() {
    if (_size == 0) {
        return;
    }
    _size--;
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

bool HHPQ::isFrontLevelEmptyAfterSort() {
    if (_lowestNonemptyLevel < _lowestUnsortedLevel)
        return _sortedLevels[_lowestNonemptyLevel]->empty();
    else {
        while (_lowestUnsortedLevel < _lowestNonemptyLevel) {
            _sortedLevels[_lowestUnsortedLevel] = new QuadHeapQueue(_levelMaxSizes[_lowestUnsortedLevel]);
            Free(_unsortedLevels[_lowestUnsortedLevel]);
            _unsortedLevels[_lowestUnsortedLevel] = 0;
            _lowestUnsortedLevel++;
        }
        _lowestUnsortedLevel++;

        _sortedLevels[_lowestNonemptyLevel] = new QuadHeapQueue(_levelMaxSizes[_lowestNonemptyLevel]);

        QItem *level = _unsortedLevels[_lowestNonemptyLevel];
        const ImgIdx levelSize = _unsortedLevelSizes[_lowestNonemptyLevel];
        QuadHeapQueue *pQ = _sortedLevels[_lowestNonemptyLevel];
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

void HHPQ::pop_queue() {
    _sortedLevels[_lowestNonemptyLevel]->pop();
    if (!_sortedLevels[_lowestNonemptyLevel]->get_cursize()) {
        do {
            _lowestNonemptyLevel++;
        } while (_lowestNonemptyLevel < _numLevels && isFrontLevelEmptyAfterSort());
    }
}