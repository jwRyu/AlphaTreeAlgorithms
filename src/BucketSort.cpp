#include <BucketSort.hpp>
#include <defines.hpp>

#define DEBUG 1

void Bucket::alloc(ImgIdx size) {
    _maxSize = size;
    _size = 0;
    _bucket = _maxSize > 0 ? (QItem *)Calloc((size_t)_maxSize * sizeof(QItem)) : nullptr;
}

Bucket::~Bucket() {
    if (_bucket)
        Free(_bucket);
}

// void Bucket::push_back(const QItem &newItem) {
// #if DEBUG
//     if (full())
//         printf("Bucket::push_back() ERROR: Push to full bucket\n");
//     if (_bucket == nullptr)
//         printf("Bucket::push_back() ERROR: Push to unallocated bucket\n");
// #endif
//     _bucket[_size++] = QItem(newItem);
// }

// void Bucket::emplace_back(const ImgIdx &index, const double &alpha) {
// #if DEBUG
//     if (full())
//         printf("Bucket::emplace_back() ERROR: Emplace to full bucket\n");
//     if (_bucket == nullptr)
//         printf("Bucket::emplace_back() ERROR: Emplace to unallocated bucket\n");
// #endif
//     _bucket[_size++] = QItem(index, alpha);
// }

void Bucket::print() const {
    for (int i = 0; i < _size; i++)
        _bucket[i].print();
    printf("\n");
}

BucketSort::BucketSort(ImgIdx *levelSizes, ImgIdx numLevels, ImgIdx totalSize)
    : _numLevels(numLevels), _totalSize(totalSize) {
    ImgIdx sizeCheck = 0;
    _buckets = (Bucket *)Malloc(_numLevels * sizeof(Bucket));
    for (ImgIdx level = 0; level < _numLevels; level++) {
        _buckets[level].alloc(levelSizes[level]);
        sizeCheck += levelSizes[level];
    }
    if (sizeCheck != _totalSize) {
        printf("BucketSort::BucketSort() Warning: size mistmatch\n");
    }
}

BucketSort::~BucketSort() {
    if (_buckets)
        Free(_buckets);
}

void BucketSort::print() const {
    printf("---------BucketSort::print()---------\n");
    printf("numBuckets = %d, sizeTotal = %d\n", _numLevels, _totalSize);
    for (int level = 0; level < _numLevels; level++) {
        const auto &bucket = _buckets[level];
        if (bucket.empty())
            continue;
        printf("bucket[%d].size() = %d: ", level, bucket.size());
        bucket.print();
    }
    printf("======================================\n");
}

void BucketSort::printPDF() const {
    printf("---------BucketSort::printPDF()---------\n");
    printf("numBuckets = %d, sizeTotal = %d\n", _numLevels, _totalSize);
    for (int level = 0; level < _numLevels; level++) {
        const auto &bucket = _buckets[level];
        if (bucket.empty())
            continue;
        printf("bucket[%d].size() = %d\n", level, bucket.size());
    }
    printf("=========================================\n");
}

void BucketSort::push(const RankItem<double> &rankItem) {
    _size++;
    _buckets[rankItem.bucketIndex].emplace_back(rankItem.dimgidx, rankItem.alpha);
}

void BucketSort::push(const ImgIdx &level, const ImgIdx &index, const double &alpha) {
    _size++;
    _buckets[level].emplace_back(index, alpha);
}

// implement logbucket sort and replace teeninga sort
void BucketSort::sort(ImgIdx *indexToRank, int32_t *rankToIndex) {
    ImgIdx *startIndices = (ImgIdx *)Calloc((_numLevels + 1) * sizeof(ImgIdx));
    startIndices[0] = 0;
    for (ImgIdx level = 0; level < _numLevels; level++)
        startIndices[level + 1] = startIndices[level] + _buckets[level].size();

#if DEBUG
    assert(_totalSize == startIndices[_numLevels]);
#endif

#pragma omp parallel for
    for (ImgIdx level = 0; level < _numLevels; level++) {
        auto &bucket = _buckets[level];
        if (bucket.empty())
            continue;
        ImgIdx rank = startIndices[level];
        // const auto startIndex = startIndices[level];
        bucket.sort();
        for (int i = 0; i < (int)bucket.size(); i++, rank++) {
            rankToIndex[rank] = bucket[i].index;
            indexToRank[bucket[i].index] = rank;
        }
    }
    Free(startIndices);
}
