#pragma once

#include <QItem.hpp>
#include <defines.hpp>

class Bucket {
  public:
    void alloc(ImgIdx size_);
    const QItem &operator[](int idx) const { return _bucket[idx]; }
    QItem &operator[](int idx) { return _bucket[idx]; }
    ImgIdx size() const { return _size; }
    bool full() const { return _size == _maxSize; }
    bool empty() const { return _size == 0; }
    void push_back(const QItem &newItem) { _bucket[_size++] = QItem(newItem); }
    void emplace_back(const ImgIdx &index, const double &alpha) { _bucket[_size++] = QItem(index, alpha); }
    void sort() { std::sort(_bucket, _bucket + _size); }
    void print() const;

    ~Bucket();

  private:
    QItem *_bucket = nullptr;
    ImgIdx _size = 0;
    ImgIdx _maxSize = 0;
};

class BucketSort {
  public:
    BucketSort(ImgIdx *levelSizes, ImgIdx numLevels, ImgIdx totalSize);
    ~BucketSort();

    void push(const ImgIdx &level, const ImgIdx &index, const double &alpha);
    void push(const RankItem<double> &rankItem);

    void sort(ImgIdx *indexToRank, int32_t *rankToIndex);

    void print() const;
    void printPDF() const;

    ImgIdx size() const { return _size; }

  private:
    ImgIdx _size = 0;
    const ImgIdx _numLevels = 0;
    const ImgIdx _totalSize = 0;

    Bucket *_buckets = nullptr;

    // std::vector<RankItem<double>> _buckets[_MAX_LEVEL + 1];
};