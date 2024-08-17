#pragma once

#include "AlphaTree.h"
#include <map>
#include <vector>

class BucketSort {
    using Bucket = std::vector<RankItem<float>>;

  public:
    static constexpr float AlPHA_MODEL_SIG = 1200.f;
    static constexpr size_t NUM_BUCKETS = 2048;

    std::vector<RankItem<float>> sort();
    std::vector<RankItem<float>> parallelSort(int numThreads = omp_get_max_threads());

    void insert(float alpha, ImgIdx idx);

    BucketSort();

  private:
    std::vector<Bucket> buckets;
    Bucket bucket0;

    int alphaToIndex(float alpha);
};