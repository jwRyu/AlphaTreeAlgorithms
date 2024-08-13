#pragma once

#include "AlphaTree.h"
#include <map>
#include <vector>

struct Bucket {
    std::vector<RankItem<float>> data;
    bool isSorted = true;
    size_t size() const { return data.size(); }
    bool empty() const { return data.empty(); }
    RankItem<float> &back() { return data.back(); }
};
using BucketArray = std::vector<Bucket>;

class BucketSort {
  private:
    static constexpr float AlPHA_MODEL_SIG = 2.0f;
    static constexpr size_t NUM_BUCKETS = 2048;

    static int alphaToIndex(float alpha);

  public:
    // static void sort(std::vector<Data> &dataAndIndex);
};