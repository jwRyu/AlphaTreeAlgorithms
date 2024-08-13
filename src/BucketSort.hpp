#pragma once

#include <map>
#include <vector>

using Data = std::pair<float, int>;
struct Bucket {
    std::vector<Data> data;
    bool isSorted = true;
    size_t size() const { return data.size(); }
    bool empty() const { return data.empty(); }
    Data &back() { return data.back(); }
};
using BucketArray = std::map<int, Bucket>;

class BucketSort {
  private:
    static constexpr float AlPHA_MODEL_SIG = 2.0f;
    static constexpr size_t NUM_BUCKETS = 2048;

    static int alphaToIndex(float alpha);

  public:
    static void sort(std::vector<Data> &dataAndIndex);
};