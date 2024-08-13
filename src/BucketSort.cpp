#include "BucketSort.hpp"
#include <cmath>
#include <cstdio>

int BucketSort::alphaToIndex(float alpha) {
    // int index = (int)((float)NUM_BUCKETS * (1.0f - exp(-(alpha / AlPHA_MODEL_SIG))));
    int index = (int)((double)NUM_BUCKETS * (1.0 - exp(-((double)alpha / 4.0))));

    index = std::min<int>(NUM_BUCKETS - 1, index);
    return index;
}

// void BucketSort::sort(std::vector<Data> &dataAndIndex) {
    // const size_t N = dataAndIndex.size();
    // BucketArray bArray;

    // for (auto &item : dataAndIndex) {
    //     {
    //         int index = alphaToIndex(item.first);
    //         if (bArray[index].empty() == false && bArray[index].back().first != item.first)
    //             bArray[index].isSorted = false;
    //         bArray[index].data.push_back(std::move(item));
    //     }
    // }
    // dataAndIndex.clear();

    // size_t preSortedPop = 0;
    // for (const auto &[bIdx, bucket] : bArray) {
    //     if (bucket.size() == 0) {
    //         printf("bucket[%d].size() = %lu\n", bIdx, bucket.size());
    //         continue;
    //     }

    //     printf("bucket[%d].size() = %lu, isSorted = %d\n", bIdx, bucket.size(), (int)bucket.isSorted);
    //     if (bucket.isSorted)
    //         preSortedPop += bucket.size();
    // }
    // printf("Presorted population ratio = %f\n", (double)preSortedPop / (double)N);
// }
