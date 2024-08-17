#include "BucketSort.hpp"
#include <cmath>
#include <cstdio>

BucketSort::BucketSort() { buckets.resize(NUM_BUCKETS); }

int BucketSort::alphaToIndex(float alpha) {
    int index = (int)((float)NUM_BUCKETS * (1.0f - exp(-(alpha / (AlPHA_MODEL_SIG * AlPHA_MODEL_SIG)))));

    index = std::min<int>(NUM_BUCKETS - 1, index);
    return index;
}

void BucketSort::insert(float alpha, ImgIdx idx) {
    if (alpha <= 0) {
        // printf("PIEP0\n");
        bucket0.emplace_back(alpha, idx);
        return;
    }

    int index = alphaToIndex(alpha);
    // printf("PIEP %d\n", index);
    buckets[index].emplace_back(alpha, idx);
    // printf("PIEPPIEP\n");
}

std::vector<RankItem<float>> BucketSort::sort() {
    std::vector<RankItem<float>> &sorted = bucket0;

    for (auto &bucket : buckets) {
        if (bucket.empty())
            continue;
        std::sort(bucket.begin(), bucket.end(),
                  [](const RankItem<float> &a, const RankItem<float> &b) { return a.alpha < b.alpha; });
        sorted.insert(sorted.end(), bucket.begin(), bucket.end());
    }

    return sorted;
}

std::vector<RankItem<float>> BucketSort::parallelSort(int numThreads) {
    numThreads = 1;
    omp_set_num_threads(numThreads);

    std::vector<RankItem<float>> &sorted = bucket0;
    // Create a vector of vectors to store sorted buckets for each thread
    std::vector<std::vector<RankItem<float>>> local_sorted_buckets(numThreads);

#pragma omp parallel for
    for (size_t i = 0; i < buckets.size(); ++i) {

        int thread_id = omp_get_thread_num();

        auto &bucket = buckets[i];
        if (bucket.empty())
            continue;

        std::sort(bucket.begin(), bucket.end(),
                  [](const RankItem<float> &a, const RankItem<float> &b) { return a.alpha < b.alpha; });

        // printf("thread_id = %d on bucket i = %d\n", thread_id, (int)i);
    }

    for (auto &bucket : buckets)
        sorted.insert(sorted.end(), bucket.begin(), bucket.end());

    return sorted;
}