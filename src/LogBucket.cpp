#include <LogBucket.hpp>

ImgIdx LogBucket::alphaToLevel(const double &alpha) {
    size_t index = (ImgIdx)(_A * log2(1.0 + alpha));
    return CLIP(index, 0, _MAX_LEVEL);
}

void LogBucket::push(ImgIdx index, double alpha) {
    ImgIdx level = alphaToLevel(alpha);
    auto &bucket = _buckets[level];
    // #pragma omp atomic
    bucket.emplace_back(index, alpha);
}

// implement logbucket sort and replace teeninga sort
void LogBucket::sort(ImgIdx *indexToRank, int32_t *rankToIndex, RankItem<double> *rankitem, ImgIdx numEdges) {
    ImgIdx startIndices[_MAX_LEVEL + 1];
    startIndices[0] = 0;
    for (int level = 0; level < _MAX_LEVEL; level++) {
        startIndices[level + 1] = startIndices[level] + _buckets[level].size();
    }
    assert(numEdges == startIndices[_MAX_LEVEL]);

    ImgIdx rank = 0;
    for (int level = 0; level <= _MAX_LEVEL; level++) {
        auto &bucket = _buckets[level];
        const auto startIndex = startIndices[level];
        std::sort(bucket.begin(), bucket.end());
        for (int i = 0; i < (int)bucket.size(); i++, rank++) {
            rankToIndex[rank] = bucket[i].dimgidx;
            indexToRank[bucket[i].dimgidx] = rank;
            rankitem[startIndex + i] = std::move(bucket[i]);
        }
    }
}
