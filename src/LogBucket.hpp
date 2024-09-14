#pragma once

#include <QItem.hpp>
#include <defines.hpp>

class LogBucket {
  public:
    static ImgIdx alphaToLevel(const double &alpha);

    void push(ImgIdx index, double alpha);

    void sort(ImgIdx *indexToRank, int32_t *rankToIndex, RankItem<double> *rankitem, ImgIdx numEdges);

  private:
    static constexpr double _A = 0;
    static constexpr ImgIdx _MAX_LEVEL = (size_t)(_A * 64.0);

    std::vector<RankItem<double>> _buckets[_MAX_LEVEL + 1];
};