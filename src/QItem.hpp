#pragma once

#include <QItem.hpp>
#include <defines.hpp>

struct QItem {
    ImgIdx index = -1;
    double alpha = std::numeric_limits<double>::infinity();

    // QItem(ImgIdx index_, Pixel alpha_) : index(index_), alpha(alpha_) {}
    QItem(ImgIdx index_, double alpha_) : index(index_), alpha(alpha_) {}
    bool operator<(const QItem &other) const { return alpha < other.alpha; }
    // bool operator<=(const QItem &other) const { return alpha <= other.alpha; }
    void print() { printf("(%d, %.2f) ", (int)index, (double)alpha); }

    static constexpr uint8_t EDGE_STANDBY = 0;
    static constexpr uint8_t EDGE_ENQUEUED = 1;
    static constexpr uint8_t EDGE_DEQUEUED = 2;
    static constexpr uint8_t EDGE_CONNECTED = 3;
    static constexpr uint8_t EDGE_REDUNDANT = 4;
    static constexpr uint8_t EDGE_ESSENTIAL = 5;
};

template <class Pixel> class RankItem {
  public:
    ImgIdx dimgidx = -1;
    Pixel alpha = 0;
    uint16_t bucketIndex = 0;

    ImgIdx get_pidx0(ImgIdx _connectivity = 4);
    ImgIdx get_pidx1(ImgIdx _width, ImgIdx _connectivity = 4);

    RankItem(ImgIdx dimgidx_, Pixel alpha_, uint16_t bucketIndex_ = 0u)
        : dimgidx(dimgidx_), alpha(alpha_), bucketIndex(bucketIndex_) {}
    bool operator<(const RankItem &other) const { return alpha < other.alpha; }
    void print() const { printf("(%d, %.2f, %d) ", (int)dimgidx, (double)alpha, (int)bucketIndex); }
};