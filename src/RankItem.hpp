#pragma once
#include <defines.hpp>

template <class Pixel> class RankItem {
  public:
    ImgIdx dimgidx = -1;
    Pixel alpha = 0;

    ImgIdx get_pidx0(ImgIdx _connectivity = 4);
    ImgIdx get_pidx1(ImgIdx _width, ImgIdx _connectivity = 4);

    RankItem(ImgIdx dimgidx_, Pixel alpha_, uint16_t bucketIndex_ = 0u) : dimgidx(dimgidx_), alpha(alpha_) {}
    bool operator<(const RankItem &other) const { return alpha < other.alpha; }
    void print() const { printf("(%d, %.2f) ", (int)dimgidx, (double)alpha); }
};