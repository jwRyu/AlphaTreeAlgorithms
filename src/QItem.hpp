#pragma once

#include <QItem.hpp>
#include <cstdio>
#include <defines.h>
#include <limits>

template <class Pixel> struct QItem {
    ImgIdx index = -1;
    Pixel alpha = std::numeric_limits<Pixel>::max();
    ImgIdx edgeIdx = -1;

    // QItem(ImgIdx index_, Pixel alpha_) : index(index_), alpha(alpha_) {}
    QItem(ImgIdx index_, Pixel alpha_, ImgIdx edgeIdx_ = -1) : index(index_), alpha(alpha_), edgeIdx(edgeIdx_) {}
    bool operator<(const QItem &other) const { return alpha < other.alpha; }
    bool operator<=(const QItem &other) const { return alpha <= other.alpha; }
    void print() { printf("(%d, %.2f) ", (int)index, (double)alpha); }
};