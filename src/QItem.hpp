#pragma once

#include <QItem.hpp>
#include <cstdio>
#include <defines.h>
#include <limits>

template <class Pixel> struct QItem {
    ImgIdx index = -1;
    Pixel alpha = std::numeric_limits<Pixel>::max();

    QItem(ImgIdx index_, Pixel alpha_) : index(index_), alpha(alpha_) {}
    bool operator<(const QItem &other) const { return alpha < other.alpha; }
    bool operator<=(const QItem &other) const { return alpha <= other.alpha; }
    void print() { printf("(%d, %.2f) ", (int)index, (double)alpha); }
};