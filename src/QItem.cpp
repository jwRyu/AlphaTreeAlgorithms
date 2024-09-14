#include <QItem.hpp>

template <class Pixel> ImgIdx RankItem<Pixel>::get_pidx0(ImgIdx _connectivity) {
    if (_connectivity == 4)
        return (this->dimgidx >> 1);
    else if (_connectivity == 8)
        return (this->dimgidx >> 2);
    else {
        return -1;
    }
}

template <class Pixel> ImgIdx RankItem<Pixel>::get_pidx1(ImgIdx _width, ImgIdx _connectivity) {
    if (_connectivity == 4)
        return (this->dimgidx >> 1) + _width + (1 - _width) * (this->dimgidx & 1);
    else if (_connectivity == 8) {
        ImgIdx neighboridx = (this->dimgidx & 2);
        return (this->dimgidx >> 2) + _width * ((ImgIdx)(neighboridx < 2) - (ImgIdx)(neighboridx == 3)) +
               (ImgIdx)(neighboridx > 0);
    } else {
        return -1;
    }
}

template class RankItem<uint8_t>;
template class RankItem<uint16_t>;
template class RankItem<uint32_t>;
template class RankItem<uint64_t>;
template class RankItem<double>;
