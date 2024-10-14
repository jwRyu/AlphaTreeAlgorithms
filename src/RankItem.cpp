#include <RankItem.hpp>

template <class Pixel> ImgIdx RankItem<Pixel>::get_pidx0(ImgIdx _connectivity) const {
    if (_connectivity == 4)
        return (this->dimgidx >> 1);
    else if (_connectivity == 8)
        return (this->dimgidx >> 2);
    else {
        return -1;
    }
}

template <class Pixel> ImgIdx RankItem<Pixel>::get_pidx1(ImgIdx _width, ImgIdx _connectivity) const {
    if (_connectivity == 4)
        return (this->dimgidx >> 1) + _width + (1 - _width) * (this->dimgidx & 1);
    else if (_connectivity == 8) {
        const ImgIdx offsets[] = {_width, _width + 1, 1, -_width + 1};
        const ImgIdx nidx = (this->dimgidx & 0b11);
        return (this->dimgidx >> 2) + offsets[nidx];
    } else {
        return -1;
    }
}

template class RankItem<uint8_t>;
template class RankItem<uint16_t>;
template class RankItem<uint32_t>;
template class RankItem<uint64_t>;
template class RankItem<float>;
template class RankItem<double>;
