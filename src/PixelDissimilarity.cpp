#include <PixelDissimilarity.hpp>

template <class Pixel> double PixelDissimilarity<Pixel>::L1(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return std::abs((double)_image[index1] - (double)_image[index2]);

    double dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        double diff = std::abs((double)_image[index1 + offset] - (double)_image[index2 + offset]);
        dist += diff;
    }
    return dist / (double)_channels;
}

template <class Pixel> double PixelDissimilarity<Pixel>::L2(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return std::abs((double)_image[index1] - (double)_image[index2]);

    double dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        double diff = _image[index1 + offset] - _image[index2 + offset];
        dist += diff * diff;
    }
    return std::sqrt(dist / (double)_channels);
}

template <class Pixel> double PixelDissimilarity<Pixel>::LInfinity(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return std::abs((double)_image[index1] - (double)_image[index2]);

    double dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        double diff = std::abs((double)_image[index1 + offset] - (double)_image[index2 + offset]);
        dist = std::max(dist, diff);
    }
    return dist;
}

template class PixelDissimilarity<_uint8>;
template class PixelDissimilarity<_uint16>;
template class PixelDissimilarity<_uint32>;
template class PixelDissimilarity<_uint64>;