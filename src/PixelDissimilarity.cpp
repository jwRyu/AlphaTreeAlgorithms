#include <PixelDissimilarity.hpp>

template <class Pixel> float PixelDissimilarity<Pixel>::L1(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return (float)std::abs((double)_image[index1] - (double)_image[index2]);

    double dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        double diff = std::abs((double)_image[index1 + offset] - (double)_image[index2 + offset]);
        dist += diff;
    }
    return (float)(dist / (double)_channels);
}

template <class Pixel> float PixelDissimilarity<Pixel>::L2(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return (float)std::abs((double)_image[index1] - (double)_image[index2]);

    double dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        double diff = _image[index1 + offset] - _image[index2 + offset];
        dist += diff * diff;
    }
    return (float)std::sqrt(dist / (double)_channels);
}

template <class Pixel> float PixelDissimilarity<Pixel>::LInfinity(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return (float)std::abs((double)_image[index1] - (double)_image[index2]);

    double dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        double diff = std::abs((double)_image[index1 + offset] - (double)_image[index2 + offset]);
        dist = std::max(dist, diff);
    }
    return (float)dist;
}

template class PixelDissimilarity<uint8_t>;
template class PixelDissimilarity<uint16_t>;
template class PixelDissimilarity<uint32_t>;
template class PixelDissimilarity<uint64_t>;