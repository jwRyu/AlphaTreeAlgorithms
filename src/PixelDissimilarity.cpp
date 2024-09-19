#include <PixelDissimilarity.hpp>

template <class Pixel> float PixelDissimilarity<Pixel>::L1(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return std::abs((float)_image[index1] - (float)_image[index2]);

    float dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        float diff = std::abs((float)_image[index1 + offset] - (float)_image[index2 + offset]);
        dist += diff;
    }
    return dist / (float)_channels;
}

template <class Pixel> float PixelDissimilarity<Pixel>::L2(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return std::abs((float)_image[index1] - (float)_image[index2]);

    float dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        float diff = _image[index1 + offset] - _image[index2 + offset];
        dist += diff * diff;
    }
    return std::sqrt(dist / (float)_channels);
}

template <class Pixel> float PixelDissimilarity<Pixel>::LInfinity(ImgIdx index1, ImgIdx index2) const {
    if (_channels == 1)
        return std::abs((float)_image[index1] - (float)_image[index2]);

    float dist = 0.0;
    for (ImgIdx ch = 0; ch < _channels; ch++) {
        ImgIdx offset = ch * _imageSize;
        float diff = std::abs((float)_image[index1 + offset] - (float)_image[index2 + offset]);
        dist = std::max(dist, diff);
    }
    return dist;
}

template class PixelDissimilarity<uint8_t>;
template class PixelDissimilarity<uint16_t>;
template class PixelDissimilarity<uint32_t>;
template class PixelDissimilarity<uint64_t>;