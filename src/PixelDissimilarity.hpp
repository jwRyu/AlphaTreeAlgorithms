#pragma once

#include <defines.h>

template <class Pixel> class PixelDissimilarity {
  public:
    // Typedef for a pointer to a dissimilarity function
    typedef double (PixelDissimilarity<Pixel>::*DissimilarityFunc)(ImgIdx index1, ImgIdx index2) const;

    // // Constructor that accepts a function pointer
    PixelDissimilarity() = default;
    PixelDissimilarity(const Pixel *image_, ImgIdx imageSize_, uint8_t channels_, DissimilarityFunc func)
        : _image(image_), _imageSize(imageSize_), _channels(channels_), dissimilarityFunc(func) {}

    // Function to compute dissimilarity using the selected method
    double computeDissimilarity(ImgIdx index1, ImgIdx index2) const {
        return (this->*dissimilarityFunc)(index1, index2);
    }

    double maximumDissmilarity() const {
        return (double)std::numeric_limits<Pixel>::max() - (double)std::numeric_limits<Pixel>::min();
    }

    double L1(ImgIdx index1, ImgIdx index2) const;
    double L2(ImgIdx index1, ImgIdx index2) const;
    double LInfinity(ImgIdx index1, ImgIdx index2) const;

  private:
    const Pixel *_image = nullptr;
    ImgIdx _imageSize = 0;
    uint8_t _channels = 0;
    DissimilarityFunc dissimilarityFunc;
};
