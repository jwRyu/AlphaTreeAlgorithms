#include "RandGenImage.hpp"
#include <random>

void RandGenImage::randomize8(_uint8 *img, int width, int height, int bit_depth, int ch) {
    _uint8 pix;
    int imgsize = width * height;
    int shamt = 8 - bit_depth;

    if (img == nullptr)
        return;

    for (int i = 0; i < imgsize * ch; i++) {
        pix = ((_uint8)(rand() & 0xff));

        img[i] = pix >> shamt;
    }
}

void RandGenImage::randomize16(_uint16 *img, int width, int height, int bit_depth, int ch) {
    _uint16 pix;
    int imgsize = width * height;
    int shamt = 16 - bit_depth;

    if (img == nullptr)
        return;

    for (_int64 i = 0; i < imgsize * ch; i++) {
        pix = ((_uint16)(rand() & 0xff) << 8);
        pix |= ((_uint16)(rand() & 0xff));

        img[i] = (_uint16)(pix >> shamt);
    }
}

void RandGenImage::randomize32(_uint32 *img, int width, int height, int bit_depth) {
    _uint32 pix;
    int imgsize = width * height;
    int shamt = 32 - bit_depth;

    if (img == nullptr)
        return;

    for (_int64 i = 0; i < imgsize; i++) {
        pix = ((_uint32)(rand() & 0xff) << 24);
        pix |= ((_uint32)(rand() & 0xff) << 16);
        pix |= ((_uint32)(rand() & 0xff) << 8);
        pix |= ((_uint32)(rand() & 0xff));
        img[i] = (_uint32)(pix >> shamt);
    }
}

void RandGenImage::randomize64(_uint64 *img, int width, int height, int bit_depth) {
    _uint64 pix;
    int imgsize = width * height;
    int shamt = 64 - bit_depth;

    if (img == nullptr)
        return;

    for (_int64 i = 0; i < imgsize; i++) {
        pix = ((_uint64)(rand() & 0xff) << 56);
        pix |= ((_uint64)(rand() & 0xff) << 48);
        pix |= ((_uint64)(rand() & 0xff) << 40);
        pix |= ((_uint64)(rand() & 0xff) << 32);
        pix |= ((_uint64)(rand() & 0xff) << 24);
        pix |= ((_uint64)(rand() & 0xff) << 16);
        pix |= ((_uint64)(rand() & 0xff) << 8);
        pix |= ((_uint64)(rand() & 0xff));
        pix = pix % 100;
        img[i] = (_uint64)(pix >> shamt);
    }
}