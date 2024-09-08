#include "RandGenImage.hpp"
#include <random>

void RandGenImage::randomize8(uint8_t *img, int width, int height, int bit_depth, int ch) {
    uint8_t pix;
    int imgsize = width * height;
    int shamt = 8 - bit_depth;

    if (img == nullptr)
        return;

    for (int i = 0; i < imgsize * ch; i++) {
        pix = ((uint8_t)(rand() & 0xff));

        img[i] = pix >> shamt;
    }
}

void RandGenImage::randomize16(uint16_t *img, int width, int height, int bit_depth, int ch) {
    uint16_t pix;
    int imgsize = width * height;
    int shamt = 16 - bit_depth;

    if (img == nullptr)
        return;

    for (int64_t i = 0; i < imgsize * ch; i++) {
        pix = ((uint16_t)(rand() & 0xff) << 8);
        pix |= ((uint16_t)(rand() & 0xff));

        img[i] = (uint16_t)(pix >> shamt);
    }
}

void RandGenImage::randomize32(uint32_t *img, int width, int height, int bit_depth) {
    uint32_t pix;
    int imgsize = width * height;
    int shamt = 32 - bit_depth;

    if (img == nullptr)
        return;

    for (int64_t i = 0; i < imgsize; i++) {
        pix = ((uint32_t)(rand() & 0xff) << 24);
        pix |= ((uint32_t)(rand() & 0xff) << 16);
        pix |= ((uint32_t)(rand() & 0xff) << 8);
        pix |= ((uint32_t)(rand() & 0xff));
        img[i] = (uint32_t)(pix >> shamt);
    }
}

void RandGenImage::randomize64(uint64_t *img, int width, int height, int bit_depth) {
    uint64_t pix;
    int imgsize = width * height;
    int shamt = 64 - bit_depth;

    if (img == nullptr)
        return;

    for (int64_t i = 0; i < imgsize; i++) {
        pix = ((uint64_t)(rand() & 0xff) << 56);
        pix |= ((uint64_t)(rand() & 0xff) << 48);
        pix |= ((uint64_t)(rand() & 0xff) << 40);
        pix |= ((uint64_t)(rand() & 0xff) << 32);
        pix |= ((uint64_t)(rand() & 0xff) << 24);
        pix |= ((uint64_t)(rand() & 0xff) << 16);
        pix |= ((uint64_t)(rand() & 0xff) << 8);
        pix |= ((uint64_t)(rand() & 0xff));
        // pix = pix % 100;
        img[i] = (uint64_t)(pix >> shamt);
    }
}