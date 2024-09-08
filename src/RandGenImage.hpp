#include <defines.hpp>

class RandGenImage {
  public:
    static void randomize8(uint8_t *img, int width, int height, int bit_depth, int ch);
    static void randomize16(uint16_t *img, int width, int height, int bit_depth, int ch);
    static void randomize32(uint32_t *img, int width, int height, int bit_depth);
    static void randomize64(uint64_t *img, int width, int height, int bit_depth);
};