#include "defines.h"

class RandGenImage {
  public:
    static void randomize8(_uint8 *img, int width, int height, int bit_depth, int ch);
    static void randomize16(_uint16 *img, int width, int height, int bit_depth, int ch);
    static void randomize32(_uint32 *img, int width, int height, int bit_depth);
    static void randomize64(_uint64 *img, int width, int height, int bit_depth);
};