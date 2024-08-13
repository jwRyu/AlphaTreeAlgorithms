
#include "pngReader.hpp"
#include <png.h>

std::tuple<std::vector<unsigned char>, int, int, int> loadPng(const std::string &filename) {
    std::vector<unsigned char> image;
    int width = 0, height = 0, channels = 0;
    FILE *fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return {image, width, height, channels};
    }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        std::cerr << "Error: Could not create png read struct" << std::endl;
        fclose(fp);
        return {image, width, height, channels};
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        std::cerr << "Error: Could not create png info struct" << std::endl;
        png_destroy_read_struct(&png, NULL, NULL);
        fclose(fp);
        return {image, width, height, channels};
    }

    if (setjmp(png_jmpbuf(png))) {
        std::cerr << "Error: setjmp failed" << std::endl;
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return {image, width, height, channels};
    }

    png_init_io(png, fp);
    png_read_info(png, info);

    width = png_get_image_width(png, info);
    height = png_get_image_height(png, info);
    channels = png_get_channels(png, info);

    png_byte bit_depth = png_get_bit_depth(png, info);
    png_byte color_type = png_get_color_type(png, info);

    if (bit_depth == 16)
        png_set_strip_16(png);

    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png);

    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(png);

    if (png_get_valid(png, info, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(png);

    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(png);

    png_read_update_info(png, info);

    png_bytep *row_pointers = (png_bytep *)malloc(sizeof(png_bytep) * height);
    for (int y = 0; y < height; y++) {
        row_pointers[y] = (png_byte *)malloc(png_get_rowbytes(png, info));
    }

    png_read_image(png, row_pointers);

    image.resize(width * height * 4);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            png_bytep px = &(row_pointers[y][x * 4]);
            image[y * width + x] = px[0];
            image[y * width + x + width * height * 1] = px[1];
            image[y * width + x + width * height * 2] = px[2];
        }
    }

    for (int y = 0; y < height; y++) {
        free(row_pointers[y]);
    }
    free(row_pointers);

    png_destroy_read_struct(&png, &info, NULL);
    fclose(fp);
    return {image, width, height, channels};
}

bool writePng(const std::string &filename, const std::vector<unsigned char> &image, int width, int height,
              int channels) {
    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        std::cerr << "Error: Could not create png write struct" << std::endl;
        fclose(fp);
        return false;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        std::cerr << "Error: Could not create png info struct" << std::endl;
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        return false;
    }

    if (setjmp(png_jmpbuf(png))) {
        std::cerr << "Error: setjmp failed" << std::endl;
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return false;
    }

    png_init_io(png, fp);

    png_set_IHDR(png, info, width, height, 8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    png_bytep row = (png_bytep)malloc(4 * width * sizeof(png_byte));
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            row[x * 4 + 0] = image[y * width + x];
            row[x * 4 + 1] = image[y * width + x + width * height * 1];
            row[x * 4 + 2] = image[y * width + x + width * height * 2];
            row[x * 4 + 3] = 255;
        }
        png_write_row(png, row);
    }
    free(row);

    png_write_end(png, NULL);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    return true;
}
