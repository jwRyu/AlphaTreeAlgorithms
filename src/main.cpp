#include "AlphaTree.h"
#include "AlphaTreeConfig.h"
#include "BucketSort.hpp"
#include "RandGenImage.hpp"
#include "defines.h"
#include "pngReader.hpp"
#include "walltime.h"
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>

// args: Filename, nchannels, numthreads, testimgsize, algorithmcode, bitdepth, tseflag
int main(int argc, char **argv) {

    auto config = alphatreeConfig.load(argc, argv);

    if (config.has_value() == false) {
        std::cerr << "Unable to open configuration file." << std::endl;
        return -1;
    }

    AlphaTreeConfig::AlphaTreeParameters params = config.value();
    auto filename = params.imageFileName;
    auto width = params.randomGenImageWidth;
    auto height = params.randomGenImageHeight;
    auto bitdepth = params.bitdepth;
    auto nch = params.nchannels;
    auto conn = params.connectivity;
    auto algCode = params.alphaTreeAlgorithmCode;
    auto nthr = params.numthreads;
    auto nitr = params.numitr;
    auto tse = params.tse;
    auto iparam1 = params.iparam1;
    auto fparam1 = params.fparam1;
    auto fparam2 = params.fparam2;

    std::cout << "File name: " << filename << std::endl;
    auto [image, w, h, c] = loadPng(filename);
    std::cout << "Image size = " << image.size() << std::endl;
    width = w;
    height = h;

    const int imgsize = w * h;
    std::vector<std::pair<float, int>> alphaAndIndex;

    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            int pixelIndex = i * width + j;

            float r = (float)image[i * width + j];
            float g = (float)image[imgsize + i * width + j];
            float b = (float)image[2 * imgsize + i * width + j];

            if (i < h - 1) {
                float rBottom = (float)image[(i + 1) * width + j];
                float gBottom = (float)image[imgsize + (i + 1) * width + j];
                float bBottom = (float)image[2 * imgsize + (i + 1) * width + j];

                float dr = r - rBottom;
                float dg = g - gBottom;
                float db = b - bBottom;
                float alpha = std::sqrt((dr * dr + dg * dg + db * db) / 3.0f);

                alphaAndIndex.emplace_back(alpha, pixelIndex * 2);
            }

            if (j < w - 1) {
                float rRight = (float)image[i * width + j + 1];
                float gRight = (float)image[imgsize + i * width + j + 1];
                float bRight = (float)image[2 * imgsize + i * width + j + 1];

                float dr = r - rRight;
                float dg = g - gRight;
                float db = b - bRight;
                float alpha = std::sqrt((dr * dr + dg * dg + db * db) / 3.0f);

                alphaAndIndex.emplace_back(alpha, pixelIndex * 2 + 1);
            }
        }
    }

    printf("alphaAndIndex.size() = %lu\n", alphaAndIndex.size());

    // for (int i = 0; i < 10000; i++) {
    //     // printf("alpha = %f / index = %d \n", alphaAndIndex[i].first, alphaAndIndex[i].second);
    //     printf("%f ", alphaAndIndex[i].first);
    // }
    // printf("\n");

    {
        auto t0 = get_wall_time();

        auto copyData = alphaAndIndex;

        for (int i = 0; i < 10; i++) {
            printf("alpha = %f / index = %d \n", copyData[i].first, copyData[i].second);
        }

        BucketSort::sort(copyData);

        auto tSort = get_wall_time() - t0;
        printf("BucketSort::sort(): %fs\n", tSort);

        for (int i = 0; i < 10; i++) {
            printf("alpha = %f / index = %d \n", copyData[i].first, copyData[i].second);
        }
    }

    {
        auto t0 = get_wall_time();

        auto copyData = alphaAndIndex;

        for (int i = 0; i < 10; i++) {
            printf("alpha = %f / index = %d \n", copyData[i].first, copyData[i].second);
        }

        std::sort(copyData.begin(), copyData.end(),
                  [](std::pair<float, int> &a, std::pair<float, int> &b) { return a.first < b.first; });

        auto tSort = get_wall_time() - t0;
        printf("std::sort(): %fs\n", tSort);

        for (int i = 0; i < 10; i++) {
            printf("alpha = %f / index = %d \n", copyData[i].first, copyData[i].second);
        }
    }

    return 0;

    printf("====================================================================================\n");
    printf("========== image [%d/%d]: imgsize = %d x %d (%d bits, %d channels) ================\n", 0, 0, (int)height,
           (int)width, (int)bitdepth, (int)nch);
    printf("====================================================================================\n");
    printf("-----------------------------------------------------------------------------------\n");
    printf("%d Running %s (%d threads)\n", (int)algCode, alphatreeConfig.getAlphaTreeAlgorithmName(algCode).c_str(),
           (int)nthr);
    printf("-----------------------------------------------------------------------------------\n");
    std::vector<double> runtimes;

    for (int itr = 0; itr < nitr; itr++) {
        double tStart = 0, tEnd = INFINITY;
        if (params.UseRandomlyGeneratedImages == true) {
            if (bitdepth > 32 && bitdepth <= 64) {
                _uint64 *image = (_uint64 *)Malloc(width * height * sizeof(_uint64));
                RandGenImage::randomize64(image, width, height, bitdepth);
                AlphaTree<_uint64> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, conn, algCode, nthr, tse, fparam1, fparam2, iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 16) {
                _uint32 *image = (_uint32 *)Malloc(width * height * sizeof(_uint32));
                RandGenImage::randomize32(image, width, height, bitdepth);
                AlphaTree<_uint32> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, conn, algCode, nthr, tse, fparam1, fparam2, iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 8) {
                _uint16 *image = (_uint16 *)Malloc(width * height * sizeof(_uint16));
                RandGenImage::randomize16(image, width, height, bitdepth, nch);
                AlphaTree<_uint16> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, conn, algCode, nthr, tse, fparam1, fparam2, iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 0) {
                _uint8 *image = (_uint8 *)Malloc(width * height * sizeof(_uint8));
                RandGenImage::randomize8(image, width, height, bitdepth, nch);
                AlphaTree<_uint8> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, conn, algCode, nthr, tse, fparam1, fparam2, iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else {
                std::cerr << "Invalid bit-depth. " << std::endl;
                return -1;
            }

        } else {
        }

        auto runtime = tEnd - tStart;
        printf("-------------------Run %d/%d: %.3f------------------\n", (int)itr + 1, nitr, runtime);
        runtimes.push_back(runtime);
    }

    if (runtimes.empty() == false) {
        double minRuntime = *std::min_element(runtimes.begin(), runtimes.end());
        double imgsize = (double)(width * height);

        printf("================== Summary ==================\n");
        for (int thridx = 0; thridx < nthr; thridx++) {
            printf("#thr = %d: %.3fMpix/s / %.3fB/pix \n", 1, (imgsize / minRuntime) * 1e-6, 0.0);
        }
        printf("=============================================\n");
    }

    return 0;
}
