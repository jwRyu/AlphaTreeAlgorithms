#include "PNGcodec.hpp"

#include "AlphaTree.h"
#include "AlphaTreeConfig.h"
#include "RandGenImage.hpp"
#include "defines.h"
#include "pgmio.h"
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

    // auto input_filename = "img03.png";
    // auto output_filename = "out.png";
    // auto [image, w, h, ch] = PNGCodec::imread("img03.png");
    // PNGCodec::imwrite(image, w, h, ch, "out.png");

    const auto &width = params.randomGenImageWidth;
    const auto &height = params.randomGenImageHeight;
    const auto &bitdepth = params.bitdepth;
    const auto &nch = params.nchannels;
    const auto &dMetric = params.dissimilarityMetric;
    const auto &conn = params.connectivity;
    const auto &algCode = params.alphaTreeAlgorithmCode;
    const auto &nthr = params.numthreads;
    const auto &nitr = params.numitr;
    const auto &tse = params.tse;
    const auto &iparam1 = params.iparam1;
    // const auto &iparam2 = params.iparam2;
    // const auto &iparam3 = params.iparam3;
    const auto &fparam1 = params.fparam1;
    const auto &fparam2 = params.fparam2;
    // const auto &fparam3 = params.fparam3;

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
                tree.BuildAlphaTree(image, height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1, fparam2,
                                    iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 16) {
                _uint32 *image = (_uint32 *)Malloc(width * height * sizeof(_uint32));
                RandGenImage::randomize32(image, width, height, bitdepth);
                AlphaTree<_uint32> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1, fparam2,
                                    iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 8) {
                _uint16 *image = (_uint16 *)Malloc(width * height * nch * sizeof(_uint16));
                RandGenImage::randomize16(image, width, height, bitdepth, nch);
                AlphaTree<_uint16> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1, fparam2,
                                    iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 0) {
                _uint8 *image = (_uint8 *)Malloc(width * height * nch * sizeof(_uint8));
                RandGenImage::randomize8(image, width, height, bitdepth, nch);
                AlphaTree<_uint8> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1, fparam2,
                                    iparam1);
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
