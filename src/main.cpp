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
    auto [image, w, h, ch] = PNGCodec::imread("img03.png");

    std::vector<RankItem<float>> rankitems;

    {
        ImgIdx width = w;
        ImgIdx height = h;
        ImgIdx imgsize = width * height;
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                ImgIdx pixelIndex = i * width + j;

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

                    rankitems.emplace_back(alpha, pixelIndex * 2);
                }

                if (j < w - 1) {
                    float rRight = (float)image[i * width + j + 1];
                    float gRight = (float)image[imgsize + i * width + j + 1];
                    float bRight = (float)image[2 * imgsize + i * width + j + 1];

                    float dr = r - rRight;
                    float dg = g - gRight;
                    float db = b - bRight;
                    float alpha = std::sqrt((dr * dr + dg * dg + db * db) / 3.0f);

                    rankitems.emplace_back(alpha, pixelIndex * 2 + 1);
                }
            }
        }
    }

    {
        auto t0 = get_wall_time();

        auto copyData = rankitems;

        for (int i = 0; i < 10; i++)
            printf("alpha = %f / index = %d \n", (double)copyData[i].alpha, (int)copyData[i].dimgidx);

        std::sort(copyData.begin(), copyData.end(),
                  [](const RankItem<float> &a, const RankItem<float> &b) { return a.alpha < b.alpha; });

        auto tSort = get_wall_time() - t0;
        printf("std::sort(): %fs\n", tSort);

        for (int i = 0; i < 10; i++)
            printf("alpha = %f / index = %d \n", (double)copyData[i].alpha, (int)copyData[i].dimgidx);
    }

    PNGCodec::imwrite(image, w, h, ch, "out.png");

    const auto &width = params.randomGenImageWidth;
    const auto &height = params.randomGenImageHeight;
    const auto &bitdepth = params.bitdepth;
    const auto &nch = params.nchannels;
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

    // if (argc == 1) {
    //     printf("args: Filename, nchannels, numthreads, testimgsize, algorithmcode, bitdepth, tseflag, outfnameheader
    //     "
    //            "connectivity float_param1 float_param2 int_param1\n");
    //     return 0;
    //     printf("No Arg - running default tests\n");

    //     int width = 4;
    //     int height = 4;

    //     _uint32 *img = getRandomizedImage<_uint32>(width * height, 32, 1);

    //     AlphaTree<_uint32> tree;

    //     // tree.BuildAlphaTree(img, height, width, 1, 4, FLOOD_HIERARHEAPQUEUE_CACHE, 1, 0, 0.0, 0.0, 0);
    //     tree.BuildAlphaTree(img, height, width, 1, 4, FLOOD_LADDERQUEUE, 1, 0, 0.0, 0.0, 32);
    //     // tree.print_tree();

    //     return 0;
    // }

    // using namespace std;

    // AlphaTreeParameters input;
    // parse_input(input, argc, argv);
    // int randomimg = (strcmp(input.name, "rand") == 0);
    // double meanrunspeed[16] =
    //     {
    //         0,
    //     },
    //        maxmemuse[16] = {
    //            0,
    //        };
    // int nthr[] = {1, 2, 4, 8, 16, 32, 48, 64, 96, 128, 192, 256, 480, 960, 1920};
    // char algname[256];
    // char outfname[128];

    // int numimg = randomimg ? 1 : 10;
    // int imgidxstart = 0;
    // int thrstart = 0;

    // // srand(time(NULL));

    // // if(!randomimg)
    // {
    //     sprintf(outfname, "%s_ch%d_nthr%d_alg%d_bit%d_tse%d_conn%d.txt", input.fnameheader, input.nchannels,
    //             input.numthreads, input.algorithmcode, input.bitdepth, input.tse, input.connectivity);
    //     printf("output file name = %s\n", outfname);
    //     fstream f(outfname);
    //     if (0 && f.good()) {
    //         char buf[128];
    //         ifstream fin(outfname);

    //         // f.seekg(0, std::ios::beg);
    //         fin.getline(buf, 128, ' ');
    //         imgidxstart = atoi(buf);
    //         fin.getline(buf, 128, '\n');
    //         thrstart = atoi(buf);
    //         int numthr = randomimg ? 1 : input.numthreads;
    //         for (int thridx = 0; thridx < numthr; thridx++) {
    //             fin.getline(buf, 128, ' ');
    //             meanrunspeed[thridx] = atof(buf);
    //             fin.getline(buf, 128, '\n');
    //             maxmemuse[thridx] = atof(buf);
    //         }
    //         printf("saved progress loaded (%d/%d)\n", imgidxstart, numimg);
    //         fin.close();
    //     } else {
    //     }
    // }

    // for (int imgidx = imgidxstart; imgidx < numimg; imgidx++) {
    //     int height, width;
    //     _uint16 *img;
    //     _uint8 *img8 = 0;

    //     height = width = input.testimgsize;
    //     int image_number = imgidx;

    //     printf("========================================================================================\n");
    //     printf("========== image [%d/%d]: imgsize = %d x %d (%d bits, %d channels) ================\n", imgidx + 1,
    //            numimg, (int)height, (int)width, (int)input.bitdepth, (int)input.nchannels);
    //     printf("========================================================================================\n");
    //     int thritr = randomimg ? 1 : input.numthreads;
    //     for (int thridx = thrstart; thridx < thritr; thridx++) {
    //         int numthreads = randomimg ? input.numthreads : nthr[thridx];
    //         alg_name(algname, (int)input.algorithmcode);
    //         printf("-----------------------------------------------------------------------------------\n");
    //         printf("%d Running %s (%d threads)\n", (int)input.algorithmcode, algname, (int)numthreads);
    //         printf("-----------------------------------------------------------------------------------\n");
    //         double minruntime = 0;

    //         for (int testrep = 0; testrep < input.numitr; testrep++) {
    //             double t = get_cpu_time();
    //             double runtime = 0.0; // = get_cpu_time() - t;
    //             // t = omp_get_wtime();

    //             if (!randomimg)
    //                 imread(img, input.name, image_number, height, width, input.nchannels,
    //                        input.bitdepth); // skip image 5 (unknown bug)

    //             if (!randomimg && input.bitdepth < 16) {
    //                 if (input.bitdepth == 8)
    //                     img8 = new _uint8[height * width];
    //                 adjust_bitdepth(img8, img, height * width, input.bitdepth);
    //                 delete[] img;
    //                 img = 0;
    //             }

    //             if (randomimg) {
    //                 img8 = 0;
    //                 img = 0;
    //                 int bitdepth = input.bitdepth;
    //                 _uint32 *img32 = 0;
    //                 _uint64 *img64 = 0;
    //                 if (bitdepth <= 8)
    //                     Randomizedimage(img8, height * width, bitdepth, 1);
    //                 else if (bitdepth <= 16)
    //                     Randomizedimage(img, height * width, bitdepth, 1);
    //                 else if (bitdepth <= 32)
    //                     Randomizedimage(img32, height * width, bitdepth, 1);
    //                 else
    //                     Randomizedimage(img64, height * width, bitdepth, 1);

    //                 t = get_cpu_time();
    //                 if (bitdepth <= 8) {
    //                     AlphaTree<_uint8> *tree = new AlphaTree<_uint8>;
    //                     tree->BuildAlphaTree(img8, height, width, input.nchannels, input.connectivity,
    //                                          input.algorithmcode, (int)numthreads, input.tse, input.fparam1,
    //                                          input.fparam2, input.iparam1);
    //                     delete tree;
    //                 } else if (bitdepth <= 16) {
    //                     AlphaTree<_uint16> *tree = new AlphaTree<_uint16>;
    //                     tree->BuildAlphaTree(img, height, width, input.nchannels, input.connectivity,
    //                                          input.algorithmcode, (int)numthreads, input.tse, input.fparam1,
    //                                          input.fparam2, input.iparam1);
    //                     delete tree;
    //                 } else if (bitdepth <= 32) {
    //                     AlphaTree<_uint32> *tree = new AlphaTree<_uint32>;
    //                     tree->BuildAlphaTree(img32, height, width, input.nchannels, input.connectivity,
    //                                          input.algorithmcode, (int)numthreads, input.tse, input.fparam1,
    //                                          input.fparam2, input.iparam1);
    //                     delete tree;
    //                 } else {
    //                     AlphaTree<_uint64> *tree = new AlphaTree<_uint64>;
    //                     tree->BuildAlphaTree(img64, height, width, input.nchannels, input.connectivity,
    //                                          input.algorithmcode, (int)numthreads, input.tse, input.fparam1,
    //                                          input.fparam2, input.iparam1);
    //                     delete tree;
    //                 }
    //                 runtime = get_cpu_time() - t;

    //                 if (img32)
    //                     delete[] img32;
    //                 if (img64)
    //                     delete[] img64;
    //             }

    //             printf("-------------------Run %d: %.3f------------------\n", (int)testrep, runtime);

    //             // cout << "Run " << testrep << ": " << runtime << endl;
    //             if (!testrep)
    //                 minruntime = runtime;
    //             minruntime = (minruntime > runtime) ? runtime : minruntime;

    //             // delete tree;
    //             if (img)
    //                 delete[] img;
    //             if (img8)
    //                 delete[] img8;
    //         }
    //         double imgsize = (double)height * (double)width;
    //         double runspeed = imgsize / (1e6 * minruntime);
    //         meanrunspeed[thridx] += runspeed;
    //         maxmemuse[thridx] += max_memuse / imgsize;

    //         if (0 && thridx < thritr - 1) {
    //             ofstream fout(outfname);
    //             fout << imgidx << ' ' << thridx + 1 << endl;
    //             for (int thridx = 0; thridx < input.numthreads; thridx++) {
    //                 fout << meanrunspeed[thridx] << ' ' << maxmemuse[thridx] << endl;
    //             }
    //             fout.close();
    //         }
    //     }

    //     thrstart = 0;
    // }
    // printf("Summary ======================================\n");
    // int numthrs = randomimg ? 1 : input.numthreads;
    // for (int thridx = 0; thridx < numthrs; thridx++) {
    //     printf("#thr = %d: %.3fMpix/s / %.3fB/pix \n", nthr[thridx], meanrunspeed[thridx] / numimg,
    //            maxmemuse[thridx] / numimg);
    // }
    // printf("==============================================\n");

    // return 0;
}
