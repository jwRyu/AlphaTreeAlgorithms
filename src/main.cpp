#include <PNGcodec.hpp>
#include <defines.hpp>

#include <AlphaTree.hpp>
#include <AlphaTreeConfig.hpp>
#include <RandGenImage.hpp>

// args: Filename, nchannels, numthreads, testimgsize, algorithmcode, bitdepth, tseflag
int main(int argc, char **argv) {
    srand(time(NULL));

    const auto configFileName = argc < 2 ? "config.txt" : std::string(argv[1]);
    alphatreeConfig.initialize(configFileName);
    auto config = alphatreeConfig.load(argc, argv);

    if (config.has_value() == false) {
        std::cerr << "Unable to open configuration file." << std::endl;
        return -1;
    }

    AlphaTreeConfig::AlphaTreeParameters params = config.value();

    const std::string filePath = params.imageFileName;

    auto [image, w, h, ch] = PNGCodec::imread(params.UseRandomlyGeneratedImages ? "RAND" : filePath);

    const bool reduceImageBitdepth = false;

    const auto &width = params.UseRandomlyGeneratedImages ? params.randomGenImageWidth : w;
    const auto &height = params.UseRandomlyGeneratedImages ? params.randomGenImageHeight : h;
    const auto &bitdepth = params.UseRandomlyGeneratedImages ? params.bitdepth : 16 * ch;
    const auto &nch = params.UseRandomlyGeneratedImages ? params.nchannels : ch;
    const auto &dMetric = params.dissimilarityMetric;
    const auto &conn = params.connectivity;
    const auto &algCode = params.alphaTreeAlgorithmCode;
    const auto &nthr = params.numthreads;
    const auto &nitr = params.numitr;
    const auto &tse = params.tse;
    const auto &iparam1 = params.iparam1;
    const auto &fparam1 = params.fparam1;
    const auto &fparam2 = params.fparam2;

    if (!params.UseRandomlyGeneratedImages)
        printf("Image file name: %s\n", params.imageFileName.c_str());
    printf("=======================================================================\n");
    printf("========== imgsize = %d x %d (%d bits, %d ch, %dN) ================\n", (int)height, (int)width,
           (int)bitdepth, (int)nch, params.connectivity);
    printf("=======================================================================\n");
    printf("-----------------------------------------------------------------------------------\n");
    printf("%d Running %s (%d threads)\n", (int)algCode, alphatreeConfig.getAlphaTreeAlgorithmName(algCode).c_str(),
           (int)nthr);
    printf("-----------------------------------------------------------------------------------\n");
    std::vector<double> runtimes;

    for (int itr = 0; itr < nitr; itr++) {
        double tStart = 0, tEnd = INFINITY;
        if (params.UseRandomlyGeneratedImages == true) {
            if (bitdepth > 32 && bitdepth <= 64) {
                uint64_t *image = (uint64_t *)Malloc(width * height * sizeof(uint64_t));
                RandGenImage::randomize64(image, width, height, bitdepth);
                AlphaTree<uint64_t> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1, fparam2,
                                    iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 16) {
                uint32_t *image = (uint32_t *)Malloc(width * height * sizeof(uint32_t));
                RandGenImage::randomize32(image, width, height, bitdepth);
                AlphaTree<uint32_t> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1, fparam2,
                                    iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 8) {
                uint16_t *image = (uint16_t *)Malloc(width * height * nch * sizeof(uint16_t));
                RandGenImage::randomize16(image, width, height, bitdepth, nch);
                AlphaTree<uint16_t> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image, height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1, fparam2,
                                    iparam1);
                tEnd = get_wall_time();
                Free(image);
            } else if (bitdepth > 0) {
                uint8_t *image = (uint8_t *)Malloc(width * height * nch * sizeof(uint8_t));
                RandGenImage::randomize8(image, width, height, bitdepth, nch);
                AlphaTree<uint8_t> tree;
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

            if (reduceImageBitdepth) {
                const size_t shamt = 16 - params.bitdepth;
                uint16_t pMin = 65535;
                uint16_t pMax = 0;

                for (auto &pixel : image) {
                    pixel = pixel >> shamt;
                    pMax = std::max(pMax, pixel);
                    pMin = std::min(pMin, pixel);
                }

                printf("reduceImageBitdepth - BitDepth = params.bitdepth = %d / DR = %d - %d\n", (int)params.bitdepth,
                       (int)pMin, (int)pMax);
            }

            uint16_t maxVal = *std::max_element(image.begin(), image.end());
            if ((int)maxVal > (int)255) {
                AlphaTree<uint16_t> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image.data(), height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1,
                                    fparam2, iparam1);

                const bool rgbFilter = true;
                if (rgbFilter) {

                    int sizeThr = 16;
                    for (int i = 1; i < 300; i++) {
                        int alphaThr = 0.1 + i * 10;
                        tree.AlphaFilter(image.data(), alphaThr, sizeThr);
                        std::string str = "out_" + std::to_string(alphaThr) + ".png";
                        PNGCodec::imwrite(image, w, h, ch, str.c_str());
                    }
                }

                tEnd = get_wall_time();
            } else {
                std::vector<uint8_t> image8(image.size());
                for (size_t i = 0; i < image.size(); i++)
                    image8[i] = (uint8_t)image[i];
                AlphaTree<uint8_t> tree;
                tStart = get_wall_time();
                tree.BuildAlphaTree(image8.data(), height, width, nch, dMetric, conn, algCode, nthr, tse, fparam1,
                                    fparam2, iparam1);
                tEnd = get_wall_time();
            }
        }

        auto runtime = tEnd - tStart;
        printf("-------------------Run %d/%d: %.3f------------------\n", (int)itr + 1, nitr, runtime);
        runtimes.push_back(runtime);
    }

    if (runtimes.empty() == false) {
        double minRuntime = *std::min_element(runtimes.begin(), runtimes.end());
        double imgsize = (double)(width * height);

        printf("================== Summary ==================\n");
        printf("Processing speed: %.3fMpix/s / Memory use %.3fB/pix \n", (imgsize / minRuntime) * 1e-6,
               (double)max_memuse / imgsize);
        printf("=============================================\n");
    }

    return 0;
}
