#include <PNGcodec.hpp>

cv::Mat PNGCodec::toCVMat(const std::vector<uint16_t> &image, int width, int height, int channels) {
    if (channels == 1) {
        cv::Mat imageCV(height, width, CV_16U);
        std::copy(image.begin(), image.end(), imageCV.begin<uint16_t>());
        return imageCV;
    } else if (channels == 3) {
        cv::Mat imageCV(height, width, CV_16UC3);
        std::vector<cv::Mat> cvChannel(3);
        cv::split(imageCV, cvChannel);

        // Blue
        std::copy(image.begin() + 2 * width * height, image.end(), cvChannel[0].begin<uint16_t>());

        // Green
        std::copy(image.begin() + 1 * width * height, image.begin() + 2 * width * height,
                  cvChannel[1].begin<uint16_t>());

        // Red
        std::copy(image.begin(), image.begin() + width * height, cvChannel[2].begin<uint16_t>());

        cv::merge(cvChannel, imageCV);
        return imageCV;
    }

    cv::Mat nullImage(height, width, CV_16U);
    return nullImage;
}

std::vector<uint16_t> PNGCodec::toImage(const cv::Mat &imageCV) {
    if (imageCV.type() == CV_16U) {
        std::vector<uint16_t> image(imageCV.cols * imageCV.rows);
        std::copy(imageCV.begin<uint16_t>(), imageCV.end<uint16_t>(), image.begin());
        return image;
    } else if (imageCV.type() == CV_16UC3) {
        cv::Mat b, g, r;

        cv::extractChannel(imageCV, b, 0);
        cv::extractChannel(imageCV, g, 1);
        cv::extractChannel(imageCV, r, 2);

        std::vector<uint16_t> image(imageCV.cols * imageCV.rows * imageCV.channels());
        std::copy(r.begin<uint16_t>(), r.end<uint16_t>(), image.begin());
        std::copy(g.begin<uint16_t>(), g.end<uint16_t>(), image.begin() + imageCV.cols * imageCV.rows);
        std::copy(b.begin<uint16_t>(), b.end<uint16_t>(), image.begin() + imageCV.cols * imageCV.rows * 2);
        return image;
    }

    return {};
}

std::tuple<std::vector<uint16_t>, int, int, int> PNGCodec::imread(const std::string &filename) {

    // // Read the image file
    cv::Mat imageCV = cv::imread(filename, cv::IMREAD_UNCHANGED);

    // Check if the image was loaded successfully
    if (imageCV.empty()) {
        std::cerr << "Error: Could not open or find the image " << filename << std::endl;
        return {{}, -1, -1, -1};
    }

    auto image = toImage(imageCV);

    return {image, imageCV.cols, imageCV.rows, imageCV.channels()};
}

bool PNGCodec::imwrite(const std::vector<uint16_t> &image, int width, int height, int channels,
                       const std::string &filename) {
    auto imageCV = toCVMat(image, width, height, channels);

    // Save the image to a new file
    if (!cv::imwrite(filename, imageCV)) {
        std::cerr << "Error: Could not save the image to " << filename << std::endl;
        return false;
    }

    return true;
}