#pragma once

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <png.h>
#include <vector>

std::tuple<std::vector<unsigned char>, int, int, int> loadPng(const std::string &filename);
bool writePng(const std::string &filename, const std::vector<unsigned char> &image, int width, int height,
              int channels);