#include "AlphaTreeConfig.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

AlphaTreeConfig::AlphaTreeConfig(const std::string &filename, const std::string &commentToken)
    : filename(filename), commentToken(commentToken) {
    AlgorithmNameToCode.clear();
    AlgorithmCodeToName.clear();

    AlgorithmNameToCode["Unionfind"] = 0;
    AlgorithmNameToCode["FloodHierQueue"] = 1;
    AlgorithmNameToCode["FloodTrieQueue"] = 2;
    AlgorithmNameToCode["FloodHeapQueue"] = 3;
    AlgorithmNameToCode["FloodHierHeapQueue"] = 4;
    AlgorithmNameToCode["FloodLadderQueue"] = 5;
    AlgorithmNameToCode["FloodHierQueueParallel"] = 6;
    AlgorithmNameToCode["HybridParallel"] = 7;
    AlgorithmNameToCode["FloodHierQueueParalle"] = 8;
    AlgorithmNameToCode["FloodHierQueueNoCache"] = 9;
    AlgorithmNameToCode["FloodTrieQueueNoCache"] = 10;
    AlgorithmNameToCode["FloodHeapQueueNoCache"] = 11;
    AlgorithmNameToCode["FloodHierQueueHypergraph"] = 12;
    AlgorithmNameToCode["FloodTrieQueueHypergraph"] = 13;
    AlgorithmNameToCode["FloodHierHeapQueueHisteq"] = 14;
    AlgorithmNameToCode["FloodHierHeapQueueNoCache"] = 15;
    AlgorithmNameToCode["FloodHeapQueueNaiveNoCache"] = 16;
    for (const auto &pair : AlgorithmNameToCode)
        AlgorithmCodeToName[pair.second] = pair.first;
}

int AlphaTreeConfig::getAlphaTreeAlgorithmCode(std::string AlphaTreeAlgorithmName) {
    return AlgorithmNameToCode[AlphaTreeAlgorithmName];
}

std::string AlphaTreeConfig::getAlphaTreeAlgorithmName(int AlphaTreeAlgorithmCode) {
    return AlgorithmCodeToName[AlphaTreeAlgorithmCode];
}

std::optional<AlphaTreeConfig::AlphaTreeParameters> AlphaTreeConfig::load(int argc, char **argv) {
    if (argc == 1) {
        std::ifstream configFile(filename);
        if (!configFile.is_open()) {
            return std::nullopt;
        }

        std::string line;
        while (std::getline(configFile, line)) {
            line = trim(line);
            if (line.empty() || line.find(commentToken) == 0) {
                continue; // Skip empty lines and comment lines
            }

            std::istringstream lineStream(line);
            std::string key;
            if (std::getline(lineStream, key, '=')) {
                std::string value;
                if (std::getline(lineStream, value)) {
                    key = trim(key);
                    value = trim(value);
                    config[key] = value;
                }
            }
        }
        configFile.close();
    } else if (argc == 19) {
        config["ImageFileName"] = std::string(argv[1]);
        config["LogFileName"] = std::string(argv[2]);
        config["UseRandomlyGeneratedImages"] = std::string(argv[3]);
        config["RandomlyGeneratedImageWidth"] = std::string(argv[4]);
        config["RandomlyGeneratedImageHeight"] = std::string(argv[5]);
        config["NumberOfChannels"] = std::string(argv[6]);
        config["NumberOfThreads"] = std::string(argv[7]);
        config["AlphaTreeAlgorithm"] = std::string(argv[8]);
        config["BitDepth"] = std::string(argv[9]);
        config["UseTreeSizeEstimation"] = std::string(argv[10]);
        config["Connectivity"] = std::string(argv[11]);
        config["NumberOfTestIterations"] = std::string(argv[12]);
        config["ParameterInteger1"] = std::string(argv[13]);
        config["ParameterInteger2"] = std::string(argv[14]);
        config["ParameterInteger3"] = std::string(argv[15]);
        config["ParameterFloat1"] = std::string(argv[16]);
        config["ParameterFloat2"] = std::string(argv[17]);
        config["ParameterFloat3"] = std::string(argv[18]);
    } else {
        std::cerr << "[AlphaTreeConfig::load] Invalid argument (requires either 0 or 18)" << std::endl;
        return std::nullopt;
    }

    AlphaTreeConfig::AlphaTreeParameters params;
    params.imageFileName = getString("ImageFileName");
    params.logFileName = getString("LogFileName");
    params.nchannels = getInteger("NumberOfChannels");
    params.numthreads = getInteger("NumberOfThreads");
    params.UseRandomlyGeneratedImages = getInteger("UseRandomlyGeneratedImages");
    params.randomGenImageWidth = getInteger("RandomlyGeneratedImageWidth");
    params.randomGenImageHeight = getInteger("RandomlyGeneratedImageHeight");
    params.alphaTreeAlgorithmCode = getAlphaTreeAlgorithmCode(getString("AlphaTreeAlgorithm"));
    params.bitdepth = getInteger("BitDepth");
    params.tse = getInteger("UseTreeSizeEstimation");
    params.connectivity = getInteger("Connectivity");
    params.numitr = getInteger("NumberOfTestIterations");
    params.iparam1 = getInteger("ParameterInteger1");
    params.iparam2 = getInteger("ParameterInteger2");
    params.iparam3 = getInteger("ParameterInteger3");
    params.fparam1 = getDouble("ParameterFloat1");
    params.fparam2 = getDouble("ParameterFloat2");
    params.fparam3 = getDouble("ParameterFloat3");

    return params;
}

std::string AlphaTreeConfig::getString(const std::string &key) const {
    auto it = config.find(key);
    if (it != config.end()) {
        return it->second;
    }
    return "";
}

int AlphaTreeConfig::getInteger(const std::string &key) const {
    auto it = config.find(key);
    if (it != config.end()) {
        if (it->second.size())
            return std::stoi(it->second);
        return -1;
    }
    return -1;
}

double AlphaTreeConfig::getDouble(const std::string &key) const {
    auto it = config.find(key);
    if (it != config.end()) {
        if (it->second.size())
            return std::stod(it->second);
        return -1.0;
    }
    return -1.0;
}

std::string AlphaTreeConfig::trim(const std::string &str) const {
    size_t first = str.find_first_not_of(' ');
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}