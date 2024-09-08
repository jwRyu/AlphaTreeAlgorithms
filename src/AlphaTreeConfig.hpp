
#ifndef ALPHATREECONFIG_H
#define ALPHATREECONFIG_H

#define CONFIGFILENAME "config.txt"

#include <defines.hpp>

class AlphaTreeConfig {
  public:
    struct AlphaTreeParameters {
        std::string imageFileName;
        std::string logFileName;
        int nchannels;
        int numthreads;
        bool UseRandomlyGeneratedImages;
        int randomGenImageWidth;
        int randomGenImageHeight;
        int alphaTreeAlgorithmCode;
        std::string dissimilarityMetric;
        int bitdepth;
        int tse;
        int connectivity;
        int numitr;
        int iparam1;
        int iparam2;
        int iparam3;
        double fparam1;
        double fparam2;
        double fparam3;
    };
    AlphaTreeConfig(const std::string &filename, const std::string &commentToken = "#");
    std::optional<AlphaTreeParameters> load(int argc, char **argv);

    int getAlphaTreeAlgorithmCode(std::string AlphaTreeAlgorithmName);
    std::string getAlphaTreeAlgorithmName(int AlphaTreeAlgorithmCode);

  private:
    std::unordered_map<std::string, int> AlgorithmNameToCode;
    std::unordered_map<int, std::string> AlgorithmCodeToName;

    std::string filename;
    std::string commentToken;
    std::map<std::string, std::string> config;

    std::string trim(const std::string &str) const;
    std::string getString(const std::string &key) const;
    int getInteger(const std::string &key) const;
    double getDouble(const std::string &key) const;
};

inline AlphaTreeConfig alphatreeConfig(CONFIGFILENAME);

#endif // CONFIGREADER_H