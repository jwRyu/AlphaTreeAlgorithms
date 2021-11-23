#pragma once

#include <string>
#include <iostream>
#include <chrono>
#include "../pmt_common.h"

NAMESPACE_PMT

class Timer
{
public:
    Timer() : stopped_(false)
    {
        start();  
    }

    Timer(std::string msg) :
        msg_(msg), stopped_(false)
    {
        start();
    }

    void start()
    {
        start_ = std::chrono::high_resolution_clock::now();
    }

    double stop()
    {
        stopped_ = true;

        std::chrono::duration<double> diff = std::chrono::high_resolution_clock::now() - start_;
 
        return diff.count();
    }

    ~Timer()
    {
        if (stopped_) return;
        
        std::cout << msg_ << " :: " << std::fixed << stop() << " s\n";
    }

private:
    std::string msg_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    bool stopped_;
};



NAMESPACE_PMT_END
