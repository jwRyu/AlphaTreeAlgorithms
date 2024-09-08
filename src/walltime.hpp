#pragma once
#include <defines.hpp>

#ifdef _WIN32
#include <Windows.h>
double get_wall_time();
double get_cpu_time();

//  Posix/Linux
#else
double get_wall_time();
double get_cpu_time();
#endif
