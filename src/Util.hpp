#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <chrono>

#include "common.hh"
#include "spdlog/spdlog.h"

using namespace std;

class Util {
    public:
        static string upper(string str) {
            std::transform(str.begin(), str.end(),str.begin(), ::toupper);
            return str;
        }
#ifdef USE_ONEAPI
        static double event_exec_time(const sycl::event &e) {
            double start_k = e.get_profiling_info<sycl::info::event_profiling::command_start>();
            double end_k = e.get_profiling_info<sycl::info::event_profiling::command_end>();
            double kernel_time = (end_k - start_k) * 1e-9; // ns to s
            return kernel_time;
        }
#endif        
};