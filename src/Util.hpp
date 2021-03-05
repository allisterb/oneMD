#pragma once

#ifdef USE_ONEAPI
// oneDPL headers should be included before standard headers
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/numeric>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/iterator>
#else
#include <execution>
#endif

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <chrono>

using namespace std;

class Util {
    public:
        static string upper(string str) {
            std::transform(str.begin(), str.end(),str.begin(), ::toupper);
            return str;
        }
        
#ifdef USE_ONEAPI
        template <typename T>
        size_t size_2d(const std::vector<std::vector<T>>& v) {
            return std::reduce(dpl::execution::par_unseq, v.begin(), v.end(), size_t(), [](auto a, auto b){ return a.size() + b.size(); });
        }
#else
        template <typename T>
        size_t size_2d(const std::vector<std::vector<T>>& v) {
            return std::reduce(std::execution::par_unseq, v.begin(), v.end(), size_t(), [](auto a, auto b){ return a.size() + b.size(); });
        }
#endif        
};