#pragma once

#include <algorithm>
#include <string>
#include <chrono>
using namespace std;

class Util {
    public:
        static string upper(string str) {
            std::transform(str.begin(), str.end(),str.begin(), ::toupper);
            return str;
        }
};