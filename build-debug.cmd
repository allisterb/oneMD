@echo off
git fetch
git pull
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cd ..
cmake --build build --config Debug