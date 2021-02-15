@echo off
echo Building oneMD in debug mode...
if not exist build\ (
    mkdir build
)
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cd ..
cmake --build build --config Debug