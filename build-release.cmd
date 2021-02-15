@echo off
echo Building oneMD in release mode...
if not exist build\ (
    mkdir build
)
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..
cmake --build build --config Release