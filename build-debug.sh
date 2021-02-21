#! /bin/bash
set -e
echo Building oneMD in debug mode...
if [ ! -d "build" ] 
then
    mkdir build
fi
cd build
cmake -DCMAKE_BUILD_TYPE=Debug $@ ..
cd ..
cmake --build build --config Debug 