#! /bin/bash
set -e
echo Building oneMD in release mode...
if [ ! -d "build" ] 
then
    mkdir build
fi
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..
cmake --build build --config Release