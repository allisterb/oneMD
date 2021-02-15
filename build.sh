#!/bin/bash
set -e
git fetch
git pull
cd build
cmake ..
make
cd ..
