#!/bin/bash
# run this script from the collidor folder to compile fcl

apt-get install cmake
apt-get install libccd-dev
apt-get install libeigen3-dev

cd ..
git clone https://github.com/flexible-collision-library/fcl.git
cd fcl
mkdir build
cd build
cmake ..
make
