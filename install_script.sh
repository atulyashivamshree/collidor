sudo apt-get install libccd-dev
git clone https://github.com/flexible-collision-library/fcl.git
cd fcl
mkdir build
cd build
cmake ..
make
cd ../../collidor
mkdir -p build
cd build
cmake ../src
make
