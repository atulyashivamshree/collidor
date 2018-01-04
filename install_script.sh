sudo apt-get install libccd-dev
git clone https://github.com/flexible-collision-library/fcl.git
cd fcl
mkdir build
cd build
cmake ..
make
cd ../..
git clone git://github.com/OctoMap/octomap.git
cd octomap
sudo apt-get install build-essential cmake doxygen libqt4-dev \
	libqt4-opengl-dev libqglviewer-dev-qt4
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
make install
cd ../../collidor
mkdir -p build
cd build
cmake ../src
make
