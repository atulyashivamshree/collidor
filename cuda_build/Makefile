CXX = nvcc

CFLAGS = -g --std=c++11 --resource-usage -Werror cross-execution-space-call -arch=sm_35 -rdc=true -lcudadevrt

SRC1 = ../src/test_triangles.cu
SRC2 = ../src/test_rss.cu
SRC3 = ../src/dist_bvh.cu
SRC4 = ../src/BVH-cuda-inl.h

FCL_PATH = /home/atulya/Documents/eecs587/fcl
INCLUDE_DIRS = -I/usr/include/eigen3 -I$(FCL_PATH)/build/include -I${FCL_PATH}/include
LIB_DIRS = -L$(FCL_PATH)/build/lib

all: dist_bvh

test: test_triangles test_rss

clean:
	rm -I *.exe

dist_bvh: $(SRC3) $(SRC4)
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC3) -o dist_bvh.exe
	./dist_bvh.exe ../CAD/sample1.bvh ../CAD/sample2.bvh

test_triangles:
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC1) -o test_triangles.exe
	./test_triangles.exe

test_rss:
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC2) -o test_rss.exe
	./test_rss.exe
