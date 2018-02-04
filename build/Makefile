# nvcc for cuda files and g++ for others
CUDA_CXX = nvcc
CXX = g++

# setting up the compiler flags
CUDA_CFLAGS = --std=c++11 --resource-usage -arch=sm_35 -rdc=true -lcudadevrt --expt-relaxed-constexpr
CFLAGS = --std=c++11

# MAIN source code for collision
DRIVER_SRC = ../src/dist_bvh.cu
CUDA_SRC = ../src/BVH-cuda-inl.h

# Testing code for verifying whether it works or not
SRC_RSS_TEST = ../src/test/test_rss.cu
SRC_BVH_TEST = ../src/test/test_bvh.cu
SRC_RSS_TEST_CPU = ../src/test/test_rss.cpp
SRC_VER_DIST = ../src/utils/distance_verification.cpp
SRC_VER_QDIST = ../src/utils/queue_verification.cpp
SRC_DISP_OBJ = ../src/utils/display_object.cpp

# Library path for FCL
FCL_PATH = /home/atulya/eecs587/fcl
INCLUDE_DIRS = -I/usr/include/eigen3 -I$(FCL_PATH)/build/include -I${FCL_PATH}/include
LIB_DIRS = -L$(FCL_PATH)/build/lib

# the default config that would be used for test_dist and test_queue
test_config = puck_puck

all: dist_bvh.exe

dist_bvh.exe: $(DRIVER_SRC) $(CUDA_SRC)
	$(CUDA_CXX) $(CUDA_CFLAGS) $(INCLUDE_DIRS) $(LIB_DIRS) $(DRIVER_SRC) -lfcl -o dist_bvh.exe

test : test_bvh.exe test_rss_gpu.exe test_rss_cpu.exe
	./test_bvh.exe ../CAD/tasks0.csv
	./test_rss_cpu.exe
	./test_rss_gpu.exe
	
test_bvh.exe: $(SRC_BVH_TEST)
	$(CUDA_CXX) $(CUDA_CFLAGS) $(INCLUDE_DIRS) $(SRC_BVH_TEST) -o test_bvh.exe

test_rss_gpu.exe: $(SRC_RSS_TEST)
	$(CUDA_CXX) $(CUDA_CFLAGS) $(INCLUDE_DIRS) $(SRC_RSS_TEST) -o test_rss_gpu.exe

test_rss_cpu.exe: $(SRC_RSS_TEST_CPU)
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(LIB_DIRS) $(SRC_RSS_TEST_CPU) -lfcl -o test_rss_cpu.exe
	

verify_dist_bvh.exe: $(SRC_VER_DIST)
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(LIB_DIRS) $(SRC_VER_DIST) -lfcl -o verify_dist_bvh.exe

verify_qdist.exe: $(SRC_VER_QDIST)
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(LIB_DIRS) $(SRC_VER_QDIST) -lfcl -o verify_qdist.exe

display_object.exe: $(SRC_DISP_OBJ)
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(LIB_DIRS) $(SRC_DISP_OBJ) -lfcl -o display_object.exe

test_dist: dist_bvh.exe verify_dist_bvh.exe
	./dist_bvh.exe ../params/$(test_config).yaml
	rm *.outp.csv
	./verify_dist_bvh.exe ../params/$(test_config).yaml

test_queue: dist_bvh.exe verify_qdist.exe
	./dist_bvh.exe ../params/$(test_config).yaml
	./verify_qdist.exe ../params/$(test_config).yaml
	rm *.outp.csv
	
clean:
	rm -f *.exe
	rm -f *.outp.csv
	rm -f *.out
