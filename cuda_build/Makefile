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

test_block_puck:
	./compute_dist_bvh ../CAD/block.obj ../CAD/puck.obj ../CAD/transforms_set1.csv ../cuda_build/block_puck.out

test_s1_s2:
	./compute_dist_bvh _ _ ../CAD/transforms_set1.csv ../cuda_build/s1s2.out

test_block_puck_queue:
	./dist_bvh.exe ../CAD/block.bvh ../CAD/puck.bvh ../CAD/transforms_set2.csv block_puck
	../build/verify_dist ../CAD/block.obj ../CAD/puck.obj ../CAD/transforms_set2.csv block_puck
	rm block_puck_*.outp.csv

test_s1_s2_queue:
	./dist_bvh.exe ../CAD/sample1.bvh ../CAD/sample2.bvh ../CAD/transforms_set1.csv s1s2
	../build/verify_dist _ _ ../CAD/transforms_set1.csv s1s2
	rm s1s2_*.outp.csv
	
clean:
	rm -I *.exe
	rm -I *.outp.csv

dist_bvh: $(SRC3) $(SRC4)
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC3) -o dist_bvh.exe

test_triangles:
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC1) -o test_triangles.exe
	./test_triangles.exe

test_rss:
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC2) -o test_rss.exe
	./test_rss.exe
