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

test:
	./dist_bvh.exe ../params/block_puck.yaml
	# rm block_puck_*.outp.csv
	cat block_puck_0.outp.csv
	# ../build/compute_dist_bvh ../params/block_puck.yaml

test_puck_puck:
	./dist_bvh.exe ../params/puck_puck_t0.yaml
	# cat puck_puck_0.outp.csv
	# echo "============"
	# cat puck_puck_1.outp.csv
	rm puck_puck_*.outp.csv
	../build/compute_dist_bvh ../params/puck_puck_t0.yaml

test_block_puck:
	./dist_bvh.exe ../params/block_puck.yaml
	rm block_puck_*.outp.csv
	../build/compute_dist_bvh ../params/block_puck.yaml

test_block_s1:
	./dist_bvh.exe ../params/block_s1_t0.yaml
	# rm block_s1_*.outp.csv
	../build/compute_dist_bvh ../params/block_s1_t0.yaml

test_s1_s2:
	./dist_bvh.exe ../params/s1s2_t0.yaml
	# rm s1s2_*.outp.csv
	../build/compute_dist_bvh ../params/s1s2_t0.yaml

test_block_puck_queue:
	./dist_bvh.exe ../params/block_puck.yaml
	../build/verify_dist ../params/block_puck.yaml
	rm block_puck_*.outp.csv

test_s1_s2_queue:
	./dist_bvh.exe ../params/s1s2.yaml
	../build/verify_dist ../params/s1s2.yaml
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
