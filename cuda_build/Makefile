CXX = nvcc

CFLAGS = -g --std=c++11 --resource-usage -Werror cross-execution-space-call -arch=sm_35 -rdc=true -lcudadevrt

SRC1 = ../src/test_triangles.cu
SRC2 = ../src/test_rss.cu

FCL_PATH = /home/atulya/Documents/eecs587/fcl
INCLUDE_DIRS = -I/usr/include/eigen3 -I$(FCL_PATH)/build/include -I${FCL_PATH}/include
LIB_DIRS = -L$(FCL_PATH)/build/lib

all: test

clean:
	rm -I *.exe

test_triangles:
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC1) -o test_triangles.exe
	./test_triangles.exe

test:
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC2) -o test_rss.exe
	./test_rss.exe
