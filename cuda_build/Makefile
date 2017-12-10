CXX = nvcc

CFLAGS = --std=c++11 --resource-usage -Werror cross-execution-space-call -arch=sm_35 -rdc=true -lcudadevrt

SRC = ../src/test_triangles.cu

FCL_PATH = /home/atulya/Documents/eecs587/fcl
INCLUDE_DIRS = -I/usr/include/eigen3 -I$(FCL_PATH)/build/include -I${FCL_PATH}/include
LIB_DIRS = -L$(FCL_PATH)/build/lib

all:
	test

clean:
	rm -I *.exe

test:
	$(CXX) $(CFLAGS) $(INCLUDE_DIRS) $(SRC) -o test_triangles.exe
	./test_triangles.exe
