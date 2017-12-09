CXX = nvcc

CFLAGS = --resource-usage -Werror cross-execution-space-call -arch=sm_35 -rdc=true -lcudadevrt

SRC = stencil_2d.cu

all:
	$(CXX) $(CFLAGS) $(SRC) -o stencil.exe

clean:
	rm -I *.exe
	rm -I *.csv
