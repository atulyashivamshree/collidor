// @author : Atulya Shivam Shree
// @description : checks triangle intersection on CUDA

#include "compile_CUDA.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "BVH-cuda-inl.h"

using namespace std;

__host__ void help();
__host__ void loadData(string file1, string file2, BVH& bvh1, BVH& bvh2);

__host__ int main(int argc, char *argv[])
{
  if(argc < 3)
  {
    help();
    exit(EXIT_FAILURE);
  }

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  cout << "Device name: " << prop.name << endl;
  
  BVH bvh1, bvh2;
  loadData(argv[1], argv[2], bvh1, bvh2);

  Matrix3 I({{1,0,0}, {0,1,0}, {0,0,1}});
  Vector3 p({0,0,0});
  Config cfg({0.4, 15, I, p, 0});

  DistanceResult result = computeDistance(&bvh1, &bvh2, cfg);

  cout << "==== RESULTS ====" << endl;
  cout << "dist: " << result.dist << " stop " << result.stop << endl;
  cout << "i1: " << result.tsk.i1 << " i2: " << result.tsk.i2 << " d: " << result.tsk.dist << endl;
  cout << "i1: " << result.tsk2.i1 << " i2: " << result.tsk2.i2 << " d: " << result.tsk2.dist << endl;
  cout << "idx: " << result.idx << " idy: " << result.idy <<  endl;
  deleteBVH(&bvh1);
  deleteBVH(&bvh2);
}

__host__ void help()
{
  cout << "USAGE : ./dist_bvh FILE1 FILE2" << endl;
}

__host__ void loadData(string file1, string file2, BVH& bvh1, BVH& bvh2)
{

  ifstream f1, f2;
  f1.open(file1.c_str());
  f2.open(file2.c_str());
  if(!f1.is_open())
  {
    cout << file1 << " could not be opened" << endl;
    exit(EXIT_FAILURE);
  }
  if(!f2.is_open())
  {
    cout << file2 << " could not be opened" << endl;
    exit(EXIT_FAILURE);
  }

  loadBVH(f1, &bvh1);
  loadBVH(f2, &bvh2);
}