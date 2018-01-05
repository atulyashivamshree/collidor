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

__host__ void help_message();
__host__ void loadData(string file1, string file2, BVH& bvh1, BVH& bvh2);
__host__ void printResult(ostream& os, const vector<DistanceResult> res);
__host__ void loadTransformations(vector<Transform3f>& transforms, const string filename);

__host__ int main(int argc, char *argv[])
{
  if(argc < 5)
  {
    help_message();
    exit(EXIT_FAILURE);
  }

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  cout << "Device name: " << prop.name << endl;
  
  BVH bvh1, bvh2;
  loadData(argv[1], argv[2], bvh1, bvh2);

  vector<Transform3f> transforms;
  loadTransformations(transforms, argv[3]);

  Config def_cfg({0.04, 15, {1,0,0, 0,1,0, 0,0,1}, {0,0,0}, 1});
  vector<DistanceResult> results = computeDistance(&bvh1, &bvh2, def_cfg, 
                              transforms,
                              string(argv[4]));

  ofstream fout;
  string outfile = argv[4];
  outfile += ".out";
  fout.open(outfile.c_str());
  printResult(fout, results);
  fout.close();

  for(const auto res : results)
  {
    cout << "==== RESULTS ====" << endl;
    cout << "dist: " << res.dist << " stop " << res.stop << endl;
    cout << "i1: " << res.tsk.i1 << " i2: " << res.tsk.i2 << " d: " << res.tsk.dist << endl;
    cout << "i1: " << res.tsk2.i1 << " i2: " << res.tsk2.i2 << " d: " << res.tsk2.dist << endl;
    cout << "idx: " << res.idx << " idy: " << res.idy <<  endl;
  }
  deleteBVH(&bvh1);
  deleteBVH(&bvh2);
}

__host__ void help_message() {
  cout << "Usage : ./compute_dist_bvh FILE1.obj FILE2.obj TRANSFORMS.csv QUEUE_OUTP_PREFIX" << endl;
  cout << " FILE1: input .obj file " << endl;
  cout << " FILE2: input .obj file " << endl;
  cout << " TRANSFORMS: input csv file containing a sample set of transformations" << endl;
  cout << " QUEUE_OUTP_PREFIX: output file for the queues" << endl;
}

__host__ void printResult(ostream& os, const vector<DistanceResult> results)
{
  for(int i = 0; i < results.size(); i++)
    os << "[" << i << "] " << results[i].dist << endl;
}

__host__ void loadTransformations(vector<Transform3f>& transforms, const string filename)
{
  ifstream fin;
  fin.open(filename.c_str());
  if(!fin.is_open())
  {
    help_message();
    exit(EXIT_FAILURE);
  }
  string header_str;
  for(int i = 0; i < 6; i++)
    fin >> header_str;

  float roll, pitch, yaw, x, y, z;

  // read in RPY, XYZ values from the file and store them as a transform
  while(fin >> roll >> pitch >> yaw >> x >> y >> z)
  {
    Eigen::Vector3f pos(x, y, z);

    Eigen::AngleAxisf rollAngle(roll, Eigen::Vector3f::UnitZ());
    Eigen::AngleAxisf yawAngle(pitch, Eigen::Vector3f::UnitY());
    Eigen::AngleAxisf pitchAngle(yaw, Eigen::Vector3f::UnitX());
    Eigen::Quaternion<float> q = rollAngle * yawAngle * pitchAngle;
    Eigen::Matrix3f R = q.matrix();
    Transform3f tf;
    tf.linear() = q.matrix();
    tf.translation() = pos;
    transforms.push_back(tf);
  }
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