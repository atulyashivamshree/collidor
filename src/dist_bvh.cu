// @author : Atulya Shivam Shree
// @description : checks triangle intersection on CUDA

#include "compile_CUDA.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

#include "BVH-cuda-inl.h"
#include "utils/parse_utils.h"

using namespace std;


__host__ void loadBVHData(BVH& bvh, const string filename);
__host__ void help_message();
__host__ void printResult(ostream& os, const vector<DistanceResult> res, const vector<float> elap_time);

__host__ int main(int argc, char *argv[])
{
  if(argc < 2)
  {
    help_message();
    exit(EXIT_FAILURE);
  }

  // Print device properties
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  cout << "Device name: " << prop.name << endl;

  // Load all configuration params from the file
  Config def_cfg({0.04, 15, {1,0,0, 0,1,0, 0,0,1}, {0,0,0}, 1});
  map<string, string> params;
  loadConfig(params, argv[1]);
  def_cfg.gamma = std::stod(params["gamma"]);
  def_cfg.enable_distance_reduction = std::stoi(params["compute_min_dist"]);
  def_cfg.max_iter = std::stoi(params["max_iter"]);
  def_cfg.max_bv_proc = std::stoi(params["max_bv_proc"]);
  
  // Load the two BVH in consideration
  BVH bvhA, bvhB;
  loadBVHData(bvhA, params["file1_bvh"]);
  loadBVHData(bvhB, params["file2_bvh"]);

  // Load the transformations
  vector<Transform3f> transforms;
  loadTransformations(transforms, params["transforms"]);

  // Compute distances and store the result
  vector<float> elap_time;
  vector<DistanceResult> results = computeDistance(&bvhA, &bvhB, def_cfg, 
                              transforms,
                              params["outp_prefix"],
                              elap_time);

  // Store the entire result in an output file
  ofstream fout;
  string outfile = params["outp_prefix"];
  outfile += ".out";
  fout.open(outfile.c_str());
  printResult(fout, results, elap_time);
  fout.close();

  // Print the final stats of the results
  for(const auto res : results)
  {
    cout << "==== RESULTS ====" << endl;
    cout << "dist: " << res.dist << " stop: " << res.stop << " num_iter: " << res.num_iter << endl;
    cout << "i1: " << res.tsk.i1 << " i2: " << res.tsk.i2 << " d: " << res.tsk.dist << endl;
    cout << "i1: " << res.tsk2.i1 << " i2: " << res.tsk2.i2 << " d: " << res.tsk2.dist << endl;
    cout << "idx: " << res.idx << " idy: " << res.idy <<  endl;
  }

  // Free up memory
  deleteBVH(&bvhA);
  deleteBVH(&bvhB);
}

__host__ void help_message() {
  cout << "Usage : ./dist_bvh.exe FILE.yaml" << endl;
  cout << "FILE.yaml : config file for the run" << endl;
}

__host__ void printResult(ostream& os, const vector<DistanceResult> results, 
                          const vector<float> elap_time)
{
  os << setprecision(7);
  for(int i = 0; i < results.size(); i++)
    os << "[" << i << "] " << results[i].dist << " " << elap_time[i] << endl;
}

HOST_PREFIX void loadBVHData(BVH& bvh, const string filename)
{
  std::ifstream fin;
  fin.open(filename.c_str());
  if(!fin.is_open())
  {
    cout << "BVH: " << filename << " could not be opened" << endl;
    exit(EXIT_FAILURE);
  }
  
  loadBVH(fin, &bvh);
}