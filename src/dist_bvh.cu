/*
 *  dist_bvh.cu
    @description : checks triangle intersection on CUDA

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#include "compile_CUDA.h"

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

#include "BVH-cuda-inl.h"
#include "utils/parse_utils.h"
#include "utils/fcl_utility.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

__host__ void help_message();
__host__ void printResult(std::ostream& os, const vector<DistanceResultGPU> res,
                          const vector<float> elap_time);

__host__ int main(int argc, char* argv[]) {
  if (argc < 2) {
    help_message();
    exit(EXIT_FAILURE);
  }

  // Print device properties
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  cout << "Device name: " << prop.name << endl;

  // Load all configuration params from the file
  Config def_cfg({0.04, 15, {1, 0, 0, 0, 1, 0, 0, 0, 1}, {0, 0, 0}, 1});
  map<string, string> params;
  loadConfig(params, argv[1]);
  def_cfg.gamma = std::stod(params["gamma"]);
  def_cfg.enable_distance_reduction = std::stoi(params["compute_min_dist"]);
  def_cfg.max_iter = std::stoi(params["max_iter"]);
  def_cfg.max_bfs_proc = std::stoi(params["max_bfs_proc"]);
  def_cfg.max_dfs_proc = std::stoi(params["max_dfs_proc"]);

  // Load the two BVH in consideration
  BVH bvhA, bvhB;
  loadOBJToBVH(params["file1"], &bvhA);
  loadOBJToBVH(params["file2"], &bvhB);

  // Load the transformations
  vector<Transform3f> transforms;
  loadTransformations(transforms, params["transforms"]);

  // Compute distances and store the result
  vector<float> elap_time;
  vector<DistanceResultGPU> results = computeDistance(
      &bvhA, &bvhB, def_cfg, transforms, params["outp_prefix"], elap_time);

  // Store the entire result in an output file
  std::ofstream fout;
  string outfile = params["outp_prefix"];
  outfile += ".out";
  fout.open(outfile.c_str());
  printResult(fout, results, elap_time);
  fout.close();

  // Print the final stats of the results
  for (const auto res : results) {
    cout << "==== RESULTS ====" << endl;
    cout << "dist: " << res.dist << " stop: " << res.stop
         << " num_iter: " << res.num_iter << endl;
    cout << "i1: " << res.tsk.i1 << " i2: " << res.tsk.i2
         << " d: " << res.tsk.dist << endl;
    cout << "i1: " << res.tsk2.i1 << " i2: " << res.tsk2.i2
         << " d: " << res.tsk2.dist << endl;
    cout << "idx: " << res.idx << " idy: " << res.idy << endl;
  }

  // Free up memory
  deleteBVH(&bvhA);
  deleteBVH(&bvhB);
}

__host__ void help_message() {
  cout << "Usage : ./dist_bvh.exe FILE.yaml" << endl;
  cout << "FILE.yaml : config file for the run" << endl;
}

__host__ void printResult(std::ostream& os,
                          const vector<DistanceResultGPU> results,
                          const vector<float> elap_time) {
  os << std::setprecision(7);
  for (int i = 0; i < results.size(); i++)
    os << "[" << i << "] " << results[i].dist << " " << elap_time[i] << endl;
}
