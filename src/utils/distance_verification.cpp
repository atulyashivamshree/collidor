/*
 *  distance_verification.cpp
    Description: Checks the working of dsitance between objects by comparing GPU
 based rresults with those obtained from the FCL algorithms on a CPU

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <sys/time.h>
#include <time.h>

#include "Eigen/Dense"
#include "fcl/geometry/bvh/detail/BVH_front.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "fcl/narrowphase/detail/traversal/distance/mesh_distance_traversal_node.h"

#include "compile_CPP.h"
#include "fcl_utility.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

using fcl::RSSf;
using fcl::Transform3;

#include "SampleObjects.h"
#include "parse_utils.h"

#define EPSILON 1e-5

template <typename S>
void test_mesh_distance(string file1, string file2,
                        vector<Transform3f> transforms, vector<S>& distance,
                        vector<float>& elap_time);

template <typename BV, typename TraversalNode>
void distance_Test_Oriented(
    const fcl::Transform3<typename BV::S>& tf, fcl::BVHModel<BV>& m1,
    fcl::BVHModel<BV>& m2, int qsize,
    fcl::test::DistanceRes<typename BV::S>& distance_result,
    vector<float>& elap_time, bool verbose);

double get_wall_time();

void validateDistances(string gpu_outfile, vector<float> distances,
                       vector<float> cpu_time);

void help_message();

int main(int argc, char* argv[]) {
  vector<Transform3f> transforms;
  vector<float> distances;
  vector<float> elap_time;

  if (argc != 2) {
    help_message();
    return -1;
  }

  map<string, string> params;
  loadConfig(params, argv[1]);

  loadTransformations(transforms, params["transforms"]);
  test_mesh_distance(params["file1"], params["file2"], transforms, distances,
                     elap_time);

  validateDistances(params["outp_prefix"] + ".out", distances, elap_time);
}

void validateDistances(string gpu_outfile, vector<float> distances,
                       vector<float> cpu_time) {
  cout << std::setprecision(7);
  std::ifstream fin;
  fin.open(gpu_outfile.c_str());
  if (!fin.is_open()) {
    cout << gpu_outfile << " could not be opened" << endl;
    exit(EXIT_FAILURE);
  }
  int count_inequality = 0;

  float tot_cpu_time = 0;
  float tot_gpu_time = 0;

  for (int i = 0; i < distances.size(); i++) {
    string tmp_id;
    float dist, del_t;
    fin >> tmp_id >> dist >> del_t;
    cout << "[" << i << "]"
         << "EXP: " << distances[i] << ", ACTUAL: " << dist
         << ", cpu dt: " << cpu_time[i] << "s"
         << ", gpu dt: " << del_t << "s";

    tot_gpu_time += del_t;
    tot_cpu_time += cpu_time[i];

    if (!approxEquals(distances[i], dist,
                      fmax(EPSILON, EPSILON * distances[i]))) {
      count_inequality++;
      cout << "DIFF";
    }
    cout << endl;
  }

  if (count_inequality == 0)
    cout << "Object Distance test : PASSED " << endl;
  else
    cout << "Object Distance test : FAILED " << endl;

  cout << "Total CPU time : " << tot_cpu_time << endl;
  cout << "Total GPU time : " << tot_gpu_time << endl;

  fin.close();
}

void help_message() {
  cout << "Usage : ./compute_dist_bvh FILE.yaml" << endl;
  cout << "FILE.yaml : config file for the run" << endl;
}

template <typename S>
void test_mesh_distance(string file1, string file2,
                        vector<Transform3f> transforms, vector<S>& distance,
                        vector<float>& elap_time) {
  // LOAD the two OBJ files or the default versions of the files
  std::vector<fcl::Vector3f> p1, p2;
  std::vector<fcl::Triangle> t1, t2;

  if (file1 == "_")
    createSampleObj1(p1, t1);
  else
    fcl::test::loadOBJFile(file1.c_str(), p1, t1);

  if (file2 == "_")
    createSampleObj2(p2, t2);
  else
    fcl::test::loadOBJFile(file2.c_str(), p2, t2);

  fcl::BVHModel<RSSf> m1;
  fcl::BVHModel<RSSf> m2;
  m1.bv_splitter.reset(
      new fcl::detail::BVSplitter<RSSf>(fcl::detail::SPLIT_METHOD_MEAN));
  m2.bv_splitter.reset(
      new fcl::detail::BVSplitter<RSSf>(fcl::detail::SPLIT_METHOD_MEAN));

  m1.beginModel();
  m1.addSubModel(p1, t1);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(p2, t2);
  m2.endModel();

  // RUN the distance computation tests on all the different transforms for the
  // two objects
  fcl::test::DistanceRes<S> res, res_now;
  for (std::size_t i = 0; i < transforms.size(); ++i) {
    distance_Test_Oriented<fcl::RSS<S>,
                           fcl::detail::MeshDistanceTraversalNodeRSS<S>>(
        transforms[i], m1, m2, 2, res, elap_time, false);
    distance.push_back(res.distance);
  }
}

template <typename BV, typename TraversalNode>
void distance_Test_Oriented(
    const fcl::Transform3<typename BV::S>& tf, fcl::BVHModel<BV>& m1,
    fcl::BVHModel<BV>& m2, int qsize,
    fcl::test::DistanceRes<typename BV::S>& distance_result,
    vector<float>& elap_time, bool verbose) {
  using S = typename BV::S;

  // START timer to measure the distance computation time
  double t_start = get_wall_time();

  fcl::DistanceResult<S> local_result;
  TraversalNode node;
  if (!initialize(node, (const fcl::BVHModel<BV>&)m1,
                  fcl::Transform3f::Identity(), (const fcl::BVHModel<BV>&)m2,
                  tf, fcl::DistanceRequest<S>(true), local_result))
    std::cout << "initialize error" << std::endl;

  node.enable_statistics = verbose;

  distance(&node, nullptr, qsize);

  // STOP timer
  elap_time.push_back(get_wall_time() - t_start);

  // points are in local coordinate, to global coordinate
  fcl::Vector3f p1 = local_result.nearest_points[0];
  fcl::Vector3f p2 = local_result.nearest_points[1];

  distance_result.distance = local_result.min_distance;
  distance_result.p1 = p1;
  distance_result.p2 = p2;

  if (verbose) {
    std::cout << "distance " << local_result.min_distance << std::endl;

    std::cout << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
    std::cout << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
    std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
  }
}

double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
