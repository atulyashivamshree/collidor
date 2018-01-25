/*
 * 	distance_verification.cpp
 *  @Desc: verifies whether the distances computed by GPU between primitive
    triangle elements is same as that computed by the fcl library
 *
 *  Created on: Dec 11, 2017
 *      Author: atulya

                Copyright (c) 2017 Atulya Shivam Shree
 *
 */

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "compile_CPP.h"
#include "../BVH.h"

#include "fcl/geometry/bvh/detail/BVH_front.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/triangle_distance.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "fcl/narrowphase/detail/traversal/distance/mesh_distance_traversal_node.h"
#include "test_fcl_utility.h"

using std::cout;
using std::endl;
using std::string;
using fcl::Transform3;
using fcl::RSSf;

#include "SampleObjects.h"
#include "parse_utils.h"

// takes in a BVH of FCL and converts to BVH of this project
void verifyDistances(string file1, string file2, string transforms_file,
                     string outp_prefix);

void help_message() {
  cout << "Usage : ./verify_dist FILE.yaml" << endl;
  cout << "FILE.yaml : config file for the run" << endl;
}

int main(int argc, char* argv[]) {
  std::vector<fcl::Vector3<float>> points;
  std::vector<fcl::Triangle> triangles;

  if (argc != 2) {
    help_message();
    return -1;
  }

  map<string, string> params;
  loadConfig(params, argv[1]);

  verifyDistances(params["file1_obj"], params["file2_obj"],
                  params["transforms"], params["outp_prefix"]);
}

void verifyDistances(string file1, string file2, string transforms_file,
                     string outp_prefix) {
  std::vector<Transform3f> transforms;
  loadTransformations(transforms, transforms_file);

  // LOAD the input .obj files using FCL library
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

  // CREATE the BVH model using default method
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

  cout << "num elements in m1 " << m1.getNumBVs() << endl;
  cout << "num elements in m2 " << m2.getNumBVs() << endl;

  // RUN the DistanceTraversal node on both trees to evaluate the distances
  fcl::DistanceResult<float> local_result;
  fcl::detail::MeshDistanceTraversalNodeRSSf node;

  for (int i = 0; i < transforms.size(); i++) {
    const Transform3f& tf = transforms[i];

    if (!initialize(node, (const fcl::BVHModel<RSSf>&)m1,
                    Transform3f::Identity(), (const fcl::BVHModel<RSSf>&)m2, tf,
                    fcl::DistanceRequest<float>(true), local_result))
      std::cout << "initialize error" << std::endl;

    // LOAD the outp.csv file
    std::ifstream if3;
    string filename = outp_prefix;
    filename += '_';
    filename += char('0' + i);
    filename += ".outp.csv";
    if3.open(filename.c_str());
    if (!if3.is_open()) {
      cout << filename << " could not be opened" << endl;
      return;
    }

    // PARSE in the outp.csv file and compare values against those obtained from
    // FCL
    cout << std::setprecision(7);

    // LOAD and evaulate each one of the bounding volumes
    cout << "========= Bounding Volumes =========" << endl;
    int num_bv, num_tri;
    string str;
    if3 >> str >> num_bv >> str >> num_tri;
    for (int i = 0; i < num_bv; i++) {
      int iA, iB;
      float dist;
      if3 >> iA >> iB >> dist;
      float correct_dist = node.BVTesting(iA, iB);

      cout << "iA: " << iA << " iB " << iB << " dist: " << dist
           << ", expected: " << correct_dist << endl;
      assert(approxEquals(correct_dist, dist, 1e-4));
    }

    cout << "============ Triangles ============" << endl;
    // LOAD and evaluate each of the triangles
    for (int i = 0; i < num_tri; i++) {
      int iA, iB;
      float dist;
      if3 >> iA >> iB >> dist;
      fcl::Vector3f temp1, temp2;
      int tA = m1.getBV(iA).primitiveId();
      int tB = m2.getBV(iB).primitiveId();
      // cout << "iA: " << iA << " iB "  << iB << " tA "
      // << tA << " tB " << tB << endl;
      float correct_dist = fcl::detail::TriangleDistancef::triDistance(
          p1[t1[tA][0]], p1[t1[tA][1]], p1[t1[tA][2]], p2[t2[tB][0]],
          p2[t2[tB][1]], p2[t2[tB][2]], tf, temp1, temp2);
      cout << "iA: " << iA << " iB " << iB << " dist: " << dist
           << ", expected: " << correct_dist << endl;
      assert(approxEquals(correct_dist, dist, 1e-4));
    }

    cout << "Transform [" << i << "] : PASSED" << endl;
    cout << endl;
  }
  cout << "All Test Transforms PASSED!" << endl;
}
