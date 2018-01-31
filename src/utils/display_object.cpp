/*
 *  fcl_cdp_conversions.cpp
    Description: peforms conversion between FCL and the custom BVH versions

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "compile_CPP.h"
#include "../BVH.h"

#include "fcl/geometry/bvh/detail/BVH_front.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "fcl_utility.h"

using std::cout;
using std::endl;
using std::vector;

#include "SampleObjects.h"

// takes in a BVH of FCL and converts to BVH of this project

void testCONV(const BVH* bvh);

int main(int argc, char* argv[]) {
  std::vector<fcl::Vector3<float>> points;
  std::vector<fcl::Triangle> triangles;

  if (argc < 2)
    createSampleObj2(points, triangles);
  else
    fcl::test::loadOBJFile(argv[1], points, triangles);

  fcl::BVHModel<fcl::RSSf> model;

  model.bv_splitter.reset(
      new fcl::detail::BVSplitter<fcl::RSSf>(fcl::detail::SPLIT_METHOD_MEAN));

  model.beginModel();
  model.addSubModel(points, triangles);
  model.endModel();

  BVH bvh_cdp;
  convertFCLToBVH(model, points, triangles, &bvh_cdp);
  printBVH(cout, &bvh_cdp);

  if (argc == 1) testCONV(&bvh_cdp);

  deleteBVH(&bvh_cdp);
}

bool approx_Equals(float a, float b, float EPSILON = 5e-7) {
  if (fabs(a - b) < EPSILON) return true;
  return false;
}

bool operator!=(const Vector3& lhs, const Vector3& rhs) {
  if (!approx_Equals(lhs.x, rhs.x)) return true;
  if (!approx_Equals(lhs.y, rhs.y)) return true;
  if (!approx_Equals(lhs.z, rhs.z)) return true;

  return false;
}

bool operator!=(const RSS& lhs, const RSS& rhs) {
  if (lhs.axis.v1 != rhs.axis.v1) return true;
  if (lhs.axis.v2 != rhs.axis.v2) return true;
  if (lhs.axis.v3 != rhs.axis.v3) return true;

  if (lhs.To != rhs.To) return true;

  if (!approx_Equals(lhs.l[0], rhs.l[0])) return true;
  if (!approx_Equals(lhs.l[1], rhs.l[1])) return true;

  if (!approx_Equals(lhs.r, rhs.r)) return true;
  if (!approx_Equals(lhs.size, rhs.size, 1e-4)) return true;

  return false;
}

bool operator==(const Triangle& lhs, const Triangle& rhs) {
  if (lhs.a != rhs.a) return false;
  if (lhs.b != rhs.b) return false;
  if (lhs.c != rhs.c) return false;

  return true;
}

bool operator==(const BV& lhs, const BV& rhs) {
  if (lhs.rss != rhs.rss) return false;
  if (lhs.id1 != rhs.id1) return false;
  if (lhs.id2 != rhs.id2) return false;
  if (lhs.idt != rhs.idt) return false;
  return true;
}

void testCONV(const BVH* bvh) {
  std::stringstream oss, iss;
  printBVH(oss, bvh);

  iss.str(oss.str());
  BVH read_bvh;
  loadBVH(iss, &read_bvh);

  assert(read_bvh.num_bv = bvh->num_bv);
  assert(read_bvh.num_tri = bvh->num_tri);
  for (int i = 0; i < bvh->num_bv; i++)
    assert(read_bvh.bv_arr[i] == bvh->bv_arr[i]);

  for (int j = 0; j < bvh->num_tri; j++)
    assert(read_bvh.tri_arr[j] == bvh->tri_arr[j]);
}
