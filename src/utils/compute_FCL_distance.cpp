/*
 * compute_FCL_distance.cpp
 *
 *  Created on: Dec 12, 2017
 *      Author: atulya
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "fcl/geometry/bvh/detail/BVH_front.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "../test_fcl_utility.h"

using namespace std;
using namespace fcl;

#include "SampleObjects.h"

template <typename S>
void test_mesh_distance(string file1, string file2,
		vector<Transform3<S>> transforms,
		vector<S>& distance);

template<typename BV, typename TraversalNode>
void distance_Test_Oriented(const fcl::Transform3<typename BV::S>& tf,
                            const std::vector<fcl::Vector3<typename BV::S>>& vertices1, const std::vector<Triangle>& triangles1,
                            const std::vector<fcl::Vector3<typename BV::S>>& vertices2, const std::vector<Triangle>& triangles2, detail::SplitMethodType split_method,
                            int qsize,
                            fcl::test::DistanceRes<typename BV::S>& distance_result,
                            bool verbose);

int main(int argc, char *argv[])
{
	vector<Transform3f> transforms;
	vector<fcl::Vector3<float>> vertices1, vertices2;
	vector<fcl::Triangle> triangles1, triangles2;
	vector<float> distances;

	string file1, file2;

	if(argc == 3)
	{
		file1 = argv[1];
		file2 = argv[2];
	}
	test_mesh_distance(file1, file2, transforms, distances);

	for(int i = 0; i < distances.size(); i++)
		cout << i << " " << distances[i] << endl;
}

template <typename S>
void test_mesh_distance(string file1, string file2,
		vector<Transform3<S>> transforms,
		vector<S>& distance)
{
  std::vector<fcl::Vector3<S>> p1, p2;
  std::vector<fcl::Triangle> t1, t2;

  if(file1 != "")
  {
	  fcl::test::loadOBJFile(file1.c_str(), p1, t1);
	  fcl::test::loadOBJFile(file2.c_str(), p2, t2);
  }
  else
  {
	  createSampleObj1(p1, t1);
	  createSampleObj2(p2, t2);
  }

  Vector3<S> position(0.0, 0.0, 0.0);
  transforms.push_back(Transform3<S>::Identity());
  transforms[0].linear() = Matrix3<S>::Identity();
  transforms[0].translation() = position;

  double dis_time = 0;

  test::DistanceRes<S> res, res_now;
  for(std::size_t i = 0; i < transforms.size(); ++i)
  {
    test::Timer timer_dist;
    timer_dist.start();
    distance_Test_Oriented<RSS<S>, fcl::detail::MeshDistanceTraversalNodeRSS<S>>(transforms[i], p1, t1, p2, t2, detail::SPLIT_METHOD_MEAN, 2, res, false);
    timer_dist.stop();
    dis_time += timer_dist.getElapsedTimeInSec();

    distance.push_back(res.distance);

  }
//  std::cout << "distance timing: " << dis_time << " sec" << std::endl;
}

template<typename BV, typename TraversalNode>
void distance_Test_Oriented(const fcl::Transform3<typename BV::S>& tf,
                            const std::vector<fcl::Vector3<typename BV::S>>& vertices1, const std::vector<Triangle>& triangles1,
                            const std::vector<fcl::Vector3<typename BV::S>>& vertices2, const std::vector<Triangle>& triangles2, detail::SplitMethodType split_method,
                            int qsize,
                            fcl::test::DistanceRes<typename BV::S>& distance_result,
                            bool verbose)
{
  using S = typename BV::S;

  fcl::BVHModel<BV> m1;
  fcl::BVHModel<BV> m2;
  m1.bv_splitter.reset(new fcl::detail::BVSplitter<BV>(split_method));
  m2.bv_splitter.reset(new fcl::detail::BVSplitter<BV>(split_method));


  m1.beginModel();
  m1.addSubModel(vertices1, triangles1);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(vertices2, triangles2);
  m2.endModel();

  fcl::DistanceResult<S> local_result;
  TraversalNode node;
  if(!initialize(node, (const fcl::BVHModel<BV>&)m1, tf,
		  (const fcl::BVHModel<BV>&)m2, fcl::Transform3<S>::Identity(),
		  fcl::DistanceRequest<S>(true), local_result))
    std::cout << "initialize error" << std::endl;

  node.enable_statistics = verbose;

  distance(&node, nullptr, qsize);

  // points are in local coordinate, to global coordinate
  Vector3<S> p1 = local_result.nearest_points[0];
  Vector3<S> p2 = local_result.nearest_points[1];


  distance_result.distance = local_result.min_distance;
  distance_result.p1 = p1;
  distance_result.p2 = p2;

  if(verbose)
  {
    std::cout << "distance " << local_result.min_distance << std::endl;

    std::cout << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
    std::cout << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
    std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
  }
}


