/*
 * fcl_cdp_conversions.cpp
 *
 *  Created on: Dec 11, 2017
 *      Author: atulya
 */

#include "../compile_CPP.h"
#include "../BVH.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>

#include "fcl/geometry/bvh/detail/BVH_front.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/triangle_distance.h"
#include "../test_fcl_utility.h"

using namespace std;
using namespace fcl;

#include "SampleObjects.h"


// takes in a BVH of FCL and converts to BVH of this project
void verifyDistances(string file1, string file2, string node_distances_file);

bool approxEquals(float a, float b, float EPSILON = 1e-5) {
	return (fabs(a - b) < EPSILON);
}

int main(int argc, char *argv[])
{
  std::vector<fcl::Vector3<float>> points;
  std::vector<fcl::Triangle> triangles;

  if(argc != 4 && argc != 2)
    cout << "Usage : ./verify_distances.exe FILE1 FILE2 FOUT_DIST" << endl;

  if(argc == 2)
	  verifyDistances("", "", argv[1]);
  else
	  verifyDistances(argv[1], argv[2], argv[3]);
}

void verifyDistances(string file1, string file2, string node_distances_file)
{
	ifstream if1, if2, if3;
	if1.open(file1.c_str());
	if(!if1.is_open())
		cout << file1 << " could not b opened" << endl;
	if2.open(file2.c_str());
	if(!if2.is_open())
		cout << file2 << " could not b opened" << endl;
	if3.open(node_distances_file.c_str());
	if(!if3.is_open())
		cout << node_distances_file << " could not b opened" << endl;

	std::vector<fcl::Vector3f> p1, p2;
	std::vector<fcl::Triangle> t1, t2;

	if(file1 == "" || file2 == "")
	{
		createSampleObj1(p1, t1);
		createSampleObj2(p2, t2);
	}
	else
	{
		fcl::test::loadOBJFile(file1.c_str(), p1, t1);
		fcl::test::loadOBJFile(file2.c_str(), p2, t2);
	}

	fcl::BVHModel<RSSf> m1;
	fcl::BVHModel<RSSf> m2;
	m1.bv_splitter.reset(new fcl::detail::BVSplitter<RSSf>(fcl::detail::SPLIT_METHOD_MEAN));
	m2.bv_splitter.reset(new fcl::detail::BVSplitter<RSSf>(fcl::detail::SPLIT_METHOD_MEAN));

	DistanceResultf local_result;

	m1.beginModel();
	m1.addSubModel(p1, t1);
	m1.endModel();

	m2.beginModel();
	m2.addSubModel(p2, t2);
	m2.endModel();

	cout << "num elements in m1 " << m1.getNumBVs() << endl;
	cout << "num elements in m2 " << m2.getNumBVs() << endl;

	detail::MeshDistanceTraversalNodeRSSf node;
	if(!initialize(node, (const BVHModel<RSSf>&)m1, Transform3f::Identity(), (const BVHModel<RSSf>&)m2, Transform3<float>::Identity(), DistanceRequestf(true), local_result))
		std::cout << "initialize error" << std::endl;

	int num_bv, num_tri;
	string str;
	if3 >> str >> num_bv >> str >> num_tri;
	for(int i = 0;i < num_bv; i++)
	{
		int iA, iB; float dist;
		if3 >> iA >> iB >> dist;
		float correct_dist = node.BVTesting(iA, iB);
//		m1.getBV(1).bv.
		cout << "iA: " << iA << " iB " << iB << " dist: " << dist <<", expected: " <<correct_dist << endl;
		assert(approxEquals(correct_dist, dist, 1e-5));

	}
	for(int i = 0;i < num_tri; i++)
	{
		int iA, iB; float dist;
		if3 >> iA >> iB >> dist;
		fcl::Vector3f temp1, temp2;
		int tA = m1.getBV(iA).primitiveId();
		int tB = m2.getBV(iB).primitiveId();
//		cout << "iA: " << iA << " iB " << iB << " tA " << tA << " tB " << tB << endl;
		float correct_dist = detail::TriangleDistancef::triDistance(p1[t1[tA][0]], p1[t1[tA][1]], p1[t1[tB][2]],
																	p2[t2[tB][0]], p2[t2[tB][1]], p2[t2[tB][2]], temp1, temp2);
		cout << "iA: " << iA << " iB " << iB << " dist: " << dist <<", expected: "<< correct_dist << endl;
		assert(approxEquals(correct_dist, dist, 1e-5));
	}

	cout << "test PASSED!" << endl;

}

