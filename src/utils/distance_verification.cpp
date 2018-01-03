/*
 * distance_verification.cpp
 *
 *  Created on: Dec 11, 2017
 *      Author: atulya
 *
 *  @Desc: verifies whether the distances computed by GPU is same as that
 *  		computed by the fcl library
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

bool approxEquals(float a, float b, float EPSILON = 1e-4) {
	return (fabs(a - b) < EPSILON);
}

void help_message() {
	cout << "Usage : ./verify_distances.exe FILE1.obj FILE2.obj FOUT_DIST.outp.csv" << endl;
	cout << " FILE1: input .obj file " << endl;
	cout << " FILE2: input .obj file " << endl;
	cout << " FOUT_DIST: the .outp.csv file. This contains the set of BVs and triangles that "
			"were computed by the GPU algorithm and their values for the distances in between  " << endl;
}

int main(int argc, char *argv[])
{
  std::vector<fcl::Vector3<float>> points;
  std::vector<fcl::Triangle> triangles;

  if(argc != 4) {
	  help_message();
	  return -1;
  }

  verifyDistances(argv[1], argv[2], argv[3]);
}

void verifyDistances(string file1, string file2, string node_distances_file)
{
	// LOAD the outp.csv file
	ifstream if3;
	if3.open(node_distances_file.c_str());
	if(!if3.is_open())
	{
		cout << node_distances_file << " could not be opened" << endl;
		return;
	}

	// LOAD the input .obj files using FCL library
	std::vector<fcl::Vector3f> p1, p2;
	std::vector<fcl::Triangle> t1, t2;

	if(file1 == "_")
		createSampleObj1(p1, t1);
	else
		fcl::test::loadOBJFile(file1.c_str(), p1, t1);

	if(file2 == "_")
		createSampleObj2(p2, t2);
	else
		fcl::test::loadOBJFile(file2.c_str(), p2, t2);

	// CREATE the BVH model using default method
	fcl::BVHModel<RSSf> m1;
	fcl::BVHModel<RSSf> m2;
	m1.bv_splitter.reset(new fcl::detail::BVSplitter<RSSf>(fcl::detail::SPLIT_METHOD_MEAN));
	m2.bv_splitter.reset(new fcl::detail::BVSplitter<RSSf>(fcl::detail::SPLIT_METHOD_MEAN));

	m1.beginModel();
	m1.addSubModel(p1, t1);
	m1.endModel();

	m2.beginModel();
	m2.addSubModel(p2, t2);
	m2.endModel();

	cout << "num elements in m1 " << m1.getNumBVs() << endl;
	cout << "num elements in m2 " << m2.getNumBVs() << endl;

	// RUN the DistanceTraversal node on both trees to evaluate the distances
	DistanceResultf local_result;
	detail::MeshDistanceTraversalNodeRSSf node;
	if(!initialize(node, (const BVHModel<RSSf>&)m1, Transform3f::Identity(), (const BVHModel<RSSf>&)m2, Transform3<float>::Identity(), DistanceRequestf(true), local_result))
		std::cout << "initialize error" << std::endl;

	// PARSE in the outp.csv file and compare values against those obtained from FCL
	cout << setprecision(7);

	// LOAD and evaulate each one of the bounding volumes
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
		assert(approxEquals(correct_dist, dist, 1e-4));

	}

	// LOAD and evaluate each of the triangles
	for(int i = 0;i < num_tri; i++)
	{
		int iA, iB; float dist;
		if3 >> iA >> iB >> dist;
		fcl::Vector3f temp1, temp2;
		int tA = m1.getBV(iA).primitiveId();
		int tB = m2.getBV(iB).primitiveId();
//		cout << "iA: " << iA << " iB "	 << iB << " tA " << tA << " tB " << tB << endl;
		float correct_dist = detail::TriangleDistancef::triDistance(p1[t1[tA][0]], p1[t1[tA][1]], p1[t1[tA][2]],
																	p2[t2[tB][0]], p2[t2[tB][1]], p2[t2[tB][2]], temp1, temp2);
		cout << "iA: " << iA << " iB " << iB << " dist: " << dist <<", expected: "<< correct_dist << endl;
		assert(approxEquals(correct_dist, dist, 1e-4));
	}

	cout << "test PASSED!" << endl;

}

