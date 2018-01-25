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

#include "fcl/geometry/bvh/detail/BVH_front.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "test_fcl_utility.h"

using namespace std;

#include "SampleObjects.h"

// takes in a BVH of FCL and converts to BVH of this project
void convertFCLToBVH(const fcl::BVHModel<fcl::RSS<float>>& bvh_fcl,
        const vector<fcl::Vector3<float>>& vertices,
        const vector<fcl::Triangle>& triangles,
        BVH *bvh);

void testCONV(const BVH *bvh);

int main(int argc, char *argv[])
{
  std::vector<fcl::Vector3<float>> points;
  std::vector<fcl::Triangle> triangles;

  if(argc < 2)
    createSampleObj1(points, triangles);
  else
    fcl::test::loadOBJFile(argv[1], points, triangles);

  fcl::BVHModel<fcl::RSSf> model;

  model.bv_splitter.reset(new fcl::detail::BVSplitter<fcl::RSSf>(
				fcl::detail::SPLIT_METHOD_MEAN));

  model.beginModel();
  model.addSubModel(points, triangles);
  model.endModel();

  // cout << "num BVs " << model.getNumBVs() << endl;
  // cout << "[model type]" << model.getNodeType() << endl;

  BVH bvh_cdp;
  convertFCLToBVH(model, points, triangles, &bvh_cdp);
  saveBVH(cout, &bvh_cdp);

  if(argc == 1)
    testCONV(&bvh_cdp);

  deleteBVH(&bvh_cdp);
}

RSS getRSSCDP(const fcl::RSSf& rss)
{
	RSS obj;
	obj.axis.v1.x = rss.axis(0,0);
	obj.axis.v1.y = rss.axis(1,0);
  obj.axis.v1.z = rss.axis(2,0);
  obj.axis.v2.x = rss.axis(0,1);
  obj.axis.v2.y = rss.axis(1,1);
  obj.axis.v2.z = rss.axis(2,1);
  obj.axis.v3.x = rss.axis(0,2);
  obj.axis.v3.y = rss.axis(1,2);
  obj.axis.v3.z = rss.axis(2,2);

	obj.To.x = rss.To(0);
  obj.To.y = rss.To(1);
  obj.To.z = rss.To(2);

	obj.l[0] = rss.l[0];
	obj.l[1] = rss.l[1];

	obj.r = rss.r;

	obj.size = rss.size();

	return obj;
}

Vector3 getVec3(const fcl::Vector3<float>& vec_fcl)
{
	return Vector3({vec_fcl(0), vec_fcl(1), vec_fcl(2)});
}

Triangle getTriangle(const vector<fcl::Vector3<float>>& vertices,
				const fcl::Triangle& tri_fcl)
{
	Triangle obj;
	obj.a = getVec3(vertices[tri_fcl[0]]);
	obj.b = getVec3(vertices[tri_fcl[1]]);
	obj.c = getVec3(vertices[tri_fcl[2]]);
	return obj;
}

void convertFCLToBVH(const fcl::BVHModel<fcl::RSS<float>>& bvh_fcl,
				const vector<fcl::Vector3<float>>& vertices,
				const vector<fcl::Triangle>& triangles,
				BVH *bvh)
{
	bvh->num_bv = bvh_fcl.getNumBVs();
	bvh->num_tri = triangles.size();
	bvh->bv_arr = new BV[bvh->num_bv];
	bvh->tri_arr = new Triangle[bvh->num_tri];

	// copy over all RSS objects and their tree heirarchy
	for(int i = 0; i < bvh_fcl.getNumBVs(); i++){
		bvh->bv_arr[i].rss = getRSSCDP(bvh_fcl.getBV(i).bv);
		bvh->bv_arr[i].id1 = bvh_fcl.getBV(i).leftChild();
		bvh->bv_arr[i].id2 = bvh_fcl.getBV(i).rightChild();
		bvh->bv_arr[i].idt = bvh_fcl.getBV(i).primitiveId();
	}

	for(int i = 0; i < triangles.size(); i++)
	{
		bvh->tri_arr[i] = getTriangle(vertices, triangles[i]);
	}
}

bool approx_Equals(float a, float b, float EPSILON = 5e-7)
{
  if(fabs(a - b) < EPSILON)
    return true;
  return false;
}

bool operator!=(const Vector3& lhs, const Vector3& rhs)
{
  if(!approx_Equals(lhs.x, rhs.x))
    return true;
  if(!approx_Equals(lhs.y, rhs.y))
    return true;
  if(!approx_Equals(lhs.z, rhs.z))
    return true;

  return false;
}

bool operator!=(const RSS& lhs, const RSS& rhs)
{
  if(lhs.axis.v1 != rhs.axis.v1)
    return true;
  if(lhs.axis.v2 != rhs.axis.v2)
    return true;
  if(lhs.axis.v3 != rhs.axis.v3)
    return true;

  if(lhs.To != rhs.To)
    return true;

  if(!approx_Equals(lhs.l[0], rhs.l[0]))
    return true;
  if(!approx_Equals(lhs.l[1], rhs.l[1]))
    return true;

  if(!approx_Equals(lhs.r, rhs.r))
    return true;
  if(!approx_Equals(lhs.size, rhs.size, 1e-4))
    return true;

  return false;
}

bool operator==(const Triangle& lhs, const Triangle& rhs)
{
  if(lhs.a != rhs.a)
    return false;
  if(lhs.b != rhs.b)
    return false;
  if(lhs.c != rhs.c)
    return false;

  return true;

}

bool operator==(const BV& lhs, const BV& rhs)
{
  if(lhs.rss != rhs.rss)
    return false;
  if(lhs.id1 != rhs.id1)
    return false;
  if(lhs.id2 != rhs.id2)
    return false;
  if(lhs.idt != rhs.idt)
    return false;
  return true;
}

void testCONV(const BVH *bvh)
{
  stringstream oss, iss;
  saveBVH(oss, bvh);

  iss.str(oss.str());
  BVH read_bvh;
  loadBVH(iss, &read_bvh);

  assert(read_bvh.num_bv = bvh->num_bv);
  assert(read_bvh.num_tri = bvh->num_tri);
  for(int i = 0; i < bvh->num_bv; i++)
    assert(read_bvh.bv_arr[i] == bvh->bv_arr[i]);

  for(int j = 0; j < bvh->num_tri; j++)
    assert(read_bvh.tri_arr[j] == bvh->tri_arr[j]);
}
