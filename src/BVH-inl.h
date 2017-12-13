/*
 * BVH-inl.h
 *
 *  Created on: Dec 11, 2017
 *      Author: atulya
 */

#include <iostream>
#include <iomanip>
#include <string>

#ifndef BVH_INL_H_
#define BVH_INL_H_

using std::endl;

// probably move them to BVH-cuda.cu
void initializeBV(BV * bv)
{
	bv->rss.axis = Matrix3({{0,0,0}, {0,0,0}, {0,0,0}});
	bv->rss.To = Vector3({0,0,0});
	bv->rss.l[0] = 0;
	bv->rss.l[1] = 0;
	bv->rss.r = 0;

	bv->id1 = -1;
	bv->id2 = -1;
	bv->idt = -1;
}

// probably move them to BVH-cuda.cu
void initializeBVH(BVH *bv)
{
	bv->bv_arr = NULL;
	bv->num_bv = 0;
	bv->tri_arr = NULL;
	bv->num_tri = 0;
}

void deleteBVH(BVH *bvh)
{
	delete[] bvh->bv_arr;
	delete[] bvh->tri_arr;
}

std::ostream& operator<<(std::ostream& os, const Triangle& tri)
{
	os << tri.a.x << " ";
	os << tri.a.y << " ";
	os << tri.a.z << " ";
	os << tri.b.x << " ";
	os << tri.b.y << " ";
	os << tri.b.z << " ";
	os << tri.c.x << " ";
	os << tri.c.y << " ";
	os << tri.c.z << endl;
	return os;
}

std::ostream& operator<<(std::ostream& os, const RSS& rss)
{
	os << rss.axis.v1.x << " ";
	os << rss.axis.v1.y << " ";
	os << rss.axis.v1.z << " ";
	os << rss.axis.v2.x << " ";
	os << rss.axis.v2.y << " ";
	os << rss.axis.v2.z << " ";
	os << rss.axis.v3.x << " ";
	os << rss.axis.v3.y << " ";
	os << rss.axis.v3.z << " ";

	os << rss.To.x << " ";
	os << rss.To.y << " ";
	os << rss.To.z << " ";

	os << rss.l[0] << " ";
	os << rss.l[1] << " ";
	os << rss.r << " ";

	os << rss.size;

	return os;
}

std::ostream& operator<<(std::ostream& os, const BV& bv)
{
	os << bv.rss << " ";
	os << bv.id1 << " " << bv.id2 << " " << bv.idt << endl;
	return os;
}

std::istream& operator>>(std::istream& is, Triangle& tri)
{
	is >> tri.a.x;
	is >> tri.a.y;
	is >> tri.a.z;
	is >> tri.b.x;
	is >> tri.b.y;
	is >> tri.b.z;
	is >> tri.c.x;
	is >> tri.c.y;
	is >> tri.c.z;
	return is;
}

std::istream& operator>>(std::istream& is, RSS& rss)
{
	is >> rss.axis.v1.x;
	is >> rss.axis.v1.y;
	is >> rss.axis.v1.z;
	is >> rss.axis.v2.x;
	is >> rss.axis.v2.y;
	is >> rss.axis.v2.z;
	is >> rss.axis.v3.x;
	is >> rss.axis.v3.y;
	is >> rss.axis.v3.z;

	is >> rss.To.x;
	is >> rss.To.y;
	is >> rss.To.z;

	is >> rss.l[0];
	is >> rss.l[1];
	is >> rss.r;

	is >> rss.size;

	return is;
}

std::istream& operator>>(std::istream& is, BV& bv)
{
	is >> bv.rss;
	is >> bv.id1 >> bv.id2 >> bv.idt;
	return is;
}
// saves the BVH object to a file
void saveBVH(std::ostream& os, const BVH* bvh){
	os << std::setprecision(7);
	os << "NUM_BV: " << bvh->num_bv << " NUM_TRI: " << bvh->num_tri << endl;
	for(int i = 0; i < bvh->num_bv; i++)
		os << bvh->bv_arr[i];

	for(int i = 0; i < bvh->num_tri; i++)
		os << bvh->tri_arr[i];
}

// loads the BVH from a file

// loads the BVH from a file where each row is a BV
void loadBVH(std::istream& is, BVH* bvh){
	std::string str1, str2;
	is >> str1 >> bvh->num_bv >> str2 >> bvh->num_tri;
	bvh->bv_arr = new BV[bvh->num_bv];
	bvh->tri_arr = new Triangle[bvh->num_tri];
	for(int i = 0; i < bvh->num_bv; i++)
		is >> bvh->bv_arr[i];

	for(int i = 0; i < bvh->num_tri; i++)
		is >> bvh->tri_arr[i];
}

#endif /* BVH_INL_H_ */
