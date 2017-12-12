/*
 * BVH.h
 *
 *  Created on: Dec 11, 2017
 *      Author: atulya
 */

#include "RSS.h"
#include <iostream>

#ifndef BVH_H_
#define BVH_H_

// stores a BV and its child information

struct id_pair
{
  int i1;
  int i2;
};

struct BV
{
  RSS rss;
  int id1;
  int id2;
  int idt;
};

struct BVH
{
  BV *bv_arr;
  int num_bv;

  Triangle *tri_arr;
  int num_tri;
};

struct BVQueue
{
  id_pair *queue;
  int num_BV;
  int max_capacity;
};

struct BVLeafQueue
{
  id_pair *queue;
  int num_tri;
  int max_capacity;
};

// probably move them to BVH-cuda.cu
void initializeBV(BV * bv);

// probably move them to BVH-cuda.cu
void initializeBVH(BVH *bv);

// deletes a BVH from memory
void deleteBVH(BVH *bv);

//// initialize the BV and Leaf Queue
//void initializeQueue(BVQueue *bvq , int num_bvq,
//		     BVLeafQueue *leafq, int num_leafq);
//
//void

// saves the BVH object to a file
void saveBVH(std::ostream& os, const BVH* bvh);

// loads the BVH from a file
void loadBVH(std::istream& is, BVH* bvh);

#include "BVH-inl.h"

#endif /* BVH_H_ */
