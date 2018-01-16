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
#define MAX_QUEUE_CAPACITY 1000000
#define MAX_DFS_SET 128
#define BFS_ROWS 12
#define BFS_COLS 32
#define MAX_BFS_BLOCKS  16

struct Task
{
  int i1;
  int i2;
  float dist;
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
  unsigned int num_bv;

  Triangle *tri_arr;
  unsigned int num_tri;
};

struct Queue
{
  Task *arr;
  unsigned int start;
  unsigned int last;
  unsigned int size;
  unsigned int max_capacity;
};

struct Set
{
  Task arr[MAX_DFS_SET];
  unsigned int size;
};

struct SetBFS
{
  Task arr[MAX_BFS_BLOCKS][BFS_COLS*BFS_ROWS*2];
  unsigned int size[MAX_BFS_BLOCKS];
};


struct Config
{
  float gamma;
  int max_iter;
  float R[3][3];
  float t[3];
  int enable_distance_reduction;  // set it to 0 to print get distance on all possible leaf elements without early termination
  int max_dfs_proc;
  int max_bfs_proc;
};

// probably move them to BVH-cuda.cu
void initializeBV(BV * bv);

// probably move them to BVH-cuda.cu
void initializeBVH(BVH *bv);

// deletes a BVH from memory
void deleteBVH(BVH *bv);

// the main function that computes distance between two BVHs
void computeDistance(const BVH* bvh1, const BVH* bvh2);

// saves the BVH object to a file
void saveBVH(std::ostream& os, const BVH* bvh);

// loads the BVH from a file
void loadBVH(std::istream& is, BVH* bvh);

#include "BVH-inl.h"

#endif /* BVH_H_ */
