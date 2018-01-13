/*
 * BVH-cuda-inl.h
 *
 *  Created on: Dec 12, 2017
 *      Author: atulya
 */

#include "BVH.h"
#include "Triangle-cuda-inl.h"
#include "RSS-cuda-inl.h"
#include "Eigen/Dense"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

#include <unistd.h>

#ifndef BVH_CUDA_INL_H_
#define BVH_CUDA_INL_H_

#define BLOCKSIZE_BV 32
#define BLOCKSIZE_LEAF 32
#define BLOCKSIZE_TASKS_ADDER 32
#define BLOCKSIZE_RED 1024
#define MAX_BLOCKS_RED  8192/(2*BLOCKSIZE_RED)

#define STOP_CONDITION_QUEUE_FULL 2
#define STOP_CONDITION_BELOW_THRESH 3
#define STOP_CONDITION_COMPLETED 1

using std::cout;
using std::endl;
using Transform3f = Eigen::Transform<float, 3, Eigen::AffineCompact>;

struct DebugVar
{
  int dfs_bv;
  int total_bv;
  int proc_bv;
  int num_tri;

  Task tsk1;
  Task tsk2;

};
struct DistanceResult
{
  float dist;
  int stop;
  int num_iter;
  Task tsk;
  Task tsk2;
  unsigned int idx;
  unsigned int idy;
};

// DECLARE DEVICE GLOBAL VARIABLES HERE
__device__ Set dfs_set;
__device__ Set dfs_extra;
__device__ Set leaf_tasks_dfs;

// REQUIRES : A valid queue
// MODIFIES : The queue structure
// EFFECTS  : initializes the start and last elements of the queue
__device__ __inline__ void initializeQueue(Queue *q)
{
  q->start = 0;
  q->last = 0;
  q->size = 0;
  q->max_capacity = MAX_QUEUE_CAPACITY;
}

// REQUIRES : A valid queue
//            ids of two elements
// MODIFIES : The queue structure
// EFFECTS  : adds a pair of ids to the queue and initializes its distance to
//            inf
__device__ __inline__ void addToQueue(Queue* bvq, const int id1, const int id2)
{
  bvq->arr[bvq->last].i1 = id1;
  bvq->arr[bvq->last].i2 = id2;
  bvq->arr[bvq->last].dist = 1e37;
  bvq->last = (bvq->last + 1)%(bvq->max_capacity);
  bvq->size = bvq->size + 1;
}

// REQUIRES : A valid queue
//            0 <= index < queue size
// MODIFIES : 
// EFFECTS  : gives the actual index of the particular element inside the queue
//            array
__device__ __inline__ int getIdFromQueue(const Queue* bvq, const int index) {
  return (index + bvq->start)%(bvq->max_capacity);
}

// REQUIRES : A valid queue
//            0 <= N < queue size
// MODIFIES : the queue data structure
// EFFECTS  : removes N elements from the queue by shifting the starting index 
//            of the queue inside the 1D array
__device__ __inline__ void removeNFromQueue(Queue* bvq, const int N) {
  bvq->start = (bvq->start + N)%(bvq->max_capacity);
  bvq->size = bvq->size - N;
}

// EFFECTS  : returns the total number of elements inside the queue
__device__ __inline__ int getSize(const Queue* bvq) {
  return bvq->size;
}

// EFFECTS  : returns if the queue is full or not
__device__ __inline__ int isFull(const Queue* bvq) {
  if(bvq->size >= bvq->max_capacity)
    return 1;
  return 0;
}

// REQUIRES : Array of tasks all of which are valid
//            starting value of tasks
//            num < 2*BLOCKDIM_RED
// MODIFIES : res
// EFFECTS  : computes the minimum distance of all tasks
__global__ void computeMin(const Task* arr, const int start, const int num,
                          float* res)
{
  // ASUMING blockDim.y == BLOCSIZE_RED
  int ty = threadIdx.y;
  int tid = blockIdx.y*2*BLOCKSIZE_RED + threadIdx.y;

  __shared__ float min_vals[2*BLOCKSIZE_RED];
  float val1 = 1e38, val2 = 1e38;

  int stepsize = BLOCKSIZE_RED;

  if(tid < num && arr[start + tid].i1 >= 0)
    min_vals[ty] = arr[start + tid].dist;
  else
    min_vals[ty] = 1e38;

  if(tid + BLOCKSIZE_RED < num && arr[start + tid + BLOCKSIZE_RED].i1 >= 0)
    min_vals[ty + BLOCKSIZE_RED] = arr[start + tid + BLOCKSIZE_RED].dist;
  else
    min_vals[ty + BLOCKSIZE_RED] = 1e38;

  __syncthreads();

  while(stepsize > 0)
  {
    if(ty < stepsize)
    {
      val1 = min_vals[ty];
      val2 = min_vals[ty + stepsize];

      min_vals[ty] = min(val1, val2);
    }

    stepsize = stepsize/2;
    __syncthreads();
  }

  if(ty == 0)
    res[blockIdx.y] = min_vals[0];

}

// REQUIRES : Array of tasks all of which are valid
//            starting value of tasks
//            num < 2*BLOCKDIM_RED
// MODIFIES : res
// EFFECTS  : computes the number of valid tasks
__global__ void countTasks(const Task* arr, const int start, const int num,
                          int* res)
{
  // ASUMING blockDim.y == BLOCSIZE_RED
  int ty = threadIdx.y;
  int tid = blockIdx.y*2*BLOCKSIZE_RED + threadIdx.y;

  __shared__ int count_vals[2*BLOCKSIZE_RED];
  int val1 = 0, val2 = 0;

  int stepsize = BLOCKSIZE_RED;

  if(tid < num && arr[start + tid].i1 >= 0)
    count_vals[ty] = 1;
  else
    count_vals[ty] = 0;

  if(tid + BLOCKSIZE_RED < num && arr[start + tid + BLOCKSIZE_RED].i1 >= 0)
    count_vals[ty + BLOCKSIZE_RED] = 1;
  else
    count_vals[ty + BLOCKSIZE_RED] = 0;

  __syncthreads();

  while(stepsize > 0)
  {
    if(ty < stepsize)
    {
      val1 = count_vals[ty];
      val2 = count_vals[ty + stepsize];

      count_vals[ty] = val1 + val2;
    }

    stepsize = stepsize/2;
    __syncthreads();
  }

  if(ty == 0)
    res[blockIdx.y] = count_vals[0];

}

// REQUIRES : Two bounding volume heirarchy bvhA, bvhB
//            Config for the computation
//            Queue of leaf elements
//            DistanceResult to update results
//            num : the number of leaf tasks that have been updated
// MODIFIES : leaf queue
//            Distance result
// EFFECTS  : given some tasks in the leaf queue it evaluates them to compute
//            the minima and updates the global values appropriately
//            num elements from leaf queue are extracted
//            The result values are updated
//            TODO : parallelize this over multiple threads
__device__ void reduceLeafTasks(const BVH* bvhA, const BVH* bvhB, const Config* cfg,
                          Queue *l_q, DistanceResult* res,
                          const int num)
{
  float d_min = res->dist;

  // first check the leaf checks to search for potential minimas
  for(int i = 0; i < num; i++)
  {
    int arr_id = getIdFromQueue(l_q, 0);
    d_min = fminf(d_min, l_q->arr[arr_id].dist);
    removeNFromQueue(l_q, 1);
  }
  res->dist = d_min;
}

// REQUIRES : Two bounding volumes bvA, bvB
// MODIFIES : new_t1, new_t2
// EFFECTS  : creates new tasks given two bounding volumes
__device__ __inline__ void getNewTasks(const BV* bvA, const BV* bvB,
                            Task* new_t1, Task* new_t2)
{
  // first over second if First is not leaf and Second is
  if(bvA->id1 > 0 && bvB->id1 < 0)
  {
    new_t1->i1 = bvA->id1;
    new_t2->i1 = bvA->id2;
  }
  else if(bvA->id1 > 0 && bvA->rss.size > bvB->rss.size)
  {
    new_t1->i1 = bvA->id1;
    new_t2->i1 = bvA->id2;
  }

  //second over first if second is not leaf and first is
  else if(bvB->id1 > 0 && bvA->id1 < 0 )
  {
    new_t1->i2 = bvB->id1;
    new_t2->i2 = bvB->id2;
  }
  else if(bvB->id1 > 0)
  {
    new_t1->i2 = bvB->id1;
    new_t2->i2 = bvB->id2;
  }
}

// REQUIRES : Two bounding volume heirarchy bvhA, bvhB
//            dfs_set
//            res->dist for deciding whether to split
// MODIFIES : dfs_set, dfs_extra, l_q
// EFFECTS  : processes tasks from dfs_set and adds them either to the same
//            set or to the dfs_xtra set if the new task is a BV
//            if the new task is leaf adds it to the leaf_tasks
//  TODO(ADV write for it to handle multiple blocks)
__global__ void addDFSTasks(const BVH* bvhA, const BVH* bvhB, const Config* cfg,
                          DistanceResult* res)
{
  __shared__ Task taskSetA[MAX_DFS_SET*2];
  __shared__ Task taskSetA_red[MAX_DFS_SET*2];
  __shared__ Task taskSetLeaf[MAX_DFS_SET*2];
  __shared__ Task taskSetLeaf_red[MAX_DFS_SET*2];

  __shared__ int set_size, set_size_max;
  __shared__ int t_size, lt_size;

  int tid = threadIdx.y;
  set_size = dfs_set.size;
  set_size_max = cfg->max_dfs_proc;
  float d_min = res->dist;

  Task new_t1, new_t2;
  new_t1.i1 = dfs_set.arr[tid].i1;
  new_t1.i2 = dfs_set.arr[tid].i2;
  new_t1.dist = 1e37;
  new_t2 = new_t1;

  if(tid < set_size)
  {
    if(dfs_set.arr[tid].dist <= d_min)
    {
      getNewTasks(&bvhA->bv_arr[new_t1.i1],
                  &bvhB->bv_arr[new_t1.i2],
                  &new_t1, &new_t2);

      // if the new elements present in tasks are leaf add them to leaf queue
      const BV *bvA = &bvhA->bv_arr[new_t1.i1];
      const BV *bvB = &bvhB->bv_arr[new_t1.i2];
    
      if(bvA->id1 < 0 && bvB->id1 < 0)
      {
        taskSetLeaf[2*tid] = new_t1;
        taskSetA[2*tid].i1 = -1;
      }
      else
      {
        taskSetLeaf[2*tid].i1 = -1;
        taskSetA[2*tid] = new_t1;
      }

      bvA = &bvhA->bv_arr[new_t2.i1];
      bvB = &bvhB->bv_arr[new_t2.i2];

      if(bvA->id1 < 0 && bvB->id1 < 0)
      {
        taskSetLeaf[2*tid + 1] = new_t2;
        taskSetA[2*tid + 1].i1 = -1;
      }
      else
      {
        taskSetLeaf[2*tid + 1].i1 = -1;
        taskSetA[2*tid + 1] = new_t2;
      }
    }
    else
    {
      //TODO(ADV OPT) remove this else and place everything before
      taskSetA[2*tid].i1 = -1;
      taskSetA[2*tid + 1].i1 = -1;
      taskSetLeaf[2*tid].i1 = -1;
      taskSetLeaf[2*tid + 1].i1 = -1;
    }

  }

  __syncthreads();

  // REDUCE THE TASKS TO A CONTIGUOUS LIST
  if(tid == 0)
  {
    t_size = 0;
    lt_size = 0;
    //TODO(OPT) avenue to further optimize
    for(int i = 0; i < 2*set_size; i++)
    {
      if(taskSetA[i].i1 >= 0)
      {
        taskSetA_red[t_size] = taskSetA[i];
        t_size++;
      }

      if(taskSetLeaf[i].i1 >= 0)
      {
        taskSetLeaf_red[lt_size] = taskSetLeaf[i];
        lt_size++;
      }
    }

    dfs_set.size = min(t_size, set_size_max);
    dfs_extra.size = max(t_size - set_size_max, 0);
    leaf_tasks_dfs.size = lt_size;
  }

  __syncthreads();

  //COPY OVER TO DFS_TASKS OR DFS_EXTRA OR LEAF_TASKS
  if(tid < dfs_set.size)
    dfs_set.arr[tid] = taskSetA_red[tid];

  if(tid < dfs_extra.size)
    dfs_extra.arr[tid] = taskSetA_red[tid + set_size_max];

  if(2*tid < lt_size)
    leaf_tasks_dfs.arr[2*tid] = taskSetLeaf_red[2*tid];

  if(2*tid + 1 < lt_size)
    leaf_tasks_dfs.arr[2*tid + 1] = taskSetLeaf_red[2*tid + 1];

  __threadfence_block();
}

// REQUIRES : Two bounding volume heirarchy bvhA, bvhB
//            bv_q
//            res->dist for deciding whether to split
// MODIFIES : bv_q, l_q
//            res
// EFFECTS  : processes tasks from the input queue and adds them to the same
//            queue if its a BV task or to the leaf task queue
//  TODO(ADV write for it to handle multiple blocks)
__global__ void addBFSTasks(const BVH* bvhA, const BVH* bvhB, const Config* cfg,
                          Queue* bv_q, Queue *l_q, DistanceResult* res,
                          const int num)
{

  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int tid = threadIdx.x * blockDim.y + threadIdx.y;
  float d_min = res->dist;

  __shared__ Task taskSetA[BFS_ROWS*BFS_COLS*2];
  __shared__ Task taskSetA_red[BFS_ROWS*BFS_COLS*2];
  __shared__ Task taskSetLeaf[BFS_ROWS*BFS_COLS*2];
  __shared__ Task taskSetLeaf_red[BFS_ROWS*BFS_COLS*2];

  __shared__ int setA_size[BFS_ROWS];
  __shared__ int setLeaf_size[BFS_ROWS];
  __shared__ int tsize, lsize;

  taskSetA[tx*2*BFS_COLS + 2*ty].i1 = -1;
  taskSetA[tx*2*BFS_COLS + 2*ty + 1].i1 = -1;
  taskSetLeaf[tx*2*BFS_COLS + 2*ty].i1 = -1;
  taskSetLeaf[tx*2*BFS_COLS + 2*ty + 1].i1 = -1;

  int start = bv_q->start;
  int last = bv_q->last;

  Task new_t1, new_t2;
  new_t1.i1 = bv_q->arr[start + tid].i1;
  new_t1.i2 = bv_q->arr[start + tid].i2;
  float bv_dist = bv_q->arr[start + tid].dist;
  new_t1.dist = 1e37;
  new_t2 = new_t1;

  if(tid < num)
  {
    if(bv_dist <= d_min)
    {
      getNewTasks(&bvhA->bv_arr[new_t1.i1],
                  &bvhB->bv_arr[new_t1.i2],
                  &new_t1, &new_t2);

      // if the new elements present in tasks are leaf add them to leaf queue
      const BV *bvA = &bvhA->bv_arr[new_t1.i1];
      const BV *bvB = &bvhB->bv_arr[new_t1.i2];
    
      if(bvA->id1 < 0 && bvB->id1 < 0)
      {
        taskSetLeaf[tx*2*BFS_COLS + 2*ty] = new_t1;
        taskSetA[tx*2*BFS_COLS + 2*ty].i1 = -1;
      }
      else
      {
        taskSetLeaf[tx*2*BFS_COLS + 2*ty].i1 = -1;
        taskSetA[tx*2*BFS_COLS + 2*ty] = new_t1;
      }

      bvA = &bvhA->bv_arr[new_t2.i1];
      bvB = &bvhB->bv_arr[new_t2.i2];

      if(bvA->id1 < 0 && bvB->id1 < 0)
      {
        taskSetLeaf[tx*2*BFS_COLS + 2*ty + 1] = new_t2;
        taskSetA[tx*2*BFS_COLS + 2*ty + 1].i1 = -1;
      }
      else
      {
        taskSetLeaf[tx*2*BFS_COLS + 2*ty + 1].i1 = -1;
        taskSetA[tx*2*BFS_COLS + 2*ty + 1] = new_t2;
      }
    }
  }
  __syncthreads();
  if(ty == 0)
  {
    setA_size[tx] = 0;
    setLeaf_size[tx] = 0;
    for(int i = 0; i < 2*BFS_COLS; i++)
    {
      if(taskSetA[tx*2*BFS_COLS + i].i1 >= 0)
      {
        taskSetA_red[tx*2*BFS_COLS + setA_size[tx]] = taskSetA[tx*2*BFS_COLS + i];
        setA_size[tx]++;
      }
      if(taskSetLeaf[tx*2*BFS_COLS + i].i1 >= 0)
      {
        taskSetLeaf_red[tx*2*BFS_COLS + setLeaf_size[tx]] = taskSetLeaf[tx*2*BFS_COLS + i];
        setLeaf_size[tx]++;
      }
    }
  }
  __syncthreads();
  if(tx == 0)
  {
    if(ty == 0)
    {
      tsize = 0;
      lsize = 0;
      // res->tsk2 = taskSetA[4];
      // res->tsk.i1 = setLeaf_size[3];
      // res->tsk.i2 = setLeaf_size[4];
      // res->tsk.dist = setLeaf_size[5];
      // res->tsk2.i1 = setA_size[3];
      // res->tsk2.i2 = setA_size[4];
      // res->tsk2.dist = setA_size[5];
    }
    
    for(int i = 0; i < BFS_ROWS; i++)
    {
      if(ty < setA_size[i])
        taskSetA[tsize + ty] = taskSetA_red[i*2*BFS_COLS + ty];

      if(ty + BFS_COLS < setA_size[i])
        taskSetA[tsize + ty + BFS_COLS] = taskSetA_red[i*2*BFS_COLS + ty + BFS_COLS];

      if(ty < setLeaf_size[i])
        taskSetLeaf[lsize + ty] = taskSetLeaf_red[i*2*BFS_COLS + ty];

      if(ty + BFS_COLS < setLeaf_size[i])
        taskSetLeaf[lsize + ty+ BFS_COLS] = taskSetLeaf_red[i*2*BFS_COLS + ty+ BFS_COLS];
      
      __syncthreads();

      if(ty == 0)
      {
        tsize += setA_size[i];
        lsize += setLeaf_size[i];
      }
    }
    if(ty == 0)
    {
      res->idx = tsize;
      res->idy = lsize;
      // res->tsk2 = taskSetA[4];
      // res->tsk.i1 = setLeaf_size[3];
      // res->tsk.i2 = setLeaf_size[4];
      // res->tsk.dist = setLeaf_size[5];
      // res->tsk2.i1 = setA_size[3];
      // res->tsk2.i2 = setA_size[4];
      // res->tsk2.dist = setA_size[5];
    }
  }

  __syncthreads();
  if(2*tid < tsize)
    bv_q->arr[last + 2*tid] = taskSetA[2*tid];
  if(2*tid + 1 < tsize)
    bv_q->arr[last + 2*tid + 1] = taskSetA[2*tid + 1];

  __syncthreads();
  if(tid == 0)
  { 
    // Assuming capacity to be really large
    bv_q->start += num;
    bv_q->last += tsize;
    bv_q->size += (tsize - num);
  }

  last = l_q->last;
  __syncthreads();

  if(2*tid < lsize)
    l_q->arr[last + 2*tid] = taskSetLeaf[2*tid];
  if(2*tid + 1 < lsize)
    l_q->arr[last + 2*tid + 1] = taskSetLeaf[2*tid + 1];

  __syncthreads();
  if(tid == 0)
  {
    l_q->last += lsize;
    l_q->size += lsize;
  }
}

// REQUIRES : dfs_set
// MODIFIES : remove from bv_q
//            add to dfs_set
// EFFECTS  : transfers (count) tasks from BV queue to dfs_set
__global__ void fillDFSSet(Queue *bv_q, int count)
{
  int tid = threadIdx.y;
  int start = bv_q->start;

  if(tid < count)
  {
    dfs_set.arr[tid] = bv_q->arr[tid + start];
  }
}

// REQUIRES : bv_q
//            dfs_extra
// MODIFIES : bv_q
//            remove from dfs_extra
// EFFECTS  : transfers tasks from dfs_extra to BV queue
__global__ void transferDFSandBFS(Queue* bv_q, int dfs_extra_size, DistanceResult* res)
{
  int tid = threadIdx.y;
  int last = bv_q->last;

  if(tid < dfs_extra_size)
  {
    bv_q->arr[last + tid] = dfs_extra.arr[tid];
  }

  __syncthreads();
}

// REQUIRES : leaf_tasks_dfs 
// MODIFIES : l_q
// EFFECTS  : transfers tasks from dfs set to leaf queue
__global__ void transferDFSandBFSLeafs(Queue* l_q, int set_size)
{
  
  int last = l_q->last;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if(j < set_size)
    l_q->arr[j + last] = leaf_tasks_dfs.arr[j];
}


// REQUIRES : Two bounding volume heirarchy bvhA, bvhB
//            Config for the computation
//            Queue of BV elements
//            DistanceResult to update results
//            num : the number of leaf tasks that have been updated
// MODIFIES : BV queue
//            Stopping conditions inside DistResult
// EFFECTS  : given some tasks in the BV queue it evaluates them to check if 
//            they add further new tasks or can be safely discarded for future
//            New tasks are either add to BV queue or Leaf queue 
//            (TODO : think of possible race conditions)
//            num elements inside the BV queue are extracted
//            Stopping conditions inside DistResult
//            TODO : parallelize this over multiple threads
__device__ void addBVTasks(const BVH* bvhA, const BVH* bvhB, const Config* cfg,
                          Queue* bv_q, Queue *l_q, 
                          DistanceResult* res,
                          int num_bfs_tasks)
{
  float d_min = res->dist;
  dim3 dimBlockTasks(1,BLOCKSIZE_TASKS_ADDER);
  dim3 dimGridTasks(1,1);

  cudaStream_t stream_dfs, stream_bfs;
  cudaStreamCreateWithFlags(&stream_dfs, cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&stream_bfs, cudaStreamNonBlocking);

  // Evaluate BV only if stopping condition has not been met
  if(d_min > cfg->gamma)
  {
    // PROCESS the DFS Set
    dimBlockTasks.y = MAX_DFS_SET;
    dimGridTasks.y = 1;
    addDFSTasks<<<dimGridTasks, dimBlockTasks, 0, stream_dfs>>>(bvhA, bvhB, cfg, res);
    // cudaDeviceSynchronize();

    // ADD BFS tasks
    dimBlockTasks.y = BFS_COLS;
    dimBlockTasks.x = BFS_ROWS;
    dimGridTasks.y = 1;
    while(num_bfs_tasks > 0)
    {
      int num_tasks_run = min(num_bfs_tasks, BFS_ROWS * BFS_COLS);
      addBFSTasks<<<dimGridTasks, dimBlockTasks, 0, stream_bfs>>>(bvhA, bvhB, cfg, 
                                                bv_q, l_q, res, num_tasks_run);
      cudaDeviceSynchronize();

      num_bfs_tasks -= num_tasks_run;
    }

    

    // ADD leaf tasks from DFS to Leaf queue
    int leaf_set_size = leaf_tasks_dfs.size;
    // TODO(run only if leaf tasks exist)
    dimBlockTasks.y = MAX_DFS_SET;
    dimBlockTasks.x = 1;
    dimGridTasks.y = (leaf_set_size - 1)/MAX_DFS_SET + 1;
    transferDFSandBFSLeafs<<<dimGridTasks, dimBlockTasks>>>(l_q, leaf_set_size);
    // Assuming capacity to be really large
    cudaDeviceSynchronize();
    l_q->last += leaf_set_size;
    l_q->size += leaf_set_size;

    // ADD extra BV tasks from DFS to BFS queue
    int dfs_extra_size = dfs_extra.size;
    // TODO(run only if leaf tasks exist)
    dimBlockTasks.y = MAX_DFS_SET;
    dimBlockTasks.x = 1;
    dimGridTasks.y = (dfs_extra_size - 1)/MAX_DFS_SET + 1;
    transferDFSandBFS<<<dimGridTasks, dimBlockTasks>>>(bv_q, dfs_extra_size, res);
    // Assuming capacity to be really large
    cudaDeviceSynchronize();
    dfs_extra.size = 0;
    bv_q->last += dfs_extra_size;
    bv_q->size += dfs_extra_size;

    // FILL up the DFS set if it is empty
    int dfs_size = dfs_set.size;
    int max_dfs_size = cfg->max_dfs_proc;
    if(dfs_size == 0 && bv_q->size > max_dfs_size)
    {
      dimBlockTasks.y = MAX_DFS_SET;
      dimBlockTasks.x = 1;
      dimGridTasks.y = (max_dfs_size - 1)/MAX_DFS_SET + 1;
      fillDFSSet<<<dimGridTasks, dimBlockTasks>>>(bv_q, max_dfs_size);
      cudaDeviceSynchronize();
      bv_q->start += max_dfs_size;
      bv_q->size -= max_dfs_size;
      dfs_set.size = max_dfs_size;
    }

  }
  else
  {
    res->stop = STOP_CONDITION_BELOW_THRESH;
  }
  cudaStreamDestroy(stream_dfs);
  cudaStreamDestroy(stream_bfs);
}

// REQUIRES : Two bounding volume heirarchy bvhA, bvhB
//            Config for the computation
//            Queue of BV elements
//            num : the number of leaf tasks that have been updated
// MODIFIES : BV queue (distance values for each task)
// EFFECTS  : Extract num BVs from the BV queue and evaluates them for pairwise 
//            distance
__global__ void processBVTasksFromSet(const BVH* bvhA, const BVH* bvhB, const Config* cfg,
                             const int num)
{
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if(j < num)
  {
    Task t = dfs_set.arr[j];
    BV bvA = bvhA->bv_arr[t.i1];
    BV bvB = bvhB->bv_arr[t.i2];

    RSSResult res;
    computeDistance(&bvA.rss, &bvB.rss, cfg->R, cfg->t, &res);

    dfs_set.arr[j].dist = res.dist;
  }
}

// REQUIRES : Two bounding volume heirarchy bvhA, bvhB
//            Config for the computation
//            Queue of BV elements
//            num : the number of leaf tasks that have been updated
// MODIFIES : BV queue (distance values for each task)
// EFFECTS  : Extract num BVs from the BV queue and evaluates them for pairwise 
//            distance
__global__ void processBVTasks(const BVH* bvhA, const BVH* bvhB, const Config* cfg,
                             Queue* bv_q, const int start,
                             const int num)
{
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int arr_id = start + j;

  if(j < num)
  {
    Task t = bv_q->arr[arr_id];
    BV bvA = bvhA->bv_arr[t.i1];
    BV bvB = bvhB->bv_arr[t.i2];

    RSSResult res;
    computeDistance(&bvA.rss, &bvB.rss, cfg->R, cfg->t, &res);

    bv_q->arr[arr_id].dist = res.dist;
  }
}

// REQUIRES : Two bounding volume heirarchy bvhA, bvhB
//            Config for the computation
//            Queue of Leaf elements
//            num : the number of leaf tasks that have been updated
// MODIFIES : Leaf queue (distance values for each task)
// EFFECTS  : Extract num triangles from Leaf queue and evaluate them for 
//            pairwise distance
__global__ void processLeafTasks(const BVH* bvhA, const BVH* bvhB, const Config* cfg,
                             Queue* leaf_q, const int num)
{
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if(j < num)
  {
    int arr_id = getIdFromQueue(leaf_q, j);
    Task t = leaf_q->arr[arr_id];
    Triangle tA = bvhA->tri_arr[bvhA->bv_arr[t.i1].idt];
    Triangle tB = bvhB->tri_arr[bvhB->bv_arr[t.i2].idt];

    TriangleResult res;
    computeDistance(&tA, &tB, cfg->R, cfg->t, &res);

    leaf_q->arr[arr_id].dist = res.dist;
  }
}

__global__ void manager(const BVH* bvhA, const BVH* bvhB, const Config* cfg,
                      Queue *bv_queue, Queue *l_queue,
                      DebugVar* dbg_var, DistanceResult *result)
{
  // start of with the roots in the queue
  dfs_set.arr[0] = Task({0, 0, 1e37});
  dfs_set.size = 1;

  initializeQueue(bv_queue);
  initializeQueue(l_queue);
  // addToQueue(bv_queue, 0, 0);

  DistanceResult res = *result;

  dim3 dimBlockBV(1, BLOCKSIZE_BV);
  dim3 dimGridBV(1,1);
  dim3 dimBlockLeaf(1, BLOCKSIZE_LEAF);
  dim3 dimGridLeaf(1,1);

  int num_dfs_tasks;
  int num_bfs_tasks;

  cudaStream_t stream_dfs, stream_bfs, stream_leaf;
  cudaStreamCreateWithFlags(&stream_dfs, cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&stream_bfs, cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&stream_leaf, cudaStreamNonBlocking);

  int i;
  for(i = 0; (i < cfg->max_iter) && (res.stop == 0); i++)
  {
    num_dfs_tasks = dfs_set.size;

    dimGridBV.y = (num_dfs_tasks - 1)/BLOCKSIZE_BV + 1;
    if(num_dfs_tasks)
    {
      processBVTasksFromSet<<<dimBlockBV, dimBlockBV, 0, stream_dfs>>>
                                          (bvhA, bvhB, cfg, num_dfs_tasks);
    }

    // PROCESS THE BV TREE NODE BFS TASKS
    num_bfs_tasks = getSize(bv_queue);
    dbg_var[i].total_bv = num_bfs_tasks;

    if(num_bfs_tasks > cfg->max_bfs_proc)
      num_bfs_tasks = cfg->max_bfs_proc;
    
    dimGridBV.y = (num_bfs_tasks - 1)/BLOCKSIZE_BV + 1;
    if(num_bfs_tasks)
    {
      int start = bv_queue->start;
      processBVTasks<<<dimGridBV, dimBlockBV, 0, stream_bfs>>>
                      (bvhA, bvhB, cfg, bv_queue, start, num_bfs_tasks);
    }
    // cudaDeviceSynchronize();

    // PROCESS THE LEAF NODE TASKS
    int num_leaf_tasks = getSize(l_queue);
    dimGridLeaf.y = (num_leaf_tasks - 1)/BLOCKSIZE_LEAF + 1;

    if(num_leaf_tasks)
    {
      processLeafTasks<<<dimGridLeaf, dimBlockLeaf, 0, stream_leaf>>>
                            (bvhA, bvhB, cfg, l_queue, num_leaf_tasks);
    }
    cudaDeviceSynchronize();

    // ADD MORE TASKS FROM BV TREE
    addBVTasks(bvhA, bvhB, cfg,
              bv_queue, l_queue, 
              result, num_bfs_tasks);
    cudaDeviceSynchronize();

    // REDUCE THE LEAF TREE DATA
    if(cfg->enable_distance_reduction)
    {
      reduceLeafTasks(bvhA, bvhB, cfg, l_queue, result, num_leaf_tasks);
      cudaDeviceSynchronize();
    }

    res = *result;

    // CHECK FOR STOPPING CONDITION
    if(num_dfs_tasks == 0 && num_bfs_tasks == 0 && num_leaf_tasks == 0)
      res.stop = STOP_CONDITION_COMPLETED;

    dbg_var[i].dfs_bv = num_dfs_tasks;
    dbg_var[i].proc_bv = num_bfs_tasks;
    dbg_var[i].num_tri = num_leaf_tasks;
    dbg_var[i].tsk1 = dfs_set.arr[0];
    dbg_var[i].tsk2 = dfs_set.arr[1];

  }

  // TODO probably remove the next 4
  BV bvA = bvhA->bv_arr[2];
  BV bvB = bvhB->bv_arr[0];

  RSSResult rss_res;
  computeDistance(&bvA.rss, &bvB.rss,cfg->R, cfg->t, &rss_res);

  *result = res;
  
  result->num_iter = i;
  cudaStreamDestroy(stream_dfs);
  cudaStreamDestroy(stream_bfs);
  cudaStreamDestroy(stream_leaf);

  // result->dist = bv_queue->arr[0].i1;
}

__host__ void printDebugInfo(DebugVar* arr, int size)
{
  cout << std::setprecision(4);
  cout << "ITER, DFS_BV, TOTAL_BV, NUM_BV_PROC, NUM_TRI " << endl;
  for(int i = 0; i < size; i++)
  {
    cout << " [" << i << "] : " << arr[i].dfs_bv << ", " << arr[i].total_bv << 
        ", " << arr[i].proc_bv << ", " << arr[i].num_tri;
    cout << "| 1: " << arr[i].tsk1.i1 << ", " << arr[i].tsk1.i2 << ", " << 
        arr[i].tsk1.dist << "| 2: " << arr[i].tsk2.i1 << ", " <<
        arr[i].tsk2.i2 << ", " << arr[i].tsk2.dist << endl;
  }
}

void initializeResult(DistanceResult& result)
{
  result.dist = 1.0e38;
  result.stop = 0;
  result.num_iter = 0;
}

__host__ void updateConfigTransform(Config& cfg, const Transform3f tf)
{
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
      cfg.R[i][j] = tf.linear()(i,j);

    cfg.t[i] = tf.translation()(i);
  }
}

__host__ void printConfig(const Config cfg)
{
  cout << "CONFIG : " << endl;
  cout << "gamma: " << cfg.gamma << ", max_iter:" << cfg.max_iter << endl;
  cout << "R: [" << cfg.R[0][0] << ", " << cfg.R[0][1] << ", " << cfg.R[0][2] <<
    "\n" << cfg.R[1][0] << ", " << cfg.R[1][1] << ", " << cfg.R[1][2] << 
    "\n" << cfg.R[2][0] << ", " << cfg.R[2][1] << ", " << cfg.R[2][2] << endl;
  cout << "T: [" << cfg.t[0] << ", " << cfg.t[1] << ", " << cfg.t[2] << endl;
  cout << "reduceLeafTasks: " << cfg.enable_distance_reduction << endl;
}

__host__ void initializeDbg(DebugVar* arr, int size)
{
  for(int i = 0; i< size; i++)
    arr[i] = DebugVar({0, 0, 0});
}

__host__ std::vector<DistanceResult> computeDistance(const BVH* bvh1, const BVH* bvh2, Config cfg,
                                        const std::vector<Transform3f> transforms,
                                        const std::string debg_queue_filename,
                                        std::vector<float>& elap_time)
{
  double t_init, t_copy1, t_run, t_copy2;

  // CREATE over the bvh onto cuda memory block
  BVH bvh1_cp = *bvh1;
  BVH bvh2_cp = *bvh2;
  BVH *d_bvh1, *d_bvh2;
  cudaMalloc(&d_bvh1, sizeof(BVH));
  cudaMalloc(&d_bvh2, sizeof(BVH));
  cudaMalloc(&bvh1_cp.bv_arr, bvh1->num_bv * sizeof(BV));
  cudaMalloc(&bvh2_cp.bv_arr, bvh2->num_bv * sizeof(BV));
  cudaMalloc(&bvh1_cp.tri_arr, bvh1->num_tri * sizeof(Triangle));
  cudaMalloc(&bvh2_cp.tri_arr, bvh2->num_tri * sizeof(Triangle));

  // CREATE THE DIFFERENT QUEUES FOR OPERATION
  Queue h_bvq({NULL, 0, 0, 0, MAX_QUEUE_CAPACITY});
  Queue *d_bvq;
  cudaMalloc(&h_bvq.arr, h_bvq.max_capacity * sizeof(Task));
  cudaMalloc(&d_bvq, sizeof(Queue));

  Queue h_leafq({NULL, 0, 0, 0, MAX_QUEUE_CAPACITY});
  Queue *d_leafq;
  cudaMalloc(&h_leafq.arr, h_leafq.max_capacity * sizeof(Task));
  cudaMalloc(&d_leafq, sizeof(Queue));

  // Set h_dfs_set({NULL, 0, MAX_DFS_SET});
  // Set *d_dfs_set;
  // cudaMalloc(&h_dfs_set.arr, h_dfs_set.max_capacity * sizeof(Task));
  // cudaMalloc(&d_dfs_set, sizeof(Set));

  //CREATE THE CONFIG VARS
  Config *d_cfg;
  cudaMalloc(&d_cfg, sizeof(Config));

  //CREATE THE OUTPUT VARIABLE
  std::vector<DistanceResult> results;
  DistanceResult h_res;
  initializeResult(h_res);
  DistanceResult* d_res;
  cudaMalloc(&d_res, sizeof(DistanceResult));

  //CREATE THE DEBUGGING VARIABLE
  DebugVar *h_dbg = new DebugVar[cfg.max_iter];
  initializeDbg(h_dbg, cfg.max_iter);
  DebugVar *d_dbg;
  cudaMalloc(&d_dbg, cfg.max_iter * sizeof(DebugVar));

  t_init = get_wall_time();
  //COPY OVER BVH
  cudaMemcpy(bvh1_cp.bv_arr, bvh1->bv_arr, bvh1->num_bv * sizeof(BV), cudaMemcpyHostToDevice);
  cudaMemcpy(bvh2_cp.bv_arr, bvh2->bv_arr, bvh2->num_bv * sizeof(BV), cudaMemcpyHostToDevice);
  cudaMemcpy(bvh1_cp.tri_arr, bvh1->tri_arr, bvh1->num_tri * sizeof(Triangle), cudaMemcpyHostToDevice);
  cudaMemcpy(bvh2_cp.tri_arr, bvh2->tri_arr, bvh2->num_tri * sizeof(Triangle), cudaMemcpyHostToDevice);
  cudaMemcpy(d_bvh1, &bvh1_cp, sizeof(BVH), cudaMemcpyHostToDevice);
  cudaMemcpy(d_bvh2, &bvh2_cp, sizeof(BVH), cudaMemcpyHostToDevice);

  //COPY OVER QUEUES
  cudaMemcpy(d_bvq, &h_bvq, sizeof(Queue), cudaMemcpyHostToDevice);
  cudaMemcpy(d_leafq, &h_leafq, sizeof(Queue), cudaMemcpyHostToDevice);
  // cudaMemcpy(d_dfs_set, &h_dfs_set, sizeof(Set), cudaMemcpyHostToDevice);

  for(int i = 0; i < transforms.size(); i++)
  {
    // update the result values and config for next iteration
    initializeDbg(h_dbg, cfg.max_iter);
    initializeResult(h_res);

    updateConfigTransform(cfg, transforms[i]); 
    printConfig(cfg);

    //COPY OVER DEBUG VARS
    cudaMemcpy(d_dbg, h_dbg, cfg.max_iter * sizeof(DebugVar), cudaMemcpyHostToDevice);

    //COPY OVER OUTPUT VARS
    cudaMemcpy(d_res, &h_res, sizeof(DistanceResult), cudaMemcpyHostToDevice);

    //COPY OVER CFG
    cudaMemcpy(d_cfg, &cfg, sizeof(Config), cudaMemcpyHostToDevice);

    t_copy1 = get_wall_time();

    // MAIN CUDA EXECUTION GOES HERE
    manager<<<1,1>>>(d_bvh1, d_bvh2, d_cfg, d_bvq, d_leafq, d_dbg, d_res);

    cudaDeviceSynchronize();
    t_run = get_wall_time();

    // COPY BACK THE QUEUE 
    cudaMemcpy(&h_bvq, d_bvq, sizeof(Queue), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_leafq, d_leafq, sizeof(Queue), cudaMemcpyDeviceToHost);

    Task *h_bv_arr, *h_leaf_arr;
    h_bv_arr = new Task[h_bvq.last];
    h_leaf_arr = new Task[h_leafq.last];
    cudaMemcpy(h_bv_arr, h_bvq.arr, h_bvq.last * sizeof(Task), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_leaf_arr, h_leafq.arr, h_leafq.last * sizeof(Task), cudaMemcpyDeviceToHost);

    // COPY BACK DEBUG VARS
    cudaMemcpy(h_dbg, d_dbg, cfg.max_iter * sizeof(DebugVar), cudaMemcpyDeviceToHost);
    // COPY BACK RESULTS
    cudaMemcpy(&h_res, d_res, sizeof(DistanceResult), cudaMemcpyDeviceToHost);

    t_copy2 = get_wall_time();

    results.push_back(h_res);
    elap_time.push_back(t_run - t_copy1);

    printDebugInfo(h_dbg, h_res.num_iter);

    cout << "======== STATS =========" << endl;
    cout << "COPY1: " << t_copy1 - t_init << "s" << endl;
    cout << "RUN: " << t_run - t_copy1 << "s" << endl;
    cout << "COPY2: " << t_copy2 - t_run << "s" << endl;

    // print out all the queue info that was gathered
    std::ofstream of;
    of << std::setprecision(7);
    std::string filename = debg_queue_filename + '_';
    filename += ('0'+i);
    filename += ".outp.csv";
    of.open(filename.c_str());

    of << "NUM_BV: " << h_bvq.last << " NUM_LEAF: " << h_leafq.last << endl;
    for(int i = 0;i < h_bvq.last; i++)
      of << h_bv_arr[i].i1 << " " << h_bv_arr[i].i2 << " " << h_bv_arr[i].dist << endl;
    for(int i = 0;i < h_leafq.last; i++)
      of << h_leaf_arr[i].i1 << " " << h_leaf_arr[i].i2 << " " << h_leaf_arr[i].dist << endl;

    delete[] h_bv_arr;
    delete[] h_leaf_arr;
  }
  
  // FREE CUDA MEMORY
  cudaFree(bvh1_cp.bv_arr);
  cudaFree(bvh2_cp.bv_arr);
  cudaFree(bvh1_cp.tri_arr);
  cudaFree(bvh2_cp.tri_arr);
  cudaFree(d_bvh1);
  cudaFree(d_bvh2);
  cudaFree(h_bvq.arr);
  cudaFree(h_leafq.arr);
  cudaFree(d_bvq);
  cudaFree(d_leafq);
  cudaFree(d_cfg);
  cudaFree(d_res);
  delete[] h_dbg;
  cudaFree(d_dbg);

  return results;
}

#endif /* BVH_CUDA_INL_H_ */
