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

#ifndef BVH_CUDA_INL_H_
#define BVH_CUDA_INL_H_

#define MAX_QUEUE_CAPACITY 1000000
#define MAX_BV_TASKS 10000
#define BLOCKSIZE_BV 32
#define BLOCKSIZE_LEAF 16

#define STOP_CONDITION_QUEUE_FULL 2
#define STOP_CONDITION_BELOW_THRESH 3
#define STOP_CONDITION_COMPLETED 1

using std::cout;
using std::endl;
using Transform3f = Eigen::Transform<float, 3, Eigen::AffineCompact>;

struct Config
{
  float gamma;
  int dbg_size;
  float R[3][3];
  float t[3];
  int enable_distance_reduction;  // set it to 0 to print get distance on all possible leaf elements without early termination
};
struct DebugVar
{
  Task tsk;
  uint idx;
  uint idy;
};
struct DistanceResult
{
  float dist;
  int stop;
  Task tsk;
  Task tsk2;
  unsigned int idx;
  unsigned int idy;
};

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
                          Queue* bv_q, Queue *l_q, DistanceResult* res,
                          const int num)
{
  BV bvA;
  BV bvB;
  float d_min = res->dist;

  // Evaluate BV only if stopping condition has not been met
  if(d_min > cfg->gamma)
  {
    // next split and add tasks if necessary
    for(int i = 0; i < num; i++)
    {
      Task new_t1, new_t2;
      int arr_id = getIdFromQueue(bv_q, 0);
      removeNFromQueue(bv_q, 1);
      new_t1.i1 = bv_q->arr[arr_id].i1;
      new_t1.i2 = bv_q->arr[arr_id].i2;
      new_t2 = new_t1;

      // task would be added only if the distance between BVs is less than d_min
      if(bv_q->arr[arr_id].dist <= d_min)
      {
        bvA = bvhA->bv_arr[bv_q->arr[arr_id].i1];
        bvB = bvhB->bv_arr[bv_q->arr[arr_id].i2];

        // first over second if First is not leaf and Second is
        if(bvA.id1 > 0 && bvB.id1 < 0)
        {
          new_t1.i1 = bvA.id1;
          new_t2.i1 = bvA.id2;
        }
        else if(bvA.id1 > 0 && bvA.rss.size > bvB.rss.size)
        {
          new_t1.i1 = bvA.id1;
          new_t2.i1 = bvA.id2;
        }

        //second over first if second is not leaf and first is
        else if(bvB.id1 > 0 && bvA.id1 < 0 )
        {
          new_t1.i2 = bvB.id1;
          new_t2.i2 = bvB.id2;
        }
        else if(bvB.id1 > 0 && bvB.rss.size > bvA.rss.size)
        {
          new_t1.i2 = bvB.id1;
          new_t2.i2 = bvB.id2;
        }
      
        // if the new elements present in tasks are leaf add them to leaf queue
        bvA = bvhA->bv_arr[new_t1.i1];
        bvB = bvhB->bv_arr[new_t1.i2];
        res->tsk.i1 = new_t1.i1;
        res->tsk.i2 = new_t1.i2;
        if(bvA.id1 < 0 && bvB.id1 < 0)
          addToQueue(l_q, new_t1.i1, new_t1.i2);
        else
          addToQueue(bv_q, new_t1.i1, new_t1.i2);

        if(isFull(bv_q))
          res->stop = STOP_CONDITION_QUEUE_FULL;

        bvA = bvhA->bv_arr[new_t2.i1];
        bvB = bvhB->bv_arr[new_t2.i2];
        res->tsk2.i1 = new_t2.i1;
        res->tsk2.i2 = new_t2.i2;

        if(bvA.id1 < 0 && bvB.id1 < 0)
          addToQueue(l_q, new_t2.i1, new_t2.i2);
        else
          addToQueue(bv_q, new_t2.i1, new_t2.i2);

        if(isFull(bv_q))
          res->stop = STOP_CONDITION_QUEUE_FULL;
      }
    }
  }
  else
  {
    res->stop = STOP_CONDITION_BELOW_THRESH;
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
                             Queue* bv_q, const int num)
{
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int arr_id = getIdFromQueue(bv_q, j);

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
  initializeQueue(bv_queue);
  initializeQueue(l_queue);
  addToQueue(bv_queue, 0, 0);

  DistanceResult res = *result;

  dim3 dimBlockBV(1, BLOCKSIZE_BV);
  dim3 dimGridBV(1,1);
  dim3 dimBlockLeaf(1, BLOCKSIZE_LEAF);
  dim3 dimGridLeaf(1,1);

  // for(int i = 0; i < cfg->dbg_size; i++)
  // {
  //   dbg_var[i].idy = i;
  //   dbg_var[i].tsk.dist = i * 3.0;
  // 4

  for(int i = 0;i < 50; i++)
  {
    // PROCESS THE BV TREE NODE TASKS
    int num_bv_tasks = getSize(bv_queue);
    if(num_bv_tasks > MAX_BV_TASKS)
      num_bv_tasks = MAX_BV_TASKS;
    
    dimGridBV.y = (num_bv_tasks - 1)/BLOCKSIZE_BV + 1;
    if(num_bv_tasks)
    {
      processBVTasks<<<dimGridBV, dimBlockBV>>>(bvhA, bvhB, cfg, bv_queue, 
                                        num_bv_tasks);
    }
    cudaDeviceSynchronize();

    // PROCESS THE LEAF NODE TASKS
    int num_leaf_tasks = getSize(l_queue);
    dimGridLeaf.y = (num_leaf_tasks - 1)/BLOCKSIZE_LEAF + 1;

    if(num_leaf_tasks)
    {
      processLeafTasks<<<dimGridLeaf, dimBlockLeaf>>>(bvhA, bvhB, cfg, l_queue, 
                                        num_leaf_tasks);
    }
    cudaDeviceSynchronize();

    // REDUCE THE LEAF TREE DATA
    if(cfg->enable_distance_reduction)
    {
      reduceLeafTasks(bvhA, bvhB, cfg, l_queue, result, num_leaf_tasks);
      cudaDeviceSynchronize();
    }

    // ADD MORE TASKS FROM BV TREE
    addBVTasks(bvhA, bvhB, cfg, bv_queue, l_queue, result, num_bv_tasks);
    cudaDeviceSynchronize();

    res = *result;

    // CHECK FOR STOPPING CONDITION
    if(num_bv_tasks == 0 && num_leaf_tasks == 0)
      res.stop = STOP_CONDITION_COMPLETED;

    dbg_var[i].idx = num_bv_tasks;
    dbg_var[i].idy = num_leaf_tasks;

  }

  BV bvA = bvhA->bv_arr[2];
  BV bvB = bvhB->bv_arr[0];

  RSSResult rss_res;
  computeDistance(&bvA.rss, &bvB.rss,cfg->R, cfg->t, &rss_res);

  for(int i = 0; i < 10; i++)
    dbg_var[i].tsk = bv_queue->arr[i];
  for(int i = 0; i < 5; i++)
    dbg_var[i + 10].tsk = l_queue->arr[i];

  *result = res;
  result->tsk.dist = rss_res.dist;
  result->idx = bv_queue->size;
  result->idy = l_queue->size;
  // result->dist = bv_queue->arr[0].i1;
}

__host__ void printDebugInfo(DebugVar* arr, int size)
{
  cout << std::setprecision(4);
  for(int i = 0; i < size; i++)
  {
    cout << "Task: " << arr[i].tsk.i1 << " " << arr[i].tsk.i2 
                    << " " << arr[i].tsk.dist << " ";
    cout << "idx: " << arr[i].idx << " " ;
    cout << "idy: " << arr[i].idy;
    cout << endl;
  }
}

void initializeResult(DistanceResult& result)
{
  result.dist = 1.0e38;
  result.stop = 0;
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

__host__ std::vector<DistanceResult> computeDistance(const BVH* bvh1, const BVH* bvh2, Config cfg,
                                        std::vector<Transform3f>& transforms,
                                        std::string debg_queue_filename)
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
  DebugVar *h_dbg = new DebugVar[cfg.dbg_size];
  for(int i = 0; i< cfg.dbg_size; i++)
    h_dbg[i] = DebugVar({{0,0,-1}, 0, 0});
  DebugVar *d_dbg;
  cudaMalloc(&d_dbg, cfg.dbg_size * sizeof(DebugVar));

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

  //COPY OVER OUTPUT VARS
  cudaMemcpy(d_res, &h_res, sizeof(DistanceResult), cudaMemcpyHostToDevice);

  //COPY OVER DEBUG VARS
  cudaMemcpy(d_dbg, h_dbg, cfg.dbg_size * sizeof(DebugVar), cudaMemcpyHostToDevice);

  for(int i = 0; i < transforms.size(); i++)
  {
    updateConfigTransform(cfg, transforms[i]); 

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
    cudaMemcpy(h_dbg, d_dbg, cfg.dbg_size * sizeof(DebugVar), cudaMemcpyDeviceToHost);
    // COPY BACK RESULTS
    cudaMemcpy(&h_res, d_res, sizeof(DistanceResult), cudaMemcpyDeviceToHost);

    t_copy2 = get_wall_time();

    results.push_back(h_res);
    printDebugInfo(h_dbg, cfg.dbg_size);

    cout << "======== STATS =========" << endl;
    cout << "COPY1: " << t_copy1 - t_init << "ms" << endl;
    cout << "RUN: " << t_run - t_copy1 << "ms" << endl;
    cout << "COPY2: " << t_copy2 - t_run << "ms" << endl;

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
