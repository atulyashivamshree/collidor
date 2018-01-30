/*
 *  test_bvh.cu
    @description : runs some basic checks on the CUDA algorithms defined in
 BVH-cuda-inl

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#include "../compile_CUDA.h"

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

#include "../BVH-cuda-inl.h"
#include "../utils/parse_utils.h"

#define EPSILON 1e-7
#define NUM_TASKS 8192
#define NUM_TESTS 10000

using std::cout;
using std::endl;
using std::string;
using std::vector;

__host__ bool approx_equal(float a, float b, float epsilon = EPSILON) {
  return (std::abs(a - b) < epsilon);
}

__host__ void loadTasks(Task *arr, int num, string filename);

__host__ void test_MinReduction(Task *);
__host__ void test_countTasks(Task *h_tasks);

__host__ int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "USAGE : ./test_bvh.exe tasks.csv" << endl;
    exit(EXIT_FAILURE);
  }

  Task h_tasks[NUM_TASKS];
  loadTasks(h_tasks, NUM_TASKS, argv[1]);

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  cout << "Device name: " << prop.name << endl;

  test_MinReduction(h_tasks);
  test_countTasks(h_tasks);

  cout << "BVH tests: COMPLETED" << endl;
}

__host__ void loadTasks(Task *arr, int num, string filename) {
  std::ifstream fin;
  fin.open(filename.c_str());
  if (!fin.is_open()) {
    cout << filename << ": could not be opened" << endl;
    exit(EXIT_FAILURE);
  }
  string tmp;
  fin >> tmp >> tmp >> tmp;

  Task t;
  int i = 0;
  while (fin >> t.i1 >> t.i2 >> t.dist) {
    arr[i] = t;
    i++;
  }

  fin.close();
}

__host__ void test_MinReduction(Task *h_tasks) {
  srand(static_cast<unsigned>(time(NULL)));

  Task *d_tasks;
  float h_min[MAX_BLOCKS_RED];
  float *d_min;

  cudaError_t err = cudaSuccess;

  err = cudaMalloc(&d_tasks, NUM_TASKS * sizeof(Task));
  if (err != cudaSuccess) {
    fprintf(stderr, "Failed to allocate device vector A (error code %s)!\n",
            cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  cudaMalloc(&d_min, MAX_BLOCKS_RED * sizeof(float));

  // RANDOM INIT of start and num
  vector<int> test_start;
  vector<int> test_num;

  test_start.push_back(0);
  test_num.push_back(8192);

  test_start.push_back(76);
  test_num.push_back(908);

  for (int i = 0; i < NUM_TESTS; i++) {
    test_start.push_back(rand() % 3000);
    test_num.push_back(rand() % (NUM_TASKS - 3000) + 1);
  }

  cudaMemcpy(d_tasks, h_tasks, NUM_TASKS * sizeof(Task),
             cudaMemcpyHostToDevice);

  dim3 dimBlockTasks(1, BLOCKSIZE_RED);
  dim3 dimGridTasks(1, 1);
  dimGridTasks.y = (NUM_TASKS - 1) / (2 * BLOCKSIZE_RED) + 1;

  int count_fail = 0;

  for (int i = 0; i < test_start.size(); i++) {
    int start = test_start[i];
    int num = test_num[i];

    float min_dist = 1e38, min_dist_gpu = 1e38;
    for (int j = 0; j < num; j++) {
      if (h_tasks[j + start].i1 >= 0)
        min_dist = min(min_dist, h_tasks[j + start].dist);
      if (min_dist < 0)
        cout << h_tasks[j + start].dist << "<0 : j " << j << ", j+start"
             << (j + start) << endl;
    }

    // computeMin<<<dimGridTasks, dimBlockTasks>>>(d_tasks, start, num, d_min);
    computeMin<<<dimGridTasks, dimBlockTasks>>>(d_tasks, start, num, d_min);

    cudaDeviceSynchronize();
    cudaMemcpy(h_min, d_min, MAX_BLOCKS_RED * sizeof(float),
               cudaMemcpyDeviceToHost);

    for (int j = 0; j < MAX_BLOCKS_RED; j++) {
      min_dist_gpu = min(min_dist_gpu, h_min[j]);

      // cout << "dist[" << j << "]: " << h_min[j] << endl;
    }

    if (!approx_equal(min_dist, min_dist_gpu)) {
      cout << "[" << i << "] EXP: " << min_dist << ", ACTUAL: " << min_dist_gpu
           << "start: " << start << ", num: " << num << endl;
      count_fail++;
    }
  }

  cudaFree(d_tasks);
  cudaFree(d_min);

  if (!count_fail)
    cout << "test Minimum computation : PASSED" << endl;
  else
    cout << "test Minimum computation : FAILED" << count_fail << "/"
         << test_start.size() << endl;
}

__host__ void test_countTasks(Task *h_tasks) {
  srand(static_cast<unsigned>(time(NULL)));

  Task *d_tasks;
  int h_count[MAX_BLOCKS_RED];
  int *d_count;

  cudaError_t err = cudaSuccess;

  err = cudaMalloc(&d_tasks, NUM_TASKS * sizeof(Task));
  if (err != cudaSuccess) {
    fprintf(stderr, "Failed to allocate device vector A (error code %s)!\n",
            cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  cudaMalloc(&d_count, MAX_BLOCKS_RED * sizeof(int));

  // RANDOM INIT of start and num
  vector<int> test_start;
  vector<int> test_num;

  test_start.push_back(0);
  test_num.push_back(8192);

  test_start.push_back(76);
  test_num.push_back(908);

  for (int i = 0; i < NUM_TESTS; i++) {
    test_start.push_back(rand() % 3000);
    test_num.push_back(rand() % (NUM_TASKS - 3000) + 1);
  }

  cudaMemcpy(d_tasks, h_tasks, NUM_TASKS * sizeof(Task),
             cudaMemcpyHostToDevice);

  dim3 dimBlockTasks(1, BLOCKSIZE_RED);
  dim3 dimGridTasks(1, 1);
  dimGridTasks.y = (NUM_TASKS - 1) / (2 * BLOCKSIZE_RED) + 1;

  int count_fail = 0;

  for (int i = 0; i < test_start.size(); i++) {
    int start = test_start[i];
    int num = test_num[i];

    float count_actual = 0, count_gpu = 0;
    for (int j = 0; j < num; j++) {
      if (h_tasks[j + start].i1 >= 0) count_actual++;
    }

    countTasks<<<dimGridTasks, dimBlockTasks>>>(d_tasks, start, num, d_count);

    cudaDeviceSynchronize();
    cudaMemcpy(h_count, d_count, MAX_BLOCKS_RED * sizeof(int),
               cudaMemcpyDeviceToHost);

    for (int j = 0; j < MAX_BLOCKS_RED; j++) {
      count_gpu += h_count[j];
    }

    if (count_actual != count_gpu) {
      cout << "[" << i << "] EXP: " << count_actual << ", ACTUAL: " << count_gpu
           << "start: " << start << ", num: " << num << endl;
      count_fail++;
    }
  }

  cudaFree(d_tasks);
  cudaFree(d_count);

  if (!count_fail)
    cout << "test Count computation : PASSED" << endl;
  else
    cout << "test Count computation : FAILED" << count_fail << "/"
         << test_start.size() << endl;
}
