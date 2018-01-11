// @author : Atulya Shivam Shree
// @description : checks triangle intersection on CUDA

#include "compile_CUDA.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>

#include "BVH-cuda-inl.h"
#include "utils/parse_utils.h"

#define EPSILON 1e-7
#define NUM_TASKS   8192
#define NUM_TESTS   20

using namespace std;

__host__ bool approx_equal(float a, float b, float epsilon = EPSILON) {
  return (std::abs(a - b) < epsilon);
}

__host__ void test_MinReduction();

__host__ int main(int argc, char *argv[])
{
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  cout << "Device name: " << prop.name << endl;
  
  test_MinReduction();
  cout << "BVH tests: PASSED" << endl;
}

__host__ void test_MinReduction()
{
  srand(static_cast<unsigned> (time(NULL)));

  Task h_tasks[NUM_TASKS];
  Task *d_tasks;
  float h_min;
  float *d_min;

  cudaError_t err = cudaSuccess;


  err = cudaMalloc(&d_tasks, NUM_TASKS*sizeof(Task));
  if (err != cudaSuccess)
  {
      fprintf(stderr, "Failed to allocate device vector A (error code %s)!\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
  }

  cudaMalloc(&d_min, sizeof(float));

  // RANDOM INIT oF TASKS
  for(int i = 0; i < NUM_TASKS; i++)
  {
    int i1 = rand()%2000;
    int i2 = rand()%2000;
    float dist = (rand()%20000 + 100)/10.0;
    h_tasks[i] = Task({i1, i2, dist});
  }

  // RANDOM INIT of start and num
  vector<int> test_start;
  vector<int> test_num;

  test_start.push_back(0);
  test_num.push_back(NUM_TASKS);

  for(int i = 0; i < 20; i++)
  {
    test_start.push_back(rand()%2000);
    test_num.push_back(rand()%(NUM_TASKS - 2000) + 1);
  }


  cudaMemcpy(d_tasks, h_tasks, NUM_TASKS*sizeof(Task), cudaMemcpyHostToDevice);

  dim3 dimBlockTasks(1, 1024);
  dim3 dimGridTasks(1, 1);
  dimGridTasks.y = (NUM_TASKS - 1)/1024 + 1;

  for(int i = 0; i < test_start.size(); i++)
  {
    int start = test_start[i];
    int num = test_num[i];

    float min_dist = 1e38;
    for(int j = 0; j < num; j++)
      min_dist = min(min_dist, h_tasks[j + start].dist);

    // computeMin<<<dimGridTasks, dimBlockTasks>>>(d_tasks, start, num, d_min);
    computeMin<<<dimGridTasks, dimBlockTasks>>>(start, num, d_min);

    cudaDeviceSynchronize();
    cudaMemcpy(&h_min, d_min, sizeof(float), cudaMemcpyDeviceToHost);
    cout << h_min;


    if(!approx_equal(min_dist, h_min))
    {
      cout << "EXP: " << min_dist << ", ACTUAL: " << h_min << endl;
      assert(false);
    }
  }

  cudaFree(d_tasks);
  cudaFree(d_min);

  cout << "test Minimum computation : PASSED" << endl;
}
