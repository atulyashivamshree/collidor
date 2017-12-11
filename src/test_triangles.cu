// @author : Atulya Shivam Shree
// @description : checks triangle intersection on CUDA

#include "compile_CUDA.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>

#include "Triangle.h"
float distTrianglesGPU(const Triangle*, const Triangle* , DistTriangleVars*);

// use the test cases on the GPU triangles function
#define DIST_TRIANGLES distTrianglesGPU
#include "Triangles_test.h"

using namespace std;

const int size_tri = sizeof(Triangle);
const int BLOCKSIZE = 32;

double get_wall_time();
// single case timing
void testSingleTiming();
// multiple instances of same in sequence
void testMultipleSerial();
// multiple instances of same in an array
void testMultipleSame();
// Check collision within an array of random triangles
void testMultipleRandom();

struct Result
{
  float dist;
};

__global__ void computeDistance(const Triangle *s1, const Triangle* s2, Result* res)
{
  __shared__ Triangle loc_s1;
  __shared__ Triangle loc_s2;
  __shared__ DistTriangleVars vars;
  loc_s1 = *s1;
  loc_s2 = *s2;
  // res->dist = 1e-6 + distTriangles(s1, s2, &vars);
  float dist = distTriangles(&loc_s1, &loc_s2, &vars);
  res->dist = dist;

}

__global__ void computeDistanceArray(
                    const Triangle *arr_s1, const Triangle* arr_s2, 
                      Result* arr_res, int n)
{
  int t_j = threadIdx.y;
  int g_j = blockIdx.y * blockDim.y + threadIdx.y;

  if(g_j < n)
  {
    __shared__ Triangle s1[BLOCKSIZE];
    __shared__ Triangle s2[BLOCKSIZE];
    __shared__ DistTriangleVars vars[BLOCKSIZE];
    s1[t_j] = arr_s1[g_j];
    s2[t_j] = arr_s2[g_j];
    // res->dist = 1e-6 + distTriangles(s1, s2, &vars);
    float dist = distTriangles(&s1[t_j], &s2[t_j], &vars[t_j]);
    arr_res[g_j].dist = dist;
  }

}

int main(int argc, char *argv[])
{
  
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  cout << "Device name: " << prop.name << endl;
  cout << "float size is " << sizeof(float) << endl;

  test_triangles_2D();
  test_triangles_3D();
  test_stress_random();

  testSingleTiming();
  testMultipleSerial();
  testMultipleSame();
  testMultipleRandom();
}

HOST_PREFIX float distTrianglesGPU(const Triangle* h_s1, const Triangle* h_s2, DistTriangleVars*)
{
  Triangle *d_s1, *d_s2;

  cudaMalloc(&d_s1, size_tri);
  cudaMalloc(&d_s2, size_tri);
  Result *d_res;
  Result h_res;
  cudaMalloc(&d_res, sizeof(Result));

  cudaMemcpy(d_s1, h_s1, size_tri, cudaMemcpyHostToDevice);
  cudaMemcpy(d_s2, h_s2, size_tri, cudaMemcpyHostToDevice);

  computeDistance<<<1, 1>>>(d_s1, d_s2, d_res);

  cudaMemcpy(&h_res, d_res, sizeof(Result), cudaMemcpyDeviceToHost);

  cudaFree(d_s1);
  cudaFree(d_s2);
  cudaFree(d_res);

  return h_res.dist;

}

void testSingleTiming()
{
  // initialize the dimension blocks for the differnet kernels
  double t_start = get_wall_time();

  // A is the input to every stage while C is the output after every stage
  Triangle h_s1, h_s2;
  Triangle *d_s1, *d_s2;

  generateRandomTriangle(&h_s1);
  generateRandomTriangle(&h_s2);

  cudaMalloc(&d_s1, size_tri);
  cudaMalloc(&d_s2, size_tri);
  Result *d_res;
  Result h_res;
  cudaMalloc(&d_res, sizeof(Result));
  
  double t_init = get_wall_time();

  cudaMemcpy(d_s1, &h_s1, size_tri, cudaMemcpyHostToDevice);
  cudaMemcpy(d_s2, &h_s2, size_tri, cudaMemcpyHostToDevice);

  computeDistance<<<1, 1>>>(d_s1, d_s2, d_res);

  cudaDeviceSynchronize();
  double t_cuda_end = get_wall_time();

  cudaMemcpy(&h_res, d_res, sizeof(Result), cudaMemcpyDeviceToHost);

  // float t_results;
  // cudaEventElapsedTime(&t_results, start, results);

  cout << setprecision(5);
  // cout << "Time for results using CUDA Event: " << t_results << "ms" <<  endl;
  cout << "Wall time (Single): " << (t_cuda_end - t_init )*1000 << "ms" <<  endl;
  cout << endl;

  cudaFree(d_s1);
  cudaFree(d_s2);
  cudaFree(d_res);

}

void testMultipleSerial()
{
  srand(static_cast<unsigned> (time(NULL)));
  // initialize the dimension blocks for the differnet kernels
  double t_start = get_wall_time();

  // A is the input to every stage while C is the output after every stage
  Triangle h_s1, h_s2;
  Triangle *d_s1, *d_s2;

  generateRandomTriangle(&h_s1);
  generateRandomTriangle(&h_s2);

  cudaMalloc(&d_s1, size_tri);
  cudaMalloc(&d_s2, size_tri);
  Result *d_res;
  Result h_res;
  cudaMalloc(&d_res, sizeof(Result));
  
  double t_init = get_wall_time();

  cudaMemcpy(d_s1, &h_s1, size_tri, cudaMemcpyHostToDevice);
  cudaMemcpy(d_s2, &h_s2, size_tri, cudaMemcpyHostToDevice);

  for(int i = 0; i < NUM_CHECK; i++)
    computeDistance<<<1, 1>>>(d_s1, d_s2, d_res);

  cudaDeviceSynchronize();
  double t_cuda_end = get_wall_time();

  cudaMemcpy(&h_res, d_res, sizeof(Result), cudaMemcpyDeviceToHost);

  // float t_results;
  // cudaEventElapsedTime(&t_results, start, results);

  cout << setprecision(5);
  cout << "Wall time multiple serial (" << NUM_CHECK << "): " << (t_cuda_end - t_init )*1000 << "ms" <<  endl;
  cout << endl;

  cudaFree(d_s1);
  cudaFree(d_s2);
  cudaFree(d_res);

}

void testMultipleSame()
{
  srand(static_cast<unsigned> (time(NULL)));
  // initialize the dimension blocks for the differnet kernels
  double t_start = get_wall_time();

  // A is the input to every stage while C is the output after every stage
  Triangle h_s1[NUM_CHECK], h_s2[NUM_CHECK];
  Triangle* d_s1, *d_s2;

  generateRandomTriangle(&h_s1[0]);
  generateRandomTriangle(&h_s2[0]);

  for(int i = 1; i < NUM_CHECK-1; i++)
  {
    h_s1[i] = h_s1[0];
    h_s2[i] = h_s2[0];
  }

  float actual_res = distTriangles_fcl(h_s1[0], h_s2[0]);
  // cout << "s1 " << h_s1[0];
  // cout << "s2 " << h_s2[0];
  // cout << "dist is " << actual_res << endl;

  // cout << "L: s1 " << h_s1[NUM_CHECK-1];
  // cout << "L: s2 " << h_s2[NUM_CHECK-1];
  // cout << "L: dist is " << distTriangles_fcl(h_s1[NUM_CHECK-1],
  //                            h_s2[NUM_CHECK-1]) << endl;

  cudaMalloc(&d_s1, NUM_CHECK*size_tri);
  cudaMalloc(&d_s2, NUM_CHECK*size_tri);
  Result *d_res;
  Result h_res[NUM_CHECK];
  cudaMalloc(&d_res, NUM_CHECK*sizeof(Result));
  
  double t_init = get_wall_time();

  cudaMemcpy(d_s1, h_s1, NUM_CHECK*size_tri, cudaMemcpyHostToDevice);
  cudaMemcpy(d_s2, h_s2, NUM_CHECK*size_tri, cudaMemcpyHostToDevice);

  int numBlocks = (NUM_CHECK - 1)/BLOCKSIZE + 1;
  dim3 dimBlock(1, BLOCKSIZE);
  dim3 dimGrid(1, numBlocks);
  computeDistanceArray<<<dimGrid, dimBlock>>>(d_s1, d_s2, d_res, NUM_CHECK);

  cudaDeviceSynchronize();
  double t_cuda_end = get_wall_time();

  cudaMemcpy(&h_res, d_res, NUM_CHECK*sizeof(Result), cudaMemcpyDeviceToHost);

  // float t_results;
  // cudaEventElapsedTime(&t_results, start, results);
  int count_correct = 0;
  for(int i = 0; i < NUM_CHECK; i++)
  {
    if(approx_equal(h_res[i].dist, actual_res))
      count_correct++;
    // cout << "actual: " << actual_res << " obtained: " << h_res[i].dist << endl;
  }

  cout << count_correct << "/" << NUM_CHECK << " are correct" << endl;

  cout << setprecision(5);
  cout << "Wall time multiple same(" << NUM_CHECK << "): " << (t_cuda_end - t_init )*1000 << "ms" <<  endl;
  cout << endl;

  cudaFree(d_s1);
  cudaFree(d_s2);
  cudaFree(d_res);

}

void testMultipleRandom()
{
  srand(static_cast<unsigned> (time(NULL)));
  // initialize the dimension blocks for the differnet kernels
  double t_start = get_wall_time();

  // A is the input to every stage while C is the output after every stage
  Triangle h_s1[NUM_CHECK], h_s2[NUM_CHECK];
  Triangle* d_s1, *d_s2;

  for(int i = 0; i < NUM_CHECK; i++)
  {
    generateRandomTriangle(&h_s1[i]);
    generateRandomTriangle(&h_s2[i]);
  }

  float actual_res = distTriangles_fcl(h_s1[0], h_s2[0]);

  cudaMalloc(&d_s1, NUM_CHECK*size_tri);
  cudaMalloc(&d_s2, NUM_CHECK*size_tri);
  Result *d_res;
  Result h_res[NUM_CHECK];
  cudaMalloc(&d_res, NUM_CHECK*sizeof(Result));
  
  double t_init = get_wall_time();

  cudaMemcpy(d_s1, h_s1, NUM_CHECK*size_tri, cudaMemcpyHostToDevice);
  cudaMemcpy(d_s2, h_s2, NUM_CHECK*size_tri, cudaMemcpyHostToDevice);

  int numBlocks = (NUM_CHECK - 1)/BLOCKSIZE + 1;
  dim3 dimBlock(1, BLOCKSIZE);
  dim3 dimGrid(1, numBlocks);
  computeDistanceArray<<<dimGrid, dimBlock>>>(d_s1, d_s2, d_res, NUM_CHECK);

  cudaDeviceSynchronize();
  double t_cuda_end = get_wall_time();

  cudaMemcpy(&h_res, d_res, NUM_CHECK*sizeof(Result), cudaMemcpyDeviceToHost);

  // float t_results;
  // cudaEventElapsedTime(&t_results, start, results);
  int count_correct = 0;
  for(int i = 0; i < NUM_CHECK; i++)
  {
    float actual_dist = distTriangles_fcl(h_s1[i], h_s2[i]);
    if(approx_equal(h_res[i].dist, actual_dist))
      count_correct++;
    // cout << "actual: " << actual_res << " obtained: " << h_res[i].dist << endl;
  }

  cout << count_correct << "/" << NUM_CHECK << " are correct" << endl;

  cout << setprecision(5);
  cout << "Wall time multiple random(" << NUM_CHECK << "): " << (t_cuda_end - t_init )*1000 << "ms" <<  endl;
  cout << endl;

  cudaFree(d_s1);
  cudaFree(d_s2);
  cudaFree(d_res);

}

HOST_PREFIX double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}