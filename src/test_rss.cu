// @author : Atulya Shivam Shree
// @description : checks triangle intersection on CUDA

#include "compile_CUDA.h"

#include <cassert>
#include <iostream>
#include <iomanip>

#include "RSS.h"
float distTrianglesGPU(const Triangle*, const Triangle* , DistTriangleVars*);

float distRSSGPU(const Matrix3* R, const Vector3* t,
            const RSS* a, const RSS* b, DistRSSVars* d);

#define DIST_TRIANGLES distTrianglesGPU
// use the test cases on the GPU triangles function
#define DIST_RSS distRSSGPU
#include "Rectangle_tests.h"

#include "RSS-cuda-inl.h"

using namespace std;

const int size_rss = sizeof(RSS);

// single case timing
void testSingleTiming();
// multiple instances of same in sequence
void testMultipleSerial();
// multiple instances of same in an array
void testMultipleSame();
// Check collision within an array of random triangles
void testMultipleRandom();

int main(int argc, char *argv[])
{
  
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  cout << "Device name: " << prop.name << endl;
  cout << "float size is " << sizeof(float) << endl;

  test_rectangles_2D();
  test_rectangles_3D();
  test_stress_random_RSS();

  print_Stats();

  testSingleTiming();
  testMultipleSerial();
  testMultipleSame();
  testMultipleRandom();
}

HOST_PREFIX float distTrianglesGPU(const Triangle* h_s1, const Triangle* h_s2, DistTriangleVars*)
{
  return 0;
}

// use this function to verify working by running it against the test cases
HOST_PREFIX float distRSSGPU(const Matrix3* R, const Vector3* t,
            const RSS* h_r1, const RSS* h_r2, DistRSSVars* d)
{
  RSS *d_r1, *d_r2;
  Matrix3 *d_R;
  Vector3 *d_t;

  cudaMalloc(&d_R, sizeof(Matrix3));
  cudaMalloc(&d_t, sizeof(Vector3));
  cudaMalloc(&d_r1, size_rss);
  cudaMalloc(&d_r2, size_rss);
  RSSResult *d_res;
  RSSResult h_res;
  cudaMalloc(&d_res, sizeof(RSSResult));

  cudaMemcpy(d_R, R, sizeof(Matrix3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_t, t, sizeof(Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_r1, h_r1, size_rss, cudaMemcpyHostToDevice);
  cudaMemcpy(d_r2, h_r2, size_rss, cudaMemcpyHostToDevice);

  computeDistance<<<1, 1>>>(d_R, d_t, d_r1, d_r2, d_res);

  cudaMemcpy(&h_res, d_res, sizeof(RSSResult), cudaMemcpyDeviceToHost);

  cudaFree(d_r1);
  cudaFree(d_r2);
  cudaFree(d_res);

  return h_res.dist;

}

void testSingleTiming()
{
  // initialize the dimension blocks for the differnet kernels
  double t_start = get_wall_time();

  // A is the input to every stage while C is the output after every stage
  RSS h_r1, h_r2;
  RSS *d_r1, *d_r2;

  //TODO change this to be more general

  h_r1 = getRandomRSS();
  h_r2 = getRandomRSS();

  Matrix3 *d_R;
  Vector3 *d_t;

  cudaMalloc(&d_R, sizeof(Matrix3));
  cudaMalloc(&d_t, sizeof(Vector3));
  cudaMalloc(&d_r1, size_rss);
  cudaMalloc(&d_r2, size_rss);
  RSSResult *d_res;
  RSSResult h_res;
  cudaMalloc(&d_res, sizeof(RSSResult));
  
  double t_init = get_wall_time();

  cudaMemcpy(d_R, &MAT_I, sizeof(Matrix3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_t, &t0, sizeof(Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_r1, &h_r1, size_rss, cudaMemcpyHostToDevice);
  cudaMemcpy(d_r2, &h_r2, size_rss, cudaMemcpyHostToDevice);

  computeDistance<<<1, 1>>>(d_R, d_t, d_r1, d_r2, d_res);

  cudaDeviceSynchronize();
  double t_cuda_end = get_wall_time();

  cudaMemcpy(&h_res, d_res, sizeof(RSSResult), cudaMemcpyDeviceToHost);

  // float t_results;
  // cudaEventElapsedTime(&t_results, start, results);

  cout << setprecision(5);
  // cout << "Time for results using CUDA Event: " << t_results << "ms" <<  endl;
  cout << "Wall time (Single): " << (t_cuda_end - t_init )*1000 << "ms" <<  endl;
  cout << endl;

  cudaFree(d_r1);
  cudaFree(d_r2);
  cudaFree(d_res);

}

void testMultipleSerial()
{
  // initialize the dimension blocks for the differnet kernels
  double t_start = get_wall_time();

  // A is the input to every stage while C is the output after every stage
  RSS h_r1, h_r2;
  RSS *d_r1, *d_r2;

  //TODO change this to be more general

  h_r1 = getRandomRSS();
  h_r2 = getRandomRSS();

  Matrix3 *d_R;
  Vector3 *d_t;

  cudaMalloc(&d_R, sizeof(Matrix3));
  cudaMalloc(&d_t, sizeof(Vector3));
  cudaMalloc(&d_r1, size_rss);
  cudaMalloc(&d_r2, size_rss);
  RSSResult *d_res;
  RSSResult h_res;
  cudaMalloc(&d_res, sizeof(RSSResult));
  
  double t_init = get_wall_time();

  cudaMemcpy(d_R, &MAT_I, sizeof(Matrix3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_t, &t0, sizeof(Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_r1, &h_r1, size_rss, cudaMemcpyHostToDevice);
  cudaMemcpy(d_r2, &h_r2, size_rss, cudaMemcpyHostToDevice);

  for(int i = 0; i < 1000; i++)
    computeDistance<<<1, 1>>>(d_R, d_t, d_r1, d_r2, d_res);

  cudaDeviceSynchronize();
  double t_cuda_end = get_wall_time();

  cudaMemcpy(&h_res, d_res, sizeof(RSSResult), cudaMemcpyDeviceToHost);

  // float t_results;
  // cudaEventElapsedTime(&t_results, start, results);

  cout << setprecision(5);
  // cout << "Time for results using CUDA Event: " << t_results << "ms" <<  endl;
  cout << "Wall time multiple serial (" << 1000 << "): " << (t_cuda_end - t_init )*1000 << "ms" <<  endl;
  cout << endl;

  cudaFree(d_r1);
  cudaFree(d_r2);
  cudaFree(d_res);

}

void testMultipleSame()
{
  srand(static_cast<unsigned> (time(NULL)));
  // initialize the dimension blocks for the differnet kernels
  double t_start = get_wall_time();

  // A is the input to every stage while C is the output after every stage
  RSS *h_r1, *h_r2;
  RSS* d_r1, *d_r2;

  h_r1 = new RSS[NUM_CHECK];
  h_r2 = new RSS[NUM_CHECK];

  if(!h_r1 || !h_r2)
  {
    cout << "Could not initialize h_r, out of memory" << endl;
    return;
  }

  h_r1[0] = getRandomRSS();
  h_r2[0] = getRandomRSS();

  for(int i = 1; i < NUM_CHECK-1; i++)
  {
    h_r1[i] = h_r1[0];
    h_r2[i] = h_r2[0];
  }

  float actual_res = distRectangles_fcl(h_r1[0], h_r2[0]);
  // cout << "s1 " << h_r1[0];
  // cout << "s2 " << h_r2[0];
  // cout << "dist is " << actual_res << endl;

  // cout << "L: s1 " << h_r1[NUM_CHECK-1];
  // cout << "L: s2 " << h_r2[NUM_CHECK-1];
  // cout << "L: dist is " << distRSSs_fcl(h_r1[NUM_CHECK-1],
  //                            h_r2[NUM_CHECK-1]) << endl;
  Matrix3 *d_R;
  Vector3 *d_t;

  cudaMalloc(&d_R, sizeof(Matrix3));
  cudaMalloc(&d_t, sizeof(Vector3));
  cudaMalloc(&d_r1, NUM_CHECK*size_rss);
  cudaMalloc(&d_r2, NUM_CHECK*size_rss);
  RSSResult *d_res;
  RSSResult *h_res;
  h_res = new RSSResult[NUM_CHECK];
  cudaMalloc(&d_res, NUM_CHECK*sizeof(RSSResult));
  
  double t_init = get_wall_time();

  cudaMemcpy(d_R, &MAT_I, sizeof(Matrix3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_t, &t0, sizeof(Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_r1, h_r1, NUM_CHECK*size_rss, cudaMemcpyHostToDevice);
  cudaMemcpy(d_r2, h_r2, NUM_CHECK*size_rss, cudaMemcpyHostToDevice);

  int numBlocks = (NUM_CHECK - 1)/BLOCKSIZE + 1;
  cout << "numBlocks is " << numBlocks << endl;
  dim3 dimBlock(1, BLOCKSIZE);
  dim3 dimGrid(1, numBlocks);
  computeDistanceArray<<<dimGrid, dimBlock>>>(d_R, d_t, d_r1, d_r2, d_res, NUM_CHECK);

  cudaDeviceSynchronize();
  double t_cuda_end = get_wall_time();

  cudaMemcpy(h_res, d_res, NUM_CHECK*sizeof(RSSResult), cudaMemcpyDeviceToHost);

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

  cudaFree(d_r1);
  cudaFree(d_r2);
  cudaFree(d_res);

  delete[] h_r1;
  delete[] h_r2;
  delete[] h_res;
}

void testMultipleRandom()
{
  srand(static_cast<unsigned> (time(NULL)));
  // initialize the dimension blocks for the differnet kernels
  double t_start = get_wall_time();

  // A is the input to every stage while C is the output after every stage
  RSS *h_r1, *h_r2;
  RSS* d_r1, *d_r2;

  h_r1 = new RSS[NUM_CHECK];
  h_r2 = new RSS[NUM_CHECK];

  if(!h_r1 || !h_r2)
  {
    cout << "Could not initialize h_r, out of memory" << endl;
    return;
  }

  cout << "Initializing the RSS ... " ;
  for(int i = 1; i < NUM_CHECK-1; i++)
  {
    h_r1[i] = getRandomRSS();
    h_r2[i] = getRandomRSS();
  }
  cout << "Done " << endl;
  // cout << "s1 " << h_r1[0];
  // cout << "s2 " << h_r2[0];
  // cout << "dist is " << actual_res << endl;

  // cout << "L: s1 " << h_r1[NUM_CHECK-1];
  // cout << "L: s2 " << h_r2[NUM_CHECK-1];
  // cout << "L: dist is " << distRSSs_fcl(h_r1[NUM_CHECK-1],
  //                            h_r2[NUM_CHECK-1]) << endl;
  Matrix3 *d_R;
  Vector3 *d_t;

  cudaMalloc(&d_R, sizeof(Matrix3));
  cudaMalloc(&d_t, sizeof(Vector3));
  cudaMalloc(&d_r1, NUM_CHECK*size_rss);
  cudaMalloc(&d_r2, NUM_CHECK*size_rss);
  RSSResult *d_res;
  RSSResult *h_res;
  h_res = new RSSResult[NUM_CHECK];
  cudaMalloc(&d_res, NUM_CHECK*sizeof(RSSResult));
  
  double t_init = get_wall_time();

  cudaMemcpy(d_R, &MAT_I, sizeof(Matrix3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_t, &t0, sizeof(Vector3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_r1, h_r1, NUM_CHECK*size_rss, cudaMemcpyHostToDevice);
  cudaMemcpy(d_r2, h_r2, NUM_CHECK*size_rss, cudaMemcpyHostToDevice);

  int numBlocks = (NUM_CHECK - 1)/BLOCKSIZE + 1;
  cout << "numBlocks is " << numBlocks << endl;
  dim3 dimBlock(1, BLOCKSIZE);
  dim3 dimGrid(1, numBlocks);
  computeDistanceArray<<<dimGrid, dimBlock>>>(d_R, d_t, d_r1, d_r2, d_res, NUM_CHECK);

  cudaDeviceSynchronize();
  double t_cuda_end = get_wall_time();

  cudaMemcpy(h_res, d_res, NUM_CHECK*sizeof(RSSResult), cudaMemcpyDeviceToHost);

  // float t_results;
  // cudaEventElapsedTime(&t_results, start, results);
  cout << "Evaluating reslts now ... " ;
  int count_correct = 0;
  for(int i = 0; i < NUM_CHECK; i++)
  {
    float actual_res = distRectangles_fcl(h_r1[i], h_r2[i]);

    if(approx_equal(h_res[i].dist, actual_res))
      count_correct++;
    // cout << "actual: " << actual_res << " obtained: " << h_res[i].dist << endl;
  }

  cout << " took " << (get_wall_time() - t_cuda_end)*1000 << "ms" << endl;

  cout << count_correct << "/" << NUM_CHECK << " are correct" << endl;
  cout << setprecision(5);
  cout << "Wall time multiple random(" << NUM_CHECK << "): " << (t_cuda_end - t_init )*1000 << "ms" <<  endl;
  cout << endl;

  cudaFree(d_r1);
  cudaFree(d_r2);
  cudaFree(d_res);

  delete[] h_r1;
  delete[] h_r2;
  delete[] h_res;
}