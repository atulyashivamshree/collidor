/*
 * Triangle-cuda-inl.h
 *
 *  Created on: Dec 12, 2017
 *      Author: atulya
 */

#ifndef TRIANGLE_CUDA_INL_H_
#define TRIANGLE_CUDA_INL_H_

const int BLOCKSIZE_TRI = 32;

__device__ void computeDistance(const Matrix3* R, const Vector3* t,
          const Triangle *s1, const Triangle* s2, TriangleResult* res)
{
  Triangle loc_s1;
  Triangle loc_s2;
  DistTriangleVars vars;
  loc_s1 = *s1;
  loc_s2 = *s2;
  // res->dist = 1e-6 + distTriangles(s1, s2, &vars);
  float dist = distTriangles(&loc_s1, &loc_s2, &vars);
  res->dist = dist;
}

__global__ void computeDistanceSingle(const Triangle *s1, const Triangle* s2, TriangleResult* res)
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
                      TriangleResult* arr_res, int n)
{
  int t_j = threadIdx.y;
  int g_j = blockIdx.y * blockDim.y + threadIdx.y;

  if(g_j < n)
  {
    __shared__ Triangle s1[BLOCKSIZE_TRI];
    __shared__ Triangle s2[BLOCKSIZE_TRI];
    __shared__ DistTriangleVars vars[BLOCKSIZE_TRI];
    s1[t_j] = arr_s1[g_j];
    s2[t_j] = arr_s2[g_j];
    // res->dist = 1e-6 + distTriangles(s1, s2, &vars);
    float dist = distTriangles(&s1[t_j], &s2[t_j], &vars[t_j]);
    arr_res[g_j].dist = dist;
  }
}




#endif /* TRIANGLE_CUDA_INL_H_ */
