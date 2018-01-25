/*
 *  Tringle-cuda-inl.h
    Description: CUDA implementations for the Triangle functions

    @author Atulya Shivam Shree
    Created on: Dec 11, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#ifndef SRC_TRIANGLE_CUDA_INL_H_
#define SRC_TRIANGLE_CUDA_INL_H_

const int BLOCKSIZE_TRI = 32;

__device__ __inline__ void deepCopy3x3(float target[3][3],
                                       const float mat[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) target[i][j] = mat[i][j];
}
__device__ __inline__ void deepCopy3x1(float target[3], const float mat[3]) {
  for (int i = 0; i < 3; i++) target[i] = mat[i];
}

__device__ void computeDistance(const Triangle* s1, const Triangle* s2,
                                const float R[3][3], const float t[3],
                                TriangleResult* res) {
  Triangle loc_s1;
  Triangle loc_s2;
  DistTriangleVars vars;
  loc_s1 = *s1;
  loc_s2 = *s2;
  // res->dist = 1e-6 + distTriangles(, &vars);
  float dist = distTriangles(&loc_s1, &loc_s2, R, t, &vars);
  res->dist = dist;
}

__global__ void computeDistanceSingle(const Triangle* s1, const Triangle* s2,
                                      TriangleResult* res) {
  __shared__ Triangle loc_s1;
  __shared__ Triangle loc_s2;
  __shared__ DistTriangleVars vars;
  float R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  float t[3] = {0, 0, 0};
  loc_s1 = *s1;
  loc_s2 = *s2;
  // res->dist = 1e-6 + distTriangles(s1, s2, &vars);
  float dist = distTriangles(&loc_s1, &loc_s2, R, t, &vars);
  res->dist = dist;
}

__global__ void computeDistanceArray(const Triangle* arr_s1,
                                     const Triangle* arr_s2,
                                     TriangleResult* arr_res, int n) {
  int t_j = threadIdx.y;
  int g_j = blockIdx.y * blockDim.y + threadIdx.y;

  if (g_j < n) {
    __shared__ Triangle s1[BLOCKSIZE_TRI];
    __shared__ Triangle s2[BLOCKSIZE_TRI];
    __shared__ DistTriangleVars vars[BLOCKSIZE_TRI];
    float R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    float t[3] = {0, 0, 0};
    s1[t_j] = arr_s1[g_j];
    s2[t_j] = arr_s2[g_j];
    // res->dist = 1e-6 + distTriangles(s1, s2, &vars);
    float dist = distTriangles(&s1[t_j], &s2[t_j], R, t, &vars[t_j]);
    arr_res[g_j].dist = dist;
  }
}

#endif  // SRC_TRIANGLE_CUDA_INL_H_
