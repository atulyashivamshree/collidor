/*
 *  RSS-cuda-inl.h
    Description: Implements functions for computing distance between two RSS

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#ifndef SRC_RSS_CUDA_INL_H_
#define SRC_RSS_CUDA_INL_H_

const int BLOCKSIZE_RECT = 32;

__device__ void computeDistance(const RSS* r1, const RSS* d2,
                                const float R[3][3], const float t[3],
                                RSSResult* res) {
  float loc_R[3][3];
  float loc_t[3];
  RSS loc_r1;
  RSS loc_d2;
  DistRSSVars vars;

  deepCopy3x3(loc_R, R);
  deepCopy3x1(loc_t, t);

  loc_r1 = *r1;
  loc_d2 = *d2;
  // res->dist = 1e-6 + distRSSs(r1, d2, &vars);
  float dist = rssDistance(loc_R, loc_t, &loc_r1, &loc_d2, &vars);
  res->dist = dist;
}


__global__ void computeDistanceArray(const Matrix3* R, const Vector3* t,
                                     const RSS* arr_r1, const RSS* arr_d2,
                                     RSSResult* arr_res, int n) {
  int t_j = threadIdx.y;
  int g_j = blockIdx.y * blockDim.y + threadIdx.y;

  __shared__ float loc_R[3][3];
  __shared__ float loc_t[3];
  __shared__ RSS r1[BLOCKSIZE_RECT];
  __shared__ RSS d2[BLOCKSIZE_RECT];
  __shared__ DistRSSVars vars[BLOCKSIZE_RECT];

  if (threadIdx.y == 0) {
    // deepCopy3x3(loc_R, R);
    // deepCopy3x1(loc_t, t);
    // loc_R = *R;
    // loc_t = *t;
  }

  if (g_j < n) {
    r1[t_j] = arr_r1[g_j];
    d2[t_j] = arr_d2[g_j];
    // res->dist = 1e-6 + distRSSs(r1, d2, &vars);
    float dist = rssDistance(loc_R, loc_t, &r1[t_j], &d2[t_j], &vars[t_j]);
    arr_res[g_j].dist = dist;
  }
}

#endif  // SRC_RSS_CUDA_INL_H_
