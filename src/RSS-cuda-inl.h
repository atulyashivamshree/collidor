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

__device__ void computeDistance(const RSS* r1, const RSS* r2,
                                const float R[3][3], const float t[3],
                                RSSResult* res) {
  float loc_R[3][3];
  float loc_t[3];
  RSS loc_r1;
  RSS loc_r2;
  DistRSSVars vars;

  deepCopy3x3(loc_R, R);
  deepCopy3x1(loc_t, t);

  loc_r1 = *r1;
  loc_r2 = *r2;
  // res->dist = 1e-6 + distRSSs(r1, r2, &vars);
  float dist = rssDistance(loc_R, loc_t, &loc_r1, &loc_r2, &vars);
  res->dist = dist;
}

__global__ void computeDistanceSingle(const RSS* r1, const RSS* r2,
                                      const Config* cfg, RSSResult* res) {
  float loc_R[3][3];
  float loc_t[3];

  deepCopy3x3(loc_R, cfg->R);
  deepCopy3x1(loc_t, cfg->t);

  computeDistance(r1, r2, loc_R, loc_t, res);
}

__global__ void computeDistanceArray(const RSS* arr_r1, const RSS* arr_r2,
                                     const Config* cfg, RSSResult* arr_res,
                                     int n) {
  int g_j = blockIdx.y * blockDim.y + threadIdx.y;

  float loc_R[3][3];
  float loc_t[3];

  deepCopy3x3(loc_R, cfg->R);
  deepCopy3x1(loc_t, cfg->t);

  if (g_j < n) {
    computeDistance(arr_r1 + g_j, arr_r2 + g_j, loc_R, loc_t, arr_res + g_j);
  }
}

#endif  // SRC_RSS_CUDA_INL_H_
