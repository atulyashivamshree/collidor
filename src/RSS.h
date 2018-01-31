/*
 *  RSS.h
    Description: Implements the RSS data structure and the utility functions for
 it

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#include "Triangle.h"

#ifndef SRC_RSS_H_
#define SRC_RSS_H_

// holds the RSS object
struct RSS {
  // Orientation of the RSS
  Matrix3 axis;

  // Orgiin of the rectangle in RSS
  Vector3 To;

  // side length of the rectangle
  float l[2];

  // radius of the sphere
  float r;

  // volume of the unit
  float size;
};

// stores temporary vals implements the RSS distance computation
struct DistRSSVars {
  DistTriangleVars dist_triangle_vars;
};

struct RSSResult {
  float id;
  float dist;
};

CUDA_PREFIX float rssDistance(const float R[3][3], const float t[3],
                              const RSS* a, const RSS* b, DistRSSVars* d);

#include "RSS-inl.h"

#endif  // SRC_RSS_H_
