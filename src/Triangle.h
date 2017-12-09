#include "vector_math.h"

#ifndef COLLIDOR_SRC_TRIANGLE
#define COLLIDOR_SRC_TRIANGLE

// holds a triangle in 3D space
struct Triangle
{
  Vector3 a;
  Vector3 b;
  Vector3 c;
};

float distTriangles(const Triangle* s1, const Triangle *s2);

#include "Triangle-inl.h"

#endif