#ifndef COLLIDOR_SRC_VECTOR_MATH_H_
#define COLLIDOR_SRC_VECTOR_MATH_H_

struct Vector3
{
  float x;
  float y;
  float z;
};

struct Matrix3
{
  Vector3 v1;
  Vector3 v2;
  Vector3 v3;
};

CUDA_PREFIX void add(const Vector3* a, const Vector3* b, Vector3* c);

CUDA_PREFIX void sub(const Vector3* a, const Vector3* b, Vector3* c);

CUDA_PREFIX float dot(const Vector3* a, const Vector3* b);

CUDA_PREFIX void cross(const Vector3* a, const Vector3* b, Vector3* c);

CUDA_PREFIX void mult(const Matrix3* a, const Vector3* b, Vector3* c);

CUDA_PREFIX void transform(const Matrix3* R, const Vector3* b, 
                                const Vector3* x, Vector3* y);

#include "vector_math-inl.h"

#endif