/*
 *  vector_math-inl.h
    Description: Implements the basic math functions

    @author Atulya Shivam Shree
    Created on: Dec 11, 2017
    Copyright (c) 2017 Atulya Shivam Shree
    TODO(atulya) : remove this file if no longer required
 */
#ifndef SRC_VECTOR_MATH_INL_H_
#define SRC_VECTOR_MATH_INL_H_

CUDA_PREFIX float dot(const Vector3* a, const Vector3* b) {
  float res = 0;
  res += a->x * b->x;
  res += a->y * b->y;
  res += a->z * b->z;
  return res;
}

#endif  // SRC_VECTOR_MATH_INL_H_
