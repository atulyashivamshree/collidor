/*
 *  vector_math-inl.h
    Description: Describes the basic data structures for math

    @author Atulya Shivam Shree
    Created on: Dec 11, 2017
    Copyright (c) 2017 Atulya Shivam Shree
    TODO(atulya) : remove this file if no longer required
 */

#ifndef SRC_VECTOR_MATH_H_
#define SRC_VECTOR_MATH_H_

struct Vector3 {
  float x;
  float y;
  float z;
};

struct Matrix3 {
  Vector3 v1;
  Vector3 v2;
  Vector3 v3;
};

#include "vector_math-inl.h"

#endif  // SRC_VECTOR_MATH_H_
