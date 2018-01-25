/*
 *  test_vector_math.cpp
    Description: Checks the basic working of vector math functions

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */
#include "../utils/compile_CPP.h"

#include <Eigen/Dense>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include "../RSS.h"

using std::cout;
using std::endl;
using Eigen::Vector3f;

const int NUM_CHECK = 100;
const float EPSILON = 1e-7;

void test_dot();
void test_add();
void test_sub();
void test_cross();
void test_mult();
void test_transform();
void test_approx();

bool approx_equal(float a, float b, float epsilon = EPSILON);

int main(int argc, char *argv[]) {
  // call all test functions
  test_approx();
  test_dot();
}

void test_approx() {
  assert(!approx_equal(1, 2));
  assert(!approx_equal(1, 1 + 1e-5));
  assert(approx_equal(1, 1 + 1e-10));
  assert(!approx_equal(1, 1 - 1e-5));
  assert(approx_equal(1, 1 - 1e-10));
  cout << "TEST Approx : PASSED" << endl;
}

void test_dot() {
  for (int i = 0; i < NUM_CHECK; i++) {
    Vector3f v1 = Vector3f::Random(3, 1);
    Vector3f v2 = Vector3f::Random(3, 1);

    Vector3 lib_v1({v1(0), v1(1), v1(2)});
    Vector3 lib_v2({v2(0), v2(1), v2(2)});

    float correct = v1.dot(v2);
    float actual = dot(&lib_v1, &lib_v2);

    assert(approx_equal(actual, correct));
  }

  cout << "TEST Dot product : PASSED" << endl;
}

bool approx_equal(float a, float b, float epsilon) {
  return (std::abs(a - b) < epsilon);
}
