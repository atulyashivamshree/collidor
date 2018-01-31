/*
 *  Triangles_test.h
    Description: Defines the basic tests for verifying that the triangle
 distance computation is working fine

    @author Atulya Shivam Shree
    Created on: Dec 10, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#include <iomanip>
#include <iostream>
#include "Eigen/Dense"
#include "../RSS.h"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/triangle_distance.h"

#ifndef SRC_TEST_TRIANGLES_TEST_H_
#define SRC_TEST_TRIANGLES_TEST_H_

using std::cout;
using std::endl;

const float SPACE_LOW = -1.0;
const float SPACE_HIGH = 1.0;
const int NUM_CHECK = 10000;
const int STRESS_CHECK = 10000;
const float EPSILON = 5e-5;

const float matI_[3][3] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
const float matR1_[3][3] = {
    sqrtf(2) / 2, sqrtf(2) / 2, 0, -sqrtf(2) / 2, sqrtf(2) / 2, 0, 0, 0, 1};

const float t0_[3] = {0, 0, 0};
const float t1_[3] = {1, 0, 0};

HOST_PREFIX std::ostream& operator<<(std::ostream& os, Triangle tri);

HOST_PREFIX Eigen::Vector3f getVector3(Vector3 v) {
  return Eigen::Vector3f(v.x, v.y, v.z);
}

HOST_PREFIX bool approx_equal(float a, float b, float epsilon = EPSILON) {
  return (std::abs(a - b) < epsilon);
}

HOST_PREFIX float distTriangles_fcl(Triangle s1, Triangle s2,
                                    const float R[3][3], const float t[3]) {
  Eigen::Vector3f p, q;
  fcl::Matrix3f R_;
  fcl::Vector3f Tl_;
  R_(0, 0) = R[0][0];
  R_(0, 1) = R[0][1];
  R_(0, 2) = R[0][2];
  R_(1, 0) = R[1][0];
  R_(1, 1) = R[1][1];
  R_(1, 2) = R[1][2];
  R_(2, 0) = R[2][0];
  R_(2, 1) = R[2][1];
  R_(2, 2) = R[2][2];

  Tl_(0) = t[0];
  Tl_(1) = t[1];
  Tl_(2) = t[2];
  return fcl::detail::TriangleDistance<float>::triDistance(
      getVector3(s1.a), getVector3(s1.b), getVector3(s1.c), getVector3(s2.a),
      getVector3(s2.b), getVector3(s2.c), R_, Tl_, p, q);
}

HOST_PREFIX void print_Stats() {
  cout << "Size of RSS is: " << sizeof(RSS) << endl;
  cout << "Size Of LineSegVars is: " << sizeof(LineSegVars) << endl;
  cout << "Size Of TriDistVars is: " << sizeof(TriDistVars) << endl;
  cout << "Size Of DistTriangleVars is: " << sizeof(DistTriangleVars) << endl;
}

HOST_PREFIX void test_triangles_2D() {
  DistTriangleVars preset_var;

  // totally inside
  Triangle p1({{0, 6, 0}, {6, 0, 0}, {0, 0, 0}});
  Triangle q1({{1, 1, 0}, {2, 1, 0}, {2, 2, 0}});
  cout << "Triangle distance inside " << distTriangles_fcl(p1, q1, matI_, t0_)
       << endl;
  assert(approx_equal(distTriangles_fcl(p1, q1, matI_, t0_), 0));
  assert(approx_equal(DIST_TRIANGLES(&p1, &q1, matI_, t0_, &preset_var), 0));

  // first distance
  Triangle p2({{0, 2, 0}, {2, 0, 0}, {0, 0, 0}});
  Triangle q2({{2, 1, 0}, {2, 4, 0}, {6, 1, 0}});
  cout << "Triangle distance away " << distTriangles_fcl(p2, q2, matI_, t0_)
       << endl;
  assert(approx_equal(distTriangles_fcl(p2, q2, matI_, t0_), sqrtf(2) / 2));
  assert(approx_equal(DIST_TRIANGLES(&p2, &q2, matI_, t0_, &preset_var),
                      sqrtf(2) / 2));

  // first intersection
  Triangle p3({{7, 3, 0}, {6, 0, 0}, {0, 0, 0}});
  Triangle q3({{1, 1, 0}, {4, 4, 0}, {6, 2, 0}});
  cout << "Triangle distance intersect "
       << distTriangles_fcl(p3, q3, matI_, t0_) << endl;
  assert(approx_equal(distTriangles_fcl(p3, q3, matI_, t0_), 0));
  assert(approx_equal(DIST_TRIANGLES(&p3, &q3, matI_, t0_, &preset_var), 0));

  Triangle p4({{0, 2, 0}, {2, 0, 0}, {0, 0, 0}});
  Triangle q4({{2, 1, 0}, {2, 4, 0}, {6, 1, 0}});
  cout << "Triangle distance away " << distTriangles_fcl(p4, q4, matI_, t1_)
       << endl;
  assert(approx_equal(distTriangles_fcl(p4, q4, matI_, t1_), sqrtf(2)));
  assert(approx_equal(DIST_TRIANGLES(&p4, &q4, matI_, t1_, &preset_var),
                      sqrtf(2)));

  Triangle p5({{0, 2, 0}, {2, 0, 0}, {0, 0, 0}});
  Triangle q5({{2, 1, 0}, {2, 4, 0}, {6, 1, 0}});
  cout << "Triangle distance away " << distTriangles_fcl(q5, p5, matR1_, t0_)
       << endl;
  cout << "Triangle distance away "
       << DIST_TRIANGLES(&q5, &p5, matR1_, t0_, &preset_var) << endl;
  assert(approx_equal(distTriangles_fcl(q5, p5, matR1_, t0_), 2 - sqrtf(2)));
  assert(approx_equal(DIST_TRIANGLES(&q5, &p5, matR1_, t0_, &preset_var),
                      2 - sqrtf(2)));

  cout << "Test Triangeles 2D : PASSED" << endl;
}

HOST_PREFIX void test_triangles_3D() {
  DistTriangleVars preset_var;
  // totally inside
  Triangle p1({{0, 2, 0}, {2, 0, 0}, {0, 0, 0}});
  Triangle q1({{-0.5, -0.5, 2}, {2.5, 2.5, 2}, {0.5, 0.5, -2}});
  cout << "Triangle distance inside " << distTriangles_fcl(p1, q1, matI_, t0_)
       << endl;
  assert(approx_equal(distTriangles_fcl(p1, q1, matI_, t0_), 0));
  assert(approx_equal(DIST_TRIANGLES(&p1, &q1, matI_, t0_, &preset_var), 0));

  // first distance
  Triangle p2({{0, 2, 0}, {2, 0, 0}, {0, 0, 0}});
  Triangle q2({{-0.5, -0.5, 2}, {2.5, 2.5, 2}, {0.5, 0.5, 0}});
  cout << "Triangle distance away " << distTriangles_fcl(p2, q2, matI_, t0_)
       << endl;
  assert(approx_equal(distTriangles_fcl(p2, q2, matI_, t0_), 0));
  assert(approx_equal(DIST_TRIANGLES(&p2, &q2, matI_, t0_, &preset_var), 0));

  // first intersection
  Triangle p3({{0, 2, 0}, {2, 0, 0}, {0, 0, 0}});
  Triangle q3({{-0.5, -0.5, 2}, {2.5, 2.5, 2}, {0.5, 0.5, 0.25}});
  cout << "Triangle distance intersect "
       << distTriangles_fcl(p3, q3, matI_, t0_) << endl;
  assert(approx_equal(distTriangles_fcl(p3, q3, matI_, t0_), 0.25));
  assert(approx_equal(DIST_TRIANGLES(&p3, &q3, matI_, t0_, &preset_var), 0.25));

  Triangle p4({{0, 2, 0}, {2, 0, 0}, {0, 0, 0}});
  Triangle q4({{4.5, 4.5, 2}, {2.5, 2.5, 2}, {2.5, 2.5, -1}});
  cout << "Triangle distance intersect "
       << distTriangles_fcl(p4, q4, matI_, t0_) << endl;
  assert(approx_equal(distTriangles_fcl(p4, q4, matI_, t0_), 1.5 * sqrtf(2)));
  assert(approx_equal(DIST_TRIANGLES(&p4, &q4, matI_, t0_, &preset_var),
                      1.5 * sqrtf(2)));

  cout << "Test Triangeles 3D : PASSED" << endl;
}

HOST_PREFIX void generateRandomTriangle(Triangle* tri) {
  float v[9];
  for (int i = 0; i < 9; i++)
    v[i] = SPACE_LOW + static_cast<float>(rand()) /
                           static_cast<float>(RAND_MAX) *
                           (SPACE_HIGH - SPACE_LOW);

  tri->a.x = v[0];
  tri->a.y = v[1];
  tri->a.z = v[2];
  tri->b.x = v[3];
  tri->b.y = v[4];
  tri->b.z = v[5];
  tri->c.x = v[6];
  tri->c.y = v[7];
  tri->c.z = v[8];
}

HOST_PREFIX void test_stress_random() {
  srand(static_cast<unsigned>(time(NULL)));

  DistTriangleVars preset_var;
  int num_failed = 0;

  for (int i = 0; i < STRESS_CHECK; i++) {
    Triangle s1, s2;
    generateRandomTriangle(&s1);
    generateRandomTriangle(&s2);

    // Eigen::Vector3f angles = Eigen::Vector3f::Random();
    // Eigen::Vector3f pos = Eigen::Vector3f::Random();

    Eigen::Vector3f angles;
    Eigen::Vector3f pos;

    Eigen::AngleAxisf rollAngle(angles(0), Eigen::Vector3f::UnitZ());
    Eigen::AngleAxisf yawAngle(angles(1), Eigen::Vector3f::UnitY());
    Eigen::AngleAxisf pitchAngle(angles(2), Eigen::Vector3f::UnitX());
    Eigen::Quaternion<float> q = rollAngle * yawAngle * pitchAngle;
    Eigen::Matrix3f R = q.matrix();

    float matR[3][3], translation[3];
    matR[0][0] = R(0, 0);
    matR[0][1] = R(0, 1);
    matR[0][2] = R(0, 2);
    matR[1][0] = R(1, 0);
    matR[1][1] = R(1, 1);
    matR[1][2] = R(1, 2);
    matR[2][0] = R(2, 0);
    matR[2][1] = R(2, 1);
    matR[2][2] = R(2, 2);

    translation[0] = pos(0);
    translation[1] = pos(1);
    translation[2] = pos(2);

    float correct = distTriangles_fcl(s1, s2, matR, translation);
    float actual = DIST_TRIANGLES(&s1, &s2, matR, translation, &preset_var);

    if (!approx_equal(actual, correct)) {
      cout << "Difference is " << fabs(actual - correct) << endl;
      num_failed++;
    }
  }

  if (num_failed) {
    cout << num_failed << "/" << STRESS_CHECK << " failed! " << endl;
    assert(false);
  }

  cout << "TEST Triangles Stress random " << endl;
}

HOST_PREFIX std::ostream& operator<<(std::ostream& os, Triangle tri) {
  os << std::setprecision(4);
  os << "[" << tri.a.x << ", " << tri.a.y << ", " << tri.a.z << "], ";
  os << "[" << tri.b.x << ", " << tri.b.y << ", " << tri.b.z << "], ";
  os << "[" << tri.c.x << ", " << tri.c.y << ", " << tri.c.z << "]" << endl;
  return os;
}

#endif  // SRC_TEST_TRIANGLES_TEST_H_
