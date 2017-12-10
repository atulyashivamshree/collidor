/*
 * Triangles_test.h
 *
 *  Created on: Dec 10, 2017
 *      Author: atulya
 */
#include <iostream>
#include "RSS.h"
#include "Eigen/Dense"
#include "fcl/narrowphase/detail/primitive_shape_algorithm/triangle_distance.h"

#ifndef TRIANGLES_TEST_H_
#define TRIANGLES_TEST_H_

using namespace Eigen;

const int NUM_CHECK = 10000;
const float EPSILON = 5e-7;

Vector3f getVector3(Vector3 v)
{
  return Vector3f(v.x, v.y, v.z);
}

bool approx_equal(float a, float b, float epsilon = EPSILON)
{
  return (std::abs(a - b) < epsilon);
}

float distTriangles_fcl(Triangle s1, Triangle s2)
{
  Vector3f p, q;
  return fcl::detail::TriangleDistance<float>::triDistance(
                      getVector3(s1.a), getVector3(s1.b), getVector3(s1.c),
                      getVector3(s2.a), getVector3(s2.b), getVector3(s2.c),
                      p, q);
}

void print_Stats()
{
  std::cout << "Size of RSS is: " << sizeof(RSS) << std::endl;
  std::cout << "Size Of LineSegVars is: " << sizeof(LineSegVars) << std::endl;
  std::cout << "Size Of TriDistVars is: " << sizeof(TriDistVars) << std::endl;
  std::cout << "Size Of DistTriangleVars is: " << sizeof(DistTriangleVars) << std::endl;
}

void test_triangles_2D()
{
  DistTriangleVars preset_var;

  // totally inside
  Triangle p1({{0,6,0}, {6,0,0}, {0,0,0}});
  Triangle q1({{1,1,0}, {2,1,0}, {2,2,0}});
  std::cout << "Triangle distance inside " <<  distTriangles_fcl(p1, q1) << std::endl;
  assert(approx_equal(distTriangles_fcl(p1, q1), 0));
  assert(approx_equal(DIST_TRIANGLES(&p1, &q1, &preset_var), 0));

  // first distance
  Triangle p2({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q2({{2,1,0}, {2,4,0}, {2,6,0}});
  std::cout << "Triangle distance away " <<  distTriangles_fcl(p2, q2) << std::endl;
  assert(approx_equal(distTriangles_fcl(p2, q2), sqrt(2)/2));
  assert(approx_equal(DIST_TRIANGLES(&p2, &q2, &preset_var), sqrt(2)/2));

  // first intersection
  Triangle p3({{7,3,0}, {6,0,0}, {0,0,0}});
  Triangle q3({{1,1,0}, {4,4,0}, {6,2,0}});
  std::cout << "Triangle distance intersect " <<  distTriangles_fcl(p3, q3) << std::endl;
  assert(approx_equal(distTriangles_fcl(p3, q3), 0));
  assert(approx_equal(DIST_TRIANGLES(&p3, &q3, &preset_var), 0));

  std::cout << "Test Triangeles 2D : PASSED" << std::endl;

}

void test_triangles_3D()
{
  DistTriangleVars preset_var;
  // totally inside
  Triangle p1({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q1({{-0.5,-0.5,2}, {2.5,2.5,2}, {0.5,0.5,-2}});
  std::cout << "Triangle distance inside " <<  distTriangles_fcl(p1, q1) << std::endl;
  assert(approx_equal(distTriangles_fcl(p1, q1), 0));
  assert(approx_equal(DIST_TRIANGLES(&p1, &q1, &preset_var), 0));

  // first distance
  Triangle p2({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q2({{-0.5,-0.5,2}, {2.5,2.5,2}, {0.5,0.5,0}});
  std::cout << "Triangle distance away " <<  distTriangles_fcl(p2, q2) << std::endl;
  assert(approx_equal(distTriangles_fcl(p2, q2), 0));
  assert(approx_equal(DIST_TRIANGLES(&p2, &q2, &preset_var), 0));

  // first intersection
  Triangle p3({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q3({{-0.5,-0.5,2}, {2.5,2.5,2}, {0.5,0.5,0.25}});
  std::cout << "Triangle distance intersect " <<  distTriangles_fcl(p3, q3) << std::endl;
  assert(approx_equal(distTriangles_fcl(p3, q3), 0.25));
  assert(approx_equal(DIST_TRIANGLES(&p3, &q3, &preset_var), 0.25));

  Triangle p4({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q4({{4.5,4.5,2}, {2.5,2.5,2}, {2.5,2.5,-1}});
  std::cout << "Triangle distance intersect " <<  distTriangles_fcl(p4, q4) << std::endl;
  assert(approx_equal(distTriangles_fcl(p4, q4), 1.5*sqrt(2)));
  assert(approx_equal(DIST_TRIANGLES(&p4, &q4, &preset_var), 1.5*sqrt(2)));

  std::cout << "Test Triangeles 3D : PASSED" << std::endl;

}

void test_stress_random()
{
  DistTriangleVars preset_var;
  int num_failed = 0;

  for(int i = 0; i < NUM_CHECK; i++)
  {
    Vector3f v1 = Vector3f::Random(3,1);
    Vector3f v2 = Vector3f::Random(3,1);
    Vector3f v3 = Vector3f::Random(3,1);
    Vector3f v4 = Vector3f::Random(3,1);
    Vector3f v5 = Vector3f::Random(3,1);
    Vector3f v6 = Vector3f::Random(3,1);

    Triangle s1({{v1(0), v1(1), v1(2)}, {v2(0), v2(1), v2(2)}, {v3(0), v3(1), v3(2)}});
    Triangle s2({{v4(0), v4(1), v4(2)}, {v5(0), v5(1), v5(2)}, {v6(0), v6(1), v6(2)}});

    float correct = distTriangles_fcl(s1, s2);
    float actual = DIST_TRIANGLES(&s1, &s2, &preset_var);

    if(!approx_equal(actual, correct))
    {
      std::cout << "Difference is " << fabs(actual - correct) << std::endl;
      num_failed++;
    }
  }

  if(num_failed)
  {
    std::cout << num_failed << "/" << NUM_CHECK << " failed! " << std::endl;
    assert(false);
  }

  std::cout << "TEST Triangles Stress random : PASSED" << std::endl;
}




#endif /* TRIANGLES_TEST_H_ */
