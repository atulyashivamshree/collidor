#define FCL_EXPORT
#include "compile_CPP.h"

#include "fcl/narrowphase/detail/primitive_shape_algorithm/triangle_distance.h"
#include "RSS.h"
#include <iostream>
#include <cstdlib>
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

const int NUM_CHECK = 10000;
const float EPSILON = 5e-7;

void test_triangles_2D();
void test_triangles_3D();
void test_stress_random();
void print_Stats();

int main(int argc, char *argv[])
{
  // call all test functions
  test_triangles_2D();
  test_triangles_3D();
  test_stress_random();

  print_Stats();
}

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
  cout << "Size of RSS is: " << sizeof(RSS) << endl;
  cout << "Size Of LineSegVars is: " << sizeof(LineSegVars) << endl;
  cout << "Size Of TriDistVars is: " << sizeof(TriDistVars) << endl;
  cout << "Size Of DistTriangleVars is: " << sizeof(DistTriangleVars) << endl;
}

void test_triangles_2D()
{
  DistTriangleVars preset_var;

  // totally inside
  Triangle p1({{0,6,0}, {6,0,0}, {0,0,0}});
  Triangle q1({{1,1,0}, {2,1,0}, {2,2,0}});
  cout << "Triangle distance inside " <<  distTriangles_fcl(p1, q1) << endl; 
  assert(approx_equal(distTriangles_fcl(p1, q1), 0));
  assert(approx_equal(distTriangles(&p1, &q1, &preset_var), 0));

  // first distance
  Triangle p2({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q2({{2,1,0}, {2,4,0}, {2,6,0}});
  cout << "Triangle distance away " <<  distTriangles_fcl(p2, q2) << endl; 
  assert(approx_equal(distTriangles_fcl(p2, q2), sqrt(2)/2));
  assert(approx_equal(distTriangles(&p2, &q2, &preset_var), sqrt(2)/2));

  // first intersection
  Triangle p3({{7,3,0}, {6,0,0}, {0,0,0}});
  Triangle q3({{1,1,0}, {4,4,0}, {6,2,0}});
  cout << "Triangle distance intersect " <<  distTriangles_fcl(p3, q3) << endl; 
  assert(approx_equal(distTriangles_fcl(p3, q3), 0));
  assert(approx_equal(distTriangles(&p3, &q3, &preset_var), 0));

  cout << "Test Triangeles 2D : PASSED" << endl;

}

void test_triangles_3D()
{
  DistTriangleVars preset_var;
  // totally inside
  Triangle p1({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q1({{-0.5,-0.5,2}, {2.5,2.5,2}, {0.5,0.5,-2}});
  cout << "Triangle distance inside " <<  distTriangles_fcl(p1, q1) << endl; 
  assert(approx_equal(distTriangles_fcl(p1, q1), 0));
  assert(approx_equal(distTriangles(&p1, &q1, &preset_var), 0));

  // first distance
  Triangle p2({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q2({{-0.5,-0.5,2}, {2.5,2.5,2}, {0.5,0.5,0}});
  cout << "Triangle distance away " <<  distTriangles_fcl(p2, q2) << endl; 
  assert(approx_equal(distTriangles_fcl(p2, q2), 0));
  assert(approx_equal(distTriangles(&p2, &q2, &preset_var), 0));

  // first intersection
  Triangle p3({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q3({{-0.5,-0.5,2}, {2.5,2.5,2}, {0.5,0.5,0.25}});
  cout << "Triangle distance intersect " <<  distTriangles_fcl(p3, q3) << endl; 
  assert(approx_equal(distTriangles_fcl(p3, q3), 0.25));
  assert(approx_equal(distTriangles(&p3, &q3, &preset_var), 0.25));

  Triangle p4({{0,2,0}, {2,0,0}, {0,0,0}});
  Triangle q4({{4.5,4.5,2}, {2.5,2.5,2}, {2.5,2.5,-1}});
  cout << "Triangle distance intersect " <<  distTriangles_fcl(p4, q4) << endl; 
  assert(approx_equal(distTriangles_fcl(p4, q4), 1.5*sqrt(2)));
  assert(approx_equal(distTriangles(&p4, &q4, &preset_var), 1.5*sqrt(2)));

  cout << "Test Triangeles 3D : PASSED" << endl;

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
    float actual = distTriangles(&s1, &s2, &preset_var);

    if(!approx_equal(actual, correct))
    {
      cout << "Difference is " << fabs(actual - correct) << endl;
      num_failed++;
    }
  }

  if(num_failed)
  {
    cout << num_failed << "/" << NUM_CHECK << " failed! " << endl;
    assert(false);
  }

  cout << "TEST Triangles Stress random : PASSED" << endl;
}
