/*
 * Rectangle_tests.h
 *
 *  Created on: Dec 10, 2017
 *      Author: atulya
 */

#include "Triangles_test.h"
#include "fcl/math/bv/RSS.h"

#ifndef RECTANGLE_TESTS_H_
#define RECTANGLE_TESTS_H_

using Eigen::Vector3f;

// TODO(@atulya) rename Matrix3 to Mat3 to avoid confusion with EIgen
float t2_[3] = {0, 0, 1.5};
const Matrix3 MAT_I = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
const Vector3 t0({0, 0, 0});
Matrix3 MAT_RXY45;
Matrix3 MAT_RXZ45;

HOST_PREFIX void initializeGlobalMAT()
{
  Vector3 v1({sqrtf(2)/2, -sqrtf(2)/2, 0});
  Vector3 v2({sqrtf(2)/2, sqrtf(2)/2, 0});

  MAT_RXY45.v1 = v1;
  MAT_RXY45.v2 = v2;
  MAT_RXY45.v3 = Vector3({0, 0, 1});

  Vector3 v3({sqrtf(2)/2, 0, sqrtf(2)/2});
  Vector3 v4({-sqrtf(2)/2, 0, sqrtf(2)/2});

  MAT_RXZ45.v1 = v3;
  MAT_RXZ45.v2 = v4;
  MAT_RXZ45.v3 = Vector3({0, 1, 0});
}

RSS getRandomRSS()
{
  float v[9];
  for(int i = 0; i < 9; i++)
      v[i] = SPACE_LOW + static_cast <float> (rand()) /
                      static_cast <float> (RAND_MAX)*(SPACE_HIGH - SPACE_LOW);
  Eigen::AngleAxisf rollAngle(v[0], Eigen::Vector3f::UnitZ());
  Eigen::AngleAxisf yawAngle(v[1], Eigen::Vector3f::UnitY());
  Eigen::AngleAxisf pitchAngle(v[2], Eigen::Vector3f::UnitX());

  Eigen::Quaternion<float> q = rollAngle * yawAngle * pitchAngle;

  Eigen::Matrix3f R = q.matrix();

  RSS obj;
  obj.axis.v1 = Vector3({R(0,0), R(1,0), R(2,0)});
  obj.axis.v2 = Vector3({R(0,1), R(1,1), R(2,1)});
  obj.axis.v3 = Vector3({R(0,2), R(1,2), R(2,2)});

  obj.To = Vector3({v[3]*100, v[4]*100, v[5]*100});
  obj.l[0] = (v[6] - SPACE_LOW) / (SPACE_HIGH - SPACE_LOW) * 100;
  obj.l[1] = (v[7] - SPACE_LOW) / (SPACE_HIGH - SPACE_LOW) * 100;
  obj.r = (v[8] - SPACE_LOW) / (SPACE_HIGH - SPACE_LOW) * 10;
  
  return obj;
}

fcl::RSS<float> getFCLRSS(RSS shp) {
  fcl::RSS<float> obj;
  obj.axis(0,0) = shp.axis.v1.x;
  obj.axis(1,0) = shp.axis.v1.y;
  obj.axis(2,0) = shp.axis.v1.z;
  obj.axis(0,1) = shp.axis.v2.x;
  obj.axis(1,1) = shp.axis.v2.y;
  obj.axis(2,1) = shp.axis.v2.z;
  obj.axis(0,2) = shp.axis.v3.x;
  obj.axis(1,2) = shp.axis.v3.y;
  obj.axis(2,2) = shp.axis.v3.z;

  obj.To(0) = shp.To.x;
  obj.To(1) = shp.To.y;
  obj.To(2) = shp.To.z;

  obj.l[0] = shp.l[0];
  obj.l[1] = shp.l[1];

  obj.r = shp.r;

  //compute volume in this or its inverse

  return obj;
}

HOST_PREFIX float distRectangles_fcl(const float R[3][3], const float t[3], RSS s1, RSS s2)
{
  fcl::RSS<float> fcl_s1, fcl_s2;
  fcl_s1 = getFCLRSS(s1);
  fcl_s2 = getFCLRSS(s2);
  Vector3f p, q;
  fcl::Transform3f tf;
  tf(0,0) = R[0][0]; tf(0,1) = R[0][1]; tf(0,2) = R[0][2]; tf(0,3) = t[0];
  tf(1,0) = R[1][0]; tf(1,1) = R[1][1]; tf(1,2) = R[1][2]; tf(1,3) = t[1];
  tf(2,0) = R[2][0]; tf(2,1) = R[2][1]; tf(2,2) = R[2][2]; tf(2,3) = t[2];
  
  return distance(tf.linear(), tf.translation(),fcl_s1 , fcl_s2, &p, &q);
}

void print_Rectangle_stats()
{
  std::cout << "Size of RSS is: " << sizeof(RSS) << std::endl;
  std::cout << "Size Of DistRSSVars is: " << sizeof(DistRSSVars) << std::endl;
}

void test_rectangles_2D()
{
  initializeGlobalMAT();

  DistRSSVars preset_var;

  // totally outside
  RSS p1, q1, q2, q3;
  p1.axis = MAT_I;
  p1.To = Vector3({0, 0, 0});
  p1.l[0] = 5;
  p1.l[1] = 3;
  p1.r = 0.0;
  q1.axis = MAT_RXY45;
  q1.To = Vector3({6.5, 1, 0});
  q1.l[0] = 5;
  q1.l[1] = 3;
  q1.r = 0.0;

  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q1), 1.5));
  std::cout << "RSS distance outside " <<  DIST_RSS(matI_, t0_, &p1, &q1, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q1, &preset_var), 1.5));

  q2.axis = MAT_RXY45;
  q2.To = Vector3({5, 1, 0});
  q2.l[0] = 5;
  q2.l[1] = 3;
  q2.r = 0.0;  
  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q2), 0));
  std::cout << "Rectangle distance on " <<  DIST_RSS(matI_, t0_, &p1, &q2, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q2, &preset_var), 0));

  // first intersection
  q3.axis = MAT_RXY45;
  q3.To = Vector3({2.5, 1, 0});
  q3.l[0] = 5;
  q3.l[1] = 3;
  q3.r = 0.0;  
  std::cout << "Rectangle distance inside " <<  DIST_RSS(matI_, t0_, &p1, &q3, &preset_var) << std::endl;
  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q3), 0));
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q3, &preset_var), 0));


  std::cout << "Rectangle distance inside " <<  DIST_RSS(matI_, t1_, &p1, &q1, &preset_var) << std::endl;
  std::cout << "Rectangle distance inside true" <<  distRectangles_fcl(matI_, t1_, p1, q1) << std::endl;
  assert(approx_equal(distRectangles_fcl(matI_, t1_, p1, q1), 2.5));
  assert(approx_equal(DIST_RSS(matI_, t1_, &p1, &q1, &preset_var), 2.5));

  std::cout << "Rectangle distance inside " <<  DIST_RSS(matI_, t2_, &p1, &q1, &preset_var) << std::endl;
    std::cout << "Rectangle distance inside true" <<  distRectangles_fcl(matI_, t2_, p1, q1) << std::endl;
    assert(approx_equal(distRectangles_fcl(matI_, t2_, p1, q1), 1.5*sqrt(2)));
    assert(approx_equal(DIST_RSS(matI_, t2_, &p1, &q1, &preset_var), 1.5*sqrt(2)));

  std::cout << "Test Rectangles 2D: PASSED" << std::endl;

}

void test_rectangles_3D()
{
  initializeGlobalMAT();

  DistRSSVars preset_var;

  // totally outside
  RSS p1, q1, q2, q3;
  p1.axis = MAT_I;
  p1.To = Vector3({0, 0, 0});
  p1.l[0] = 5;
  p1.l[1] = 3;
  p1.r = 0.0;
  q1.axis = MAT_RXZ45;
  q1.To = Vector3({2.5, 1, 1.23});
  q1.l[0] = 5;
  q1.l[1] = 3;
  q1.r = 0.0;

  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q1), 1.23));
  std::cout << "Rectangle distance outside " <<  DIST_RSS(matI_, t0_, &p1, &q1, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q1, &preset_var), 1.23));

  q2.axis = MAT_RXZ45;
  q2.To = Vector3({2.5, 1, 0});
  q2.l[0] = 5;
  q2.l[1] = 3;
  q2.r = 0.0;
  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q2), 0));
  std::cout << "Rectangle distance on " <<  DIST_RSS(matI_, t0_, &p1, &q2, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q2, &preset_var), 0));

  // first intersection
  q3.axis = MAT_RXZ45;
  q3.To = Vector3({2.5, 1, -0.5});
  q3.l[0] = 5;
  q3.l[1] = 3;
  q3.r = 0.0;

  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q3), 0));
  std::cout << "Rectangle distance intersection " <<  DIST_RSS(matI_, t0_, &p1, &q3, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q3, &preset_var), 0));

  std::cout << "Test Rectangles 3D : PASSED" << std::endl;
}

void test_RSS_2D()
{
  initializeGlobalMAT();

  DistRSSVars preset_var;

  // totally outside
  RSS p1, q1, q2, q3;
  p1.axis = MAT_I;
  p1.To = Vector3({0, 0, 0});
  p1.l[0] = 5;
  p1.l[1] = 3;
  p1.r = 0.5;
  q1.axis = MAT_RXY45;
  q1.To = Vector3({6.5, 1, 0});
  q1.l[0] = 5;
  q1.l[1] = 3;
  q1.r = 0.5;

  std::cout << "RSS distance outside " <<  distRectangles_fcl(matI_, t0_, p1, q1) << std::endl;
  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q1), 0.5));
  std::cout << "RSS distance outisede " <<  DIST_RSS(matI_, t0_, &p1, &q1, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q1, &preset_var), 0.5));

  q2.axis = MAT_RXY45;
  q2.To = Vector3({6, 1, 0});
  q2.l[0] = 5;
  q2.l[1] = 3;
  q2.r = 0.5;  
  std::cout << "RSS distance on " <<  distRectangles_fcl(matI_, t0_, p1, q2) << std::endl;
  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q2), 0));
  std::cout << "RSS distance on " <<  DIST_RSS(matI_, t0_, &p1, &q2, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q2, &preset_var), 0));

  // first intersection
  q3.axis = MAT_RXY45;
  q3.To = Vector3({2.5, 1, 0});
  q3.l[0] = 5;
  q3.l[1] = 3;
  q3.r = 0.5;  
  std::cout << "RSS distance intersection " <<  distRectangles_fcl(matI_, t0_, p1, q3) << std::endl;
  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q3), 0));
  std::cout << "RSS distance intersection " <<  DIST_RSS(matI_, t0_, &p1, &q3, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q3, &preset_var), 0));

  std::cout << "Test RSS 2D: PASSED" << std::endl;

}

void test_RSS_3D()
{
  initializeGlobalMAT();

  DistRSSVars preset_var;

  // totally outside
  RSS p1, q1, q2, q3;
  p1.axis = MAT_I;
  p1.To = Vector3({0, 0, 0});
  p1.l[0] = 5;
  p1.l[1] = 3;
  p1.r = 0.5;
  q1.axis = MAT_RXZ45;
  q1.To = Vector3({2.5, 1, 1.23});
  q1.l[0] = 5;
  q1.l[1] = 3;
  q1.r = 0.5;

  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q1), 0.23));
  std::cout << "RSS distance outside " <<  DIST_RSS(matI_, t0_, &p1, &q1, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q1, &preset_var), 0.23));

  q2.axis = MAT_RXZ45;
  q2.To = Vector3({2.5, 1, 1});
  q2.l[0] = 5;
  q2.l[1] = 3;
  q2.r = 0.5;
  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q2), 0));
  std::cout << "RSS distance on " <<  DIST_RSS(matI_, t0_, &p1, &q2, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q2, &preset_var), 0));

  // first intersection
  q3.axis = MAT_RXZ45;
  q3.To = Vector3({2.5, 1, -0.5});
  q3.l[0] = 5;
  q3.l[1] = 3;
  q3.r = 0.5;
  assert(approx_equal(distRectangles_fcl(matI_, t0_, p1, q3), 0));
  std::cout << "RSS distance intersection " <<  DIST_RSS(matI_, t0_, &p1, &q3, &preset_var) << std::endl;
  assert(approx_equal(DIST_RSS(matI_, t0_, &p1, &q3, &preset_var), 0));

  std::cout << "Test RSS 3D : PASSED" << std::endl;
}


void test_stress_random_RSS()
{
  srand(static_cast<unsigned> (time(NULL)));

  DistRSSVars preset_var;
  int num_failed = 0;

  for(int i = 0; i < STRESS_CHECK; i++)
  {

    RSS r1 = getRandomRSS();
    RSS r2 = getRandomRSS();

    Eigen::Vector3f angles = Eigen::Vector3f::Random();
    Eigen::Vector3f pos = Eigen::Vector3f::Random();

    Eigen::AngleAxisf rollAngle(angles(0), Eigen::Vector3f::UnitZ());
    Eigen::AngleAxisf yawAngle(angles(1), Eigen::Vector3f::UnitY());
    Eigen::AngleAxisf pitchAngle(angles(2), Eigen::Vector3f::UnitX());
    Eigen::Quaternion<float> q = rollAngle * yawAngle * pitchAngle;
    Eigen::Matrix3f R = q.matrix();

    float matR[3][3], translation[3];
    matR[0][0] = R(0,0); matR[0][1] = R(0,1); matR[0][2] = R(0,2);
    matR[1][0] = R(1,0); matR[1][1] = R(1,1); matR[1][2] = R(1,2);
    matR[2][0] = R(2,0); matR[2][1] = R(2,1); matR[2][2] = R(2,2);

    translation[0] = pos(0);
    translation[1] = pos(1);
    translation[2] = pos(2);

   float correct = distRectangles_fcl(matR, translation, r1, r2);
//   matR[0][0] = 1;
   float actual = DIST_RSS(matR, translation, &r1, &r2, &preset_var);

   if(!approx_equal(actual, correct))
   {
     std::cout << "Difference is " << fabs(actual - correct) << std::endl;
     num_failed++;
   }
   // cout << "actual " << actual << "; correct" << correct << endl;
  }

  std::cout << num_failed << "/" << STRESS_CHECK << " failed! " << std::endl;
  std::cout << "TEST Rectangles Stress random : PASSED" << std::endl;
}



#endif /* RECTANGLE_TESTS_H_ */
