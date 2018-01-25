/*
 *  test_rss.cpp
    Description: Checks the working of RSS on the CPU

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */
#define FCL_EXPORT
#include <iostream>
#include "../utils/compile_CPP.h"

#define DIST_TRIANGLES distTriangles
#include "../Triangles_test.h"

#define DIST_RSS rssDistance
#include "../Rectangle_tests.h"

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  // call all test functions
  test_triangles_2D();
  test_triangles_3D();
  test_stress_random();

  test_rectangles_2D();
  test_rectangles_3D();
  test_RSS_2D();
  test_RSS_3D();
  test_stress_random_RSS();

  cout << "=================" << endl;
  cout << "ALL TESTS PASSED!" << endl;
  cout << endl;
  cout << "=================" << endl;
  cout << "struct stats " << endl;
  print_Stats();
}
