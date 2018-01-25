#define FCL_EXPORT
#include "../compile_CPP.h"
#include <iostream>

#define DIST_TRIANGLES distTriangles
#include "../Triangles_test.h"

#define DIST_RSS rssDistance
#include "../Rectangle_tests.h"

using namespace std;

int main(int argc, char *argv[])
{
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
