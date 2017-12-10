#define FCL_EXPORT
#include "compile_CPP.h"

#define DIST_TRIANGLES distTriangles
#include "Triangles_test.h"

int main(int argc, char *argv[])
{
  // call all test functions
  test_triangles_2D();
  test_triangles_3D();
  test_stress_random();

  print_Stats();
}
