#include "vector_math.h"
#include "Triangle.h"

#ifndef COLLIDOR_SRC_RSS
#define COLLIDOR_SRC_RSS

// holds the RSS object
struct RSS
{
	// Orientation of the RSS
	Matrix3 axis;

	// Orgiin of the rectangle in RSS
	Vector3 To;

	// side length of the rectangle
	float l[2];

	// radius of the sphere
	float r;

  // volume of the unit
  float vol;
};

// stores temporary vals implements the RSS distance computation
struct DistRSSVars
{
	float vert_a[4][3];
  float vert_b[4][3];

  float d[4];

  Triangle a1, a2, b1, b2;

	DistTriangleVars dist_triangle_vars;
};

CUDA_PREFIX float rssDistance(const Matrix3* R, const Vector3* t,
						const RSS* a, const RSS* b, DistRSSVars* d);

#include "RSS-inl.h"

#endif
