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
	double l[2];

	// radius of the sphere
	double r;
};

// stores temporary vals implements the RSS distance computation
struct RSS_dist
{

};

// stores temporary vals for triangle distance computation
struct Triangle_dist
{

};

CUDA_PREFIX float rectDistance(const RSS* a, const RSS* b, RSS_dist* d);

CUDA_PREFIX float rssDistance(const RSS* a, const RSS* b, RSS_dist* d);

CUDA_PREFIX float triangleDistance(const Triangle* a, const Triangle* b,
                                                      Triangle_dist* d);

#include "RSS-inl.h"

#endif