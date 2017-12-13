#include "vector_math.h"

#ifndef COLLIDOR_SRC_TRIANGLE
#define COLLIDOR_SRC_TRIANGLE

// holds a triangle in 3D space
struct Triangle
{
  Vector3 a;
  Vector3 b;
  Vector3 c;
};

// stores temporary vals for triangle distance computation

struct LineSegVars
{
  float T[3], A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;

  float TMP[3];

  float t,u;

  float denom;
};

struct TriangleResult
{
  float dist;
};

struct TriDistVars
{
  float Sv[3][3];
  float Tv[3][3];

  float VEC[3];

  float V[3];
  float Z[3];
  float minP[3], minQ[3], mindd;
  int shown_disjoint;

  float dd;

  float a;
  float b;

  float Sn[3], Snl;

  float Tp[3];

  int point;

  float Tn[3], Tnl;

  float Sp[3];

  float p;

  LineSegVars line_seg_vars;
};

struct DistTriangleVars
{
  float S[3][3];
  float T[3][3];

  float X[3];
  float Y[3];

  TriDistVars tri_dist_vars;
};

CUDA_PREFIX INLINE_PREFIX float distTriangles(const Triangle* s1, const Triangle *s2,
                        DistTriangleVars* p_var);

CUDA_PREFIX INLINE_PREFIX float distTriangles(const Triangle* s1, const Triangle *s2,
						const float R[3][3], const float t[3],
                        DistTriangleVars* p_var);

#include "Triangle-inl.h"

#endif
