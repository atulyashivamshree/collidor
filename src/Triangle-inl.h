/*
 *  Triangle-inl.h
    Description: Implements the functions for distance between triangles

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */
#include "LineSeg-inl.h"

#ifndef SRC_TRIANGLE_INL_H_
#define SRC_TRIANGLE_INL_H_

CUDA_PREFIX float TriDist(float P[3], float Q[3], const float S[3][3],
                          const float T[3][3], TriDistVars* p_var) {
  // Compute vectors along the 6 sides

  VmV(p_var->Sv[0], S[1], S[0]);
  VmV(p_var->Sv[1], S[2], S[1]);
  VmV(p_var->Sv[2], S[0], S[2]);

  VmV(p_var->Tv[0], T[1], T[0]);
  VmV(p_var->Tv[1], T[2], T[1]);
  VmV(p_var->Tv[2], T[0], T[2]);

  // For each edge pair, the vector connecting the closest p_var->points
  // of the edges defines a slab (parallel planes at head and tail
  // enclose the slab). If we can show that the off-edge vertex of
  // each triangle is outside of the slab, then the closest points
  // of the edges are the closest points for the triangles.
  // Even if these tests fail, it may be helpful to know the closest
  // points found, and whether the triangles were shown disjoint

  p_var->shown_disjoint = 0;

  p_var->mindd = VdistV2(S[0], T[0]) + 1;  // Set first minimum safely high

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      // Find closest points on edges i & j, plus the
      // vector (and distance squared) between these points

      SegPoints(p_var->VEC, P, Q, S[i], p_var->Sv[i], T[j], p_var->Tv[j],
                &p_var->line_seg_vars);

      VmV(p_var->V, Q, P);
      p_var->dd = VdotV(p_var->V, p_var->V);

      // Verify this closest point pair only if the distance
      // squared is less than the minimum found thus far.

      if (p_var->dd <= p_var->mindd) {
        VcV(p_var->minP, P);
        VcV(p_var->minQ, Q);
        p_var->mindd = p_var->dd;

        VmV(p_var->Z, S[(i + 2) % 3], P);
        p_var->a = VdotV(p_var->Z, p_var->VEC);
        VmV(p_var->Z, T[(j + 2) % 3], Q);
        p_var->b = VdotV(p_var->Z, p_var->VEC);

        if ((p_var->a <= 0) && (p_var->b >= 0)) return sqrt(p_var->dd);

        p_var->p = VdotV(p_var->V, p_var->VEC);

        if (p_var->a < 0) p_var->a = 0;
        if (p_var->b > 0) p_var->b = 0;
        if ((p_var->p - p_var->a + p_var->b) > 0) p_var->shown_disjoint = 1;
      }
    }
  }

  // No edge pairs contained the closest points.
  // either:
  // 1. one of the closest points is a vertex, and the
  //    other point is interior to a face.
  // 2. the triangles are overlapping.
  // 3. an edge of one triangle is parallel to the other's face. If
  //    cases 1 and 2 are not true, then the closest points from the 9
  //    edge pairs checks above can be taken as closest points for the
  //    triangles.
  // 4. possibly, the triangles were degenerate.  When the
  //    triangle points are nearly colinear or coincident, one
  //    of above tests might fail even though the edges tested
  //    contain the closest points.

  // First check for case 1

  VcrossV(p_var->Sn, p_var->Sv[0],
          p_var->Sv[1]);  // Compute normal to S triangle
  p_var->Snl =
      VdotV(p_var->Sn, p_var->Sn);  // Compute square of length of normal

  // If cross product is long enough,

  if (p_var->Snl > 1e-15) {
    // Get projection lengths of T points

    VmV(p_var->V, S[0], T[0]);
    p_var->Tp[0] = VdotV(p_var->V, p_var->Sn);

    VmV(p_var->V, S[0], T[1]);
    p_var->Tp[1] = VdotV(p_var->V, p_var->Sn);

    VmV(p_var->V, S[0], T[2]);
    p_var->Tp[2] = VdotV(p_var->V, p_var->Sn);

    // If p_var->Sn is a separating direction,
    // find point with smallest projection

    p_var->point = -1;
    if ((p_var->Tp[0] > 0) && (p_var->Tp[1] > 0) && (p_var->Tp[2] > 0)) {
      if (p_var->Tp[0] < p_var->Tp[1])
        p_var->point = 0;
      else
        p_var->point = 1;
      if (p_var->Tp[2] < p_var->Tp[p_var->point]) p_var->point = 2;
    } else if ((p_var->Tp[0] < 0) && (p_var->Tp[1] < 0) && (p_var->Tp[2] < 0)) {
      if (p_var->Tp[0] > p_var->Tp[1])
        p_var->point = 0;
      else
        p_var->point = 1;
      if (p_var->Tp[2] > p_var->Tp[p_var->point]) p_var->point = 2;
    }

    // If p_var->Sn is a separating direction,

    if (p_var->point >= 0) {
      p_var->shown_disjoint = 1;

      // Test whether the p_var->point found, when projected onto the
      // other triangle, lies within the face.

      VmV(p_var->V, T[p_var->point], S[0]);
      VcrossV(p_var->Z, p_var->Sn, p_var->Sv[0]);
      if (VdotV(p_var->V, p_var->Z) > 0) {
        VmV(p_var->V, T[p_var->point], S[1]);
        VcrossV(p_var->Z, p_var->Sn, p_var->Sv[1]);
        if (VdotV(p_var->V, p_var->Z) > 0) {
          VmV(p_var->V, T[p_var->point], S[2]);
          VcrossV(p_var->Z, p_var->Sn, p_var->Sv[2]);
          if (VdotV(p_var->V, p_var->Z) > 0) {
            // T[p_var->point] passed the test - it's a closest p_var->point for
            // the T triangle; the other p_var->point is on the face of S

            VpVxS(P, T[p_var->point], p_var->Sn,
                  p_var->Tp[p_var->point] / p_var->Snl);
            VcV(Q, T[p_var->point]);
            return sqrt(VdistV2(P, Q));
          }
        }
      }
    }
  }

  VcrossV(p_var->Tn, p_var->Tv[0], p_var->Tv[1]);
  p_var->Tnl = VdotV(p_var->Tn, p_var->Tn);

  if (p_var->Tnl > 1e-15) {
    VmV(p_var->V, T[0], S[0]);
    p_var->Sp[0] = VdotV(p_var->V, p_var->Tn);

    VmV(p_var->V, T[0], S[1]);
    p_var->Sp[1] = VdotV(p_var->V, p_var->Tn);

    VmV(p_var->V, T[0], S[2]);
    p_var->Sp[2] = VdotV(p_var->V, p_var->Tn);

    p_var->point = -1;
    if ((p_var->Sp[0] > 0) && (p_var->Sp[1] > 0) && (p_var->Sp[2] > 0)) {
      if (p_var->Sp[0] < p_var->Sp[1])
        p_var->point = 0;
      else
        p_var->point = 1;
      if (p_var->Sp[2] < p_var->Sp[p_var->point]) p_var->point = 2;
    } else if ((p_var->Sp[0] < 0) && (p_var->Sp[1] < 0) && (p_var->Sp[2] < 0)) {
      if (p_var->Sp[0] > p_var->Sp[1])
        p_var->point = 0;
      else
        p_var->point = 1;
      if (p_var->Sp[2] > p_var->Sp[p_var->point]) p_var->point = 2;
    }

    if (p_var->point >= 0) {
      p_var->shown_disjoint = 1;

      VmV(p_var->V, S[p_var->point], T[0]);
      VcrossV(p_var->Z, p_var->Tn, p_var->Tv[0]);
      if (VdotV(p_var->V, p_var->Z) > 0) {
        VmV(p_var->V, S[p_var->point], T[1]);
        VcrossV(p_var->Z, p_var->Tn, p_var->Tv[1]);
        if (VdotV(p_var->V, p_var->Z) > 0) {
          VmV(p_var->V, S[p_var->point], T[2]);
          VcrossV(p_var->Z, p_var->Tn, p_var->Tv[2]);
          if (VdotV(p_var->V, p_var->Z) > 0) {
            VcV(P, S[p_var->point]);
            VpVxS(Q, S[p_var->point], p_var->Tn,
                  p_var->Sp[p_var->point] / p_var->Tnl);
            return sqrt(VdistV2(P, Q));
          }
        }
      }
    }
  }

  // Case 1 can't be shown.
  // If one of these tests showed the triangles disjoint,
  // we assume case 3 or 4, otherwise we conclude case 2,
  // that the triangles overlap.

  if (p_var->shown_disjoint) {
    VcV(P, p_var->minP);
    VcV(Q, p_var->minQ);
    return sqrt(p_var->mindd);
  } else {
    return 0;
  }
}

CUDA_PREFIX float distTriangles(const Triangle* s1, const Triangle* s2,
                                const float R[3][3], const float t[3],
                                DistTriangleVars* p_var) {
  p_var->S[0][0] = s1->a.x;
  p_var->S[0][1] = s1->a.y;
  p_var->S[0][2] = s1->a.z;
  p_var->S[1][0] = s1->b.x;
  p_var->S[1][1] = s1->b.y;
  p_var->S[1][2] = s1->b.z;
  p_var->S[2][0] = s1->c.x;
  p_var->S[2][1] = s1->c.y;
  p_var->S[2][2] = s1->c.z;

  float s2_a[3] = {s2->a.x, s2->a.y, s2->a.z};
  float s2_b[3] = {s2->b.x, s2->b.y, s2->b.z};
  float s2_c[3] = {s2->c.x, s2->c.y, s2->c.z};

  MxV(p_var->T[0], R, s2_a);
  VpV(p_var->T[0], p_var->T[0], t);
  MxV(p_var->T[1], R, s2_b);
  VpV(p_var->T[1], p_var->T[1], t);
  MxV(p_var->T[2], R, s2_c);
  VpV(p_var->T[2], p_var->T[2], t);

  return TriDist(p_var->X, p_var->Y, p_var->S, p_var->T, &p_var->tri_dist_vars);
}
CUDA_PREFIX float distTriangles(const Triangle* s1, const Triangle* s2,
                                DistTriangleVars* p_var) {
  p_var->S[0][0] = s1->a.x;
  p_var->S[0][1] = s1->a.y;
  p_var->S[0][2] = s1->a.z;
  p_var->S[1][0] = s1->b.x;
  p_var->S[1][1] = s1->b.y;
  p_var->S[1][2] = s1->b.z;
  p_var->S[2][0] = s1->c.x;
  p_var->S[2][1] = s1->c.y;
  p_var->S[2][2] = s1->c.z;

  p_var->T[0][0] = s2->a.x;
  p_var->T[0][1] = s2->a.y;
  p_var->T[0][2] = s2->a.z;
  p_var->T[1][0] = s2->b.x;
  p_var->T[1][1] = s2->b.y;
  p_var->T[1][2] = s2->b.z;
  p_var->T[2][0] = s2->c.x;
  p_var->T[2][1] = s2->c.y;
  p_var->T[2][2] = s2->c.z;

  return TriDist(p_var->X, p_var->Y, p_var->S, p_var->T, &p_var->tri_dist_vars);
}

#endif  // SRC_TRIANGLE_INL_H_
