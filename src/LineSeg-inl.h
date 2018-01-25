/*
 *  BVH.h
    Description: Describes the basic structure of a BVH

    @author Atulya Shivam Shree
    Created on: Dec 11, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#ifndef SRC_LINESEG_INL_H_
#define SRC_LINESEG_INL_H_

CUDA_PREFIX INLINE_PREFIX void initV(float V[3]) {
  V[0] = 0;
  V[1] = 0;
  V[2] = 0;
}

CUDA_PREFIX INLINE_PREFIX void VcV(float Vr[3], const float V[3]) {
  Vr[0] = V[0];
  Vr[1] = V[1];
  Vr[2] = V[2];
}

CUDA_PREFIX INLINE_PREFIX void VmV(float Vr[3], const float V1[3],
                                   const float V2[3]) {
  Vr[0] = V1[0] - V2[0];
  Vr[1] = V1[1] - V2[1];
  Vr[2] = V1[2] - V2[2];
}

CUDA_PREFIX INLINE_PREFIX void VcrossV(float Vr[3], const float V1[3],
                                       const float V2[3]) {
  Vr[0] = V1[1] * V2[2] - V1[2] * V2[1];
  Vr[1] = V1[2] * V2[0] - V1[0] * V2[2];
  Vr[2] = V1[0] * V2[1] - V1[1] * V2[0];
}

CUDA_PREFIX INLINE_PREFIX float Vlength(float V[3]) {
  return sqrtf(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
}

CUDA_PREFIX INLINE_PREFIX void Vnormalize(float V[3]) {
  float d = (float)1.0 / sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
  V[0] *= d;
  V[1] *= d;
  V[2] *= d;
}

CUDA_PREFIX INLINE_PREFIX float VdotV(const float V1[3], const float V2[3]) {
  return (V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2]);
}

CUDA_PREFIX INLINE_PREFIX float VdistV2(const float V1[3], const float V2[3]) {
  return ((V1[0] - V2[0]) * (V1[0] - V2[0]) +
          (V1[1] - V2[1]) * (V1[1] - V2[1]) +
          (V1[2] - V2[2]) * (V1[2] - V2[2]));
}

CUDA_PREFIX INLINE_PREFIX void VpV(float Vr[3], const float V1[3],
                                   const float V2[3]) {
  Vr[0] = V1[0] + V2[0];
  Vr[1] = V1[1] + V2[1];
  Vr[2] = V1[2] + V2[2];
}

CUDA_PREFIX INLINE_PREFIX void VpVxS(float Vr[3], const float V1[3],
                                     const float V2[3], float s) {
  Vr[0] = V1[0] + V2[0] * s;
  Vr[1] = V1[1] + V2[1] * s;
  Vr[2] = V1[2] + V2[2] * s;
}

CUDA_PREFIX INLINE_PREFIX void VxS(float Vr[3], const float V[3], float s) {
  Vr[0] = V[0] * s;
  Vr[1] = V[1] * s;
  Vr[2] = V[2] * s;
}

CUDA_PREFIX INLINE_PREFIX void MxV(float Vr[3], const float M[3][3],
                                   const float V[3]) {
  Vr[0] = M[0][0] * V[0] + M[0][1] * V[1] + M[0][2] * V[2];
  Vr[1] = M[1][0] * V[0] + M[1][1] * V[1] + M[1][2] * V[2];
  Vr[2] = M[2][0] * V[0] + M[2][1] * V[1] + M[2][2] * V[2];
}

//--------------------------------------------------------------------------
// SegPoints()
//
// Returns closest points between an segment pair.
// Implemented from an algorithm described in
//
// Vladimir J. Lumelsky,
// On fast computation of distance between line segments.
// In Information Processing Letters, no. 21, pages 55-61, 1985.
//--------------------------------------------------------------------------

CUDA_PREFIX void SegPoints(float VEC[3], float X[3],
                           float Y[3],  // closest points
                           const float P[3],
                           const float A[3],  // seg 1 origin, vector
                           const float Q[3],
                           const float B[3],  // seg 2 origin, vector
                           LineSegVars* p_var) {
  VmV(p_var->T, Q, P);
  p_var->A_dot_A = VdotV(A, A);
  p_var->B_dot_B = VdotV(B, B);
  p_var->A_dot_B = VdotV(A, B);
  p_var->A_dot_T = VdotV(A, p_var->T);
  p_var->B_dot_T = VdotV(B, p_var->T);

  // t parameterizes ray P,A
  // u parameterizes ray Q,B

  // compute t for the closest point on ray P,A to
  // ray Q,B

  p_var->denom =
      p_var->A_dot_A * p_var->B_dot_B - p_var->A_dot_B * p_var->A_dot_B;

  p_var->t =
      (p_var->A_dot_T * p_var->B_dot_B - p_var->B_dot_T * p_var->A_dot_B) /
      p_var->denom;

  // clamp result so t is on the segment P,A

  if ((p_var->t < 0) || isnan(p_var->t))
    p_var->t = 0;
  else if (p_var->t > 1)
    p_var->t = 1;

  // find u for poinp_var->t on ray Q,B closest to point at t

  p_var->u = (p_var->t * p_var->A_dot_B - p_var->B_dot_T) / p_var->B_dot_B;

  // if u is on segment Q,B, t and u correspond to
  // closest points, otherwise, clamp u, recompute and
  // clamp t

  if ((p_var->u <= 0) || isnan(p_var->u)) {
    VcV(Y, Q);

    p_var->t = p_var->A_dot_T / p_var->A_dot_A;

    if ((p_var->t <= 0) || isnan(p_var->t)) {
      VcV(X, P);
      VmV(VEC, Q, P);
    } else if (p_var->t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Q, X);
    } else {
      VpVxS(X, P, A, p_var->t);
      VcrossV(p_var->TMP, p_var->T, A);
      VcrossV(VEC, A, p_var->TMP);
    }
  } else if (p_var->u >= 1) {
    VpV(Y, Q, B);

    p_var->t = (p_var->A_dot_B + p_var->A_dot_T) / p_var->A_dot_A;

    if ((p_var->t <= 0) || isnan(p_var->t)) {
      VcV(X, P);
      VmV(VEC, Y, P);
    } else if (p_var->t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Y, X);
    } else {
      VpVxS(X, P, A, p_var->t);
      VmV(p_var->T, Y, P);
      VcrossV(p_var->TMP, p_var->T, A);
      VcrossV(VEC, A, p_var->TMP);
    }
  } else {
    VpVxS(Y, Q, B, p_var->u);

    if ((p_var->t <= 0) || isnan(p_var->t)) {
      VcV(X, P);
      VcrossV(p_var->TMP, p_var->T, B);
      VcrossV(VEC, B, p_var->TMP);
    } else if (p_var->t >= 1) {
      VpV(X, P, A);
      VmV(p_var->T, Q, X);
      VcrossV(p_var->TMP, p_var->T, B);
      VcrossV(VEC, B, p_var->TMP);
    } else {
      VpVxS(X, P, A, p_var->t);
      VcrossV(VEC, A, B);
      if (VdotV(VEC, p_var->T) < 0) {
        VxS(VEC, VEC, -1);
      }
    }
  }
}

#endif  // SRC_LINESEG_INL_H_
