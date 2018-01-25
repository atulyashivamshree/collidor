/*
 *  BVH.h
    Description: Describes the basic structure of a BVH

    @author Atulya Shivam Shree
    Created on: Dec 11, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#ifndef SRC_RSS_INL_H_
#define SRC_RSS_INL_H_

CUDA_PREFIX void transformPoints(const float R[3][3], const float t[3],
                                 const RSS* a, float rot_v[4][3]) {
  float n1[3] = {a->axis.v1.x, a->axis.v1.y, a->axis.v1.z};
  float n2[3] = {a->axis.v2.x, a->axis.v2.y, a->axis.v2.z};
  float l[2] = {a->l[0], a->l[1]};
  float v0[3] = {a->To.x, a->To.y, a->To.z};

  // rot v0  = R*v0 + t;
  MxV(rot_v[0], R, v0);
  VpV(rot_v[0], rot_v[0], t);

  // v1 = v0 + l1*n1;
  // rot_v[1] = R*v1 + t = R*v0 + l1*R*n1 + t = rot_v[0] + l1*R*n1
  MxV(rot_v[1], R, n1);
  VxS(rot_v[1], rot_v[1], l[0]);
  VpV(rot_v[1], rot_v[1], rot_v[0]);

  // same for v2
  MxV(rot_v[2], R, n2);
  VxS(rot_v[2], rot_v[2], l[1]);
  VpV(rot_v[2], rot_v[2], rot_v[0]);

  // v3 = (v2-v0) + (v1-v0) + v0
  VcV(rot_v[3], rot_v[1]);
  VpV(rot_v[3], rot_v[3], rot_v[2]);
  VmV(rot_v[3], rot_v[3], rot_v[0]);
}

CUDA_PREFIX void getPoints(const RSS* a, float v[4][3]) {
  float n1[3] = {a->axis.v1.x, a->axis.v1.y, a->axis.v1.z};
  float n2[3] = {a->axis.v2.x, a->axis.v2.y, a->axis.v2.z};
  float l[2] = {a->l[0], a->l[1]};
  float v0[3] = {a->To.x, a->To.y, a->To.z};

  VcV(v[0], v0);

  // v1 = v0 + l1*n1
  VxS(v[1], n1, l[0]);
  VpV(v[1], v[1], v[0]);

  VxS(v[2], n2, l[1]);
  VpV(v[2], v[2], v[0]);

  // v3 = (v2-v0) + (v1-v0) + v0
  VcV(v[3], v[1]);
  VpV(v[3], v[3], v[2]);
  VmV(v[3], v[3], v[0]);
}

CUDA_PREFIX INLINE_PREFIX void initTriangle(Triangle* tri, const float v0[3],
                                            const float v1[3],
                                            const float v2[3]) {
  tri->a.x = v0[0];
  tri->a.y = v0[1];
  tri->a.z = v0[2];
  tri->b.x = v1[0];
  tri->b.y = v1[1];
  tri->b.z = v1[2];
  tri->c.x = v2[0];
  tri->c.y = v2[1];
  tri->c.z = v2[2];
}

CUDA_PREFIX float rssDistance(const float R[3][3], const float t[3],
                              const RSS* a, const RSS* b, DistRSSVars* p_var) {
  float vert_a[4][3];
  float vert_b[4][3];

  transformPoints(R, t, b, vert_b);
  getPoints(a, vert_a);

  float d[4];
  Triangle a1, a2, b1, b2;
  initTriangle(&a1, vert_a[0], vert_a[1], vert_a[3]);
  initTriangle(&a2, vert_a[0], vert_a[2], vert_a[3]);
  initTriangle(&b1, vert_b[0], vert_b[1], vert_b[3]);
  initTriangle(&b2, vert_b[0], vert_b[2], vert_b[3]);

  d[0] = distTriangles(&a1, &b1, &p_var->dist_triangle_vars);
  d[1] = distTriangles(&a1, &b2, &p_var->dist_triangle_vars);
  d[2] = distTriangles(&a2, &b1, &p_var->dist_triangle_vars);
  d[3] = distTriangles(&a2, &b2, &p_var->dist_triangle_vars);

  float min_d = fminf(d[0], d[1]);
  min_d = fminf(min_d, d[2]);
  min_d = fminf(min_d, d[3]);
  min_d = min_d - (a->r + b->r);
  min_d = fmaxf(min_d, 0);

  return min_d;
}

#endif  // SRC_RSS_INL_H_
