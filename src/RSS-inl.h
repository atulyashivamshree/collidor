CUDA_PREFIX void transformPoints(const Matrix3* R, const Vector3* t, 
                                      const RSS *a, float rot_v[4][3])
{
  float matR[3][3];
  matR[0][0] = R->v1.x;
  matR[1][0] = R->v1.y;
  matR[2][0] = R->v1.z;
  matR[0][1] = R->v2.x;
  matR[1][1] = R->v2.y;
  matR[2][1] = R->v2.z;
  matR[0][2] = R->v3.x;
  matR[1][2] = R->v3.y;
  matR[2][2] = R->v3.z;
  float vecT[3] = {t->x, t->y, t->z};

  float n1[3] = {a->axis.v1.x, a->axis.v1.y, a->axis.v1.z};
  float n2[3] = {a->axis.v2.x, a->axis.v2.y, a->axis.v2.z};
  float l[2] = {a->l[0], a->l[1]};
  float v0[3] = {a->To.x, a->To.y, a->To.z};

  // rot v0  = R*v0 + t;
  MxV(rot_v[0], matR, v0);
  VpV(rot_v[0], rot_v[0], vecT);

  // v1 = v0 + l1*n1;
  // rot_v[1] = R*v1 + t = R*v0 + l1*R*n1 + t = rot_v[0] + l1*R*n1
  MxV(rot_v[1], matR, n1);
  VxS(rot_v[1], rot_v[1], l[0]);
  VpV(rot_v[1], rot_v[1], rot_v[0]);

  // same for v2
  MxV(rot_v[2], matR, n2);
  VxS(rot_v[2], rot_v[2], l[1]);
  VpV(rot_v[2], rot_v[2], rot_v[0]);

  // v3 = (v2-v0) + (v1-v0) + v0
  VcV(rot_v[3], rot_v[1]);
  VpV(rot_v[3], rot_v[3], rot_v[2]);
  VmV(rot_v[3], rot_v[3], rot_v[0]);
}

CUDA_PREFIX void getPoints(const RSS *a, float v[4][3])
{
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

CUDA_PREFIX INLINE_PREFIX void initTriangle(Triangle* tri, 
            const float v0[3], const float v1[3], const float v2[3])
{
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

CUDA_PREFIX float rssDistance(const Matrix3* R, const Vector3* t,
						const RSS* a, const RSS* b, DistRSSVars* p_var)
{
  

  transformPoints(R, t, a, p_var->vert_a);
  getPoints(b, p_var->vert_b);

  initTriangle(&p_var->a1, p_var->vert_a[0], p_var->vert_a[1], p_var->vert_a[3]);
  initTriangle(&p_var->a2, p_var->vert_a[0], p_var->vert_a[2], p_var->vert_a[3]);
  initTriangle(&p_var->b1, p_var->vert_b[0], p_var->vert_b[1], p_var->vert_b[3]);
  initTriangle(&p_var->b2, p_var->vert_b[0], p_var->vert_b[2], p_var->vert_b[3]);

  p_var->d[0] = distTriangles(&p_var->a1, &p_var->b1, &p_var->dist_triangle_vars);
  p_var->d[1] = distTriangles(&p_var->a1, &p_var->b2, &p_var->dist_triangle_vars);
  p_var->d[2] = distTriangles(&p_var->a2, &p_var->b1, &p_var->dist_triangle_vars);
  p_var->d[3] = distTriangles(&p_var->a2, &p_var->b2, &p_var->dist_triangle_vars);

  float min_d = fminf(p_var->d[0], p_var->d[1]);
  min_d = fminf(min_d, p_var->d[2]);
  min_d = fminf(min_d, p_var->d[3]);
  min_d = min_d - (a->r + b->r);
  min_d = fmaxf(min_d, 0);

	return min_d;
}
