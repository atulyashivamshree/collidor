CUDA_PREFIX void add(const Vector3* a, const Vector3* b, Vector3* c)
{

}

CUDA_PREFIX void sub(const Vector3* a, const Vector3* b, Vector3* c)
{

}

CUDA_PREFIX float dot(const Vector3* a, const Vector3* b)
{
  float res = 0;
  res += a->x * b->x;
  res += a->y * b->y;
  res += a->z * b->z;
  return res;
}

CUDA_PREFIX void cross(const Vector3* a, const Vector3* b, Vector3* c)
{

}

CUDA_PREFIX void mult(const Matrix3* a, const Vector3* b, Vector3* c)
{

}

CUDA_PREFIX void transform(const Matrix3* R, const Vector3* b, 
                                const Vector3* x, Vector3* y)
{

}