//------------------------------------------------------------------------------
// File : vec3fv.hpp
//------------------------------------------------------------------------------
// GLVU : Copyright 1997 - 2002 
//        The University of North Carolina at Chapel Hill
//------------------------------------------------------------------------------
// Permission to use, copy, modify, distribute and sell this software and its 
// documentation for any purpose is hereby granted without fee, provided that 
// the above copyright notice appear in all copies and that both that copyright 
// notice and this permission notice appear in supporting documentation. 
// Binaries may be compiled with this software without any royalties or 
// restrictions. 
//
// The University of North Carolina at Chapel Hill makes no representations 
// about the suitability of this software for any purpose. It is provided 
// "as is" without express or implied warranty.

//============================================================================
// vec3fv.hpp
//============================================================================

#ifndef _VEC3FV_H_
#define _VEC3FV_H_

#include <stdio.h>
#include <math.h>

inline void Set3fv(float v[3], float x, float y, float z)
{
  v[0]=x;
  v[1]=y;
  v[2]=z;
}

inline void Copy3fv(float A[3], const float B[3]) // A=B
{
  A[0]=B[0];
  A[1]=B[1];
  A[2]=B[2];
}

inline void ScalarMult3fv(float c[3], const float a[3], float s) // c=a*s
{
  c[0] = a[0] * s;
  c[1] = a[1] * s;
  c[2] = a[2] * s;
}

inline void ScalarDiv3fv(float v[3], float s)
{
  v[0] /= s;
  v[1] /= s;
  v[2] /= s;
}

inline void Add3fv(float c[3], const float a[3], const float b[3]) // c = a + b
{
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

inline void Subtract3fv(float c[3], const float a[3], const float b[3]) // c = a - b
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

inline void Negate3fv(float a[3], const float b[3])  // a = -b
{
  a[0] = -b[0];
  a[1] = -b[1];
  a[2] = -b[2];
}

inline float Length3fv(const float v[3])
{
  return( (float)sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) );
}

inline void Normalize3fv(float v[3])
{
  float l = Length3fv(v);
  v[0] /= l;
  v[1] /= l;
  v[2] /= l;
}

inline float DotProd3fv(const float a[3], const float b[3])
{
  return( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

inline void CrossProd3fv(float* C, const float* A, const float* B) // C = A X B
{ 
  Set3fv(C, A[1]*B[2]-A[2]*B[1], A[2]*B[0]-A[0]*B[2], A[0]*B[1]-A[1]*B[0]);
}

#endif
