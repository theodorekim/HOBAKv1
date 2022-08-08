//------------------------------------------------------------------------------
// File : mat16fv.hpp
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
// mat16fv.hpp : opengl-style float[16] matrix routines.
//----------------------------------------------------------------------------
// $Id: mat16fv.hpp,v 1.5 2002/03/13 08:21:44 harrism Exp $
//============================================================================

#ifndef _MAT16FV_H_
#define _MAT16FV_H_

#include <stdio.h>
#include <math.h>
#include "glvuVec3fv.h"

inline float* Copy16fv(float* A, const float* B) // A=B
{
  A[0]=B[0];   A[1]=B[1];   A[2]=B[2];   A[3]=B[3]; 
  A[4]=B[4];   A[5]=B[5];   A[6]=B[6];   A[7]=B[7]; 
  A[8]=B[8];   A[9]=B[9];   A[10]=B[10]; A[11]=B[11]; 
  A[12]=B[12]; A[13]=B[13]; A[14]=B[14]; A[15]=B[15]; 
  return A;
}

inline float* Mult16fv(float* C, const float* A, const float* B) // C=A*B
{
  float tC[16];

  tC[0]  = A[0]*B[0] + A[4]*B[1] + A[8]*B[2] + A[12]*B[3];
  tC[1]  = A[1]*B[0] + A[5]*B[1] + A[9]*B[2] + A[13]*B[3];
  tC[2]  = A[2]*B[0] + A[6]*B[1] + A[10]*B[2] + A[14]*B[3];
  tC[3]  = A[3]*B[0] + A[7]*B[1] + A[11]*B[2] + A[15]*B[3];

  tC[4]  = A[0]*B[4] + A[4]*B[5] + A[8]*B[6] + A[12]*B[7];
  tC[5]  = A[1]*B[4] + A[5]*B[5] + A[9]*B[6] + A[13]*B[7];
  tC[6]  = A[2]*B[4] + A[6]*B[5] + A[10]*B[6] + A[14]*B[7];
  tC[7]  = A[3]*B[4] + A[7]*B[5] + A[11]*B[6] + A[15]*B[7];

  tC[8]  = A[0]*B[8] + A[4]*B[9] + A[8]*B[10] + A[12]*B[11];
  tC[9]  = A[1]*B[8] + A[5]*B[9] + A[9]*B[10] + A[13]*B[11];
  tC[10] = A[2]*B[8] + A[6]*B[9] + A[10]*B[10] + A[14]*B[11];
  tC[11] = A[3]*B[8] + A[7]*B[9] + A[11]*B[10] + A[15]*B[11];

  tC[12] = A[0]*B[12] + A[4]*B[13] + A[8]*B[14] + A[12]*B[15];
  tC[13] = A[1]*B[12] + A[5]*B[13] + A[9]*B[14] + A[13]*B[15];
  tC[14] = A[2]*B[12] + A[6]*B[13] + A[10]*B[14] + A[14]*B[15];
  tC[15] = A[3]*B[12] + A[7]*B[13] + A[11]*B[14] + A[15]*B[15];

  Copy16fv(C,tC);

  return(C);
}

inline float* Mult16fv3fv(float *NewV, const float* M, const float *V)
{
  NewV[0] = M[0]*V[0] + M[4]*V[1] + M[8]*V[2];
  NewV[1] = M[1]*V[0] + M[5]*V[1] + M[9]*V[2];
  NewV[2] = M[2]*V[0] + M[6]*V[1] + M[10]*V[2];
  return(NewV);
}

inline float* Mult16fv3fvPerspDiv(float *NewV, const float* M, const float *V)
{
  float W = M[3]*V[0] + M[7]*V[1] + M[11]*V[2] + M[15];
  NewV[0] = (M[0]*V[0] + M[4]*V[1] + M[8]*V[2] + M[12]) / W;
  NewV[1] = (M[1]*V[0] + M[5]*V[1] + M[9]*V[2] + M[13]) / W;
  NewV[2] = (M[2]*V[0] + M[6]*V[1] + M[10]*V[2] + M[14]) / W;
  return(NewV);
}

inline float* Mult16fv4fv(float *NewV, const float* M, const float *V)
{
  NewV[0] = M[0]*V[0] + M[4]*V[1] + M[8]*V[2] + M[12]*V[3];
  NewV[1] = M[1]*V[0] + M[5]*V[1] + M[9]*V[2] + M[13]*V[3];
  NewV[2] = M[2]*V[0] + M[6]*V[1] + M[10]*V[2] + M[14]*V[3];
  NewV[3] = M[3]*V[0] + M[7]*V[1] + M[11]*V[2] + M[15]*V[3];
  return(NewV);
}

inline float* Identity16fv(float* M)
{
  M[0]=M[5]=M[10]=M[15]=1;
  M[1]=M[2]=M[3]=M[4]=M[6]=M[7]=M[8]=M[9]=M[11]=M[12]=M[13]=M[14]=0;
  return(M);
}

inline float* Transpose16fv(float* M)
{
  #define SWAP(a,b,t) (t)=(a);(a)=(b);(b)=(t);
  float t;
  SWAP(M[1],M[4],t);
  SWAP(M[2],M[8],t);
  SWAP(M[6],M[9],t);
  SWAP(M[3],M[12],t);
  SWAP(M[7],M[13],t);
  SWAP(M[11],M[14],t);
  return(M);
}

inline float* Rotate16fv(float *M, float DegAng, const float Axis[3])
{
  float RadAng = DegAng * 0.0174532f;
  float ca=(float)cos(RadAng),
        sa=(float)sin(RadAng);
  if (Axis[0]==1 && Axis[1]==0 && Axis[2]==0)  // ABOUT X-AXIS
  {
   M[0]=1; M[4]=0;  M[8]=0;   M[12]=0;
   M[1]=0; M[5]=ca; M[9]=-sa; M[13]=0;
   M[2]=0; M[6]=sa; M[10]=ca; M[14]=0;
   M[3]=0; M[7]=0;  M[11]=0;  M[15]=1;
  }
  if (Axis[0]==0 && Axis[1]==1 && Axis[2]==0)  // ABOUT Y-AXIS
  {
   M[0]=ca;  M[4]=0; M[8]=sa;  M[12]=0;
   M[1]=0;   M[5]=1; M[9]=0;   M[13]=0;
   M[2]=-sa; M[6]=0; M[10]=ca; M[14]=0;
   M[3]=0;   M[7]=0; M[11]=0;  M[15]=1;
  }
  if (Axis[0]==0 && Axis[1]==0 && Axis[2]==1)  // ABOUT Z-AXIS
  {
   M[0]=ca; M[4]=-sa; M[8]=0;  M[12]=0;
   M[1]=sa; M[5]=ca;  M[9]=0;  M[13]=0;
   M[2]=0;  M[6]=0;   M[10]=1; M[14]=0;
   M[3]=0;  M[7]=0;   M[11]=0; M[15]=1;
  }
  else                                      // ARBITRARY AXIS
  {
   float l = Axis[0]*Axis[0]+Axis[1]*Axis[1]+Axis[2]*Axis[2];
   float x, y, z;
   x=Axis[0],y=Axis[1],z=Axis[2];
   if (l > 1.0001f || (l < 0.9999f && l!=0))
   {
     // needs normalization
     l=1.0f/(float)sqrt(l);
     x*=l; y*=l; z*=l;
   }
   float x2=x*x, y2=y*y, z2=z*z;
   M[0]=x2+ca*(1-x2); M[4]=(x*y)+ca*(-x*y)+sa*(-z); M[8]=(x*z)+ca*(-x*z)+sa*y;
   M[1]=(x*y)+ca*(-x*y)+sa*z; M[5]=y2+ca*(1-y2); M[9]=(y*z)+ca*(-y*z)+sa*(-x);
   M[2]=(x*z)+ca*(-x*z)+sa*(-y); M[6]=(y*z)+ca*(-y*z)+sa*x; M[10]=z2+ca*(1-z2);
   M[12]=M[13]=M[14]=M[3]=M[7]=M[11]=0;
   M[15]=1;
  }
  return(M);
}

inline float* invRotate16fv(float *M, float DegAng, const float Axis[3])
{
  Rotate16fv(M,DegAng,Axis);
  Transpose16fv(M);
  return(M);
}

inline float* Scale16fv(float* M, float sx, float sy, float sz)
{
  M[0]=sx; M[4]=0;  M[8]=0;   M[12]=0;
  M[1]=0;  M[5]=sy; M[9]=0;   M[13]=0;
  M[2]=0;  M[6]=0;  M[10]=sz; M[14]=0;
  M[3]=0;  M[7]=0;  M[11]=0;  M[15]=1;
  return(M);
}

inline float* invScale16fv(float* M, float sx, float sy, float sz)
{
  M[0]=1/sx; M[4]=0;    M[8]=0;     M[12]=0;
  M[1]=0;    M[5]=1/sy; M[9]=0;     M[13]=0;
  M[2]=0;    M[6]=0;    M[10]=1/sz; M[14]=0;
  M[3]=0;    M[7]=0;    M[11]=0;    M[15]=1;
  return(M);
}

inline float* Translate16fv(float* M, float tx, float ty, float tz)
{
  M[0]=1; M[4]=0;  M[8]=0;  M[12]=tx;
  M[1]=0; M[5]=1;  M[9]=0;  M[13]=ty;
  M[2]=0; M[6]=0;  M[10]=1; M[14]=tz;
  M[3]=0; M[7]=0;  M[11]=0; M[15]=1;
  return(M);
}

inline float* invTranslate16fv(float* M, float tx, float ty, float tz)
{
  M[0]=1; M[4]=0;  M[8]=0;  M[12]=-tx;
  M[1]=0; M[5]=1;  M[9]=0;  M[13]=-ty;
  M[2]=0; M[6]=0;  M[10]=1; M[14]=-tz;
  M[3]=0; M[7]=0;  M[11]=0; M[15]=1;
  return(M);
}

inline float* LookAt(float* M,
              const float Eye[3], 
              const float LookAtPt[3],
              const float ViewUp[3])
{
  float X[3], Y[3], Z[3];
  Subtract3fv(Z,Eye,LookAtPt);  Normalize3fv(Z);
  CrossProd3fv(X,ViewUp,Z);     Normalize3fv(X);
  CrossProd3fv(Y,Z,X);          Normalize3fv(Y);
  M[0]=X[0];  M[4]=X[1];  M[8]=X[2];  M[12]=-DotProd3fv(X,Eye);  // TRANS->ROT
  M[1]=Y[0];  M[5]=Y[1];  M[9]=Y[2];  M[13]=-DotProd3fv(Y,Eye);
  M[2]=Z[0];  M[6]=Z[1];  M[10]=Z[2]; M[14]=-DotProd3fv(Z,Eye);
  M[3]=0;     M[7]=0;     M[11]=0;    M[15]=1;
  return(M);
}

inline float* invLookAt(float* M, 
                 const float Eye[3], 
                 const float LookAtPt[3], 
                 const float ViewUp[3])
{
  float X[3], Y[3], Z[3];
  Subtract3fv(Z,Eye,LookAtPt);  Normalize3fv(Z);
  CrossProd3fv(X,ViewUp,Z);     Normalize3fv(X);
  CrossProd3fv(Y,Z,X);          Normalize3fv(Y);
  M[0]=X[0];  M[4]=Y[0];  M[8]=Z[0];  M[12]=Eye[0];  // ROT->TRANS
  M[1]=X[1];  M[5]=Y[1];  M[9]=Z[1];  M[13]=Eye[1];
  M[2]=X[2];  M[6]=Y[2];  M[10]=Z[2]; M[14]=Eye[2];
  M[3]=0;     M[7]=0;     M[11]=0;    M[15]=1;
  return(M);
}

inline float* Frustum16fv(float* M, float l, float r, float b, float t, 
                   float n, float f)
{
  M[0]=(2*n)/(r-l); M[4]=0;           M[8]=(r+l)/(r-l);   M[12]=0;
  M[1]=0;           M[5]=(2*n)/(t-b); M[9]=(t+b)/(t-b);   M[13]=0;
  M[2]=0;           M[6]=0;           M[10]=-(f+n)/(f-n); M[14]=(-2*f*n)/(f-n);
  M[3]=0;           M[7]=0;           M[11]=-1;           M[15]=0;
  return(M);  
}

inline float* invFrustum16fv(float* M, float l, float r, float b, float t, 
                      float n, float f)
{
  M[0]=(r-l)/(2*n); M[4]=0;           M[8]=0;               M[12]=(r+l)/(2*n);
  M[1]=0;           M[5]=(t-b)/(2*n); M[9]=0;               M[13]=(t+b)/(2*n);
  M[2]=0;           M[6]=0;           M[10]=0;              M[14]=-1;
  M[3]=0;           M[7]=0;           M[11]=-(f-n)/(2*f*n); M[15]=(f+n)/(2*f*n);
  return(M);  
}

inline float* Perspective(float* M, float Yfov, float Aspect, 
                   float Ndist, float Fdist)
{
  Yfov *= 0.0174532f;  // CONVERT TO RADIANS
  float wT=(float)tan(Yfov*0.5f)*Ndist, wB=-wT;
  float wR=wT*Aspect, wL=-wR;
  Frustum16fv(M,wL,wR,wB,wT,Ndist,Fdist);
  return(M);
}

inline float* invPerspective(float* M, float Yfov, float Aspect, 
                      float Ndist, float Fdist)
{
  Yfov *= 0.0174532f;  // CONVERT TO RADIANS
  float wT=(float)tan(Yfov*0.5f)*Ndist, wB=-wT;
  float wR=wT*Aspect, wL=-wR;
  invFrustum16fv(M,wL,wR,wB,wT,Ndist,Fdist);
  return(M);
}

inline float* Viewing16fv(
  float* M,
  const float X[3], const float Y[3], const float Z[3], const float O[3])
{
  M[0]=X[0];  M[4]=X[1];   M[8]=X[2]; M[12]=-DotProd3fv(X,O);
  M[1]=Y[0];  M[5]=Y[1];   M[9]=Y[2]; M[13]=-DotProd3fv(Y,O);
  M[2]=Z[0];  M[6]=Z[1];  M[10]=Z[2]; M[14]=-DotProd3fv(Z,O);
  M[3]=0;    M[7]=0;    M[11]=0;   M[15]=1;
  return(M);
}

// THE INVERSE OF Viewing16fv,
// THIS TAKES A VIEW MATRIX AND RETURNS VIEWING AXES.  
// MATRIX ASSUMED TO BE ORTHONORMAL
inline void Viewing2CoordFrame16fv(
  const float *M, float X[3], float Y[3], float Z[3], float O[3])
{
  X[0]=M[0];  X[1]=M[4];   X[2]=M[8]; O[0]=-DotProd3fv(M,M+12);
  Y[0]=M[1];  Y[1]=M[5];   Y[2]=M[9]; O[1]=-DotProd3fv(M+4,M+12);
  Z[0]=M[2];  Z[1]=M[6];  Z[2]=M[10]; O[2]=-DotProd3fv(M+8,M+12);
};

inline float* invViewing16fv(
  float* M, 
  const float X[3], const float Y[3], const float Z[3], const float O[3])
{  
  M[0]=X[0];  M[4]=Y[0];   M[8]=Z[0]; M[12]=O[0];
  M[1]=X[1];  M[5]=Y[1];   M[9]=Z[1]; M[13]=O[1];
  M[2]=X[2];  M[6]=Y[2];  M[10]=Z[2]; M[14]=O[2];
  M[3]=0;    M[7]=0;    M[11]=0;   M[15]=1;
  return(M);
}

inline float* Viewport16fv(float* M, int WW, int WH)
{
  float WW2=(float)WW*0.5f, WH2=(float)WH*0.5f;
  M[0]=WW2;  M[4]=0;     M[8]=0;     M[12]=WW2;
  M[1]=0;    M[5]=WH2;   M[9]=0;     M[13]=WH2;
  M[2]=0;    M[6]=0;     M[10]=0.5f; M[14]=0.5f;
  M[3]=0;    M[7]=0;     M[11]=0;    M[15]=1;
  return(M);
}

inline float* invViewport16fv(float* M, int WW, int WH)
{
  float WW2=2.0f/(float)WW, WH2=2.0f/(float)WH;
  M[0]=WW2;  M[4]=0;     M[8]=0;    M[12]=-1.0;
  M[1]=0;    M[5]=WH2;   M[9]=0;    M[13]=-1.0;
  M[2]=0;    M[6]=0;     M[10]=2.0; M[14]=-1.0;
  M[3]=0;    M[7]=0;     M[11]=0;   M[15]=1;
  return(M);
}

//--------------------------------------------------------------------------
// Given the coefficient [A B C D] of a plane in the implicit form
// Ax+By+Cz+D=0 (see plane.hpp), this routine generates a reflection matrix
// that will "reflect" all points/vectors about the given plane.
// NOTE: the plane is assumed to be normalized: normal vector (A,B,C) is
// normalized (unit-length), and D is the negative distance from the origin
// to the plane along the normal.
//--------------------------------------------------------------------------
inline float* PlanarReflection16fv(float M[16], const float P[4])
{
  float AA=P[0]*P[0], AB=P[0]*P[1], AC=P[0]*P[2], AD=P[0]*P[3],
        BB=P[1]*P[1], BC=P[1]*P[2], BD=P[1]*P[3],
        CC=P[2]*P[2], CD=P[2]*P[3];
  M[0]=1-2*AA;  M[4]=-2*AB;   M[8]=-2*AC;   M[12]=-2*AD;
  M[1]=-2*AB;   M[5]=1-2*BB;  M[9]=-2*BC;   M[13]=-2*BD;
  M[2]=-2*AC;   M[6]=-2*BC;   M[10]=1-2*CC; M[14]=-2*CD;
  M[3]=0;       M[7]=0;       M[11]=0;      M[15]=1;
  return(M);
}

//--------------------------------------------------------------------------
// Returns a matrix that will xform a point from the given object-space
// coordinate frame to world space
//   A composite matrix is formed as follows:
//   C = Translation * Rotation * Scaling (using column vectors/pre-multiplication)
//   WorldPt = C * ModelPt
// The corresponding inverse Xform is also provided (GetWorld2ObjXform)
//--------------------------------------------------------------------------
inline float* Obj2WorldXform16fv(
  float *M, 
  const float X[3], const float Y[3], const float Z[3], 
  const float O[3], float Scale)
{
  float sX[3], sY[3], sZ[3]; // CREATE SCALED VERSION OF ROT AXES
  ScalarMult3fv(sX,X,Scale);
  ScalarMult3fv(sY,Y,Scale);
  ScalarMult3fv(sZ,Z,Scale);
  M[0]=sX[0]; M[4]=sY[0]; M[8]=sZ[0];  M[12]=O[0];
  M[1]=sX[1]; M[5]=sY[1]; M[9]=sZ[1];  M[13]=O[1];
  M[2]=sX[2]; M[6]=sY[2]; M[10]=sZ[2]; M[14]=O[2];
  M[3]=0;     M[7]=0;     M[11]=0;     M[15]=1;
  return(M);
}

inline float* World2ObjXform16fv(
  float *M,
  const float X[3], const float Y[3], const float Z[3], 
  const float O[3], float Scale)
{
  if (Scale<=0) { printf("Too small scale!\n"); Scale=1; }
  float invScale = 1/Scale;
  float sX[3], sY[3], sZ[3];
  ScalarMult3fv(sX,X,invScale);
  ScalarMult3fv(sY,Y,invScale);
  ScalarMult3fv(sZ,Z,invScale);
  M[0]=sX[0]; M[4]=sX[1]; M[8]=sX[2];  M[12]=-DotProd3fv(O,sX);
  M[1]=sY[0]; M[5]=sY[1]; M[9]=sY[2];  M[13]=-DotProd3fv(O,sY);
  M[2]=sZ[0]; M[6]=sZ[1]; M[10]=sZ[2]; M[14]=-DotProd3fv(O,sZ);
  M[3]=0;     M[7]=0;     M[11]=0;     M[15]=1;
  return(M);
}

// ONLY TRANSLATES, ROTATES, AND SCALES ARE ALLOWED
inline float XformCoordFrame16fv(
  const float *M, float X[3], float Y[3], float Z[3], float O[3])
{
  Set3fv(X, X[0]*M[0] + X[1]*M[4] + X[2]*M[8],
             X[0]*M[1] + X[1]*M[5] + X[2]*M[9],
             X[0]*M[2] + X[1]*M[6] + X[2]*M[10] );
  Set3fv(Y, Y[0]*M[0] + Y[1]*M[4] + Y[2]*M[8],
             Y[0]*M[1] + Y[1]*M[5] + Y[2]*M[9],
             Y[0]*M[2] + Y[1]*M[6] + Y[2]*M[10] );
  Set3fv(Z, Z[0]*M[0] + Z[1]*M[4] + Z[2]*M[8],
             Z[0]*M[1] + Z[1]*M[5] + Z[2]*M[9],
             Z[0]*M[2] + Z[1]*M[6] + Z[2]*M[10] );
  Set3fv(O, O[0]*M[0] + O[1]*M[4] + O[2]*M[8] + M[12],
             O[0]*M[1] + O[1]*M[5] + O[2]*M[9] + M[13],
             O[0]*M[2] + O[1]*M[6] + O[2]*M[10] + M[14] );

  // MUST RENORMALIZE AXES TO FIND THE UNIFORM SCALE
  float Scale = Length3fv(X);
  ScalarDiv3fv(X,Scale);
  ScalarDiv3fv(Y,Scale);
  ScalarDiv3fv(Z,Scale);

  // RETURN UNIFORM SCALING OF AXES (how much coordinate frame was scaled)
  return(Scale);
};


//----------------------------------------------------------------------------
// Given a complete definition for a particular view (viewing, projection,
// and viewport), returns the COMPOSITE xform matrix that takes a point in the 
// world space to a screen space (pixel) point. The inverse is also provided.
//----------------------------------------------------------------------------
inline float* Screen2WorldXform16fv(
  float* M, 
  const float X[3], const float Y[3], const float Z[3],  // VIEWING AXES
  const float O[3], // VIEWING ORIGIN
  float l, float r, float b, float t, float n, float f, // PROJECTION
  int WW, int WH) // VIEWPORT
{
  // C = InverseModelview * InverseProjection * InverseViewport
  float N[16];
  invViewing16fv(M,X,Y,Z,O);
  invFrustum16fv(N,l,r,b,t,n,f);
  Mult16fv(M,M,N);  // M=M*N;
  invViewport16fv(N,WW,WH);
  Mult16fv(M,M,N);  // M=M*N;
  return(M);
}

inline float* World2ScreenXform16fv(
  float* M,
  const float X[3], const float Y[3], const float Z[3], // VIEWING AXES
  const float O[3], // VIEWING ORIGIN
  float l, float r, float b, float t, float n, float f, // PROJECTION
  int WW, int WH) // VIEWPORT
{
  // C = Viewport * Projection * Modelview
  float N[16];
  Viewport16fv(M,WW,WH);
  Frustum16fv(N,l,r,b,t,n,f);
  Mult16fv(M,M,N);  // M=M*N;
  Viewing16fv(N,X,Y,Z,O);
  Mult16fv(M,M,N);  // M=M*N;
  return(M);
}

#endif
