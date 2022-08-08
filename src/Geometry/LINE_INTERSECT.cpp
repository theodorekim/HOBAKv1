/*
This file is part of HOBAK.

Permission is hereby granted to use this software solely for non-commercial applications
and purposes including academic or industrial research, evaluation and not-for-profit media
production.
 
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES, YALE, OR
THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/* Copyright (C) Graham Rhodes, 2001. 
 * All rights reserved worldwide.
 *
 * This software is provided "as is" without express or implied
 * warranties. You may freely copy and compile this source into
 * applications you distribute provided that the copyright text
 * below is included in the resulting source code, for example:
 * "Portions Copyright (C) Graham Rhodes, 2001"
 */
/**************************************************************************************
|
|           File: lineintersect_utils.cpp
|
|        Purpose: Implementation of line segment intersection utility functions
|
|     Book Title: Game Programming Gems II
|
|  Chapter Title: Fast, Robust Intersection of 3D Line Segments
|
|         Author: Graham Rhodes
|
|      Revisions: 05-Apr-2001 - GSR. Original.
|
**************************************************************************************/
#include <math.h>
#include "LINE_INTERSECT.h"
#include <assert.h>

// uncomment the following line to have the code check intermediate results
//#define CHECK_ANSWERS

// uncomment the following line to use Cramer's rule instead of Gaussian elimination
//#define USE_CRAMERS_RULE

#define FMAX(a,b) ((a) > (b) ? (a) : (b))
#define FMIN(a,b) ((a) > (b) ? (b) : (a))
#define FABS(a) ((a) < 0.0f ? -(a) : (a))
#define OUT_OF_RANGE(a) ((a) < 0.0f || (a) > 1.f)

// pragma to get rid of math.h inline function removal warnings.
//#pragma warning(disable:4514)

#include <iostream>
void IntersectLineSegments(const VECTOR3& a0, const VECTOR3& a1,
                           const VECTOR3& b0, const VECTOR3& b1,
                           VECTOR3& aPoint, VECTOR3& bPoint)
{
  VECTOR3 midpoint, normal;
  bool intersect;
  IntersectLineSegments(a0[0], a0[1], a0[2], a1[0], a1[1], a1[2],
                        b0[0], b0[1], b0[2], b1[0], b1[1], b1[2],
                        false, 1e-8,
                        aPoint[0], aPoint[1], aPoint[2],
                        bPoint[0], bPoint[1], bPoint[2],
                        midpoint[0], midpoint[1], midpoint[2],
                        normal[0], normal[1], normal[2], intersect);

  //using namespace std;
  //cout << " Midpoint: " << endl << midpoint << endl;
}

/**************************************************************************
|
|     Method: IntersectLineSegments
|
|    Purpose: Find the nearest point between two finite length line segments
|             or two infinite lines in 3-dimensional space. The function calculates
|             the point on each line/line segment that is closest to the other
|             line/line segment, the midpoint between the nearest points, and
|             the vector between these two points. If the two nearest points
|             are close within a tolerance, a flag is set indicating the lines
|             have a "true" intersection.
|
| Parameters: Input:
|             ------
|             A1x, A1y, A1z   - Coordinates of first defining point of line/segment A
|             A2x, A2y, A2z   - Coordinates of second defining point of line/segment A
|             B1x, B1y, B1z   - Coordinates of first defining point of line/segment B
|             B2x, B2y, B2z   - Coordinates of second defining point of line/segment B
|             infinite_lines  - set to true if lines are to be treated as infinite
|             epsilon         - tolerance value to be used to check for degenerate
|                               and parallel lines, and to check for true intersection.
|
|             Output:
|             -------
|             PointOnSegAx,   - Coordinates of the point on segment A that are nearest
|             PointOnSegAy,     to segment B. This corresponds to point C in the text.
|             PointOnSegAz
|             PointOnSegBx,   - Coordinates of the point on segment B that are nearest
|             PointOnSegBy,     to segment A. This corresponds to point D in the text.
|             PointOnSegBz
|             NearestPointX,  - Midpoint between the two nearest points. This can be
|             NearestPointY,    treated as *the* intersection point if nearest points
|             NearestPointZ     are sufficiently close. This corresponds to point P
|                               in the text.
|             NearestVectorX, - Vector between the nearest point on A to the nearest
|                               point on segment B. This vector is normal to both
|                               lines if the lines are infinite, but is not guaranteed
|                               to be normal to both lines if both lines are finite
|                               length.
|           true_intersection - true if the nearest points are close within a small
|                               tolerance.
**************************************************************************/
void IntersectLineSegments(const REAL A1x, const REAL A1y, const REAL A1z,
                           const REAL A2x, const REAL A2y, const REAL A2z,
                           const REAL B1x, const REAL B1y, const REAL B1z,
                           const REAL B2x, const REAL B2y, const REAL B2z,
                           bool infinite_lines, REAL epsilon, REAL &PointOnSegAx,
                           REAL &PointOnSegAy, REAL &PointOnSegAz, REAL &PointOnSegBx,
                           REAL &PointOnSegBy, REAL &PointOnSegBz, REAL &NearestPointX,
                           REAL &NearestPointY, REAL &NearestPointZ, REAL &NearestVectorX,
                           REAL &NearestVectorY, REAL &NearestVectorZ, bool &true_intersection)
{
  REAL temp = 0.f;
  REAL epsilon_squared = epsilon * epsilon;

  // Compute parameters from Equations (1) and (2) in the text
  REAL Lax = A2x - A1x;
  REAL Lay = A2y - A1y;
  REAL Laz = A2z - A1z;
  REAL Lbx = B2x - B1x;
  REAL Lby = B2y - B1y;
  REAL Lbz = B2z - B1z;
  // From Equation (15)
  REAL L11 =  (Lax * Lax) + (Lay * Lay) + (Laz * Laz);
  REAL L22 =  (Lbx * Lbx) + (Lby * Lby) + (Lbz * Lbz);

  // Line/Segment A is degenerate ---- Special Case #1
  if (L11 < epsilon_squared)
  {
    PointOnSegAx = A1x;
    PointOnSegAy = A1y;
    PointOnSegAz = A1z;
    FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz, A1x, A1y, A1z,
                                  infinite_lines, epsilon, PointOnSegBx, PointOnSegBy,
                                  PointOnSegBz, temp);
  }
  // Line/Segment B is degenerate ---- Special Case #1
  else if (L22 < epsilon_squared)
  {
    PointOnSegBx = B1x;
    PointOnSegBy = B1y;
    PointOnSegBz = B1z;
    FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, B1x, B1y, B1z,
                                  infinite_lines, epsilon, PointOnSegAx, PointOnSegAy,
                                  PointOnSegAz, temp);
  }
  // Neither line/segment is degenerate
  else
  {
    // Compute more parameters from Equation (3) in the text.
    REAL ABx = B1x - A1x;
    REAL ABy = B1y - A1y;
    REAL ABz = B1z - A1z;

    // and from Equation (15).
    REAL L12 = -(Lax * Lbx) - (Lay * Lby) - (Laz * Lbz);

    REAL DetL = L11 * L22 - L12 * L12;
    // Lines/Segments A and B are parallel ---- special case #2.
    if (FABS(DetL) < epsilon)
    {
      FindNearestPointOfParallelLineSegments(A1x, A1y, A1z, A2x, A2y, A2z,
                                             Lax, Lay, Laz,
                                             B1x, B1y, B1z, B2x, B2y, B2z,
                                             Lbx, Lby, Lbz,
                                             infinite_lines, epsilon,
                                             PointOnSegAx, PointOnSegAy, PointOnSegAz,
                                             PointOnSegBx, PointOnSegBy, PointOnSegBz);
    }
    // The general case
    else
    {
      // from Equation (15)
      REAL ra = Lax * ABx + Lay * ABy + Laz * ABz;
      REAL rb = -Lbx * ABx - Lby * ABy - Lbz * ABz;

      REAL t = (L11 * rb - ra * L12)/DetL; // Equation (12)

#ifdef USE_CRAMERS_RULE
      REAL s = (L22 * ra - rb * L12)/DetL;
#else
      REAL s = (ra-L12*t)/L11;             // Equation (13)
#endif // USE_CRAMERS_RULE

#ifdef CHECK_ANSWERS
      REAL check_ra = s*L11 + t*L12;
      REAL check_rb = s*L12 + t*L22;
      assert(FABS(check_ra-ra) < epsilon);
      assert(FABS(check_rb-rb) < epsilon);
#endif // CHECK_ANSWERS

      // if we are dealing with infinite lines or if parameters s and t both
      // lie in the range [0,1] then just compute the points using Equations
      // (1) and (2) from the text.
      PointOnSegAx = (A1x + s * Lax);
      PointOnSegAy = (A1y + s * Lay);
      PointOnSegAz = (A1z + s * Laz);
      PointOnSegBx = (B1x + t * Lbx);
      PointOnSegBy = (B1y + t * Lby);
      PointOnSegBz = (B1z + t * Lbz);
      // otherwise, at least one of s and t is outside of [0,1] and we have to
      // handle this case.
      if (false == infinite_lines && (OUT_OF_RANGE(s) || OUT_OF_RANGE(t)))
      {
        AdjustNearestPoints(A1x, A1y, A1z, Lax, Lay, Laz,
                            B1x, B1y, B1z, Lbx, Lby, Lbz,
                            epsilon, s, t,
                            PointOnSegAx, PointOnSegAy, PointOnSegAz,
                            PointOnSegBx, PointOnSegBy, PointOnSegBz);
      }
    }
  }

  NearestPointX = 0.5f * (PointOnSegAx + PointOnSegBx);
  NearestPointY = 0.5f * (PointOnSegAy + PointOnSegBy);
  NearestPointZ = 0.5f * (PointOnSegAz + PointOnSegBz);

  NearestVectorX = PointOnSegBx - PointOnSegAx;
  NearestVectorY = PointOnSegBy - PointOnSegAy;
  NearestVectorZ = PointOnSegBz - PointOnSegAz;

  // optional check to indicate if the lines truly intersect
  true_intersection = (FABS(NearestVectorX) +
                       FABS(NearestVectorY) +
                       FABS(NearestVectorZ)) < epsilon ? true : false;
}

/**************************************************************************
|
|     Method: FindNearestPointOnLineSegment
|
|    Purpose: Given a line (segment) and a point in 3-dimensional space,
|             find the point on the line (segment) that is closest to the
|             point.
|
| Parameters: Input:
|             ------
|             A1x, A1y, A1z   - Coordinates of first defining point of the line/segment
|             Lx, Ly, Lz      - Vector from (A1x, A1y, A1z) to the second defining point
|                               of the line/segment.
|             Bx, By, Bz      - Coordinates of the point
|             infinite_lines  - set to true if lines are to be treated as infinite
|             epsilon_squared - tolerance value to be used to check for degenerate
|                               and parallel lines, and to check for true intersection.
|
|             Output:
|             -------
|             NearestPointX,  - Point on line/segment that is closest to (Bx, By, Bz)
|             NearestPointY,
|             NearestPointZ
|             parameter       - Parametric coordinate of the nearest point along the
|                               line/segment. parameter = 0 at (A1x, A1y, A1z) and
|                               parameter = 1 at the second defining point of the line/
|                               segmetn
**************************************************************************/
void FindNearestPointOnLineSegment(const REAL A1x, const REAL A1y, const REAL A1z,
                                   const REAL Lx, const REAL Ly, const REAL Lz,
                                   const REAL Bx, const REAL By, const REAL Bz,
                                   bool infinite_line, REAL epsilon_squared, REAL &NearestPointX,
                                   REAL &NearestPointY, REAL &NearestPointZ,
                                   REAL &parameter)
{
  // Line/Segment is degenerate --- special case #1
  REAL D = Lx * Lx + Ly * Ly + Lz * Lz;
  if (D < epsilon_squared)
  {
    NearestPointX = A1x;
    NearestPointY = A1y;
    NearestPointZ = A1z;
    return;
  }

  REAL ABx = Bx - A1x;
  REAL ABy = By - A1y;
  REAL ABz = Bz - A1z;

  // parameter is computed from Equation (20).
  parameter = (Lx * ABx + Ly * ABy + Lz * ABz) / D;

  if (false == infinite_line) parameter = FMAX(0.0f, FMIN(1.0f, parameter));

  NearestPointX = A1x + parameter * Lx;
  NearestPointY = A1y + parameter * Ly;
  NearestPointZ = A1z + parameter * Lz;
  return;
}

/**************************************************************************
|
|     Method: FindNearestPointOfParallelLineSegments
|
|    Purpose: Given two lines (segments) that are known to be parallel, find
|             a representative point on each that is nearest to the other. If
|             the lines are considered to be finite then it is possible that there
|             is one true point on each line that is nearest to the other. This
|             code properly handles this case.
|
|             This is the most difficult line intersection case to handle, since
|             there is potentially a family, or locus of points on each line/segment
|             that are nearest to the other.
| Parameters: Input:
|             ------
|             A1x, A1y, A1z   - Coordinates of first defining point of line/segment A
|             A2x, A2y, A2z   - Coordinates of second defining point of line/segment A
|             Lax, Lay, Laz   - Vector from (A1x, A1y, A1z) to the (A2x, A2y, A2z).
|             B1x, B1y, B1z   - Coordinates of first defining point of line/segment B
|             B2x, B2y, B2z   - Coordinates of second defining point of line/segment B
|             Lbx, Lby, Lbz   - Vector from (B1x, B1y, B1z) to the (B2x, B2y, B2z).
|             infinite_lines  - set to true if lines are to be treated as infinite
|             epsilon_squared - tolerance value to be used to check for degenerate
|                               and parallel lines, and to check for true intersection.
|
|             Output:
|             -------
|             PointOnSegAx,   - Coordinates of the point on segment A that are nearest
|             PointOnSegAy,     to segment B. This corresponds to point C in the text.
|             PointOnSegAz
|             PointOnSegBx,   - Coordinates of the point on segment B that are nearest
|             PointOnSegBy,     to segment A. This corresponds to point D in the text.
|             PointOnSegBz

**************************************************************************/
void FindNearestPointOfParallelLineSegments(REAL A1x, REAL A1y, REAL A1z,
                                            REAL A2x, REAL A2y, REAL A2z,
                                            REAL Lax, REAL Lay, REAL Laz,
                                            REAL B1x, REAL B1y, REAL B1z,
                                            REAL B2x, REAL B2y, REAL B2z,
                                            REAL Lbx, REAL Lby, REAL Lbz,
                                            bool infinite_lines, REAL epsilon_squared,
                                            REAL &PointOnSegAx, REAL &PointOnSegAy, REAL &PointOnSegAz,
                                            REAL &PointOnSegBx, REAL &PointOnSegBy, REAL &PointOnSegBz)
{
  REAL s[2] = {0, 0};
  REAL temp;
  FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, B1x, B1y, B1z,
                                true, epsilon_squared, PointOnSegAx, PointOnSegAy, PointOnSegAz, s[0]);
  if (true == infinite_lines)
  {
    PointOnSegBx = B1x;
    PointOnSegBy = B1y;
    PointOnSegBz = B1z;
  }
  else
  {
    REAL tp[3];
    FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, B2x, B2y, B2z,
                                  true, epsilon_squared, tp[0], tp[1], tp[2], s[1]);
    if (s[0] < 0.f && s[1] < 0.f)
    {
      PointOnSegAx = A1x;
      PointOnSegAy = A1y;
      PointOnSegAz = A1z;
      if (s[0] < s[1])
      {
        PointOnSegBx = B2x;
        PointOnSegBy = B2y;
        PointOnSegBz = B2z;
      }
      else
      {
        PointOnSegBx = B1x;
        PointOnSegBy = B1y;
        PointOnSegBz = B1z;
      }
    }
    else if (s[0] > 1.f && s[1] > 1.f)
    {
      PointOnSegAx = A2x;
      PointOnSegAy = A2y;
      PointOnSegAz = A2z;
      if (s[0] < s[1])
      {
        PointOnSegBx = B1x;
        PointOnSegBy = B1y;
        PointOnSegBz = B1z;
      }
      else
      {
        PointOnSegBx = B2x;
        PointOnSegBy = B2y;
        PointOnSegBz = B2z;
      }
    }
    else
    {
      temp = 0.5f*(FMAX(0.0f, FMIN(1.0f, s[0])) + FMAX(0.0f, FMIN(1.0f, s[1])));
      PointOnSegAx = (A1x + temp * Lax);
      PointOnSegAy = (A1y + temp * Lay);
      PointOnSegAz = (A1z + temp * Laz);
      FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz,
                                    PointOnSegAx, PointOnSegAy, PointOnSegAz, true,
                                    epsilon_squared, PointOnSegBx, PointOnSegBy, PointOnSegBz, temp);
    }
  }
}

/**************************************************************************
|
|     Method: AdjustNearestPoints
|
|    Purpose: Given nearest point information for two infinite lines, adjust
|             to model finite line segments.
|
| Parameters: Input:
|             ------
|             A1x, A1y, A1z   - Coordinates of first defining point of line/segment A
|             Lax, Lay, Laz   - Vector from (A1x, A1y, A1z) to the (A2x, A2y, A2z).
|             B1x, B1y, B1z   - Coordinates of first defining point of line/segment B
|             Lbx, Lby, Lbz   - Vector from (B1x, B1y, B1z) to the (B2x, B2y, B2z).
|             epsilon_squared - tolerance value to be used to check for degenerate
|                               and parallel lines, and to check for true intersection.
|             s               - parameter representing nearest point on infinite line A
|             t               - parameter representing nearest point on infinite line B
|
|             Output:
|             -------
|             PointOnSegAx,   - Coordinates of the point on segment A that are nearest
|             PointOnSegAy,     to segment B. This corresponds to point C in the text.
|             PointOnSegAz
|             PointOnSegBx,   - Coordinates of the point on segment B that are nearest
|             PointOnSegBy,     to segment A. This corresponds to point D in the text.
|             PointOnSegBz
**************************************************************************/
void AdjustNearestPoints(REAL A1x, REAL A1y, REAL A1z,
                         REAL Lax, REAL Lay, REAL Laz,
                         REAL B1x, REAL B1y, REAL B1z,
                         REAL Lbx, REAL Lby, REAL Lbz,
                         REAL epsilon_squared, REAL s, REAL t,
                         REAL &PointOnSegAx, REAL &PointOnSegAy, REAL &PointOnSegAz,
                         REAL &PointOnSegBx, REAL &PointOnSegBy, REAL &PointOnSegBz)
{
  // handle the case where both parameter s and t are out of range
  if (OUT_OF_RANGE(s) && OUT_OF_RANGE(t))
  {
    s = FMAX(0.0f, FMIN(1.0f, s));
    PointOnSegAx = (A1x + s * Lax);
    PointOnSegAy = (A1y + s * Lay);
    PointOnSegAz = (A1z + s * Laz);
    FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz, PointOnSegAx,
                                  PointOnSegAy, PointOnSegAz, true, epsilon_squared,
                                  PointOnSegBx, PointOnSegBy, PointOnSegBz, t);
    if (OUT_OF_RANGE(t))
    {
      t = FMAX(0.0f, FMIN(1.0f, t));
      PointOnSegBx = (B1x + t * Lbx);
      PointOnSegBy = (B1y + t * Lby);
      PointOnSegBz = (B1z + t * Lbz);
      FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, PointOnSegBx,
                                    PointOnSegBy, PointOnSegBz, false, epsilon_squared,
                                    PointOnSegAx, PointOnSegAy, PointOnSegAz, s);
      FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz, PointOnSegAx,
                                    PointOnSegAy, PointOnSegAz, false, epsilon_squared,
                                    PointOnSegBx, PointOnSegBy, PointOnSegBz, t);
    }
  }
  // otherwise, handle the case where the parameter for only one segment is
  // out of range
  else if (OUT_OF_RANGE(s))
  {
    s = FMAX(0.0f, FMIN(1.0f, s));
    PointOnSegAx = (A1x + s * Lax);
    PointOnSegAy = (A1y + s * Lay);
    PointOnSegAz = (A1z + s * Laz);
    FindNearestPointOnLineSegment(B1x, B1y, B1z, Lbx, Lby, Lbz, PointOnSegAx,
                                  PointOnSegAy, PointOnSegAz, false, epsilon_squared,
                                  PointOnSegBx, PointOnSegBy, PointOnSegBz, t);
  }
  else if (OUT_OF_RANGE(t))
  {
    t = FMAX(0.0f, FMIN(1.0f, t));
    PointOnSegBx = (B1x + t * Lbx);
    PointOnSegBy = (B1y + t * Lby);
    PointOnSegBz = (B1z + t * Lbz);
    FindNearestPointOnLineSegment(A1x, A1y, A1z, Lax, Lay, Laz, PointOnSegBx,
                                  PointOnSegBy, PointOnSegBz, false, epsilon_squared,
                                  PointOnSegAx, PointOnSegAy, PointOnSegAz, s);
  }
  else
  {
    assert(0);
  }
}
