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
|           File: lineintersect_utils.h
|
|        Purpose: Function prototypes for line segment intersection utility functions
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
#ifndef _lineintersect_utils_h
#define _lineintersect_utils_h

#include "SETTINGS.h"

void IntersectLineSegments(const VECTOR3& a0, const VECTOR3& a1,
                           const VECTOR3& b0, const VECTOR3& b1,
                           VECTOR3& aPoint, VECTOR3& bPoint);

void IntersectLineSegments(const REAL A1x, const REAL A1y, const REAL A1z,
                           const REAL A2x, const REAL A2y, const REAL A2z,
                           const REAL B1x, const REAL B1y, const REAL B1z,
                           const REAL B2x, const REAL B2y, const REAL B2z,
                           bool infinite_lines, REAL epsilon, REAL &PointOnSegAx,
                           REAL &PointOnSegAy, REAL &PointOnSegAz, REAL &PointOnSegBx,
                           REAL &PointOnSegBy, REAL &PointOnSegBz, REAL &NearestPointX,
                           REAL &NearestPointY, REAL &NearestPointZ, REAL &NearestVectorX,
                           REAL &NearestVectorY, REAL &NearestVectorZ, bool &true_intersection);

void FindNearestPointOnLineSegment(const REAL A1x, const REAL A1y, const REAL A1z,
                                   const REAL Lx, const REAL Ly, const REAL Lz,
                                   const REAL Bx, const REAL By, const REAL Bz,
                                   bool infinite_line, REAL epsilon_squared, REAL &NearestPointX,
                                   REAL &NearestPointY, REAL &NearestPointZ,
                                   REAL &parameter);

void FindNearestPointOfParallelLineSegments(REAL A1x, REAL A1y, REAL A1z,
                                            REAL A2x, REAL A2y, REAL A2z,
                                            REAL Lax, REAL Lay, REAL Laz,
                                            REAL B1x, REAL B1y, REAL B1z,
                                            REAL B2x, REAL B2y, REAL B2z,
                                            REAL Lbx, REAL Lby, REAL Lbz,
                                            bool infinite_lines, REAL epsilon_squared,
                                            REAL &PointOnSegAx, REAL &PointOnSegAy, REAL &PointOnSegAz,
                                            REAL &PointOnSegBx, REAL &PointOnSegBy, REAL &PointOnSegBz);

void AdjustNearestPoints(REAL A1x, REAL A1y, REAL A1z,
                         REAL Lax, REAL Lay, REAL Laz,
                         REAL B1x, REAL B1y, REAL B1z,
                         REAL Lbx, REAL Lby, REAL Lbz,
                         REAL epsilon_squared, REAL s, REAL t,
                         REAL &PointOnSegAx, REAL &PointOnSegAy, REAL &PointOnSegAz,
                         REAL &PointOnSegBx, REAL &PointOnSegBy, REAL &PointOnSegBz);

#endif // _lineintersect_utils_h
