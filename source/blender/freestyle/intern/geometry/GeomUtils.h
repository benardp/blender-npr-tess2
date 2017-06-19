/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __GEOMUTILS_H__
#define __GEOMUTILS_H__

/** \file blender/freestyle/intern/geometry/GeomUtils.h
 *  \ingroup freestyle
 *  \brief Various tools for geometry
 *  \author Stephane Grabli
 *  \date 12/04/2002
 */

#include <vector>

#include "Geom.h"

#include "../system/FreestyleConfig.h"

using namespace std;

namespace Freestyle {

using namespace Geometry;

namespace GeomUtils {

//
// Templated procedures
//
/////////////////////////////////////////////////////////////////////////////

/*! Computes the distance from a point P to a segment AB */
template<class T>
real distPointSegment(const T& P, const T& A, const T& B)
{
    T AB, AP, BP;
    AB = B - A;
    AP = P - A;
    BP = P - B;

    real c1(AB * AP);
    if (c1 <= 0)
        return AP.norm();

    real c2(AB * AB);
    if (c2 <= c1)
        return BP.norm();

    real b = c1 / c2;
    T Pb, PPb;
    Pb = A + b * AB;
    PPb = P - Pb;

    return PPb.norm();
}

//
// Non-templated procedures
//
/////////////////////////////////////////////////////////////////////////////
typedef enum {
    DONT_INTERSECT,
    DO_INTERSECT,
    COLINEAR,
    COINCIDENT,
} intersection_test;

intersection_test intersect2dSeg2dSeg(const Vec2r& p1, const Vec2r& p2, // first segment
                                      const Vec2r& p3, const Vec2r& p4, // second segment
                                      Vec2r& res);                      // found intersection point

intersection_test intersect2dLine2dLine(const Vec2r& p1, const Vec2r& p2, // first segment
                                        const Vec2r& p3, const Vec2r& p4, // second segment
                                        Vec2r& res);                      // found intersection point

intersection_test intersect2dSeg2dSegParametric(const Vec2r& p1, const Vec2r& p2, // first segment
                                                const Vec2r& p3, const Vec2r& p4, // second segment
                                                real& t,                          // I = P1 + t * P1P2)
                                                real& u,                          // I = P3 + u * P3P4
                                                real epsilon = M_EPSILON);

/*! check whether a 2D segment intersect a 2D region or not */
bool intersect2dSeg2dArea(const Vec2r& min, const Vec2r& max, const Vec2r& A, const Vec2r& B);

/*! check whether a 2D segment is included in a 2D region or not */
bool include2dSeg2dArea(const Vec2r& min, const Vec2r& max, const Vec2r& A, const Vec2r& B);

/*! Box-triangle overlap test, adapted from Tomas Akenine-Möller code */
bool overlapTriangleBox(Vec3r& boxcenter, Vec3r& boxhalfsize, Vec3r triverts[3]);

/*! Fast, Minimum Storage Ray-Triangle Intersection, adapted from Tomas Möller and Ben Trumbore code. */
bool intersectRayTriangle(const Vec3r& orig, const Vec3r& dir, const Vec3r& v0, const Vec3r& v1, const Vec3r& v2,
                          real& t,                         // I = orig + t * dir
                          real& u, real& v,                // I = (1 - u - v) * v0 + u * v1 + v * v2
                          const real epsilon = M_EPSILON); // the epsilon to use

/*! Intersection between plane and ray adapted from Graphics Gems, Didier Badouel */
intersection_test intersectRayPlane(const Vec3r& orig, const Vec3r& dir, // ray origin and direction
                                    // plane's normal and offset (plane = { P / P.N + d = 0 })
                                    const Vec3r& norm, const real d,
                                    real& t,                         // I = orig + t * dir
                                    const real epsilon = M_EPSILON); // the epsilon to use

/*! Intersection Ray-Bounding box (axis aligned).
 *  Adapted from Williams et al, "An Efficient Robust Ray-Box Intersection Algorithm", JGT 10:1 (2005), pp. 49-54.
 */
bool intersectRayBBox(const Vec3r& orig, const Vec3r& dir,      // ray origin and direction
                      const Vec3r& boxMin, const Vec3r& boxMax, // the bbox
                      // the interval in which at least on of the intersections must happen
                      real t0, real t1,
                      real& tmin,                               // Imin = orig + tmin * dir is the first intersection
                      real& tmax,                               // Imax = orig + tmax * dir is the second intersection
                      real epsilon = M_EPSILON);                // the epsilon to use

/*! Checks whether 3D point P lies inside or outside of the triangle ABC */
bool includePointTriangle(const Vec3r& P, const Vec3r& A, const Vec3r& B, const Vec3r& C);

void transformVertex(const Vec3r& vert, const Matrix44r& matrix, Vec3r& res);

void transformVertices(const vector<Vec3r>& vertices, const Matrix44r& trans, vector<Vec3r>& res);

Vec3r rotateVector(const Matrix44r& mat, const Vec3r& v);

//
// Coordinates systems changing procedures
//
/////////////////////////////////////////////////////////////////////////////

/*! From world to image
 *  p
 *    point's coordinates expressed in world coordinates system
 *  q
 *    vector in which the result will be stored
 *  model_view_matrix
 *    The model view matrix expressed in line major order (OpenGL
 *    matrices are column major ordered)
 *  projection_matrix
 *    The projection matrix expressed in line major order (OpenGL
 *    matrices are column major ordered)
 *  viewport
 *    The viewport: x,y coordinates followed by width and height (OpenGL like viewport)
 */
void fromWorldToImage(const Vec3r& p, Vec3r& q, const real model_view_matrix[4][4], const real projection_matrix[4][4],
const int viewport[4]);

/*! From world to image
 *  p
 *    point's coordinates expressed in world coordinates system
 *  q
 *    vector in which the result will be stored
 *  transform
 *    The transformation matrix (gathering model view and projection),
 *    expressed in line major order (OpenGL matrices are column major ordered)
 *  viewport
 *    The viewport: x,y coordinates followed by width and height (OpenGL like viewport)
 */
void fromWorldToImage(const Vec3r& p, Vec3r& q, const real transform[4][4], const int viewport[4]);

/*! Projects from world coordinates to camera coordinates 
 *  Returns the point's coordinates expressed in the camera's
 *  coordinates system.
 *  p
 *    point's coordinates expressed in world coordinates system
 *  q
 *    vector in which the result will be stored
 *  model_view_matrix
 *    The model view matrix expressed in line major order (OpenGL
 *    matrices are column major ordered)
 */
void fromWorldToCamera(const Vec3r& p, Vec3r& q, const real model_view_matrix[4][4]);

/*! Projects from World Coordinates to retina coordinates
 *  Returns the point's coordinates expressed in Retina system.
 *  p
 *    point's coordinates expressed in camera system
 *  q
 *    vector in which the result will be stored
 *  projection_matrix
 *    The projection matrix expressed in line major order (OpenGL
 *    matrices are column major ordered)
 */
void fromCameraToRetina(const Vec3r& p, Vec3r& q, const real projection_matrix[4][4]);

/*! From retina to image.
 *  Returns the coordinates expressed in Image coordinates system.
 *  p
 *    point's coordinates expressed in retina system
 *  q
 *    vector in which the result will be stored
 *  viewport
 *    The viewport: x,y coordinates followed by width and height (OpenGL like viewport).
 */
void fromRetinaToImage(const Vec3r& p, Vec3r& q, const int viewport[4]);

/*! From image to retina
 *  p
 *    point's coordinates expressed in image system
 *  q
 *    vector in which the result will be stored
 *  viewport
 *    The viewport: x,y coordinates followed by width and height (OpenGL like viewport).
 */
void fromImageToRetina(const Vec3r& p, Vec3r& q, const int viewport[4]);

/*! computes the coordinates of q in the camera coordinates system, 
 *  using the known z coordinates of the 3D point.
 *  That means that this method does not inverse any matrices,
 *  it only computes X and Y from x,y and Z)
 *  p
 *    point's coordinates expressed in retina system
 *  q
 *    vector in which the result will be stored
 *  projection_matrix
 *    The projection matrix expressed in line major order (OpenGL
 *    matrices are column major ordered)
 */
void fromRetinaToCamera(const Vec3r& p, Vec3r& q, real z, const real projection_matrix[4][4]);

/*! Projects from camera coordinates to world coordinates
 *  Returns the point's coordinates expressed in the world's
 *  coordinates system.
 *  p
 *    point's coordinates expressed in the camera coordinates system
 *  q
 *    vector in which the result will be stored
 *  model_view_matrix
 *    The model view matrix expressed in line major order (OpenGL
 *    matrices are column major ordered)
 */
void fromCameraToWorld(const Vec3r& p, Vec3r& q, const real model_view_matrix[4][4]);


inline Vec3r Bary2Point(const Vec3r A,const Vec3r B,const Vec3r C,const Vec3r P)
{
    return A*P.x() + B*P.y() + C*P.z();
}

inline real segmentParam(const Vec2r A, const Vec2r B, const Vec2r P)
{
    Vec2r PA = P-A;
    Vec2r BA = B-A;
    return (PA * BA) / (BA*BA);
}

inline real segmentParam(const Vec3r A, const Vec3r B, const Vec3r P)
{
    Vec3r PA = P-A;
    Vec3r BA = B-A;
    return (PA * BA) / (BA*BA);
}

/*
*
*  Three-dimensional Triangle-Triangle Intersection
*
*/

#define CROSS(dest,v1,v2) \
    dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
    dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
    dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; \
    dest[1]=v1[1]-v2[1]; \
    dest[2]=v1[2]-v2[2];

#define SCALAR(dest,alpha,v) dest[0] = alpha * v[0]; \
    dest[1] = alpha * v[1]; \
    dest[2] = alpha * v[2];

/*
   This macro is called when the triangles surely intersect
   It constructs the segment of intersection of the two triangles
   if they are not coplanar.
*/

#define CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) { \
    SUB(v1,q1,p1) \
    SUB(v2,r2,p1) \
    CROSS(N,v1,v2) \
    SUB(v,p2,p1) \
    if (DOT(v,N) > 0.0f) {\
    SUB(v1,r1,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) <= 0.0f) { \
    SUB(v2,q2,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) > 0.0f) { \
    SUB(v1,p1,p2) \
    SUB(v2,p1,r1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p1,v1) \
    SUB(v1,p2,p1) \
    SUB(v2,p2,r2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p2,v1) \
    return 1; \
} else { \
    SUB(v1,p2,p1) \
    SUB(v2,p2,q2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p2,v1) \
    SUB(v1,p2,p1) \
    SUB(v2,p2,r2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p2,v1) \
    return 1; \
} \
} else { \
    return 0; \
} \
} else { \
    SUB(v2,q2,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) < 0.0f) { \
    return 0; \
} else { \
    SUB(v1,r1,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) >= 0.0f) { \
    SUB(v1,p1,p2) \
    SUB(v2,p1,r1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p1,v1) \
    SUB(v1,p1,p2) \
    SUB(v2,p1,q1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p1,v1) \
    return 1; \
} else { \
    SUB(v1,p2,p1) \
    SUB(v2,p2,q2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p2,v1) \
    SUB(v1,p1,p2) \
    SUB(v2,p1,q1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p1,v1) \
    return 1; \
}}}}

#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))

#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
    if (ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
    if (ORIENT_2D(P1,P2,Q1) > 0.0f) {\
    if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1; \
    else return 0;} else {\
    if (ORIENT_2D(P1,P2,R1) >= 0.0f)\
    if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1; \
    else return 0;\
    else return 0;}\
    else \
    if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
    if (ORIENT_2D(R2,Q2,R1) <= 0.0f)\
    if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1; \
    else return 0;\
    else return 0;\
    else return 0;\
    else\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) \
    if (ORIENT_2D(Q1,R1,R2) >= 0.0f)\
    if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;\
    else return 0;\
    else \
    if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
    if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1; \
    else return 0; }\
    else return 0; \
    else  return 0; \
};

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
    if (ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
    if (ORIENT_2D(P1,Q1,R2) >= 0.0f) return 1; \
    else return 0;} else { \
    if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
    if (ORIENT_2D(R1,P1,P2) >= 0.0f) return 1; else return 0;} \
    else return 0; } \
} else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,R1) >= 0.0f) {\
    if (ORIENT_2D(P1,R1,R2) >= 0.0f) return 1;  \
    else {\
    if (ORIENT_2D(Q1,R1,R2) >= 0.0f) return 1; else return 0;}}\
    else  return 0; }\
    else return 0; }}


inline int ccw_tri_tri_intersection_2d(real p1[2], real q1[2], real r1[2],
real p2[2], real q2[2], real r2[2]) {
    if ( ORIENT_2D(p2,q2,p1) >= 0.0f ) {
        if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
            if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) return 1;
            else INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
        } else {
            if ( ORIENT_2D(r2,p2,p1) >= 0.0f )
                INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
                        else INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
    else {
        if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
            if ( ORIENT_2D(r2,p2,p1) >= 0.0f )
                INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
                        else  INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
        else INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
}

inline int tri_tri_overlap_test_2d(real p1[2], real q1[2], real r1[2],
real p2[2], real q2[2], real r2[2]) {
    if ( ORIENT_2D(p1,q1,r1) < 0.0f )
        if ( ORIENT_2D(p2,q2,r2) < 0.0f )
            return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
        else
            return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
    else
        if ( ORIENT_2D(p2,q2,r2) < 0.0f )
            return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
        else
            return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);
}

inline int coplanar_tri_tri3d(real p1[3], real q1[3], real r1[3],
real p2[3], real q2[3], real r2[3],
real normal_1[3], real normal_2[3]){

    real P1[2],Q1[2],R1[2];
    real P2[2],Q2[2],R2[2];

    real n_x, n_y, n_z;

    n_x = ((normal_1[0]<0)?-normal_1[0]:normal_1[0]);
    n_y = ((normal_1[1]<0)?-normal_1[1]:normal_1[1]);
    n_z = ((normal_1[2]<0)?-normal_1[2]:normal_1[2]);


    /* Projection of the triangles in 3D onto 2D such that the area of
     the projection is maximized. */


    if (( n_x > n_z ) && ( n_x >= n_y )) {
        // Project onto plane YZ

        P1[0] = q1[2]; P1[1] = q1[1];
        Q1[0] = p1[2]; Q1[1] = p1[1];
        R1[0] = r1[2]; R1[1] = r1[1];

        P2[0] = q2[2]; P2[1] = q2[1];
        Q2[0] = p2[2]; Q2[1] = p2[1];
        R2[0] = r2[2]; R2[1] = r2[1];

    } else if (( n_y > n_z ) && ( n_y >= n_x )) {
        // Project onto plane XZ

        P1[0] = q1[0]; P1[1] = q1[2];
        Q1[0] = p1[0]; Q1[1] = p1[2];
        R1[0] = r1[0]; R1[1] = r1[2];

        P2[0] = q2[0]; P2[1] = q2[2];
        Q2[0] = p2[0]; Q2[1] = p2[2];
        R2[0] = r2[0]; R2[1] = r2[2];

    } else {
        // Project onto plane XY

        P1[0] = p1[0]; P1[1] = p1[1];
        Q1[0] = q1[0]; Q1[1] = q1[1];
        R1[0] = r1[0]; R1[1] = r1[1];

        P2[0] = p2[0]; P2[1] = p2[1];
        Q2[0] = q2[0]; Q2[1] = q2[1];
        R2[0] = r2[0]; R2[1] = r2[1];
    }

    return tri_tri_overlap_test_2d(P1,Q1,R1,P2,Q2,R2);
}

#define TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
    if (dp2 > 0.0f) { \
    if (dq2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2) \
    else if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) }\
    else if (dp2 < 0.0f) { \
    if (dq2 < 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
} else { \
    if (dq2 < 0.0f) { \
    if (dr2 >= 0.0f)  CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2)\
} \
    else if (dq2 > 0.0f) { \
    if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
    else  CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
} \
    else  { \
    if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2)\
    else { \
    *coplanar = 1; \
    return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);\
} \
}} }


/*
   The following version computes the segment of intersection of the
   two triangles if it exists.
   coplanar returns whether the triangles are coplanar
   source and target are the endpoints of the line segment of intersection
*/

inline int tri_tri_intersection_test_3d(real p1[3], real q1[3], real r1[3],
real p2[3], real q2[3], real r2[3],
int * coplanar, real source[3], real target[3])
{
    real dp1, dq1, dr1, dp2, dq2, dr2;
    real v1[3], v2[3], v[3];
    real N1[3], N2[3], N[3];
    real alpha;

    // Compute distance signs  of p1, q1 and r1
    // to the plane of triangle(p2,q2,r2)


    SUB(v1,p2,r2)
            SUB(v2,q2,r2)
            CROSS(N2,v1,v2)

            SUB(v1,p1,r2)
            dp1 = DOT(v1,N2);
    SUB(v1,q1,r2)
            dq1 = DOT(v1,N2);
    SUB(v1,r1,r2)
            dr1 = DOT(v1,N2);

    if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0;

    // Compute distance signs  of p2, q2 and r2
    // to the plane of triangle(p1,q1,r1)


    SUB(v1,q1,p1)
            SUB(v2,r1,p1)
            CROSS(N1,v1,v2)

            SUB(v1,p2,r1)
            dp2 = DOT(v1,N1);
    SUB(v1,q2,r1)
            dq2 = DOT(v1,N1);
    SUB(v1,r2,r1)
            dr2 = DOT(v1,N1);

    if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

    // Permutation in a canonical form of T1's vertices

    if (dp1 > 0.0f) {
        if (dq1 > 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
                else if (dr1 > 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)

                else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    } else if (dp1 < 0.0f) {
        if (dq1 < 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
                else if (dr1 < 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
                else TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
    } else {
        if (dq1 < 0.0f) {
            if (dr1 >= 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
                    else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
        }
        else if (dq1 > 0.0f) {
            if (dr1 > 0.0f) TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
                    else TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
        }
        else  {
            if (dr1 > 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
                    else if (dr1 < 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
                    else {
                // triangles are co-planar

                *coplanar = 1;
                return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);
            }
        }
    }
}

} // end of namespace GeomUtils

} /* namespace Freestyle */

#endif // __GEOMUTILS_H__
