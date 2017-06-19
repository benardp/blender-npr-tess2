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

/** \file blender/freestyle/intern/view_map/ViewMapBuilder.cpp
 *  \ingroup freestyle
 *  \brief Class to build silhouette edges from a Winged-Edge structure
 *  \author Stephane Grabli
 *  \date 25/03/2002
 */

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <sstream>

#include "FRS_freestyle.h"

#include "BoxGrid.h"
#include "CulledOccluderSource.h"
#include "HeuristicGridDensityProviderFactory.h"
#include "OccluderSource.h"
#include "SphericalGrid.h"
#include "ViewMapBuilder.h"

#include "../geometry/GridHelpers.h"
#include "../geometry/GeomUtils.h"

#include "../winged_edge/WFillGrid.h"

#include "carve/exact.hpp"

#include "BKE_global.h"

namespace Freestyle {

// XXX Grmll... G is used as template's typename parameter :/
static const Global &_global = G;

#define LOGGING 0

using namespace std;

template <typename G, typename I>
static void findOccludee(FEdge *fe, G& /*grid*/, I& occluders, real epsilon, WFace **oaWFace,
                         Vec3r& u, Vec3r& A, Vec3r& origin, Vec3r& edge, vector<WVertex*>& faceVertices)
{
    WFace *face = NULL;
    if (fe->isSmooth()) {
        FEdgeSmooth *fes = dynamic_cast<FEdgeSmooth*>(fe);
        face = (WFace *)fes->face();
    }
    WFace *oface;
    bool skipFace;

    WVertex::incoming_edge_iterator ie;

    *oaWFace = NULL;
    if (((fe)->getNature() & Nature::SILHOUETTE) || ((fe)->getNature() & Nature::BORDER)) {
        // we cast a ray from A in the same direction but looking behind
        Vec3r v(-u[0], -u[1], -u[2]);
        bool noIntersection = true;
        real mint = FLT_MAX;

        for (occluders.initAfterTarget(); occluders.validAfterTarget(); occluders.nextOccludee()) {
#if LOGGING
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "\t\tEvaluating intersection for occludee " << occluders.getWFace() << " and ray " << A <<
                        " * " << u << endl;
            }
#endif
            oface = occluders.getWFace();
            Polygon3r *p = occluders.getCameraSpacePolygon();
            real d = -((p->getVertices())[0] * p->getNormal());
            real t, t_u, t_v;

            if (0 != face) {
                skipFace = false;

                if (face == oface)
                    continue;

                if (faceVertices.empty())
                    continue;

                for (vector<WVertex*>::iterator fv = faceVertices.begin(), fvend = faceVertices.end();
                     fv != fvend;
                     ++fv)
                {
                    if ((*fv)->isBoundary())
                        continue;
                    WVertex::incoming_edge_iterator iebegin = (*fv)->incoming_edges_begin();
                    WVertex::incoming_edge_iterator ieend = (*fv)->incoming_edges_end();
                    for (ie = iebegin; ie != ieend; ++ie) {
                        if ((*ie) == 0)
                            continue;

                        WFace *sface = (*ie)->GetbFace();
                        if (sface == oface) {
                            skipFace = true;
                            break;
                        }
                    }
                    if (skipFace)
                        break;
                }
                if (skipFace)
                    continue;
            }
            else {
                // check whether the edge and the polygon plane are coincident:
                //-------------------------------------------------------------
                //first let us compute the plane equation.
                if (GeomUtils::COINCIDENT == GeomUtils::intersectRayPlane(origin, edge, p->getNormal(), d, t, epsilon))
                {
#if LOGGING
                    if (_global.debug & G_DEBUG_FREESTYLE) {
                        cout << "\t\tRejecting occluder for target coincidence." << endl;
                    }
#endif
                    continue;
                }
            }

            if (p->rayIntersect(A, v, t, t_u, t_v)) {
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\t\tRay " << A << " * " << v << " intersects at time " << t << endl;
                    cout << "\t\t(v * normal) == " << (v * p->getNormal()) << " for normal " << p->getNormal() << endl;
                }
#endif
                if (fabs(v * p->getNormal()) > 0.0001) {
                    if ((t > 0.0)) { // && (t<1.0))
                        if (t < mint) {
                            *oaWFace = oface;
                            mint = t;
                            noIntersection = false;
                            fe->setOccludeeIntersection(Vec3r(A + t * v));
#if LOGGING
                            if (_global.debug & G_DEBUG_FREESTYLE) {
                                cout << "\t\tIs occludee" << endl;
                            }
#endif
                        }
                    }
                }

                occluders.reportDepth(A, v, t);
            }
        }

        if (noIntersection)
            *oaWFace = NULL;
    }
}

template <typename G, typename I>
static void findOccludee(FEdge *fe, G& grid, real epsilon, ViewEdge * /*ve*/, WFace **oaFace)
{
    Vec3r A;
    Vec3r edge;
    Vec3r origin;
    A = Vec3r(((fe)->vertexA()->point3D() + (fe)->vertexB()->point3D()) / 2.0);
    edge = Vec3r((fe)->vertexB()->point3D() - (fe)->vertexA()->point3D());
    origin = Vec3r((fe)->vertexA()->point3D());
    Vec3r u;
    if (grid.orthographicProjection()) {
        u = Vec3r(0.0, 0.0, grid.viewpoint().z() - A.z());
    }
    else {
        u = Vec3r(grid.viewpoint() - A);
    }
    u.normalize();

    vector<WVertex*> faceVertices;

    WFace *face = NULL;
    if (fe->isSmooth()) {
        FEdgeSmooth *fes = dynamic_cast<FEdgeSmooth*>(fe);
        face = (WFace *)fes->face();
    }

    if (face) {
        face->RetrieveVertexList(faceVertices);
    }

    I occluders(grid, A, epsilon);
    findOccludee<G, I>(fe, grid, occluders, epsilon, oaFace, u, A, origin, edge, faceVertices);
}

// computeVisibility takes a pointer to foundOccluders, instead of using a reference,
// so that computeVeryFastVisibility can skip the AddOccluders step with minimal overhead.
template <typename G, typename I>
static int computeVisibility(ViewMap *viewMap, FEdge *fe, G& grid, real epsilon, ViewEdge * /*ve*/, WFace **oaWFace,
                             set<ViewShape*> *foundOccluders)
{
    int qi = 0;

    Vec3r center;
    Vec3r edge;
    Vec3r origin;

    center = fe->center3d();
    edge = Vec3r(fe->vertexB()->point3D() - fe->vertexA()->point3D());
    origin = Vec3r(fe->vertexA()->point3D());

    Vec3r vp;
    if (grid.orthographicProjection()) {
        vp = Vec3r(center.x(), center.y(), grid.viewpoint().z());
    }
    else {
        vp = Vec3r(grid.viewpoint());
    }
    Vec3r u(vp - center);
    real raylength = u.norm();
    u.normalize();

    WFace *face = NULL;
    if (fe->isSmooth()) {
        FEdgeSmooth *fes = dynamic_cast<FEdgeSmooth*>(fe);
        face = (WFace *)fes->face();
    }
    vector<WVertex*> faceVertices;
    WVertex::incoming_edge_iterator ie;

    WFace *oface;
    bool skipFace;

    if (face)
        face->RetrieveVertexList(faceVertices);

    I occluders(grid, center, epsilon);

    for (occluders.initBeforeTarget(); occluders.validBeforeTarget(); occluders.nextOccluder()) {
        // If we're dealing with an exact silhouette, check whether we must take care of this occluder of not.
        // (Indeed, we don't consider the occluders that share at least one vertex with the face containing this edge).
        //-----------
        oface = occluders.getWFace();
        Polygon3r *p = occluders.getCameraSpacePolygon();
        real t, t_u, t_v;
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\t\tEvaluating intersection for occluder " << (p->getVertices())[0] << (p->getVertices())[1] <<
                    (p->getVertices())[2] << endl << "\t\t\tand ray " << vp << " * " << u << " (center " << center <<
                    ")" << endl;
        }
#endif

#if LOGGING
        Vec3r v(vp - center);
        real rl = v.norm();
        v.normalize();
        vector<Vec3r> points;
        // Iterate over vertices, storing projections in points
        for (vector<WOEdge*>::const_iterator woe = oface->getEdgeList().begin(), woend = oface->getEdgeList().end();
             woe != woend;
             woe++)
        {
            points.push_back(Vec3r((*woe)->GetaVertex()->GetVertex()));
        }
        Polygon3r p1(points, oface->GetNormal());
        Vec3r v1((p1.getVertices())[0]);
        real d = -(v1 * p->getNormal());
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\t\tp:  " << (p->getVertices())[0] << (p->getVertices())[1] << (p->getVertices())[2] << ", norm: " <<
                    p->getNormal() << endl;
            cout << "\t\tp1: " << (p1.getVertices())[0] << (p1.getVertices())[1] << (p1.getVertices())[2] << ", norm: " <<
                    p1.getNormal() << endl;
        }
#else
        real d = -((p->getVertices())[0] * p->getNormal());
#endif

        if (face) {
#if LOGGING
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "\t\tDetermining face adjacency...";
            }
#endif
            skipFace = false;

            if (face == oface) {
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "  Rejecting occluder for face concurrency." << endl;
                }
#endif
                continue;
            }


            for (vector<WVertex*>::iterator fv = faceVertices.begin(), fvend = faceVertices.end(); fv != fvend; ++fv) {
                if ((*fv)->isBoundary())
                    continue;

                WVertex::incoming_edge_iterator iebegin = (*fv)->incoming_edges_begin();
                WVertex::incoming_edge_iterator ieend = (*fv)->incoming_edges_end();
                for (ie = iebegin; ie != ieend; ++ie) {
                    if ((*ie) == 0)
                        continue;

                    WFace *sface = (*ie)->GetbFace();
                    //WFace *sfacea = (*ie)->GetaFace();
                    //if ((sface == oface) || (sfacea == oface))
                    if (sface == oface) {
                        skipFace = true;
                        break;
                    }
                }
                if (skipFace)
                    break;
            }
            if (skipFace) {
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "  Rejecting occluder for face adjacency." << endl;
                }
#endif
                continue;
            }
        }
        else {
            // check whether the edge and the polygon plane are coincident:
            //-------------------------------------------------------------
            //first let us compute the plane equation.
            if (GeomUtils::COINCIDENT == GeomUtils::intersectRayPlane(origin, edge, p->getNormal(), d, t, epsilon)) {
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\t\tRejecting occluder for target coincidence." << endl;
                }
#endif
                continue;
            }
        }

#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            real x;
            if (p1.rayIntersect(center, v, x, t_u, t_v)) {
                cout << "\t\tRay should intersect at time " << (rl - x) << ". Center: " << center << ", V: " << v <<
                        ", RL: " << rl << ", T:" << x << endl;
            }
            else {
                cout << "\t\tRay should not intersect. Center: " << center << ", V: " << v <<  ", RL: " << rl << endl;
            }
        }
#endif

        if (p->rayIntersect(center, u, t, t_u, t_v)) {
#if LOGGING
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "\t\tRay " << center << " * " << u << " intersects at time " << t << " (raylength is " <<
                        raylength << ")" << endl;
                cout << "\t\t(u * normal) == " << (u * p->getNormal()) << " for normal " << p->getNormal() << endl;
            }
#endif
            if (fabs(u * p->getNormal()) > 0.0001) {
                if ((t > 0.0) && (t < raylength)) {
#if LOGGING
                    if (_global.debug & G_DEBUG_FREESTYLE) {
                        cout << "\t\tIs occluder" << endl;
                    }
#endif
                    if ( foundOccluders != NULL ) {
                        ViewShape *vshape = viewMap->viewShape(oface->GetVertex(0)->shape()->GetId());
                        foundOccluders->insert(vshape);
                    }
                    ++qi;

                    if (! grid.enableQI())
                        break;
                }

                occluders.reportDepth(center, u, t);
            }
        }
    }

    // Find occludee
    findOccludee<G, I>(fe, grid, occluders, epsilon, oaWFace, u, center, origin, edge, faceVertices);

    return qi;
}

// computeCumulativeVisibility returns the lowest x such that the majority of FEdges have QI <= x
//
// This was probably the original intention of the "normal" algorithm on which computeDetailedVisibility is based.
// But because the "normal" algorithm chooses the most popular QI, without considering any other values, a ViewEdge
// with FEdges having QIs of 0, 21, 22, 23, 24 and 25 will end up having a total QI of 0, even though most of the
// FEdges are heavily occluded. computeCumulativeVisibility will treat this case as a QI of 22 because 3 out of
// 6 occluders have QI <= 22.

template <typename G, typename I>
static void computeCumulativeVisibility(ViewMap *ioViewMap, G& grid, real epsilon, RenderMonitor *iRenderMonitor)
{
    vector<ViewEdge*>& vedges = ioViewMap->ViewEdges();

    FEdge *fe, *festart;
    int nSamples = 0;
    vector<WFace*> wFaces;
    WFace *wFace = NULL;
    unsigned cnt = 0;
    unsigned cntStep = (unsigned)ceil(0.01f * vedges.size());
    unsigned tmpQI = 0;
    unsigned qiClasses[256];
    unsigned maxIndex, maxCard;
    unsigned qiMajority;
    for (vector<ViewEdge*>::iterator ve = vedges.begin(), veend = vedges.end(); ve != veend; ve++) {
        if (iRenderMonitor) {
            if (iRenderMonitor->testBreak())
                break;
            if (cnt % cntStep == 0) {
                stringstream ss;
                ss << "Freestyle: Visibility computations " << (100 * cnt / vedges.size()) << "%";
                iRenderMonitor->setInfo(ss.str());
                iRenderMonitor->progress((float)cnt / vedges.size());
            }
            cnt++;
        }
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Processing ViewEdge " << (*ve)->getId() << endl;
        }
#endif
        // Find an edge to test
        if (!(*ve)->isInImage()) {
            // This view edge has been proscenium culled
            (*ve)->setQI(255);
            (*ve)->setaShape(0);
#if LOGGING
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "\tCulled." << endl;
            }
#endif
            continue;
        }

        // Test edge
        festart = (*ve)->fedgeA();
        fe = (*ve)->fedgeA();
        qiMajority = 0;
        do {
            if (fe != NULL && fe->isInImage()) {
                qiMajority++;
            }
            fe = fe->nextEdge();
        } while (fe && fe != festart);

        if (qiMajority == 0) {
            // There are no occludable FEdges on this ViewEdge
            // This should be impossible.
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "View Edge in viewport without occludable FEdges: " << (*ve)->getId() << endl;
            }
            // We can recover from this error:
            // Treat this edge as fully visible with no occludee
            (*ve)->setQI(0);
            (*ve)->setaShape(0);
            continue;
        }
        else {
            ++qiMajority;
            qiMajority >>= 1;
        }
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tqiMajority: " << qiMajority << endl;
        }
#endif

        tmpQI = 0;
        maxIndex = 0;
        maxCard = 0;
        nSamples = 0;
        memset(qiClasses, 0, 256 * sizeof(*qiClasses));
        set<ViewShape*> foundOccluders;

        fe = (*ve)->fedgeA();
        do {
            if (!fe || !fe->isInImage()) {
                fe = fe->nextEdge();
                continue;
            }
            if ((maxCard < qiMajority)) {
                //ARB: change &wFace to wFace and use reference in called function
                tmpQI = computeVisibility<G, I>(ioViewMap, fe, grid, epsilon, *ve, &wFace, &foundOccluders);
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFEdge: visibility " << tmpQI << endl;
                }
#endif

                //ARB: This is an error condition, not an alert condition.
                // Some sort of recovery or abort is necessary.
                if (tmpQI >= 256) {
                    cerr << "Warning: too many occluding levels" << endl;
                    //ARB: Wild guess: instead of aborting or corrupting memory, treat as tmpQI == 255
                    tmpQI = 255;
                }

                if (++qiClasses[tmpQI] > maxCard) {
                    maxCard = qiClasses[tmpQI];
                    maxIndex = tmpQI;
                }
            }
            else {
                //ARB: FindOccludee is redundant if ComputeRayCastingVisibility has been called
                //ARB: change &wFace to wFace and use reference in called function
                findOccludee<G, I>(fe, grid, epsilon, *ve, &wFace);
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFEdge: occludee only (" << (wFace != NULL ? "found" : "not found") << ")" << endl;
                }
#endif
            }

            // Store test results
            if (wFace) {
                vector<Vec3r> vertices;
                for (int i = 0, numEdges = wFace->numberOfEdges(); i < numEdges; ++i) {
                    vertices.push_back(Vec3r(wFace->GetVertex(i)->GetVertex()));
                }
                Polygon3r poly(vertices, wFace->GetNormal());
                poly.userdata = (void *)wFace;
                fe->setaFace(poly);
                wFaces.push_back(wFace);
                fe->setOccludeeEmpty(false);
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFound occludee" << endl;
                }
#endif
            }
            else {
                fe->setOccludeeEmpty(true);
            }

            ++nSamples;
            fe = fe->nextEdge();
        } while ((maxCard < qiMajority) && (fe) && (fe != festart));

#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tFinished with " << nSamples << " samples, maxCard = " << maxCard << endl;
        }
#endif

        // ViewEdge
        // qi --
        // Find the minimum value that is >= the majority of the QI
        for (unsigned count = 0, i = 0; i < 256; ++i) {
            count += qiClasses[i];
            if (count >= qiMajority) {
                (*ve)->setQI(i);
                break;
            }
        }
        // occluders --
        // I would rather not have to go through the effort of creating this set and then copying out its contents.
        // Is there a reason why ViewEdge::_Occluders cannot be converted to a set<>?
        for (set<ViewShape*>::iterator o = foundOccluders.begin(), oend = foundOccluders.end(); o != oend; ++o) {
            (*ve)->AddOccluder((*o));
        }
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tConclusion: QI = " << maxIndex << ", " << (*ve)->occluders_size() << " occluders." << endl;
        }
#else
        (void)maxIndex;
#endif
        // occludee --
        if (!wFaces.empty()) {
            if (wFaces.size() <= (float)nSamples / 2.0f) {
                (*ve)->setaShape(0);
            }
            else {
                ViewShape *vshape = ioViewMap->viewShape((*wFaces.begin())->GetVertex(0)->shape()->GetId());
                (*ve)->setaShape(vshape);
            }
        }

        wFaces.clear();
    }
    if (iRenderMonitor && vedges.size()) {
        stringstream ss;
        ss << "Freestyle: Visibility computations " << (100 * cnt / vedges.size()) << "%";
        iRenderMonitor->setInfo(ss.str());
        iRenderMonitor->progress((float)cnt / vedges.size());
    }
}

template <typename G, typename I>
static void computeDetailedVisibility(ViewMap *ioViewMap, G& grid, real epsilon, RenderMonitor *iRenderMonitor)
{
    vector<ViewEdge*>& vedges = ioViewMap->ViewEdges();

    FEdge *fe, *festart;
    int nSamples = 0;
    vector<WFace*> wFaces;
    WFace *wFace = NULL;
    unsigned tmpQI = 0;
    unsigned qiClasses[256];
    unsigned maxIndex, maxCard;
    unsigned qiMajority;
    for (vector<ViewEdge*>::iterator ve = vedges.begin(), veend = vedges.end(); ve != veend; ve++) {
        if (iRenderMonitor && iRenderMonitor->testBreak())
            break;
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Processing ViewEdge " << (*ve)->getId() << endl;
        }
#endif
        // Find an edge to test
        if (!(*ve)->isInImage()) {
            // This view edge has been proscenium culled
            (*ve)->setQI(255);
            (*ve)->setaShape(0);
#if LOGGING
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "\tCulled." << endl;
            }
#endif
            continue;
        }

        // Test edge
        festart = (*ve)->fedgeA();
        fe = (*ve)->fedgeA();
        qiMajority = 0;
        do {
            if (fe != NULL && fe->isInImage()) {
                qiMajority++;
            }
            fe = fe->nextEdge();
        } while (fe && fe != festart);

        if (qiMajority == 0) {
            // There are no occludable FEdges on this ViewEdge
            // This should be impossible.
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "View Edge in viewport without occludable FEdges: " << (*ve)->getId() << endl;
            }
            // We can recover from this error:
            // Treat this edge as fully visible with no occludee
            (*ve)->setQI(0);
            (*ve)->setaShape(0);
            continue;
        }
        else {
            ++qiMajority;
            qiMajority >>= 1;
        }
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tqiMajority: " << qiMajority << endl;
        }
#endif

        tmpQI = 0;
        maxIndex = 0;
        maxCard = 0;
        nSamples = 0;
        memset(qiClasses, 0, 256 * sizeof(*qiClasses));
        set<ViewShape*> foundOccluders;

        fe = (*ve)->fedgeA();
        do {
            if (fe == NULL || ! fe->isInImage()) {
                fe = fe->nextEdge();
                continue;
            }
            if ((maxCard < qiMajority)) {
                //ARB: change &wFace to wFace and use reference in called function
                tmpQI = computeVisibility<G, I>(ioViewMap, fe, grid, epsilon, *ve, &wFace, &foundOccluders);
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFEdge: visibility " << tmpQI << endl;
                }
#endif

                //ARB: This is an error condition, not an alert condition.
                // Some sort of recovery or abort is necessary.
                if (tmpQI >= 256) {
                    cerr << "Warning: too many occluding levels" << endl;
                    //ARB: Wild guess: instead of aborting or corrupting memory, treat as tmpQI == 255
                    tmpQI = 255;
                }

                if (++qiClasses[tmpQI] > maxCard) {
                    maxCard = qiClasses[tmpQI];
                    maxIndex = tmpQI;
                }
            }
            else {
                //ARB: FindOccludee is redundant if ComputeRayCastingVisibility has been called
                //ARB: change &wFace to wFace and use reference in called function
                findOccludee<G, I>(fe, grid, epsilon, *ve, &wFace);
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFEdge: occludee only (" << (wFace != NULL ? "found" : "not found") << ")" << endl;
                }
#endif
            }

            // Store test results
            if (wFace) {
                vector<Vec3r> vertices;
                for (int i = 0, numEdges = wFace->numberOfEdges(); i < numEdges; ++i) {
                    vertices.push_back(Vec3r(wFace->GetVertex(i)->GetVertex()));
                }
                Polygon3r poly(vertices, wFace->GetNormal());
                poly.userdata = (void *)wFace;
                fe->setaFace(poly);
                wFaces.push_back(wFace);
                fe->setOccludeeEmpty(false);
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFound occludee" << endl;
                }
#endif
            }
            else {
                fe->setOccludeeEmpty(true);
            }

            ++nSamples;
            fe = fe->nextEdge();
        } while ((maxCard < qiMajority) && (fe) && (fe != festart));

#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tFinished with " << nSamples << " samples, maxCard = " << maxCard << endl;
        }
#endif

        // ViewEdge
        // qi --
        (*ve)->setQI(maxIndex);
        // occluders --
        // I would rather not have to go through the effort of creating this this set and then copying out its contents.
        // Is there a reason why ViewEdge::_Occluders cannot be converted to a set<>?
        for (set<ViewShape*>::iterator o = foundOccluders.begin(), oend = foundOccluders.end(); o != oend; ++o) {
            (*ve)->AddOccluder((*o));
        }
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tConclusion: QI = " << maxIndex << ", " << (*ve)->occluders_size() << " occluders." << endl;
        }
#endif
        // occludee --
        if (!wFaces.empty()) {
            if (wFaces.size() <= (float)nSamples / 2.0f) {
                (*ve)->setaShape(0);
            }
            else {
                ViewShape *vshape = ioViewMap->viewShape((*wFaces.begin())->GetVertex(0)->shape()->GetId());
                (*ve)->setaShape(vshape);
            }
        }

        wFaces.clear();
    }
}

template <typename G, typename I>
static void computeFastVisibility(ViewMap *ioViewMap, G& grid, real epsilon)
{
    vector<ViewEdge*>& vedges = ioViewMap->ViewEdges();

    FEdge *fe, *festart;
    unsigned nSamples = 0;
    vector<WFace*> wFaces;
    WFace *wFace = NULL;
    unsigned tmpQI = 0;
    unsigned qiClasses[256];
    unsigned maxIndex, maxCard;
    unsigned qiMajority;
    bool even_test;
    for (vector<ViewEdge*>::iterator ve = vedges.begin(), veend = vedges.end(); ve != veend; ve++) {
        // Find an edge to test
        if (!(*ve)->isInImage()) {
            // This view edge has been proscenium culled
            (*ve)->setQI(255);
            (*ve)->setaShape(0);
            continue;
        }

        // Test edge
        festart = (*ve)->fedgeA();
        fe = (*ve)->fedgeA();

        even_test = true;
        qiMajority = 0;
        do {
            if (even_test && fe && fe->isInImage()) {
                qiMajority++;
                even_test = !even_test;
            }
            fe = fe->nextEdge();
        } while (fe && fe != festart);

        if (qiMajority == 0 ) {
            // There are no occludable FEdges on this ViewEdge
            // This should be impossible.
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "View Edge in viewport without occludable FEdges: " << (*ve)->getId() << endl;
            }
            // We can recover from this error:
            // Treat this edge as fully visible with no occludee
            (*ve)->setQI(0);
            (*ve)->setaShape(0);
            continue;
        }
        else {
            ++qiMajority;
            qiMajority >>= 1;
        }

        even_test = true;
        maxIndex = 0;
        maxCard = 0;
        nSamples = 0;
        memset(qiClasses, 0, 256 * sizeof(*qiClasses));
        set<ViewShape*> foundOccluders;

        fe = (*ve)->fedgeA();
        do {
            if (!fe || !fe->isInImage()) {
                fe = fe->nextEdge();
                continue;
            }
            if (even_test) {
                if ((maxCard < qiMajority)) {
                    //ARB: change &wFace to wFace and use reference in called function
                    tmpQI = computeVisibility<G, I>(ioViewMap, fe, grid, epsilon, *ve, &wFace, &foundOccluders);

                    //ARB: This is an error condition, not an alert condition.
                    // Some sort of recovery or abort is necessary.
                    if (tmpQI >= 256) {
                        cerr << "Warning: too many occluding levels" << endl;
                        //ARB: Wild guess: instead of aborting or corrupting memory, treat as tmpQI == 255
                        tmpQI = 255;
                    }

                    if (++qiClasses[tmpQI] > maxCard) {
                        maxCard = qiClasses[tmpQI];
                        maxIndex = tmpQI;
                    }
                }
                else {
                    //ARB: FindOccludee is redundant if ComputeRayCastingVisibility has been called
                    //ARB: change &wFace to wFace and use reference in called function
                    findOccludee<G, I>(fe, grid, epsilon, *ve, &wFace);
                }

                if (wFace) {
                    vector<Vec3r> vertices;
                    for (int i = 0, numEdges = wFace->numberOfEdges(); i < numEdges; ++i) {
                        vertices.push_back(Vec3r(wFace->GetVertex(i)->GetVertex()));
                    }
                    Polygon3r poly(vertices, wFace->GetNormal());
                    poly.userdata = (void *)wFace;
                    fe->setaFace(poly);
                    wFaces.push_back(wFace);
                }
                ++nSamples;
            }

            even_test = ! even_test;
            fe = fe->nextEdge();
        } while ((maxCard < qiMajority) && (fe) && (fe != festart));

        // qi --
        (*ve)->setQI(maxIndex);

        // occluders --
        for (set<ViewShape*>::iterator o = foundOccluders.begin(), oend = foundOccluders.end(); o != oend; ++o) {
            (*ve)->AddOccluder((*o));
        }

        // occludee --
        if (!wFaces.empty()) {
            if (wFaces.size() < nSamples / 2) {
                (*ve)->setaShape(0);
            }
            else {
                ViewShape *vshape = ioViewMap->viewShape((*wFaces.begin())->GetVertex(0)->shape()->GetId());
                (*ve)->setaShape(vshape);
            }
        }

        wFaces.clear();
    }
}

template <typename G, typename I>
static void computeVeryFastVisibility(ViewMap *ioViewMap, G& grid, real epsilon)
{
    vector<ViewEdge*>& vedges = ioViewMap->ViewEdges();

    FEdge *fe;
    unsigned qi = 0;
    WFace *wFace = 0;

    for (vector<ViewEdge*>::iterator ve = vedges.begin(), veend = vedges.end(); ve != veend; ve++) {
        // Find an edge to test
        if (!(*ve)->isInImage()) {
            // This view edge has been proscenium culled
            (*ve)->setQI(255);
            (*ve)->setaShape(0);
            continue;
        }
        fe = (*ve)->fedgeA();
        // Find a FEdge inside the occluder proscenium to test for visibility
        FEdge *festart = fe;
        while (fe && !fe->isInImage() && fe != festart) {
            fe = fe->nextEdge();
        }

        // Test edge
        if (!fe || !fe->isInImage()) {
            // There are no occludable FEdges on this ViewEdge
            // This should be impossible.
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "View Edge in viewport without occludable FEdges: " << (*ve)->getId() << endl;
            }
            // We can recover from this error:
            // Treat this edge as fully visible with no occludee
            qi = 0;
            wFace = NULL;
        }
        else {
            qi = computeVisibility<G, I>(ioViewMap, fe, grid, epsilon, *ve, &wFace, NULL);
        }

        // Store test results
        if (wFace) {
            vector<Vec3r> vertices;
            for (int i = 0, numEdges = wFace->numberOfEdges(); i < numEdges; ++i) {
                vertices.push_back(Vec3r(wFace->GetVertex(i)->GetVertex()));
            }
            Polygon3r poly(vertices, wFace->GetNormal());
            poly.userdata = (void *)wFace;
            fe->setaFace(poly);  // This works because setaFace *copies* the polygon
            ViewShape *vshape = ioViewMap->viewShape(wFace->GetVertex(0)->shape()->GetId());
            (*ve)->setaShape(vshape);
        }
        else {
            (*ve)->setaShape(0);
        }
        (*ve)->setQI(qi);
    }
}

void ViewMapBuilder::BuildGrid(WingedEdge& we, const BBox<Vec3r>& bbox, unsigned int sceneNumFaces)
{
    _Grid->clear();
    Vec3r size;
    for (unsigned int i = 0; i < 3; i++) {
        size[i] = fabs(bbox.getMax()[i] - bbox.getMin()[i]);
        // let make the grid 1/10 bigger to avoid numerical errors while computing triangles/cells intersections.
        size[i] += size[i] / 10.0;
        if (size[i] == 0) {
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "Warning: the bbox size is 0 in dimension " << i << endl;
            }
        }
    }
    _Grid->configure(Vec3r(bbox.getMin() - size / 20.0), size, sceneNumFaces);

    // Fill in the grid:
    WFillGrid fillGridRenderer(_Grid, &we);
    fillGridRenderer.fillGrid();

    // DEBUG
    _Grid->displayDebug();
}

ViewMap *ViewMapBuilder::BuildViewMap(WingedEdge& we, visibility_algo iAlgo, real epsilon,
                                      const BBox<Vec3r>& bbox, unsigned int sceneNumFaces)
{
    _ViewMap = new ViewMap;
    _currentId = 1;
    _currentFId = 0;
    _currentSVertexId = 0;

    // Builds initial view edges
    computeInitialViewEdges(we);

    // Detects cusps
    computeCusps(_ViewMap);

    // Compute intersections
    ComputeIntersections(_ViewMap, sweep_line, epsilon);

    // Compute visibility
    ComputeEdgesVisibility(_ViewMap, we, bbox, sceneNumFaces, iAlgo, epsilon);

    PropagateVisibilty(_ViewMap);

    CheckVisibilityCoherency(_ViewMap);

    return _ViewMap;
}

static inline real distance2D(const Vec3r & point, const real origin[2])
{
    return ::hypot((point[0] - origin[0]), (point[1] - origin[1]));
}

static inline bool crossesProscenium(real proscenium[4], FEdge *fe)
{
    Vec2r min(proscenium[0], proscenium[2]);
    Vec2r max(proscenium[1], proscenium[3]);
    Vec2r A(fe->vertexA()->getProjectedX(), fe->vertexA()->getProjectedY());
    Vec2r B(fe->vertexB()->getProjectedX(), fe->vertexB()->getProjectedY());

    return GeomUtils::intersect2dSeg2dArea (min, max, A, B);
}

static inline bool insideProscenium(real proscenium[4], const Vec3r& point)
{
    return !(point[0] < proscenium[0] || point[0] > proscenium[1] ||
            point[1] < proscenium[2] || point[1] > proscenium[3]);
}

void ViewMapBuilder::CullViewEdges(ViewMap *ioViewMap, real viewProscenium[4], real occluderProscenium[4],
bool extensiveFEdgeSearch)
{
    // Cull view edges by marking them as non-displayable.
    // This avoids the complications of trying to delete edges from the ViewMap.

    // Non-displayable view edges will be skipped over during visibility calculation.

    // View edges will be culled according to their position w.r.t. the viewport proscenium (viewport + 5% border,
    // or some such).

    // Get proscenium boundary for culling
    GridHelpers::getDefaultViewProscenium(viewProscenium);
    real prosceniumOrigin[2];
    prosceniumOrigin[0] = (viewProscenium[1] - viewProscenium[0]) / 2.0;
    prosceniumOrigin[1] = (viewProscenium[3] - viewProscenium[2]) / 2.0;
    if (_global.debug & G_DEBUG_FREESTYLE) {
        cout << "Proscenium culling:" << endl;
        cout << "Proscenium: [" << viewProscenium[0] << ", " << viewProscenium[1] << ", " << viewProscenium[2] <<
                ", " << viewProscenium[3] << "]"<< endl;
        cout << "Origin: [" << prosceniumOrigin[0] << ", " << prosceniumOrigin[1] << "]"<< endl;
    }

    // A separate occluder proscenium will also be maintained, starting out the same as the viewport proscenium, and
    // expanding as necessary so that it encompasses the center point of at least one feature edge in each retained view
    // edge.
    // The occluder proscenium will be used later to cull occluding triangles before they are inserted into the Grid.
    // The occluder proscenium starts out the same size as the view proscenium
    GridHelpers::getDefaultViewProscenium(occluderProscenium);

    // N.B. Freestyle is inconsistent in its use of ViewMap::viewedges_container and vector<ViewEdge*>::iterator.
    // Probably all occurences of vector<ViewEdge*>::iterator should be replaced ViewMap::viewedges_container
    // throughout the code.
    // For each view edge
    ViewMap::viewedges_container::iterator ve, veend;

    for (ve = ioViewMap->ViewEdges().begin(), veend = ioViewMap->ViewEdges().end(); ve != veend; ve++) {
        // Overview:
        //    Search for a visible feature edge
        //    If none: mark view edge as non-displayable
        //    Otherwise:
        //        Find a feature edge with center point inside occluder proscenium.
        //        If none exists, find the feature edge with center point closest to viewport origin.
        //            Expand occluder proscenium to enclose center point.

        // For each feature edge, while bestOccluderTarget not found and view edge not visibile
        bool bestOccluderTargetFound = false;
        FEdge *bestOccluderTarget = NULL;
        real bestOccluderDistance = 0.0;
        FEdge *festart = (*ve)->fedgeA();
        FEdge *fe = festart;
        // All ViewEdges start culled
        (*ve)->setIsInImage(false);

        // For simple visibility calculation: mark a feature edge that is known to have a center point inside the
        // occluder proscenium. Cull all other feature edges.
        do {
            // All FEdges start culled
            fe->setIsInImage(false);

            // Look for the visible edge that can most easily be included in the occluder proscenium.
            if (!bestOccluderTargetFound) {
                // If center point is inside occluder proscenium,
                if (insideProscenium(occluderProscenium, fe->center2d())) {
                    // Use this feature edge for visibility deterimination
                    fe->setIsInImage(true);
                    // Mark bestOccluderTarget as found
                    bestOccluderTargetFound = true;
                    bestOccluderTarget = fe;
                }
                else {
                    real d = distance2D(fe->center2d(), prosceniumOrigin);
                    // If center point is closer to viewport origin than current target
                    if (bestOccluderTarget == NULL || d < bestOccluderDistance) {
                        // Then store as bestOccluderTarget
                        bestOccluderDistance = d;
                        bestOccluderTarget = fe;
                    }
                }
            }

            // If feature edge crosses the view proscenium
            if (!(*ve)->isInImage() && crossesProscenium(viewProscenium, fe)) {
                // Then the view edge will be included in the image
                (*ve)->setIsInImage(true);
            }
            fe = fe->nextEdge();
        } while (fe && fe != festart && !(bestOccluderTargetFound && (*ve)->isInImage()));

        // Either we have run out of FEdges, or we already have the one edge we need to determine visibility
        // Cull all remaining edges.
        while (fe && fe != festart) {
            fe->setIsInImage(false);
            fe = fe->nextEdge();
        }

        // If bestOccluderTarget was not found inside the occluder proscenium, we need to expand the occluder
        // proscenium to include it.
        if ((*ve)->isInImage() && bestOccluderTarget != NULL && !bestOccluderTargetFound) {
            // Expand occluder proscenium to enclose bestOccluderTarget
            Vec3r point = bestOccluderTarget->center2d();
            if (point[0] < occluderProscenium[0]) {
                occluderProscenium[0] = point[0];
            }
            else if (point[0] > occluderProscenium[1]) {
                occluderProscenium[1] = point[0];
            }
            if (point[1] < occluderProscenium[2]) {
                occluderProscenium[2] = point[1];
            }
            else if (point[1] > occluderProscenium[3]) {
                occluderProscenium[3] = point[1];
            }
            // Use bestOccluderTarget for visibility determination
            bestOccluderTarget->setIsInImage(true);
        }
    }

    // We are done calculating the occluder proscenium.
    // Expand the occluder proscenium by an epsilon to avoid rounding errors.
    const real epsilon = 1.0e-6;
    occluderProscenium[0] -= epsilon;
    occluderProscenium[1] += epsilon;
    occluderProscenium[2] -= epsilon;
    occluderProscenium[3] += epsilon;

    // For "Normal" or "Fast" style visibility computation only:

    // For more detailed visibility calculation, make a second pass through the view map, marking all feature edges
    // with center points inside the final occluder proscenium. All of these feature edges can be considered during
    // visibility calculation.

    // So far we have only found one FEdge per ViewEdge.  The "Normal" and "Fast" styles of visibility computation
    // want to consider many FEdges for each ViewEdge.
    // Here we re-scan the view map to find any usable FEdges that we skipped on the first pass, or that have become
    // usable because the occluder proscenium has been expanded since the edge was visited on the first pass.
    if (extensiveFEdgeSearch) {
        // For each view edge,
        for (ve = ioViewMap->ViewEdges().begin(), veend = ioViewMap->ViewEdges().end(); ve != veend; ve++) {
            if (!(*ve)->isInImage()) {
                continue;
            }
            // For each feature edge,
            FEdge *festart = (*ve)->fedgeA();
            FEdge *fe = festart;
            do {
                // If not (already) visible and center point inside occluder proscenium,
                if (!fe->isInImage() && insideProscenium(occluderProscenium, fe->center2d())) {
                    // Use the feature edge for visibility determination
                    fe->setIsInImage(true);
                }
                fe = fe->nextEdge();
            } while (fe && fe != festart);
        }
    }
}

void ViewMapBuilder::computeInitialViewEdges(WingedEdge& we)
{
    vector<WShape*> wshapes = we.getWShapes();
    SShape *psShape;

    for (vector<WShape*>::const_iterator it = wshapes.begin(); it != wshapes.end(); it++) {
        if (_pRenderMonitor && _pRenderMonitor->testBreak())
            break;

        // create the embedding
        psShape = new SShape;
        psShape->setId((*it)->GetId());
        psShape->setName((*it)->getName());
        psShape->setFrsMaterials((*it)->frs_materials()); // FIXME

        // create the view shape
        ViewShape *vshape = new ViewShape(psShape);
        // add this view shape to the view map:
        _ViewMap->AddViewShape(vshape);

        // we want to number the view edges in a unique way for the while scene.
        _pViewEdgeBuilder->setCurrentViewId(_currentId);
        // we want to number the feature edges in a unique way for the while scene.
        _pViewEdgeBuilder->setCurrentFId(_currentFId);
        // we want to number the SVertex in a unique way for the while scene.
        _pViewEdgeBuilder->setCurrentSVertexId(_currentFId);
        _pViewEdgeBuilder->BuildViewEdges(dynamic_cast<WXShape*>(*it), vshape, _ViewMap->ViewEdges(),
                                          _ViewMap->ViewVertices(), _ViewMap->FEdges(), _ViewMap->SVertices());

        _currentId = _pViewEdgeBuilder->currentViewId() + 1;
        _currentFId = _pViewEdgeBuilder->currentFId() + 1;
        _currentSVertexId = _pViewEdgeBuilder->currentSVertexId() + 1;

        psShape->ComputeBBox();
    }
}


WFace * GetNearFace(WEdge * edge, Vec3r viewpoint, real * discriminant = NULL)
{
    // find the opp vertices
    WVertex * oppA = NULL;
    for(int j=0;j<3;j++)
    {
        WVertex * vert = edge->GetaFace()->GetVertex(j);
        if (vert != edge->GetaVertex() && vert != edge->GetbVertex())
        {
            oppA = vert;
            break;
        }
    }
    WVertex * oppB = NULL;
    for(int j=0;j<3;j++)
    {
        WVertex * vert = edge->GetbFace()->GetVertex(j);
        if (vert != edge->GetaVertex() && vert != edge->GetbVertex())
        {
            oppB = vert;
            break;
        }
    }

    Vec3r vA = edge->GetaVertex()->GetVertex();
    Vec3r vB = edge->GetbVertex()->GetVertex();
    Vec3r midpoint = (vA+vB)/2;
    Vec3r edgeDir = vA-vB;

    edgeDir = edgeDir/edgeDir.norm();

    Vec3r vv = viewpoint - midpoint;

    // unit vectors in the plane of each face (?), perpendicular to the edge
    Vec3r edgePerpA =   oppA->GetVertex() - midpoint;
    edgePerpA = edgePerpA -  (edgePerpA*edgeDir) * edgeDir;
    edgePerpA = edgePerpA/edgePerpA.norm();
    Vec3r edgePerpB =   oppB->GetVertex() - midpoint;
    edgePerpB = edgePerpB - (edgePerpB * edgeDir) * edgeDir;
    edgePerpB = edgePerpB/edgePerpB.norm();

    // unit vectors in the plane of each face, parallel to the view vector to the midpoint and with the same orientation
    Vec3r adir = vv - (vv * edge->GetaFace()->GetNormal()) * edge->GetaFace()->GetNormal();
    Vec3r bdir = vv - (vv * edge->GetbFace()->GetNormal()) * edge->GetbFace()->GetNormal();

    if (discriminant != NULL)
        *discriminant = (adir -bdir)* vv / vv.norm();

    if (adir * edgePerpA < 0)
        adir = adir/(-adir.norm());
    else
        adir = adir/adir.norm();
    if (bdir * edgePerpB < 0)
        bdir = bdir/(-bdir.norm());
    else
        bdir = bdir/bdir.norm();

    return ((adir - bdir) * vv > 0)  ? edge->GetaFace() : edge->GetbFace();
}

// determine which side of the triangle PQR point A lies on.
// positive for one side, negative for the other
// magnitude is approximately six times the volume of the tetrahedron (PQRA)
inline real orient3d(Vec3r P, Vec3r Q, Vec3r R, Vec3r A)
{
    real p[3] = { P[0], P[1], P[2] };
    real q[3] = { Q[0], Q[1], Q[2] };
    real r[3] = { R[0], R[1], R[2] };
    real a[3] = { A[0], A[1], A[2] };

    return carve::exact::orient3dexact(p,q,r,a);
}

inline bool sameSide(Vec3r P, Vec3r Q, Vec3r R, Vec3r A1, Vec3r A2)
{
    real orient1 = orient3d(P,Q,R,A1);
    real orient2 = orient3d(P,Q,R,A2);

    if (orient1 == 0 || orient2 == 0)
        return false;

    return (orient1 > 0) == (orient2 > 0);
}

inline bool sameSide(WFace * face, Vec3r A1, Vec3r A2)
{
    return sameSide(face->GetVertex(0)->GetVertex(), face->GetVertex(1)->GetVertex(), face->GetVertex(2)->GetVertex(), A1, A2);
}

// is the sharp edge "fesh" that ends at the vertex "wv" occluded by another face in the one ring of
// the vertex?
// result: 1: overlap exists, 0: no overlap, -1: can't tell due to inconsistency
int OneRingOcclusion(SVertex * sv, WVertex * wv, FEdgeSharp * fesh, set<WXFace*> & oneRing,
                     Vec3r & viewpoint, bool useConsistency)
{
    int result = 0;

    // iterate over the one-ring of wv and check if any face occludes the edge fesh
    for(set<WXFace*>::iterator fit = oneRing.begin(); fit != oneRing.end(); ++fit)
    {
        WXFace * face = *fit;

        // find the vertex's index in this triangle
        int vind = face->GetIndex(wv);

        // check if the face is back-facing
        if (!face->front(useConsistency))
            continue;

        // check if the face is inconsistent
        if (useConsistency && !face->consistent())
        {
            result = -1;
            continue;
        }

        // check if the face is adjacent to the edge
        if (fesh->edge()->GetaFace() == face || fesh->edge()->GetbFace() == face)
            continue;

        // check for image-space overlap.  this is reduces to three clipping tests

        SVertex * endpoint = (fesh->vertexA() == sv ? fesh->vertexB() : fesh->vertexA());
        Vec3r endPos = endpoint->getPoint3D();
        Vec3r vPos = wv->GetVertex();
        Vec3r v1Pos = face->GetVertex( (vind+1)%3)->GetVertex();
        Vec3r v2Pos = face->GetVertex( (vind+2)%3)->GetVertex();

        if (!sameSide(face, endPos, viewpoint) &&
                sameSide(viewpoint, vPos, v1Pos, endPos, v2Pos) &&
                sameSide(viewpoint, vPos, v2Pos, endPos, v1Pos))
            return 1;
    }
    return result;
}

// check if the sharp edge is adjacent to a vertex, and occluded by the one-ring of that vertex
// result: 1: overlap exists, 0: no overlap, -1: can't tell due to inconsistency
int OneRingOcclusion(FEdgeSharp * fesh, Vec3r & viewpoint, bool useConsistency)
{
    int result = 0;

    for(int i=0;i<2;i++)
    {
        SVertex * sv = (i == 0 ? fesh->vertexA() : fesh->vertexB());
        WVertex * wv = sv->sourceVertex();
        if (wv == NULL)
            continue;

        set<WXFace*> oneRing;

        // collect the one-ring
        for(vector<WEdge*>::iterator wit = wv->GetEdges().begin(); wit != wv->GetEdges().end(); ++wit)
        {
            if ( (*wit)->GetaFace() != NULL)
                oneRing.insert( (WXFace*)(*wit)->GetaFace());
            if ( (*wit)->GetbFace() != NULL)
                oneRing.insert( (WXFace*)(*wit)->GetbFace());
        }

        int r = OneRingOcclusion(sv, wv, fesh, oneRing, viewpoint, useConsistency);

        if (r == 1)
            return 1;

        if (r == -1)
            result = -1;
    }

    return result;
}

// check if the vertex sv is a cusp  based on overlaps
// result: 1: cusp, 0: not a cusp, -1: can't tell due to inconsistency
int IsCusp(SVertex * sv, WVertex * wv, FEdgeSharp * fesh0, FEdgeSharp * fesh1, set<WXFace*> & oneRing,
           Vec3r & viewpoint, bool useConsistency)
{
    int r1 = OneRingOcclusion(sv, wv, fesh0, oneRing, viewpoint, useConsistency);

    if (r1 == -1)
        return -1;

    if (r1 == 1)
        return 1;

    int r2 = OneRingOcclusion(sv, wv, fesh1, oneRing, viewpoint, useConsistency);

    if (r2 == -1)
        return -1;

    if (r2 == 1)
        return 1;

    return 0;
}

bool DiscreteRadialSign(FEdgeSmooth * fesh, Vec3r & viewpoint)
{
    Vec3r A = fesh->vertexA()->point3d();
    Vec3r B = fesh->vertexB()->point3d();
    Vec3r AB = B-A;
    AB.normalize();
    Vec3r m((A+B)/2.0);
    Vec3r crossP(AB^fesh->normal());
    crossP.normalize();
    Vec3r viewvector(m-viewpoint);
    viewvector.normalize();

    //---- point crossP in the positive direction ----
    WFace * face = (WFace*)fesh->face();
    real maxNdotVMag = 0;
    int vind = -1;
    for(int i=0;i<3;i++)
    {
        real ndotv = ((WXVertex*)face->GetVertex(i))->ndotv();
        if (fabs(ndotv) > maxNdotVMag)
        {
            vind = i;
            maxNdotVMag = fabs(ndotv);
        }
    }
    if (vind != -1)
    {
        real ndotv = ((WXVertex*)face->GetVertex(vind))->ndotv();
        Vec3r v = face->GetVertex(vind)->GetVertex() - m;
        v.normalize();
        if ( (v * crossP > 0) != (ndotv > 0))
            crossP = crossP * -1;
    }

    return crossP * viewvector > 0;
}

void ViewMapBuilder::computeCusps(ViewMap *ioViewMap)
{
#if 0
    vector<ViewVertex*> newVVertices;
    vector<ViewEdge*> newVEdges;
    ViewMap::viewedges_container& vedges = ioViewMap->ViewEdges();
    ViewMap::viewedges_container::iterator ve = vedges.begin(), veend = vedges.end();
    for (; ve != veend; ++ve) {
        if (_pRenderMonitor && _pRenderMonitor->testBreak())
            break;
        if ((!((*ve)->getNature() & Nature::SILHOUETTE)) || (!((*ve)->fedgeA()->isSmooth())))
            continue;
        FEdge *fe = (*ve)->fedgeA();
        FEdge *fefirst = fe;
        bool first = true;
        bool positive = true;
        do {
            FEdgeSmooth *fes = dynamic_cast<FEdgeSmooth*>(fe);
            Vec3r A((fes)->vertexA()->point3d());
            Vec3r B((fes)->vertexB()->point3d());
            Vec3r AB(B - A);
            AB.normalize();
            Vec3r m((A + B) / 2.0);
            Vec3r crossP(AB ^ (fes)->normal());
            crossP.normalize();
            Vec3r viewvector;
            if (_orthographicProjection) {
                viewvector = Vec3r(0.0, 0.0, m.z() - _viewpoint.z());
            }
            else {
                viewvector = Vec3r(m - _viewpoint);
            }
            viewvector.normalize();
            if (first) {
                if (((crossP) * (viewvector)) > 0)
                    positive = true;
                else
                    positive = false;
                first = false;
            }
            // If we're in a positive part, we need a stronger negative value to change
            NonTVertex *cusp = NULL;
            if (positive) {
                if (((crossP) * (viewvector)) < -0.1) {
                    // state changes
                    positive = false;
                    // creates and insert cusp
                    cusp = dynamic_cast<NonTVertex*>(ioViewMap->InsertViewVertex(fes->vertexA(), newVEdges));
                    if (cusp)
                        cusp->setNature(cusp->getNature() | Nature::CUSP);
                }
            }
            else {
                // If we're in a negative part, we need a stronger negative value to change
                if (((crossP) * (viewvector)) > 0.1) {
                    positive = true;
                    cusp = dynamic_cast<NonTVertex*>(ioViewMap->InsertViewVertex(fes->vertexA(), newVEdges));
                    if (cusp)
                        cusp->setNature(cusp->getNature() | Nature::CUSP);
                }
            }
            fe = fe->nextEdge();
        } while (fe && fe != fefirst);
    }
    for (ve = newVEdges.begin(), veend = newVEdges.end(); ve != veend; ++ve) {
        (*ve)->viewShape()->AddEdge(*ve);
        vedges.push_back(*ve);
    }
#else

    vector<ViewEdge*> newVEdges;

    // --------------------- cusps for smooth silhouettes ---------------------

    for(vector<SVertex*>::iterator sit = ioViewMap->SVertices().begin(); sit != ioViewMap->SVertices().end(); ++sit)
    {
        SVertex * sv = *sit;

        // only need to insert cusps in the interiors of ViewEdges
        if (sv->viewvertex() != NULL)
            continue;

        // retreive the two adjacent FEdges

        FEdge * fes[2];

        fes[0] = sv->fedge();
        fes[1] = fes[0]->vertexB() == sv ? fes[0]->nextEdge() : fes[0]->previousEdge();

        if (fes[0] == NULL || fes[1] == NULL)
            continue;

        // only consider smooth silhouettes and boundaries
        if (!(fes[0]->getNature() & Nature::SILHOUETTE) || !fes[0]->isSmooth())
            continue;

        FEdgeSmooth * fesh[2] = { dynamic_cast<FEdgeSmooth*>(fes[0]), dynamic_cast<FEdgeSmooth*>(fes[1]) };

        // sign (viewvec * (AB ^ normal))

        bool frontFacing[2];
        for(int i=0;i<2;i++)
            frontFacing[i] = DiscreteRadialSign(fesh[i], _viewpoint);

        if (frontFacing[0] != frontFacing[1])
        {
            // this point is a cusp
            ViewVertex* cusp = ioViewMap->InsertViewVertex(sv, newVEdges);
            if (cusp != NULL)
            {
                cusp->setNature(cusp->getNature()|Nature::CUSP);
            }
        }
    }
    ViewMap::viewedges_container& vedges = ioViewMap->ViewEdges();

    for(ViewMap::viewedges_container::iterator ve=newVEdges.begin(), veend=newVEdges.end(); ve!=veend; ++ve){
        (*ve)->viewShape()->AddEdge(*ve);
        vedges.push_back(*ve);
    }


    // --------------------- cusps for mesh boundaries and silhouettes ---------------------

    for(vector<SVertex*>::iterator sit = ioViewMap->SVertices().begin(); sit != ioViewMap->SVertices().end(); ++sit)
    {
        SVertex * sv = *sit;

        // only need to insert cusps in the interiors of ViewEdges
        if (sv->viewvertex() != NULL)
            continue;

        // retreive the two adjacent FEdges
        FEdge * fes[2];

        fes[0] = sv->fedge();
        fes[1] = fes[0]->vertexB() == sv ? fes[0]->nextEdge() : fes[0]->previousEdge();

        if (fes[0] == NULL || fes[1] == NULL)
            continue;

        // only consider mesh silhouettes and boundaries.
        if (!(fes[0]->getNature() & (Nature::SILHOUETTE | Nature::BORDER)) || fes[0]->isSmooth())
            continue;

        FEdgeSharp * fesh[2] = { dynamic_cast<FEdgeSharp*>(fes[0]), dynamic_cast<FEdgeSharp*>(fes[1]) };

        // find the WVertex shared by both edges
        WVertex * wv = WEdge::CommonVertex(fesh[0]->edge(), fesh[1]->edge());
        //assert(wv == sv->GetSourceVertex()); // another way to get the same vertex

        set<WXFace*> oneRing;

        // collect the one-ring
        for(vector<WEdge*>::iterator wit = wv->GetEdges().begin(); wit != wv->GetEdges().end(); ++wit)
        {
            if ( (*wit)->GetaFace() != NULL)
                oneRing.insert( (WXFace*)(*wit)->GetaFace());
            if ( (*wit)->GetbFace() != NULL)
                oneRing.insert( (WXFace*)(*wit)->GetbFace());
        }

        // check if there are inconsistencies or overlaps in the neighborhood
        int cuspTest = IsCusp(sv, wv, fesh[0], fesh[1], oneRing, _viewpoint, _useConsistency);
        bool isCusp;

        if (fesh[0]->getNature() & Nature::SILHOUETTE)
        {
            // use Markosian's algorithm for silhouettes with inconsistency, because it can kinda handle inconsistency
            //
            // Markosian et al. define a silhouette edge as front-facing if the adjacent face nearer to the camera is front-facing.
            // A vertex is a cusp if it connects a front-facing to a back-facing edge.

            bool facing[2];

            for(int i=0;i<2;i++)
            {
                WFace * nearerFace = GetNearFace(fesh[i]->edge(),_viewpoint);

                facing[i] = ((WXFace*)nearerFace)->front(_useConsistency);
            }
            isCusp = (facing[0] != facing[1]);
        }
        else
            if (cuspTest == 1)
                isCusp = true;
            else
                if (cuspTest == 0)
                    isCusp = false;
                else
                    isCusp = true; // being conservative

        if (isCusp)
        {
            // insert the cusp and split the ViewEdge
            ViewShape * shape = fes[0]->viewedge()->viewShape();

            NonTVertex * ntv = new NonTVertex(sv);

            ioViewMap->AddViewVertex(ntv);
            shape->AddVertex(ntv);
            if (cuspTest != -1)
                ntv->setNature( ntv->getNature() | Nature::CUSP);
            else
                ntv->setNature( ntv->getNature() | Nature::CUSP | Nature::AMBIG_CUSP);

            ViewEdge * vEdge = fes[0]->viewedge();
            ViewEdge * newViewEdge = NULL;

            shape->SplitEdge(ntv, vEdge, newViewEdge, false);

            if (newViewEdge != vEdge)
                ioViewMap->AddViewEdge(newViewEdge);
        }
    }
#endif
}

void ViewMapBuilder::ComputeCumulativeVisibility(ViewMap *ioViewMap, WingedEdge& we, const BBox<Vec3r>& bbox,
                                                 real epsilon, bool cull, GridDensityProviderFactory& factory)
{
    AutoPtr<GridHelpers::Transform> transform;
    AutoPtr<OccluderSource> source;

    if (_orthographicProjection) {
        transform.reset(new BoxGrid::Transform);
    }
    else {
        transform.reset(new SphericalGrid::Transform);
    }

    if (cull) {
        source.reset(new CulledOccluderSource(*transform, we, *ioViewMap, true));
    }
    else {
        source.reset(new OccluderSource(*transform, we));
    }

    AutoPtr<GridDensityProvider> density(factory.newGridDensityProvider(*source, bbox, *transform));

    if (_orthographicProjection) {
        BoxGrid grid(*source, *density, ioViewMap, _viewpoint, _EnableQI);
        computeCumulativeVisibility<BoxGrid, BoxGrid::Iterator>(ioViewMap, grid, epsilon, _pRenderMonitor);
    }
    else {
        SphericalGrid grid(*source, *density, ioViewMap, _viewpoint, _EnableQI);
        computeCumulativeVisibility<SphericalGrid, SphericalGrid::Iterator>(ioViewMap, grid, epsilon, _pRenderMonitor);
    }
}

void ViewMapBuilder::ComputeDetailedVisibility(ViewMap *ioViewMap, WingedEdge& we, const BBox<Vec3r>& bbox,
                                               real epsilon, bool cull, GridDensityProviderFactory& factory)
{
    AutoPtr<GridHelpers::Transform> transform;
    AutoPtr<OccluderSource> source;

    if (_orthographicProjection) {
        transform.reset(new BoxGrid::Transform);
    }
    else {
        transform.reset(new SphericalGrid::Transform);
    }

    if (cull) {
        source.reset(new CulledOccluderSource(*transform, we, *ioViewMap, true));
    }
    else {
        source.reset(new OccluderSource(*transform, we));
    }

    AutoPtr<GridDensityProvider> density(factory.newGridDensityProvider(*source, bbox, *transform));

    if (_orthographicProjection) {
        BoxGrid grid(*source, *density, ioViewMap, _viewpoint, _EnableQI);
        computeDetailedVisibility<BoxGrid, BoxGrid::Iterator>(ioViewMap, grid, epsilon, _pRenderMonitor);
    }
    else {
        SphericalGrid grid(*source, *density, ioViewMap, _viewpoint, _EnableQI);
        computeDetailedVisibility<SphericalGrid, SphericalGrid::Iterator>(ioViewMap, grid, epsilon, _pRenderMonitor);
    }
}

void ViewMapBuilder::ComputeEdgesVisibility(ViewMap *ioViewMap, WingedEdge& we, const BBox<Vec3r>& bbox,
                                            unsigned int sceneNumFaces, visibility_algo iAlgo, real epsilon)
{
#if 0
    iAlgo = ray_casting; // for testing algorithms equivalence
#endif
    switch (iAlgo) {
    case ray_casting:
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Using ordinary ray casting" << endl;
        }
        BuildGrid(we, bbox, sceneNumFaces);
        ComputeRayCastingVisibility(ioViewMap, epsilon);
        break;
    case ray_casting_fast:
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Using fast ray casting" << endl;
        }
        BuildGrid(we, bbox, sceneNumFaces);
        ComputeFastRayCastingVisibility(ioViewMap, epsilon);
        break;
    case ray_casting_very_fast:
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Using very fast ray casting" << endl;
        }
        BuildGrid(we, bbox, sceneNumFaces);
        ComputeVeryFastRayCastingVisibility(ioViewMap, epsilon);
        break;
    case ray_casting_culled_adaptive_traditional:
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Using culled adaptive grid with heuristic density and traditional QI calculation" << endl;
        }
        try {
        HeuristicGridDensityProviderFactory factory(0.5f, sceneNumFaces);
        ComputeDetailedVisibility(ioViewMap, we, bbox, epsilon, true, factory);
    }
        catch (...) {
        // Last resort catch to make sure RAII semantics hold for OptimizedGrid. Can be replaced with
        // try...catch block around main() if the program as a whole is converted to RAII

        // This is the little-mentioned caveat of RAII: RAII does not work unless destructors are always
        // called, but destructors are only called if all exceptions are caught (or std::terminate() is
        // replaced).

        // We don't actually handle the exception here, so re-throw it now that our destructors have had a
        // chance to run.
        throw;
    }
        break;
    case ray_casting_adaptive_traditional:
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Using unculled adaptive grid with heuristic density and traditional QI calculation" << endl;
        }
        try {
        HeuristicGridDensityProviderFactory factory(0.5f, sceneNumFaces);
        ComputeDetailedVisibility(ioViewMap, we, bbox, epsilon, false, factory);
    }
        catch (...) {
        throw;
    }
        break;
    case ray_casting_culled_adaptive_cumulative:
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Using culled adaptive grid with heuristic density and cumulative QI calculation" << endl;
        }
        try {
        HeuristicGridDensityProviderFactory factory(0.5f, sceneNumFaces);
        ComputeCumulativeVisibility(ioViewMap, we, bbox, epsilon, true, factory);
    }
        catch (...) {
        throw;
    }
        break;
    case ray_casting_adaptive_cumulative:
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Using unculled adaptive grid with heuristic density and cumulative QI calculation" << endl;
        }
        try {
        HeuristicGridDensityProviderFactory factory(0.5f, sceneNumFaces);
        ComputeCumulativeVisibility(ioViewMap, we, bbox, epsilon, false, factory);
    }
        catch (...) {
        throw;
    }
        break;
    default:
        break;
    }
}

real ArcLength2D(ViewEdge * ve, real maxLength)
{
    real arcLength = 0;

    FEdge * fe = ve->fedgeA();

    do
    {
        arcLength += fe->getLength2D();
        fe = fe->nextEdge();
    } while (fe != NULL && fe != ve->fedgeA() && arcLength < maxLength);

    return arcLength;
}

void ViewMapBuilder::PropagateVisibilty(ViewMap *ioViewMap) {
    printf("Propagating visibility to ambiguous edges\n");

    // gather all chains that do not have any visibility votes
    set<ViewEdge*> ambiguousEdges;

    // mark chains with short image-space extent as ambiguous
    for(vector<ViewEdge*>::iterator vit = ioViewMap->ViewEdges().begin(); vit != ioViewMap->ViewEdges().end(); ++vit){
        if (_graftThreshold > 0)
            if (ArcLength2D(*vit, 2*_graftThreshold) < _graftThreshold)
                (*vit)->markAmbiguous();

        // Spurious cusps heuristic
        NonTVertex * vertA = dynamic_cast<NonTVertex*>((*vit)->A());
        NonTVertex * vertB = dynamic_cast<NonTVertex*>((*vit)->B());
        if (vertA != NULL && vertB != NULL)
            if((vertA->getNature() & Nature::CUSP) && (vertB->getNature() & Nature::CUSP)
                    && ((*vit)->fedgeA()->nextEdge() == NULL)){
                (*vit)->markAmbiguous();
            }
    }

    for(vector<ViewEdge*>::iterator vit = ioViewMap->ViewEdges().begin(); vit != ioViewMap->ViewEdges().end(); ++vit)
        if ( (*vit)->ambiguousVisibility())
            ambiguousEdges.insert(*vit);

    bool changed = true;
    while (ambiguousEdges.size() > 0 && changed)
    {
        changed = false;

        // copy over the current list of ambiguous edges so that we can iterate over it
        vector<ViewEdge*> edgeList;
        for(set<ViewEdge*>::iterator it = ambiguousEdges.begin(); it!= ambiguousEdges.end(); ++it)
            edgeList.push_back(*it);
        //      copy(ambiguousEdges.begin(), ambiguousEdges.end(), edgeList.begin());

        for(vector<ViewEdge*>::iterator vit = edgeList.begin(); vit != edgeList.end(); ++vit)
        {
            // see if we can resolve the ambiguity
            ViewEdge * edge = *vit;

            if (!edge->ambiguousVisibility())
                continue;

            // Spurious cusps heuristic
            NonTVertex * vertA = dynamic_cast<NonTVertex*>(edge->A());
            NonTVertex * vertB = dynamic_cast<NonTVertex*>(edge->B());
            if (vertA != NULL && vertB != NULL)
                if(vertA->getNature() & Nature::CUSP && vertB->getNature() & Nature::CUSP){
                    ViewEdge * mateA = NULL;
                    ViewEdge * mateB = NULL;
                    if(vertA->viewedges().size() == 2)
                        mateA = (vertA->viewedges()[0].first != edge ? vertA->viewedges()[0].first : vertA->viewedges()[1].first);
                    if(vertB->viewedges().size() == 2)
                        mateB = (vertB->viewedges()[0].first != edge ? vertB->viewedges()[0].first : vertB->viewedges()[1].first);
                    if(mateA != NULL && mateB != NULL)
                        if(mateA->qi() == 0 && mateB->qi() == 0){
                            edge->setQI(0);
                            ambiguousEdges.erase(ambiguousEdges.find(edge));
                            edge->fixAmbiguous();
                            changed = true;
                            break;
                        }
                }

            // Fix visibility of tiny bits
            real arcLength = ArcLength2D(edge,FLT_MAX);
            if(arcLength<=0.01){
                TVertex * vertA = dynamic_cast<TVertex*>(edge->A());
                TVertex * vertB = dynamic_cast<TVertex*>(edge->B());
                if (vertA != NULL && vertB != NULL){
                    ViewEdge * mateA = vertA->mate(edge);
                    ViewEdge * mateB = vertB->mate(edge);
                    if(mateA && mateB && mateA->qi()==0 && mateB->qi()==0){
                        edge->setQI(0);
                        ambiguousEdges.erase(ambiguousEdges.find(edge));
                        edge->fixAmbiguous();
                        changed = true;
                        break;
                    }
                }
            }

            ViewVertex * v[2] = {edge->A(), edge->B() };

            for(int i=0;i<2;i++)
            {
                ViewEdge * mate = NULL;

                if (v[i] == NULL)
                    continue;

                TVertex * tvert = dynamic_cast<TVertex*>(v[i]);
                if (tvert != NULL)
                {
                    if (tvert->frontEdgeA().first == edge || tvert->frontEdgeB().first == edge)
                    {
                        mate = tvert->mate(edge);

                        if (mate == NULL || mate->ambiguousVisibility())
                            continue;
                    }
                }
                else
                {
                    NonTVertex * ntv = (NonTVertex*)v[i];
                    assert(ntv);

                    if (ntv->viewedges().size() == 2 && (!(ntv->getNature() & Nature::CUSP) || edge->getNature() == Nature::BORDER ))
                        mate = (ntv->viewedges()[0].first != edge ? ntv->viewedges()[0].first : ntv->viewedges()[1].first);
                }

                if (mate != NULL && mate->qi() == 0)
                {
                    printf("%08X found visible mate %08X\n",edge,mate);
                    edge->setQI(0);
                    ambiguousEdges.erase(ambiguousEdges.find(edge));
                    edge->fixAmbiguous();
                    changed = true;
                    break;
                }
            }
        }
    }
}

void ViewMapBuilder::CheckVisibilityCoherency(ViewMap *ioViewMap)
{
    printf("Check coherence of visibility\n");
    for(vector<ViewVertex*>::iterator vit = _ViewMap->ViewVertices().begin(); vit != _ViewMap->ViewVertices().end(); ++vit)
    {
        ViewVertex * vVertex = *vit;
        TVertex* tVertex = dynamic_cast<TVertex*>(vVertex);
        if(tVertex!=NULL && tVertex->backEdgeA().first && tVertex->backEdgeB().first &&
                tVertex->frontEdgeA().first && tVertex->frontEdgeB().first)
        {
            // Check if visibility at T-junctions involving CUSPs is correct (>=2)
            // If not, make the ViewEdge connected to the front CUSP visible

            int visibleFont = 0, visibleBack = 0;;
            if(tVertex->backEdgeA().first->qi() == 0)
                visibleBack++;
            if(tVertex->frontEdgeA().first->qi() == 0)
                visibleFont++;
            if(tVertex->backEdgeB().first->qi() == 0)
                visibleBack++;
            if(tVertex->frontEdgeB().first->qi() == 0)
                visibleFont++;

            if((visibleFont + visibleBack)==1){
                printf("INCOHERENT VISIBILITY\n");
                char str[200];

                ViewEdge* mate, *edge;
                NonTVertex* ntVertex = NULL;
                ViewEdge* nextE = NULL;
                ViewVertex* nextV = NULL;
                if(visibleFont==1){
                    edge = tVertex->backEdgeB().first;
                    if(edge->B()->getNature() & Nature::CUSP){
                        ntVertex = dynamic_cast<NonTVertex*>(edge->B());
                        assert(ntVertex->viewedges().size() == 2);
                        mate = (ntVertex->viewedges()[0].first != edge ? ntVertex->viewedges()[0].first : ntVertex->viewedges()[1].first);
                        if(mate->qi() == 0){
                            sprintf(str, "FIX 1\n");
                        }else{
                            ntVertex = NULL;
                        }
                    }else if(edge->B()->getNature() & Nature::T_VERTEX){
                        TVertex* tv = dynamic_cast<TVertex*>(edge->B());
                        if(tv->numEdges()<4)
                            continue;
                        if(tv->backEdgeA().first == edge){
                            nextE = tv->backEdgeB().first;
                            nextV = nextE->B();
                        }else if(tv->backEdgeB().first == edge){
                            nextE = tv->backEdgeA().first;
                            nextV = nextE->A();
                        }else if(tv->frontEdgeB().first == edge){
                            nextE = tv->frontEdgeA().first;
                            nextV = nextE->A();
                        }else if(tv->frontEdgeA().first == edge){
                            nextE = tv->frontEdgeB().first;
                            nextV = nextE->B();
                        }
                        if(nextV->getNature() & Nature::CUSP){
                            ntVertex = dynamic_cast<NonTVertex*>(nextV);
                            assert(ntVertex->viewedges().size() == 2);
                            mate = (ntVertex->viewedges()[0].first != edge ? ntVertex->viewedges()[0].first : ntVertex->viewedges()[1].first);
                            if(mate->qi() == 0){
                                sprintf(str, "FIX 1.1\n");
                            }else{
                                ntVertex = NULL;
                                nextE = NULL;
                            }
                        }
                    }
                    if(!ntVertex){
                        edge = tVertex->frontEdgeA().first;
                        if(edge->A()->getNature() & Nature::CUSP){
                            ntVertex = dynamic_cast<NonTVertex*>(edge->A());
                            assert(ntVertex->viewedges().size() == 2);
                            mate = (ntVertex->viewedges()[0].first != edge ? ntVertex->viewedges()[0].first : ntVertex->viewedges()[1].first);
                            if(mate->qi() == 0){
                                sprintf(str, "FIX 2\n");
                            }else{
                                ntVertex = NULL;
                            }
                        }else if(edge->A()->getNature() & Nature::T_VERTEX){
                            TVertex* tv = dynamic_cast<TVertex*>(edge->A());
                            if(tv->numEdges()<4)
                                continue;
                            if(tv->backEdgeA().first == edge){
                                nextE = tv->backEdgeB().first;
                                nextV = nextE->B();
                            }else if(tv->backEdgeB().first == edge){
                                nextE = tv->backEdgeA().first;
                                nextV = nextE->A();
                            }else if(tv->frontEdgeB().first == edge){
                                nextE = tv->frontEdgeA().first;
                                nextV = nextE->A();
                            }else if(tv->frontEdgeA().first == edge){
                                nextE = tv->frontEdgeB().first;
                                nextV = nextE->B();
                            }
                            if(nextV->getNature() & Nature::CUSP){
                                ntVertex = dynamic_cast<NonTVertex*>(nextV);
                                assert(ntVertex->viewedges().size() == 2);
                                mate = (ntVertex->viewedges()[0].first != edge ? ntVertex->viewedges()[0].first : ntVertex->viewedges()[1].first);
                                if(mate->qi() == 0){
                                    sprintf(str, "FIX 2.1\n");
                                }else{
                                    ntVertex = NULL;
                                    nextE = NULL;
                                }
                            }
                        }
                    }
                }else if(visibleBack==1){
                    edge = tVertex->frontEdgeB().first;
                    if(edge->B()->getNature() & Nature::CUSP){
                        ntVertex = dynamic_cast<NonTVertex*>(edge->B());
                        assert(ntVertex->viewedges().size() == 2);
                        mate = (ntVertex->viewedges()[0].first != edge ? ntVertex->viewedges()[0].first : ntVertex->viewedges()[1].first);
                        if(mate->qi() == 0){
                            sprintf(str, "FIX 3\n");
                        }else{
                            ntVertex = NULL;
                        }
                    }else if(edge->B()->getNature() & Nature::T_VERTEX){
                        TVertex* tv = dynamic_cast<TVertex*>(edge->B());
                        if(tv->numEdges()<4)
                            continue;
                        if(tv->backEdgeA().first == edge){
                            nextE = tv->backEdgeB().first;
                            nextV = nextE->B();
                        }else if(tv->backEdgeB().first == edge){
                            nextE = tv->backEdgeA().first;
                            nextV = nextE->A();
                        }else if(tv->frontEdgeB().first == edge){
                            nextE = tv->frontEdgeA().first;
                            nextV = nextE->A();
                        }else if(tv->frontEdgeA().first == edge){
                            nextE = tv->frontEdgeB().first;
                            nextV = nextE->B();
                        }
                        if(nextV->getNature() & Nature::CUSP){
                            ntVertex = dynamic_cast<NonTVertex*>(nextV);
                            assert(ntVertex->viewedges().size() == 2);
                            mate = (ntVertex->viewedges()[0].first != edge ? ntVertex->viewedges()[0].first : ntVertex->viewedges()[1].first);
                            if(mate->qi() == 0){
                                sprintf(str, "FIX 3.1\n");
                            }else{
                                ntVertex = NULL;
                                nextE = NULL;
                            }
                        }
                    }
                    if(!ntVertex){
                        edge = tVertex->frontEdgeA().first;
                        if(edge->A()->getNature() & Nature::CUSP){
                            ntVertex = dynamic_cast<NonTVertex*>(edge->A());
                            assert(ntVertex->viewedges().size() == 2);
                            mate = (ntVertex->viewedges()[0].first != edge ? ntVertex->viewedges()[0].first : ntVertex->viewedges()[1].first);
                            if(mate->qi() == 0){
                                sprintf(str, "FIX 4\n");
                            }else{
                                ntVertex = NULL;
                            }
                        }else if(edge->A()->getNature() & Nature::T_VERTEX){
                            TVertex* tv = dynamic_cast<TVertex*>(edge->A());
                            if(tv->numEdges()<4)
                                continue;
                            if(tv->backEdgeA().first == edge){
                                nextE = tv->backEdgeB().first;
                                nextV = nextE->B();
                            }else if(tv->backEdgeB().first == edge){
                                nextE = tv->backEdgeA().first;
                                nextV = nextE->A();
                            }else if(tv->frontEdgeB().first == edge){
                                nextE = tv->frontEdgeA().first;
                                nextV = nextE->A();
                            }else if(tv->frontEdgeA().first == edge){
                                nextE = tv->frontEdgeB().first;
                                nextV = nextE->B();
                            }
                            if(nextV->getNature() & Nature::CUSP){
                                ntVertex = dynamic_cast<NonTVertex*>(nextV);
                                assert(ntVertex->viewedges().size() == 2);
                                mate = (ntVertex->viewedges()[0].first != edge ? ntVertex->viewedges()[0].first : ntVertex->viewedges()[1].first);
                                if(mate->qi() == 0){
                                    sprintf(str, "FIX 1.4\n");
                                }else{
                                    ntVertex = NULL;
                                    nextE = NULL;
                                }
                            }
                        }
                    }
                }
                if(ntVertex){
                    assert(ntVertex->viewedges().size() == 2);
                    assert(mate->qi() == 0);
                    real arcLength = ArcLength2D(edge,FLT_MAX);
                    if(arcLength<=2.0){
                        printf("%s",str);
                        edge->setQI(0);
                        edge->markInconsistent(true);
                    }
                    if(nextE != NULL){
                        nextE->setQI(0);
                        nextE->markInconsistent(true);
                    }
                }
            }
        }else{
            // If 2 CUSPs are directly (no image space intersection) connected by an invisible edge with arclength<=1px, change it to visible
            NonTVertex* ntVertex = dynamic_cast<NonTVertex*>(vVertex);
            if(ntVertex && ntVertex->getNature() & Nature::CUSP){
                assert(ntVertex->viewedges().size()==2);
                for(int i=0; i<2; i++){
                    ViewEdge* edge = ntVertex->viewedges()[i].first;
                    ViewVertex* otherVVertex = (edge->A() != vVertex) ? edge->A() : edge->B();
                    NonTVertex* otherNTVertex = dynamic_cast<NonTVertex*>(otherVVertex);
                    if(otherNTVertex != NULL && otherNTVertex->getNature() & Nature::CUSP){
                        assert(otherNTVertex->viewedges().size()==2);

                        real arcLength = ArcLength2D(edge,FLT_MAX);
                        if(arcLength<=1.0){
                            ViewEdge* otherEdge = (otherNTVertex->viewedges()[0].first != edge ? otherNTVertex->viewedges()[0].first : otherNTVertex->viewedges()[1].first);
                            if(edge->qi()!=0 && ntVertex->viewedges()[(i+1)%2].first->qi()==0 && otherEdge->qi()==0){
                                edge->setQI(0);
                                edge->markInconsistent();
                            }
                        }
                    }
                }
            }
        }
    }

    if (_cuspTrimThreshold > 0)
    {
        printf("Hiding small loops and dead ends\n");
        bool changed1 = false;
        bool changed2 = false;
        do
        {
            changed1 = HideSmallBits(_ViewMap);
            changed2 = HideSmallLoopsAndDeadEnds(_ViewMap);
        } while (changed1 || changed2);
    }
}

inline int NumVisibleEdgesTV(TVertex * tv)
{
    int vis = 0;

    for(int i=0;i<tv->numEdges();i++)
        if ( tv->getEdge(i)->first->qi() == 0)
            vis ++;
    return vis;
}

inline int NumVisibleEdgesNTV(NonTVertex * ntv)
{
    int vis =0;
    for(int i=0;i<ntv->viewedges().size();i++)
        if (ntv->viewedges()[i].first->qi() == 0)
            vis ++;

    return vis;
}

inline int NumVisibleEdges(ViewVertex * vv)
{
    TVertex * tv = dynamic_cast<TVertex*>(vv);

    if (tv != NULL)
        return NumVisibleEdgesTV(tv);

    assert(dynamic_cast<NonTVertex*>(vv) != NULL);
    return NumVisibleEdgesNTV((NonTVertex*)vv);
}

bool ViewMapBuilder::HideSmallBits(ViewMap * ioViewMap)
{
    bool changed = false;

    for(vector<ViewVertex*>::iterator vit = ioViewMap->ViewVertices().begin(); vit != ioViewMap->ViewVertices().end(); ++vit)
    {
        ViewVertex * startVertex = *vit;
        TVertex* startTVertex = dynamic_cast<TVertex*>(startVertex);
        if(!startTVertex)
            continue;

        int numVisibleEdges = NumVisibleEdges(startVertex);
        if(startTVertex->numEdges()==4 && numVisibleEdges == 3){
            ViewEdge* startEdge = NULL;
            ViewVertex* endVertex = NULL;
            ViewEdge* otherEdges[2];
            if(startTVertex->backEdgeA().first->qi()==0 && startTVertex->backEdgeB().first->qi()!=0){
                startEdge = startTVertex->backEdgeA().first;
                endVertex = startEdge->A();
                otherEdges[0] = startTVertex->frontEdgeA().first;
                otherEdges[1] = startTVertex->frontEdgeB().first;
            }else if(startTVertex->backEdgeA().first->qi()!=0 && startTVertex->backEdgeB().first->qi()==0){
                startEdge = startTVertex->backEdgeB().first;
                endVertex = startEdge->B();
                otherEdges[0] = startTVertex->frontEdgeA().first;
                otherEdges[1] = startTVertex->frontEdgeB().first;
            }else if(startTVertex->frontEdgeA().first->qi()==0 && startTVertex->frontEdgeB().first->qi()!=0){
                startEdge = startTVertex->frontEdgeA().first;
                endVertex = startEdge->A();
                otherEdges[0] = startTVertex->backEdgeA().first;
                otherEdges[1] = startTVertex->backEdgeB().first;
            }else if(startTVertex->frontEdgeA().first->qi()!=0 && startTVertex->frontEdgeB().first->qi()==0){
                startEdge = startTVertex->frontEdgeB().first;
                endVertex = startEdge->B();
                otherEdges[0] = startTVertex->backEdgeA().first;
                otherEdges[1] = startTVertex->backEdgeB().first;
            }
            assert(startEdge && endVertex);
            if(startVertex == endVertex)
                continue;

            real arclength = ArcLength2D(startEdge,2*_cuspTrimThreshold);
            if(arclength < _cuspTrimThreshold){
                TVertex* endTVertex = dynamic_cast<TVertex*>(endVertex);

                if(!endTVertex){
                    // NonTVertex case
                    NonTVertex* endNTVertex = dynamic_cast<NonTVertex*>(endVertex);
                    if(endNTVertex->viewedges().size()==3){
                        bool found = false;
                        ViewEdge* thirdEdge = NULL;
                        for(int i=0;i<3;i++){
                            ViewEdge* e = endNTVertex->viewedges()[i].first;
                            found = found || (e == otherEdges[0]) || (e == otherEdges[1]);
                            if(e!=otherEdges[0] && e!=otherEdges[1] && e!=startEdge)
                                thirdEdge = e;
                        }
                        if(found && thirdEdge && thirdEdge->qi()==0){
                            printf("HIDE SMALL BITS %08X\n",startEdge);
                            startEdge->setQI(1339);
                            changed = true;
                            continue;
                        }
                    }
                    continue;
                }

                // TVertex case
                numVisibleEdges = NumVisibleEdges(endVertex);

                ViewEdge* midEdge = NULL;
                TVertex* midVertex = NULL;
                if(endTVertex->numEdges()==4 && (numVisibleEdges == 2 || numVisibleEdges == 4)){
                    midVertex = endTVertex;
                    if(endTVertex->frontEdgeA().first == startEdge && endTVertex->frontEdgeB().first->qi() == 0){
                        midEdge = startEdge;
                        startEdge = endTVertex->frontEdgeB().first;
                        endVertex = startEdge->B();
                    }else if(endTVertex->frontEdgeB().first == startEdge && endTVertex->frontEdgeA().first->qi() == 0){
                        midEdge = startEdge;
                        startEdge = endTVertex->frontEdgeA().first;
                        endVertex = startEdge->A();
                    }else if(endTVertex->backEdgeA().first == startEdge && endTVertex->backEdgeB().first->qi() == 0){
                        midEdge = startEdge;
                        startEdge = endTVertex->backEdgeB().first;
                        endVertex = startEdge->B();
                    }else if(endTVertex->backEdgeB().first == startEdge && endTVertex->backEdgeA().first->qi() == 0){
                        midEdge = startEdge;
                        startEdge = endTVertex->backEdgeA().first;
                        endVertex = startEdge->A();
                    }
                }
                if(midEdge){
                    endTVertex = dynamic_cast<TVertex*>(endVertex);
                    if(!endTVertex)
                        continue;
                    numVisibleEdges = NumVisibleEdges(endVertex);
                }

                if(endTVertex->numEdges()==4 && numVisibleEdges == 3){
                    bool valid = false;

                    //   o-o----o
                    //   |      |   Detect these kind of cases (one invisible TVertex allowed)x
                    // --o---o--o---
                    if((endTVertex->backEdgeA().first==startEdge && endTVertex->backEdgeB().first->qi()!=0)
                            || (endTVertex->backEdgeB().first==startEdge && endTVertex->backEdgeA().first->qi()!=0)){
                        valid = (endTVertex->frontEdgeA().first == otherEdges[0]) || (endTVertex->frontEdgeA().first == otherEdges[1])
                                || (endTVertex->frontEdgeB().first == otherEdges[0]) || (endTVertex->frontEdgeB().first == otherEdges[1]);
                        if(!valid){
                            TVertex* midV = dynamic_cast<TVertex*>(endTVertex->frontEdgeA().first->A());
                            if(midV && (midV == otherEdges[0]->A() || midV == otherEdges[1]->B())){
                                if(midV == midVertex){
                                    valid = true;
                                }else if(midV->frontEdgeA().first == endTVertex->frontEdgeA().first ||
                                         midV->frontEdgeB().first == endTVertex->frontEdgeA().first){
                                    valid = ((!midV->backEdgeA().first || midV->backEdgeA().first->qi()!=0)
                                             && (!midV->backEdgeB().first || midV->backEdgeB().first->qi()!=0));
                                }else{
                                    valid = ((!midV->frontEdgeA().first || midV->frontEdgeA().first->qi()!=0)
                                             && (!midV->frontEdgeB().first || midV->frontEdgeB().first->qi()!=0));
                                }
                            }
                        }
                        if(!valid){
                            TVertex* midV = dynamic_cast<TVertex*>(endTVertex->frontEdgeB().first->B());
                            if(midV && (midV == otherEdges[0]->A() || midV == otherEdges[1]->B())){
                                if(midV == midVertex){
                                    valid = true;
                                }else if(midV->frontEdgeA().first == endTVertex->frontEdgeB().first ||
                                         midV->frontEdgeB().first == endTVertex->frontEdgeB().first){
                                    valid = ((!midV->backEdgeA().first || midV->backEdgeA().first->qi()!=0)
                                             && (!midV->backEdgeB().first || midV->backEdgeB().first->qi()!=0));
                                }else{
                                    valid = ((!midV->frontEdgeA().first || midV->frontEdgeA().first->qi()!=0)
                                             && (!midV->frontEdgeB().first || midV->frontEdgeB().first->qi()!=0));
                                }
                            }
                        }
                    }
                    if(!valid && (endTVertex->frontEdgeA().first==startEdge && endTVertex->frontEdgeB().first->qi()!=0)
                            || (endTVertex->frontEdgeB().first==startEdge && endTVertex->frontEdgeA().first->qi()!=0)){
                        valid = (endTVertex->backEdgeA().first == otherEdges[0]) || (endTVertex->backEdgeA().first == otherEdges[1])
                                || (endTVertex->backEdgeB().first == otherEdges[0]) || (endTVertex->backEdgeB().first == otherEdges[1]);
                        if(!valid){
                            TVertex* midV = dynamic_cast<TVertex*>(endTVertex->backEdgeA().first->A());
                            if(midV && (midV == otherEdges[0]->A() || midV == otherEdges[1]->B())){
                                if(midV == midVertex){
                                    valid = true;
                                }else if(midV->frontEdgeA().first == endTVertex->backEdgeA().first ||
                                         midV->frontEdgeB().first == endTVertex->backEdgeA().first){
                                    valid = ((!midV->backEdgeA().first || midV->backEdgeA().first->qi()!=0)
                                             && (!midV->backEdgeB().first || midV->backEdgeB().first->qi()!=0));
                                }else{
                                    valid = ((!midV->frontEdgeA().first || midV->frontEdgeA().first->qi()!=0)
                                             && (!midV->frontEdgeB().first || midV->frontEdgeB().first->qi()!=0));
                                }
                            }
                        }
                        if(!valid){
                            TVertex* midV = dynamic_cast<TVertex*>(endTVertex->backEdgeB().first->B());
                            if(midV && (midV == otherEdges[0]->A() || midV == otherEdges[1]->B())){
                                if(midV == midVertex){
                                    valid = true;
                                }else if(midV->frontEdgeA().first == endTVertex->backEdgeB().first ||
                                         midV->frontEdgeB().first == endTVertex->backEdgeB().first){
                                    valid = ((!midV->backEdgeA().first || midV->backEdgeA().first->qi()!=0)
                                             && (!midV->backEdgeB().first || midV->backEdgeB().first->qi()!=0));
                                }else{
                                    valid = ((!midV->frontEdgeA().first || midV->frontEdgeA().first->qi()!=0)
                                             && (!midV->frontEdgeB().first || midV->frontEdgeB().first->qi()!=0));
                                }
                            }
                        }
                    }

                    //   o-o---o--
                    //   |     |   Detect these kind of cases (one invisible TVertex allowed)
                    // --o---o-o
                    if(!valid && (endTVertex->backEdgeA().first==startEdge || endTVertex->backEdgeB().first==startEdge)){
                        if(endTVertex->frontEdgeA().first->qi()==0 && endTVertex->frontEdgeB().first->qi()!=0)
                            valid = (endTVertex->frontEdgeA().first == otherEdges[0]) || (endTVertex->frontEdgeA().first == otherEdges[1]);
                        else if(endTVertex->frontEdgeA().first->qi()!=0 && endTVertex->frontEdgeB().first->qi()==0)
                            valid = (endTVertex->frontEdgeB().first == otherEdges[0]) || (endTVertex->frontEdgeB().first == otherEdges[1]);
                        if(!valid){
                            TVertex* midV = dynamic_cast<TVertex*>(endTVertex->frontEdgeA().first->A());
                            if(midV && (midV == otherEdges[0]->A() || midV == otherEdges[1]->B())){
                                if(midV == midVertex){
                                    valid = true;
                                }else if(midV->frontEdgeA().first == endTVertex->frontEdgeA().first ||
                                         midV->frontEdgeB().first == endTVertex->frontEdgeA().first){
                                    valid = ((!midV->backEdgeA().first || midV->backEdgeA().first->qi()!=0)
                                             && (!midV->backEdgeB().first || midV->backEdgeB().first->qi()!=0));
                                }else{
                                    valid = ((!midV->frontEdgeA().first || midV->frontEdgeA().first->qi()!=0)
                                             && (!midV->frontEdgeB().first || midV->frontEdgeB().first->qi()!=0));
                                }
                            }
                        }
                        if(!valid){
                            TVertex* midV = dynamic_cast<TVertex*>(endTVertex->frontEdgeB().first->B());
                            if(midV && (midV == otherEdges[0]->A() || midV == otherEdges[1]->B())){
                                if(midV == midVertex){
                                    valid = true;
                                }else if(midV->frontEdgeA().first == endTVertex->frontEdgeB().first ||
                                         midV->frontEdgeB().first == endTVertex->frontEdgeB().first){
                                    valid = ((!midV->backEdgeA().first || midV->backEdgeA().first->qi()!=0)
                                             && (!midV->backEdgeB().first || midV->backEdgeB().first->qi()!=0));
                                }else{
                                    valid = ((!midV->frontEdgeA().first || midV->frontEdgeA().first->qi()!=0)
                                             && (!midV->frontEdgeB().first || midV->frontEdgeB().first->qi()!=0));
                                }
                            }
                        }
                    }
                    if(!valid && (endTVertex->frontEdgeA().first==startEdge || endTVertex->frontEdgeB().first==startEdge)){
                        if(endTVertex->backEdgeA().first->qi()==0 && endTVertex->backEdgeB().first->qi()!=0)
                            valid = (endTVertex->backEdgeA().first == otherEdges[0]) || (endTVertex->backEdgeA().first == otherEdges[1]);
                        else if(endTVertex->backEdgeA().first->qi()!=0 && endTVertex->backEdgeB().first->qi()==0)
                            valid = (endTVertex->backEdgeB().first == otherEdges[0]) || (endTVertex->backEdgeB().first == otherEdges[1]);
                        if(!valid){
                            TVertex* midV = dynamic_cast<TVertex*>(endTVertex->backEdgeA().first->A());
                            if(midV && (midV == otherEdges[0]->A() || midV == otherEdges[1]->B())){
                                if(midV == midVertex){
                                    valid = true;
                                }else if(midV->frontEdgeA().first == endTVertex->backEdgeA().first ||
                                         midV->frontEdgeB().first == endTVertex->backEdgeA().first){
                                    valid = ((!midV->backEdgeA().first || midV->backEdgeA().first->qi()!=0)
                                             && (!midV->backEdgeB().first || midV->backEdgeB().first->qi()!=0));
                                }else{
                                    valid = ((!midV->frontEdgeA().first || midV->frontEdgeA().first->qi()!=0)
                                             && (!midV->frontEdgeB().first || midV->frontEdgeB().first->qi()!=0));
                                }
                            }
                        }
                        if(!valid){
                            TVertex* midV = dynamic_cast<TVertex*>(endTVertex->backEdgeB().first->B());
                            if(midV && (midV == otherEdges[0]->A() || midV == otherEdges[1]->B())){
                                if(midV == midVertex){
                                    valid = true;
                                }else if(midV->frontEdgeA().first == endTVertex->backEdgeB().first ||
                                         midV->frontEdgeB().first == endTVertex->backEdgeB().first){
                                    valid = ((!midV->backEdgeA().first || midV->backEdgeA().first->qi()!=0)
                                             && (!midV->backEdgeB().first || midV->backEdgeB().first->qi()!=0));
                                }else{
                                    valid = ((!midV->frontEdgeA().first || midV->frontEdgeA().first->qi()!=0)
                                             && (!midV->frontEdgeB().first || midV->frontEdgeB().first->qi()!=0));
                                }
                            }
                        }
                    }

                    if(valid){
                        printf("HIDE SMALL BITS %08X\n",startEdge);
                        startEdge->setQI(1339);
                        if(midEdge){
                            printf("HIDE SMALL BITS %08X\n",midEdge);
                            midEdge->setQI(1339);
                        }
                        changed = true;
                        continue;
                    }
                }
            }
        }
    }
    return changed;
}

inline int NumEdges(ViewVertex * vv)
{
    TVertex * tv = dynamic_cast<TVertex*>(vv);
    if (tv != NULL)
        return tv->numEdges();
    return  ((NonTVertex*)vv)->viewedges().size();
}

inline ViewEdge * GetEdge(ViewVertex * vv, int i)
{
    TVertex * tv = dynamic_cast<TVertex*>(vv);
    if (tv != NULL)
        return tv->getEdge(i)->first;
    return ((NonTVertex*)vv)->viewedges()[i].first;
}

ViewEdge * AdvanceAlongVertex2(ViewVertex * nextVertex, ViewEdge * lastEdge)
{
    NonTVertex * ntv = dynamic_cast<NonTVertex*>(nextVertex);
    if (ntv != NULL)
    {
        for(int i=0;i<ntv->viewedges().size();i++)
            if (ntv->viewedges()[i].first != lastEdge && ntv->viewedges()[i].first->qi() == 0)
                return ntv->viewedges()[i].first;

        return NULL;
    }

    TVertex * tv = (TVertex*)nextVertex;
    for(int i=0;i<tv->numEdges();i++)
        if (tv->getEdge(i)->first != lastEdge && tv->getEdge(i)->first->qi() == 0)
            return tv->getEdge(i)->first;

    return NULL;
}

// let:
//  * dead-end be a vertex with 1 visible outgoing edge
//  * junction be a vertex > 2 visible outgoing edges

// hide curve segments that:
// * connect a junction to a dead-end,
// * connect a dead-end to a dead-end
// * connect any vertex to itself

// ONLY along segments that
// * contain no other junctions,
// * no invisible edges, and
// * are within the length threshold
//
// some of the code below only looks at t-vertices, not general junctions
bool ViewMapBuilder::HideSmallLoopsAndDeadEnds(ViewMap * ioViewMap)
{
    bool changed = false;

    for(vector<ViewVertex*>::iterator vit = ioViewMap->ViewVertices().begin(); vit != ioViewMap->ViewVertices().end(); ++vit)
    {
        ViewVertex * startVertex = *vit;

        for(int i=0;i<NumEdges(startVertex);i++)
        {
            ViewEdge * startEdge = GetEdge(startVertex, i);

            if (startEdge->qi() != 0)
                continue;

            // hide single-edge loop (not sure if this case needs special treatment)
            if (startEdge->A() == startVertex && startEdge->B() == startVertex && ArcLength2D(startEdge, 2*_cuspTrimThreshold))
            {
                startEdge->setQI(1338);
                changed = true;
                continue;
            }

            // trace from the start edge
            set<ViewEdge*> visitedEdges;
            ViewVertex * lastVertex = startVertex;
            ViewEdge * currentEdge = startEdge;
            real arcLength = 0;

            bool foundSegmentToDelete = false;

            arcLength += ArcLength2D(currentEdge, 2*_cuspTrimThreshold);
            while(arcLength < _cuspTrimThreshold)
            {
                if (currentEdge->qi() != 0)
                    break;

                visitedEdges.insert(currentEdge);

                ViewVertex * nextVertex = currentEdge->A() != lastVertex ? currentEdge->A() : currentEdge->B();

                assert(nextVertex != NULL);

                if (nextVertex == startVertex) // check if we found a small loop
                {
                    foundSegmentToDelete = true;
                    break;
                }

                int numVis = NumVisibleEdges(nextVertex);
                assert(numVis > 0);
                if (numVis == 1) // check if we found a dead-end
                {
                    if (NumVisibleEdges(startVertex) != 2) // delete dead-end to dead-end, or dead-end to junction
                        foundSegmentToDelete = true;
                    break;
                }

                if (numVis > 2) // check if found a junction
                {
                    if (NumVisibleEdges(startVertex) == 1)  // delete junction to dead-end
                        foundSegmentToDelete = true;
                    break;
                }

                // continue along the next viewedge
                currentEdge = AdvanceAlongVertex2(nextVertex, currentEdge);
                assert(currentEdge != NULL);

                if (visitedEdges.find(currentEdge) != visitedEdges.end())
                {
                    break;
                }

                lastVertex = nextVertex;
                arcLength += ArcLength2D(currentEdge, 2*_cuspTrimThreshold);
            }

            if (foundSegmentToDelete)
            {
                for(set<ViewEdge*>::iterator it = visitedEdges.begin(); it != visitedEdges.end(); ++it)
                {
                    if ( (*it)->qi() == 0)
                        changed = true;
                    (*it)->setQI(1337);
                }
            }
        }
    }

    return changed;
}

static const unsigned gProgressBarMaxSteps = 10;
static const unsigned gProgressBarMinSize = 2000;

void ViewMapBuilder::ComputeRayCastingVisibility(ViewMap *ioViewMap, real epsilon)
{
    vector<ViewEdge*>& vedges = ioViewMap->ViewEdges();
    bool progressBarDisplay = false;
    unsigned progressBarStep = 0;
    unsigned vEdgesSize = vedges.size();
    unsigned fEdgesSize = ioViewMap->FEdges().size();

    if (_pProgressBar != NULL && fEdgesSize > gProgressBarMinSize) {
        unsigned progressBarSteps = min(gProgressBarMaxSteps, vEdgesSize);
        progressBarStep = vEdgesSize / progressBarSteps;
        _pProgressBar->reset();
        _pProgressBar->setLabelText("Computing Ray casting Visibility");
        _pProgressBar->setTotalSteps(progressBarSteps);
        _pProgressBar->setProgress(0);
        progressBarDisplay = true;
    }

    unsigned counter = progressBarStep;
    FEdge *fe, *festart;
    int nSamples = 0;
    vector<Polygon3r*> aFaces;
    Polygon3r *aFace = NULL;
    unsigned tmpQI = 0;
    unsigned qiClasses[256];
    unsigned maxIndex, maxCard;
    unsigned qiMajority;
    static unsigned timestamp = 1;
    for (vector<ViewEdge*>::iterator ve = vedges.begin(), veend = vedges.end(); ve != veend; ve++) {
        int visVotes = 0;
        int invisVotes = 0;

        if (_pRenderMonitor && _pRenderMonitor->testBreak())
            break;
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "Processing ViewEdge " << (*ve)->getId() << endl;
        }
#endif
        festart = (*ve)->fedgeA();
        fe = (*ve)->fedgeA();
        qiMajority = 1;
        do {
            qiMajority++;
            fe = fe->nextEdge();
        } while (fe && fe != festart);
        qiMajority >>= 1;
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tqiMajority: " << qiMajority << endl;
        }
#endif

        tmpQI = 0;
        maxIndex = 0;
        maxCard = 0;
        nSamples = 0;
        fe = (*ve)->fedgeA();
        memset(qiClasses, 0, 256 * sizeof(*qiClasses));
        set<ViewShape*> occluders;
        do {
            if ((maxCard < qiMajority)) {
                tmpQI = ComputeRayCastingVisibility(fe, _Grid, epsilon, occluders, &aFace, timestamp++);

#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFEdge: visibility " << tmpQI << endl;
                }
#endif
                //ARB: This is an error condition, not an alert condition.
                // Some sort of recovery or abort is necessary.
                if (tmpQI >= 256) {
                    cerr << "Warning: too many occluding levels" << endl;
                    //ARB: Wild guess: instead of aborting or corrupting memory, treat as tmpQI == 255
                    tmpQI = 255;
                }

                if (++qiClasses[tmpQI] > maxCard) {
                    maxCard = qiClasses[tmpQI];
                    maxIndex = tmpQI;
                }

                if (tmpQI == 0)
                    visVotes ++;
                else
                    invisVotes ++;
            }
            else {
                //ARB: FindOccludee is redundant if ComputeRayCastingVisibility has been called
                FindOccludee(fe, _Grid, epsilon, &aFace, timestamp++);
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFEdge: occludee only (" << (aFace != NULL ? "found" : "not found") << ")" << endl;
                }
#endif
            }

            if (aFace) {
                fe->setaFace(*aFace);
                aFaces.push_back(aFace);
                fe->setOccludeeEmpty(false);
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\tFound occludee" << endl;
                }
#endif
            }
            else {
                //ARB: We are arbitrarily using the last observed value for occludee (almost always the value observed
                //     for the edge before festart). Is that meaningful?
                // ...in fact, _occludeeEmpty seems to be unused.
                fe->setOccludeeEmpty(true);
            }

            ++nSamples;
            fe = fe->nextEdge();
        } while ((maxCard < qiMajority) && (fe) && (fe != festart));
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tFinished with " << nSamples << " samples, maxCard = " << maxCard << endl;
        }
#endif

        // I don't care about estimating QI.
        if (visVotes > invisVotes)
            (*ve)->setQI(0);
        else
            (*ve)->setQI(100);

        if (invisVotes > 0 && visVotes > 0)
            (*ve)->markInconsistent();

        if (invisVotes == 0 && visVotes == 0){
            (*ve)->markAmbiguous();
        }

        // ViewEdge
        // qi --
        (*ve)->setQI(maxIndex);
        // occluders --
        for (set<ViewShape*>::iterator o = occluders.begin(), oend = occluders.end(); o != oend; ++o)
            (*ve)->AddOccluder((*o));
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\tConclusion: QI = " << maxIndex << ", " << (*ve)->occluders_size() << " occluders." << endl;
        }
#endif
        // occludee --
        if (!aFaces.empty()) {
            if (aFaces.size() <= (float)nSamples / 2.0f) {
                (*ve)->setaShape(0);
            }
            else {
                vector<Polygon3r*>::iterator p = aFaces.begin();
                WFace *wface = (WFace *)((*p)->userdata);
                ViewShape *vshape = ioViewMap->viewShape(wface->GetVertex(0)->shape()->GetId());
                ++p;
                (*ve)->setaShape(vshape);
            }
        }

        if (progressBarDisplay) {
            counter--;
            if (counter <= 0) {
                counter = progressBarStep;
                _pProgressBar->setProgress(_pProgressBar->getProgress() + 1);
            }
        }
        aFaces.clear();
    }
}

void ViewMapBuilder::ComputeFastRayCastingVisibility(ViewMap *ioViewMap, real epsilon)
{
    vector<ViewEdge*>& vedges = ioViewMap->ViewEdges();
    bool progressBarDisplay = false;
    unsigned progressBarStep = 0;
    unsigned vEdgesSize = vedges.size();
    unsigned fEdgesSize = ioViewMap->FEdges().size();

    if (_pProgressBar != NULL && fEdgesSize > gProgressBarMinSize) {
        unsigned progressBarSteps = min(gProgressBarMaxSteps, vEdgesSize);
        progressBarStep = vEdgesSize / progressBarSteps;
        _pProgressBar->reset();
        _pProgressBar->setLabelText("Computing Ray casting Visibility");
        _pProgressBar->setTotalSteps(progressBarSteps);
        _pProgressBar->setProgress(0);
        progressBarDisplay = true;
    }

    unsigned counter = progressBarStep;
    FEdge *fe, *festart;
    unsigned nSamples = 0;
    vector<Polygon3r*> aFaces;
    Polygon3r *aFace = NULL;
    unsigned tmpQI = 0;
    unsigned qiClasses[256];
    unsigned maxIndex, maxCard;
    unsigned qiMajority;
    static unsigned timestamp = 1;
    bool even_test;
    for (vector<ViewEdge*>::iterator ve = vedges.begin(), veend = vedges.end(); ve != veend; ve++) {
        if (_pRenderMonitor && _pRenderMonitor->testBreak())
            break;

        festart = (*ve)->fedgeA();
        fe = (*ve)->fedgeA();
        qiMajority = 1;
        do {
            qiMajority++;
            fe = fe->nextEdge();
        } while (fe && fe != festart);
        if (qiMajority >= 4)
            qiMajority >>= 2;
        else
            qiMajority = 1;

        set<ViewShape*> occluders;

        even_test = true;
        maxIndex = 0;
        maxCard = 0;
        nSamples = 0;
        memset(qiClasses, 0, 256 * sizeof(*qiClasses));
        fe = (*ve)->fedgeA();
        do {
            if (even_test) {
                if ((maxCard < qiMajority)) {
                    tmpQI = ComputeRayCastingVisibility(fe, _Grid, epsilon, occluders, &aFace, timestamp++);

                    //ARB: This is an error condition, not an alert condition.
                    // Some sort of recovery or abort is necessary.
                    if (tmpQI >= 256) {
                        cerr << "Warning: too many occluding levels" << endl;
                        //ARB: Wild guess: instead of aborting or corrupting memory, treat as tmpQI == 255
                        tmpQI = 255;
                    }

                    if (++qiClasses[tmpQI] > maxCard) {
                        maxCard = qiClasses[tmpQI];
                        maxIndex = tmpQI;
                    }
                }
                else {
                    //ARB: FindOccludee is redundant if ComputeRayCastingVisibility has been called
                    FindOccludee(fe, _Grid, epsilon, &aFace, timestamp++);
                }

                if (aFace) {
                    fe->setaFace(*aFace);
                    aFaces.push_back(aFace);
                }
                ++nSamples;
                even_test = false;
            }
            else {
                even_test = true;
            }
            fe = fe->nextEdge();
        } while ((maxCard < qiMajority) && (fe) && (fe != festart));

        (*ve)->setQI(maxIndex);

        if (!aFaces.empty()) {
            if (aFaces.size() < nSamples / 2) {
                (*ve)->setaShape(0);
            }
            else {
                vector<Polygon3r*>::iterator p = aFaces.begin();
                WFace *wface = (WFace *)((*p)->userdata);
                ViewShape *vshape = ioViewMap->viewShape(wface->GetVertex(0)->shape()->GetId());
                ++p;
#if 0
                for (; p != pend; ++p) {
                    WFace *f = (WFace*)((*p)->userdata);
                    ViewShape *vs = ioViewMap->viewShape(f->GetVertex(0)->shape()->GetId());
                    if (vs != vshape) {
                        sameShape = false;
                        break;
                    }
                }
                if (sameShape)
#endif
                    (*ve)->setaShape(vshape);
            }
        }

        //(*ve)->setaFace(aFace);

        if (progressBarDisplay) {
            counter--;
            if (counter <= 0) {
                counter = progressBarStep;
                _pProgressBar->setProgress(_pProgressBar->getProgress() + 1);
            }
        }
        aFaces.clear();
    }
}

void ViewMapBuilder::ComputeVeryFastRayCastingVisibility(ViewMap *ioViewMap, real epsilon)
{
    vector<ViewEdge*>& vedges = ioViewMap->ViewEdges();
    bool progressBarDisplay = false;
    unsigned progressBarStep = 0;
    unsigned vEdgesSize = vedges.size();
    unsigned fEdgesSize = ioViewMap->FEdges().size();

    if (_pProgressBar != NULL && fEdgesSize > gProgressBarMinSize) {
        unsigned progressBarSteps = min(gProgressBarMaxSteps, vEdgesSize);
        progressBarStep = vEdgesSize / progressBarSteps;
        _pProgressBar->reset();
        _pProgressBar->setLabelText("Computing Ray casting Visibility");
        _pProgressBar->setTotalSteps(progressBarSteps);
        _pProgressBar->setProgress(0);
        progressBarDisplay = true;
    }

    unsigned counter = progressBarStep;
    FEdge *fe;
    unsigned qi = 0;
    Polygon3r *aFace = NULL;
    static unsigned timestamp = 1;
    for (vector<ViewEdge*>::iterator ve = vedges.begin(), veend = vedges.end(); ve != veend; ve++) {
        if (_pRenderMonitor && _pRenderMonitor->testBreak())
            break;

        set<ViewShape*> occluders;

        fe = (*ve)->fedgeA();
        qi = ComputeRayCastingVisibility(fe, _Grid, epsilon, occluders, &aFace, timestamp++);
        if (aFace) {
            fe->setaFace(*aFace);
            WFace *wface = (WFace *)(aFace->userdata);
            ViewShape *vshape = ioViewMap->viewShape(wface->GetVertex(0)->shape()->GetId());
            (*ve)->setaShape(vshape);
        }
        else {
            (*ve)->setaShape(0);
        }

        (*ve)->setQI(qi);

        if (progressBarDisplay) {
            counter--;
            if (counter <= 0) {
                counter = progressBarStep;
                _pProgressBar->setProgress(_pProgressBar->getProgress() + 1);
            }
        }
    }
}

void ViewMapBuilder::FindOccludee(FEdge *fe, Grid *iGrid, real epsilon, Polygon3r **oaPolygon, unsigned timestamp,
                                  Vec3r& u, Vec3r& A, Vec3r& origin, Vec3r& edge, vector<WVertex*>& faceVertices)
{
    WFace *face = NULL;
    if (fe->isSmooth()) {
        FEdgeSmooth *fes = dynamic_cast<FEdgeSmooth*>(fe);
        face = (WFace *)fes->face();
    }
    OccludersSet occluders;
    WFace *oface;
    bool skipFace;

    WVertex::incoming_edge_iterator ie;
    OccludersSet::iterator p, pend;

    *oaPolygon = NULL;
    if (((fe)->getNature() & Nature::SILHOUETTE) || ((fe)->getNature() & Nature::BORDER)) {
        occluders.clear();
        // we cast a ray from A in the same direction but looking behind
        Vec3r v(-u[0], -u[1], -u[2]);
        iGrid->castInfiniteRay(A, v, occluders, timestamp);

        bool noIntersection = true;
        real mint = FLT_MAX;
        // we met some occluders, let us fill the aShape field with the first intersected occluder
        for (p = occluders.begin(), pend = occluders.end(); p != pend; p++) {
            // check whether the edge and the polygon plane are coincident:
            //-------------------------------------------------------------
            //first let us compute the plane equation.
            oface = (WFace *)(*p)->userdata;
            Vec3r v1(((*p)->getVertices())[0]);
            Vec3r normal((*p)->getNormal());
            real d = -(v1 * normal);
            real t, t_u, t_v;

            if (face) {
                skipFace = false;

                if (face == oface)
                    continue;

                if (faceVertices.empty())
                    continue;

                for (vector<WVertex*>::iterator fv = faceVertices.begin(), fvend = faceVertices.end();
                     fv != fvend;
                     ++fv)
                {
                    if ((*fv)->isBoundary())
                        continue;
                    WVertex::incoming_edge_iterator iebegin = (*fv)->incoming_edges_begin();
                    WVertex::incoming_edge_iterator ieend = (*fv)->incoming_edges_end();
                    for (ie = iebegin; ie != ieend; ++ie) {
                        if ((*ie) == 0)
                            continue;

                        WFace *sface = (*ie)->GetbFace();
                        if (sface == oface) {
                            skipFace = true;
                            break;
                        }
                    }
                    if (skipFace)
                        break;
                }
                if (skipFace)
                    continue;
            }
            else {
                if (GeomUtils::COINCIDENT == GeomUtils::intersectRayPlane(origin, edge, normal, d, t, epsilon))
                    continue;
            }
            if ((*p)->rayIntersect(A, v, t, t_u, t_v)) {
                if (fabs(v * normal) > 0.0001) {
                    if (t > 0.0) { // && t < 1.0) {
                        if (t < mint) {
                            *oaPolygon = (*p);
                            mint = t;
                            noIntersection = false;
                            fe->setOccludeeIntersection(Vec3r(A + t * v));
                        }
                    }
                }
            }
        }

        if (noIntersection)
            *oaPolygon = NULL;
    }
}

void ViewMapBuilder::FindOccludee(FEdge *fe, Grid *iGrid, real epsilon, Polygon3r **oaPolygon, unsigned timestamp)
{
    OccludersSet occluders;

    Vec3r A;
    Vec3r edge;
    Vec3r origin;
    A = Vec3r(((fe)->vertexA()->point3D() + (fe)->vertexB()->point3D()) / 2.0);
    edge = Vec3r((fe)->vertexB()->point3D() - (fe)->vertexA()->point3D());
    origin = Vec3r((fe)->vertexA()->point3D());
    Vec3r u;
    if (_orthographicProjection) {
        u = Vec3r(0.0, 0.0, _viewpoint.z() - A.z());
    }
    else {
        u = Vec3r(_viewpoint - A);
    }
    u.normalize();
    if (A < iGrid->getOrigin())
        cerr << "Warning: point is out of the grid for fedge " << fe->getId().getFirst() << "-" <<
                fe->getId().getSecond() << endl;

    vector<WVertex*> faceVertices;

    WFace *face = NULL;
    if (fe->isSmooth()) {
        FEdgeSmooth *fes = dynamic_cast<FEdgeSmooth*>(fe);
        face = (WFace *)fes->face();
    }
    if (face)
        face->RetrieveVertexList(faceVertices);

    return FindOccludee(fe, iGrid, epsilon, oaPolygon, timestamp, u, A, origin, edge, faceVertices);
}

int ViewMapBuilder::ComputeRayCastingVisibility(FEdge *fe, Grid *iGrid, real epsilon, set<ViewShape*>& oOccluders,
                                                Polygon3r **oaPolygon, unsigned timestamp)
{
    OccludersSet occluders;
    int qi = 0;

    Vec3r center;
    Vec3r edge;
    Vec3r origin;

    WXFace * face1 = dynamic_cast<WXFace*>(fe->getFace1());
    WXFace * face2 = dynamic_cast<WXFace*>(fe->getFace2());

//    if (fe->isSmooth())
//    {
//        // check if this edge is on a back-face
//        if (!(fe->getNature() & Nature::SILHOUETTE))
//        {
//            // may have two source faces
//            if (face1 != NULL && !face1->front(_useConsistency))
//                return 101;
//            if (face2 != NULL && !face2->front(_useConsistency))
//                return 101;
//        }

//        if (face1 != NULL && _useConsistency && !face1->consistent())
//            return -1;

//        if (face2 != NULL && _useConsistency && !face2->consistent())
//            return -1;
//    }
//    else
//    {

//        // markosian test
//        if (fe->getNature() & Nature::SILHOUETTE)
//        {
//            assert(dynamic_cast<FEdgeSharp*>(fe) != NULL);
//            WEdge * edge = dynamic_cast<FEdgeSharp*>(fe)->edge();
//            assert(edge != NULL);
//            real discriminant;
//            WXFace * nearFace = (WXFace*)GetNearFace(edge, _viewpoint, &discriminant);

//            if ((discriminant > 1 || discriminant < -1) && !nearFace->front(false))
//                return 123;
//        }
//    }

    center = fe->center3d();
    edge = Vec3r(fe->vertexB()->point3D() - fe->vertexA()->point3D());
    origin = Vec3r(fe->vertexA()->point3D());
    // Is the edge outside the view frustum ?
    Vec3r gridOrigin(iGrid->getOrigin());
    Vec3r gridExtremity(iGrid->getOrigin() + iGrid->gridSize());

    if ((center.x() < gridOrigin.x()) || (center.y() < gridOrigin.y()) || (center.z() < gridOrigin.z()) ||
            (center.x() > gridExtremity.x()) || (center.y() > gridExtremity.y()) || (center.z() > gridExtremity.z()))
    {
        cerr << "Warning: point is out of the grid for fedge " << fe->getId() << endl;
        //return 0;
    }

#if 0
    Vec3r A(fe->vertexA()->point2d());
    Vec3r B(fe->vertexB()->point2d());
    int viewport[4];
    SilhouetteGeomEngine::retrieveViewport(viewport);
    if ((A.x() < viewport[0]) || (A.x() > viewport[2]) || (A.y() < viewport[1]) || (A.y() > viewport[3]) ||
            (B.x() < viewport[0]) || (B.x() > viewport[2]) || (B.y() < viewport[1]) || (B.y() > viewport[3])) {
        cerr << "Warning: point is out of the grid for fedge " << fe->getId() << endl;
        //return 0;
    }
#endif

    Vec3r vp;
    if (_orthographicProjection) {
        vp = Vec3r(center.x(), center.y(), _viewpoint.z());
    }
    else {
        vp = Vec3r(_viewpoint);
    }
    Vec3r u(vp - center);
    real raylength = u.norm();
    u.normalize();
#if 0
    if (_global.debug & G_DEBUG_FREESTYLE) {
        cout << "grid origin " << iGrid->getOrigin().x() << "," << iGrid->getOrigin().y() << ","
             << iGrid->getOrigin().z() << endl;
        cout << "center " << center.x() << "," << center.y() << "," << center.z() << endl;
    }
#endif

    iGrid->castRay(center, vp, occluders, timestamp);

    WFace *face = NULL;
    if (fe->isSmooth()) {
        FEdgeSmooth *fes = dynamic_cast<FEdgeSmooth *>(fe);
        face = (WFace *)fes->face();
    }
    vector<WVertex *> faceVertices;
    WVertex::incoming_edge_iterator ie;

    WFace *oface;
    bool skipFace;
    OccludersSet::iterator p, pend;
    if (face)
        face->RetrieveVertexList(faceVertices);

    for (p = occluders.begin(), pend = occluders.end(); p != pend; p++) {
        // If we're dealing with an exact silhouette, check whether we must take care of this occluder of not.
        // (Indeed, we don't consider the occluders that share at least one vertex with the face containing this edge).
        //-----------
        oface = (WFace *)(*p)->userdata;
#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\t\tEvaluating intersection for occluder " << ((*p)->getVertices())[0] <<
                    ((*p)->getVertices())[1] << ((*p)->getVertices())[2] << endl << "\t\t\tand ray " << vp <<
                    " * " << u << " (center " << center << ")" << endl;
        }
#endif
        Vec3r v1(((*p)->getVertices())[0]);
        Vec3r normal((*p)->getNormal());
        real d = -(v1 * normal);
        real t, t_u, t_v;

#if LOGGING
        if (_global.debug & G_DEBUG_FREESTYLE) {
            cout << "\t\tp:  " << ((*p)->getVertices())[0] << ((*p)->getVertices())[1] << ((*p)->getVertices())[2] <<
                    ", norm: " << (*p)->getNormal() << endl;
        }
#endif

        if (face) {
#if LOGGING
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "\t\tDetermining face adjacency...";
            }
#endif
            skipFace = false;

            if (face == oface) {
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "  Rejecting occluder for face concurrency." << endl;
                }
#endif
                continue;
            }

            for (vector<WVertex*>::iterator fv = faceVertices.begin(), fvend = faceVertices.end();
                 fv != fvend;
                 ++fv)
            {
                if ((*fv)->isBoundary())
                    continue;

                WVertex::incoming_edge_iterator iebegin = (*fv)->incoming_edges_begin();
                WVertex::incoming_edge_iterator ieend = (*fv)->incoming_edges_end();
                for (ie = iebegin; ie != ieend; ++ie) {
                    if ((*ie) == 0)
                        continue;

                    WFace *sface = (*ie)->GetbFace();
                    //WFace *sfacea = (*ie)->GetaFace();
                    //if ((sface == oface) || (sfacea == oface)) {
                    if (sface == oface) {
                        skipFace = true;
                        break;
                    }
                }
                if (skipFace)
                    break;
            }
            if (skipFace) {
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "  Rejecting occluder for face adjacency." << endl;
                }
#endif
                continue;
            }
        }
        else {
            // check whether the edge and the polygon plane are coincident:
            //-------------------------------------------------------------
            //first let us compute the plane equation.

            if (GeomUtils::COINCIDENT == GeomUtils::intersectRayPlane(origin, edge, normal, d, t, epsilon)) {
#if LOGGING
                if (_global.debug & G_DEBUG_FREESTYLE) {
                    cout << "\t\tRejecting occluder for target coincidence." << endl;
                }
#endif
                continue;
            }
        }

        if ((*p)->rayIntersect(center, u, t, t_u, t_v)) {
#if LOGGING
            if (_global.debug & G_DEBUG_FREESTYLE) {
                cout << "\t\tRay " << vp << " * " << u << " intersects at time " << t << " (raylength is " <<
                        raylength << ")" << endl;
                cout << "\t\t(u * normal) == " << (u * normal) << " for normal " << normal << endl;
            }
#endif
            if (fabs(u * normal) > 0.0001) {
                if ((t>0.0) && (t<raylength)) {
#if LOGGING
                    if (_global.debug & G_DEBUG_FREESTYLE) {
                        cout << "\t\tIs occluder" << endl;
                    }
#endif
                    WFace *f = (WFace *)((*p)->userdata);
                    ViewShape *vshape = _ViewMap->viewShape(f->GetVertex(0)->shape()->GetId());
                    oOccluders.insert(vshape);
                    ++qi;
                    if (!_EnableQI)
                        break;
                }
            }
        }
    }

    // Find occludee
    FindOccludee(fe, iGrid, epsilon, oaPolygon, timestamp, u, center, edge, origin, faceVertices);

    return qi;
}

void ViewMapBuilder::ComputeIntersections(ViewMap *ioViewMap, intersection_algo iAlgo, real epsilon)
{
    switch (iAlgo) {
    case sweep_line:
        ComputeSweepLineIntersections(ioViewMap, epsilon);
        break;
    default:
        break;
    }
#if 0
    if (_global.debug & G_DEBUG_FREESTYLE) {
        ViewMap::viewvertices_container& vvertices = ioViewMap->ViewVertices();
        for (ViewMap::viewvertices_container::iterator vv = vvertices.begin(), vvend = vvertices.end();
             vv != vvend; ++vv)
        {
            if ((*vv)->getNature() == Nature::T_VERTEX) {
                TVertex *tvertex = (TVertex *)(*vv);
                cout << "TVertex " << tvertex->getId() << " has :" << endl;
                cout << "FrontEdgeA: " << tvertex->frontEdgeA().first << endl;
                cout << "FrontEdgeB: " << tvertex->frontEdgeB().first << endl;
                cout << "BackEdgeA: " << tvertex->backEdgeA().first << endl;
                cout << "BackEdgeB: " << tvertex->backEdgeB().first << endl << endl;
            }
        }
    }
#endif
}

struct less_SVertex2D : public binary_function<SVertex *, SVertex *, bool>
{
    real epsilon;

    less_SVertex2D(real eps) : binary_function<SVertex *, SVertex *, bool>()
    {
        epsilon = eps;
    }

    bool operator()(SVertex *x, SVertex *y)
    {
        Vec3r A = x->point2D();
        Vec3r B = y->point2D();
        for (unsigned int i = 0; i < 3; i++) {
            if ((fabs(A[i] - B[i])) < epsilon)
                continue;
            if (A[i] < B[i])
                return true;
            if (A[i] > B[i])
                return false;
        }
        return false;
    }
};

typedef Segment<FEdge *, Vec3r> segment;
typedef Intersection<segment> intersection;

struct less_Intersection : public binary_function<intersection *, intersection *, bool>
{
    segment *edge;

    less_Intersection(segment *iEdge) : binary_function<intersection *, intersection *, bool>()
    {
        edge = iEdge;
    }

    bool operator()(intersection *x, intersection *y)
    {
        real tx = x->getParameter(edge);
        real ty = y->getParameter(edge);
        if (tx > ty)
            return true;
        return false;
    }
};

struct silhouette_binary_rule : public binary_rule<segment, segment>
{
    silhouette_binary_rule() : binary_rule<segment, segment>() {}

    virtual bool operator()(segment& s1, segment& s2)
    {
        FEdge *f1 = s1.edge();
        FEdge *f2 = s2.edge();

        if ((!(((f1)->getNature() & Nature::SILHOUETTE) || ((f1)->getNature() & Nature::BORDER))) &&
                (!(((f2)->getNature() & Nature::SILHOUETTE) || ((f2)->getNature() & Nature::BORDER))))
        {
            return false;
        }

        return true;
    }
};

void ViewMapBuilder::ComputeSweepLineIntersections(ViewMap *ioViewMap, real epsilon)
{
    vector<SVertex *>& svertices = ioViewMap->SVertices();
    bool progressBarDisplay = false;
    unsigned sVerticesSize = svertices.size();
    unsigned fEdgesSize = ioViewMap->FEdges().size();
#if 0
    if (_global.debug & G_DEBUG_FREESTYLE) {
        ViewMap::fedges_container& fedges = ioViewMap->FEdges();
        for (ViewMap::fedges_container::const_iterator f = fedges.begin(), end = fedges.end(); f != end; ++f) {
            cout << (*f)->aMaterialIndex() << "-" << (*f)->bMaterialIndex() << endl;
        }
    }
#endif
    unsigned progressBarStep = 0;

    if (_pProgressBar != NULL && fEdgesSize > gProgressBarMinSize) {
        unsigned progressBarSteps = min(gProgressBarMaxSteps, sVerticesSize);
        progressBarStep = sVerticesSize / progressBarSteps;
        _pProgressBar->reset();
        _pProgressBar->setLabelText("Computing Sweep Line Intersections");
        _pProgressBar->setTotalSteps(progressBarSteps);
        _pProgressBar->setProgress(0);
        progressBarDisplay = true;
    }

    unsigned counter = progressBarStep;

    sort(svertices.begin(), svertices.end(), less_SVertex2D(epsilon));

    SweepLine<FEdge *, Vec3r> SL;

    vector<FEdge *>& ioEdges = ioViewMap->FEdges();

    vector<segment*> segments;

    vector<FEdge*>::iterator fe, fend;

    for (fe = ioEdges.begin(), fend = ioEdges.end(); fe != fend; fe++) {
        segment *s = new segment((*fe), (*fe)->vertexA()->point2D(), (*fe)->vertexB()->point2D());
        (*fe)->userdata = s;
        segments.push_back(s);
    }

    vector<segment*> vsegments;
    for (vector<SVertex*>::iterator sv = svertices.begin(), svend = svertices.end(); sv != svend; sv++) {
        if (_pRenderMonitor && _pRenderMonitor->testBreak())
            break;

        const vector<FEdge*>& vedges = (*sv)->fedges();

        for (vector<FEdge *>::const_iterator sve = vedges.begin(), sveend = vedges.end(); sve != sveend; sve++) {
            vsegments.push_back((segment *)((*sve)->userdata));
        }

        Vec3r evt((*sv)->point2D());
        silhouette_binary_rule sbr;
        SL.process(evt, vsegments, sbr, epsilon);

        if (progressBarDisplay) {
            counter--;
            if (counter <= 0) {
                counter = progressBarStep;
                _pProgressBar->setProgress(_pProgressBar->getProgress() + 1);
            }
        }
        vsegments.clear();
    }

    if (_pRenderMonitor && _pRenderMonitor->testBreak()) {
        // delete segments
        if (!segments.empty()) {
            vector<segment*>::iterator s, send;
            for (s = segments.begin(), send = segments.end(); s != send; s++) {
                delete *s;
            }
        }
        return;
    }

    // list containing the new edges resulting from splitting operations.
    vector<FEdge*> newEdges;

    // retrieve the intersected edges:
    vector<segment*>& iedges = SL.intersectedEdges();
    // retrieve the intersections:
    vector<intersection*>& intersections = SL.intersections();

    int id = 0;
    // create a view vertex for each intersection and linked this one with the intersection object
    vector<intersection*>::iterator i, iend;
    for (i = intersections.begin(), iend = intersections.end(); i != iend; i++) {
        FEdge *fA = (*i)->EdgeA->edge();
        FEdge *fB = (*i)->EdgeB->edge();

        Vec3r A1 = fA->vertexA()->point3D();
        Vec3r A2 = fA->vertexB()->point3D();
        Vec3r B1 = fB->vertexA()->point3D();
        Vec3r B2 = fB->vertexB()->point3D();

        Vec3r a1 = fA->vertexA()->point2D();
        Vec3r a2 = fA->vertexB()->point2D();
        Vec3r b1 = fB->vertexA()->point2D();
        Vec3r b2 = fB->vertexB()->point2D();

        real ta = (*i)->tA;
        real tb = (*i)->tB;

        if ((ta < -epsilon) || (ta > 1 + epsilon))
            cerr << "Warning: 2D intersection out of range for edge " << fA->vertexA()->getId() << " - " <<
                    fA->vertexB()->getId() << endl;

        if ((tb < -epsilon) || (tb > 1 + epsilon))
            cerr << "Warning: 2D intersection out of range for edge " << fB->vertexA()->getId() << " - " <<
                    fB->vertexB()->getId() << endl;

        real Ta = SilhouetteGeomEngine::ImageToWorldParameter(fA, ta);
        real Tb = SilhouetteGeomEngine::ImageToWorldParameter(fB, tb);

        if ((Ta < -epsilon) || (Ta > 1 + epsilon))
            cerr << "Warning: 3D intersection out of range for edge " << fA->vertexA()->getId() << " - " <<
                    fA->vertexB()->getId() << endl;

        if ((Tb < -epsilon) || (Tb > 1 + epsilon))
            cerr << "Warning: 3D intersection out of range for edge " << fB->vertexA()->getId() << " - " <<
                    fB->vertexB()->getId() << endl;

#if 0
        if (_global.debug & G_DEBUG_FREESTYLE) {
            if ((Ta < -epsilon) || (Ta > 1 + epsilon) || (Tb < -epsilon) || (Tb > 1 + epsilon)) {
                printf("ta %.12e\n", ta);
                printf("tb %.12e\n", tb);
                printf("a1 %e, %e -- a2 %e, %e\n", a1[0], a1[1], a2[0], a2[1]);
                printf("b1 %e, %e -- b2 %e, %e\n", b1[0], b1[1], b2[0], b2[1]);
                //printf("line([%e, %e], [%e, %e]);\n", a1[0], a2[0], a1[1], a2[1]);
                //printf("line([%e, %e], [%e, %e]);\n", b1[0], b2[0], b1[1], b2[1]);
                if ((Ta < -epsilon) || (Ta > 1 + epsilon))
                    printf("Ta %.12e\n", Ta);
                if ((Tb < -epsilon) || (Tb > 1 + epsilon))
                    printf("Tb %.12e\n", Tb);
                printf("A1 %e, %e, %e -- A2 %e, %e, %e\n", A1[0], A1[1], A1[2], A2[0], A2[1], A2[2]);
                printf("B1 %e, %e, %e -- B2 %e, %e, %e\n", B1[0], B1[1], B1[2], B2[0], B2[1], B2[2]);
            }
        }
#endif

        TVertex *tvertex = ioViewMap->CreateTVertex(Vec3r(A1 + Ta * (A2 - A1)), Vec3r(a1 + ta * (a2 - a1)), fA,
                                                    Vec3r(B1 + Tb * (B2 - B1)), Vec3r(b1 + tb * (b2 - b1)), fB, id);

        (*i)->userdata = tvertex;
        ++id;
    }



    // -------------------------- do smooth-sharp intersections on the surface --------

    // because we haven't yet done any splits, there should be at most one sharp edge per mesh edge
    map<WEdge*, FEdgeSharp*>  edgemap;

    // create edgemap
    for(vector<FEdge*>::iterator feit = ioViewMap->FEdges().begin();
        feit != ioViewMap->FEdges().end(); ++feit)
        if (! (*feit)->isSmooth())
        {
            assert(dynamic_cast<FEdgeSharp*>(*feit) != NULL);
            FEdgeSharp * fes = (FEdgeSharp*)*feit;

            assert(edgemap.find(fes->edge()) == edgemap.end());
            edgemap.insert(pair<WEdge*,FEdgeSharp*>(fes->edge(), fes));
        }

    // find SVertex-fedge pairs, such that
    //  * the SVertex lies on a mesh edge (and thus is from a smooth curve)
    //  * the fedge is sharp
    //  * the two come from the same edge

    vector<SVertex*> newSVertices;

    for(vector<SVertex*>::iterator svit = ioViewMap->SVertices().begin();
        svit != ioViewMap->SVertices().end(); ++ svit)
    {
        SVertex * svSmooth = *svit;

        WEdge * meshEdge = svSmooth->sourceEdge();

        if (meshEdge == NULL) // only smooth curves have sources
            continue;

        map<WEdge*,FEdgeSharp*>::iterator it = edgemap.find(meshEdge);

        if (it == edgemap.end())
            continue;

        assert(svSmooth->viewvertex() == NULL || dynamic_cast<NonTVertex*>(svSmooth->viewvertex())!=NULL);

        FEdgeSharp * feSharp = (*it).second;

        Vec3r P3d = svSmooth->getPoint3D();
        Vec2f P2d2 = svSmooth->getPoint2D();
        Vec3r P2d(P2d2.x(), P2d2.y(), svSmooth->getProjectedZ());

        real t_sm = -1;
        real t_sh = GeomUtils::segmentParam(feSharp->vertexA()->getPoint2D(),
                                            feSharp->vertexB()->getPoint2D(), P2d2);

        segment * segsm = NULL;
        segment * segsh = (segment*)feSharp->userdata;
        Intersection<segment> * inter = new Intersection<segment>(segsm, t_sm, segsh, t_sh);

        intersections.push_back(inter);
        iedges.push_back(segsh);
        segsh->AddIntersection(inter);

        SVertex * svSharp = feSharp->shape()->CreateSVertex(P3d,P2d,feSharp->vertexA()->getId());
        svSharp->setSourceEdge(svSmooth->sourceEdge());
        // smooth vertex is put in "front" position

        ViewVertex * oldVV = svSmooth->viewvertex();


        TVertex * tvertex = new TVertex(svSmooth, svSharp);

        tvertex->setId(id);
        tvertex->setSameFace(true);

        ioViewMap->AddViewVertex(tvertex);
        newSVertices.push_back(svSharp);

        FEdge * feSmooth = svSmooth->fedges()[0]; // any edge on the smooth chain

        if (oldVV != NULL)
        {
            // if we're at a boundary, replace the viewvertex
            assert(dynamic_cast<NonTVertex*>(oldVV) != NULL);
            NonTVertex * ntv = (NonTVertex*)oldVV;
            assert(ntv->viewedges().size() == 1);
            assert(ntv->viewedges()[0].first == feSmooth->viewedge());

            ioViewMap->RemoveVertex(ntv);

            ViewEdge * viewedge = feSmooth->viewedge();
            assert (viewedge->A() == ntv || viewedge->B() == ntv);
            if (viewedge->A() == ntv)
            {
                viewedge->setA(tvertex);
                tvertex->setFrontEdgeB(viewedge, false);
            }
            else
            {
                viewedge->setB(tvertex);
                tvertex->setFrontEdgeA(viewedge, true);
            }
        }
        else
        {
            // split the smooth chain into two chains
            ViewEdge * newVEdge;
            feSmooth->viewedge()->viewShape()->SplitEdge(tvertex, feSmooth->viewedge(), newVEdge, true );
            if (newVEdge != feSmooth->viewedge())
                ioViewMap->AddViewEdge(newVEdge);
        }

        feSharp->viewedge()->viewShape()->AddVertex(tvertex);
        if (feSharp->viewedge()->viewShape() != feSmooth->viewedge()->viewShape())
            feSmooth->viewedge()->viewShape()->AddVertex(tvertex);

        svSmooth->setViewVertex(tvertex);
        svSharp->setViewVertex(tvertex);

        inter->userdata = tvertex;
        ++id;
    }

    // do this outside the above loop to avoid messing it up
    for(vector<SVertex*>::iterator sit = newSVertices.begin(); sit != newSVertices.end(); ++sit)
        ioViewMap->AddSVertex(*sit);


    edgemap.clear();


    progressBarStep = 0;

    if (progressBarDisplay) {
        unsigned iEdgesSize = iedges.size();
        unsigned progressBarSteps = min(gProgressBarMaxSteps, iEdgesSize);
        progressBarStep = iEdgesSize / progressBarSteps;
        _pProgressBar->reset();
        _pProgressBar->setLabelText("Splitting intersected edges");
        _pProgressBar->setTotalSteps(progressBarSteps);
        _pProgressBar->setProgress(0);
    }

    counter = progressBarStep;

    // -------- for each edge that has at least one intersection, find its intersections and split it ----
    // TODO: this should only be for edges that are not begin intersected by a sharp edge
    // e.g., when intersecting a mesh silhouette with a surface intersection, different code should split the silhouette

    vector<TVertex*> edgeVVertices;
    vector<ViewEdge*> newVEdges;
    vector<segment*>::iterator s, send;
    for (s = iedges.begin(), send = iedges.end(); s != send; s++) {
        edgeVVertices.clear();
        newEdges.clear();
        newVEdges.clear();

        FEdge *fedge = (*s)->edge();
        ViewEdge *vEdge = fedge->viewedge();
        ViewShape *shape = vEdge->viewShape();

        vector<intersection*>& eIntersections = (*s)->intersections();
        // we first need to sort these intersections from farther to closer to A
        sort(eIntersections.begin(), eIntersections.end(), less_Intersection(*s));
        for (i = eIntersections.begin(), iend = eIntersections.end(); i != iend; i++)
            edgeVVertices.push_back((TVertex *)(*i)->userdata);

        shape->SplitEdge(fedge, edgeVVertices, ioViewMap->FEdges(), ioViewMap->ViewEdges());

        if (progressBarDisplay) {
            counter--;
            if (counter <= 0) {
                counter = progressBarStep;
                _pProgressBar->setProgress(_pProgressBar->getProgress() + 1);
            }
        }
    }

    // reset userdata:
    for (fe = ioEdges.begin(), fend = ioEdges.end(); fe != fend; fe++)
        (*fe)->userdata = NULL;

    // delete segments
    if (!segments.empty()) {
        for (s = segments.begin(), send = segments.end(); s != send; s++) {
            delete *s;
        }
    }
}

} /* namespace Freestyle */
