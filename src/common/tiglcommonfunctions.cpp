/*
* Copyright (C) 2007-2013 German Aerospace Center (DLR/SC)
*
* Created: 2013-06-01 Martin Siggel <Martin.Siggel@dlr.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

/**
 * @brief This file implements some commonly used functions/algorithms on opencascade
 * data structures.
 */

#include <limits>
#include "tiglcommonfunctions.h"

#include "CTiglError.h"
#include "CNamedShape.h"
#include "boolean_operations/CBooleanOperTools.h"
#include "boolean_operations/BRepSewingToBRepBuilderShapeAdapter.h"
#include "ListPNamedShape.h"
#include "CNamedShape.h"
#include "PNamedShape.h"

#include "Geom_Curve.hxx"
#include "Geom_Surface.hxx"
#include "BRep_Tool.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Edge.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopTools_HSequenceOfShape.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "GeomAdaptor_Curve.hxx"
#include "BRepAdaptor_CompCurve.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GCPnts_AbscissaPoint.hxx"
#include "BRep_Builder.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "BRepTools.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"
#include "BRepFill.hxx"
#include "BRepExtrema_ExtCF.hxx"
#include "BRepAlgoAPI_Section.hxx"
#include "ShapeExtend_WireData.hxx"
#include "ShapeAnalysis_WireOrder.hxx"
#include "ShapeAnalysis_Edge.hxx"

#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dAPI_InterCurveCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomConvert.hxx>

#include "ShapeAnalysis_FreeBounds.hxx"

#include <list>
#include <algorithm>
#include <cassert>

namespace
{
    struct IsSame
    {
        IsSame(double tolerance) : _tol(tolerance) {}
        bool operator() (double first, double second)
        {
            return (fabs(first-second)<_tol);
        }

        double _tol;
    };
} // anonymous namespace

// calculates a wire's circumference
Standard_Real GetWireLength(const TopoDS_Wire& wire)
{
    BRepAdaptor_CompCurve aCompoundCurve(wire, Standard_True);
    return GCPnts_AbscissaPoint::Length( aCompoundCurve );
}

Standard_Real GetEdgeLength(const TopoDS_Edge &edge)
{
    Standard_Real umin, umax;
    Handle_Geom_Curve curve = BRep_Tool::Curve(edge, umin, umax);
    GeomAdaptor_Curve adaptorCurve(curve, umin, umax);
    Standard_Real length = GCPnts_AbscissaPoint::Length(adaptorCurve, umin, umax);
    return length;
}

unsigned int GetNumberOfEdges(const TopoDS_Shape& shape)
{
    TopExp_Explorer edgeExpl(shape, TopAbs_EDGE);
    unsigned int iEdges = 0;
    for (; edgeExpl.More(); edgeExpl.Next()) {
        iEdges++;
    }

    return iEdges;
}

unsigned int GetNumberOfFaces(const TopoDS_Shape& shape)
{
    TopExp_Explorer faceExpl(shape, TopAbs_FACE);
    unsigned int iFaces = 0;
    for (; faceExpl.More(); faceExpl.Next()) {
        iFaces++;
    }

    return iFaces;
}

unsigned int GetNumberOfSubshapes(const TopoDS_Shape &shape)
{
    if (shape.ShapeType() == TopAbs_COMPOUND) {
        unsigned int n = 0;
        for (TopoDS_Iterator anIter(shape); anIter.More(); anIter.Next()) {
            n++;
        }
        return n;
    }
    else {
        return 0;
    }
}

// Gets a point on the wire line in dependence of a parameter alpha with
// 0.0 <= alpha <= 1.0. For alpha = 0.0 this is the line starting point,
// for alpha = 1.0 the last point on the intersection line.
gp_Pnt WireGetPoint(const TopoDS_Wire& wire, double alpha)
{
    gp_Pnt point;
    gp_Vec tangent;
    WireGetPointTangent(wire, alpha, point, tangent);
    return point;
}

void WireGetPointTangent(const TopoDS_Wire& wire, double alpha, gp_Pnt& point, gp_Vec& tangent)
{
    if (alpha < 0.0 || alpha > 1.0) {
        throw tigl::CTiglError("Error: Parameter alpha not in the range 0.0 <= alpha <= 1.0 in WireGetPointTangent", TIGL_ERROR);
    }
    // ETA 3D point
    BRepAdaptor_CompCurve aCompoundCurve(wire, Standard_True);

    Standard_Real len =  GCPnts_AbscissaPoint::Length( aCompoundCurve );
    GCPnts_AbscissaPoint algo(aCompoundCurve, len*alpha, aCompoundCurve.FirstParameter());
    if (algo.IsDone()) {
        double par = algo.Parameter();
        aCompoundCurve.D1( par, point, tangent );
        // normalize tangent to length of the curve
        tangent = len*tangent/tangent.Magnitude();
    }
    else {
        throw tigl::CTiglError("WireGetPointTangent: Cannot compute point on curve.", TIGL_MATH_ERROR);
    }
}

gp_Pnt EdgeGetPoint(const TopoDS_Edge& edge, double alpha)
{
    gp_Pnt point;
    gp_Vec tangent;
    EdgeGetPointTangent(edge, alpha, point, tangent);
    return point;
}

void EdgeGetPointTangent(const TopoDS_Edge& edge, double alpha, gp_Pnt& point, gp_Vec& tangent)
{
    if (alpha < 0.0 || alpha > 1.0) {
        throw tigl::CTiglError("Error: Parameter alpha not in the range 0.0 <= alpha <= 1.0 in EdgeGetPointTangent", TIGL_ERROR);
    }
    // ETA 3D point
    Standard_Real umin, umax;
    Handle_Geom_Curve curve = BRep_Tool::Curve(edge, umin, umax);
    GeomAdaptor_Curve adaptorCurve(curve, umin, umax);
    Standard_Real len =  GCPnts_AbscissaPoint::Length( adaptorCurve, umin, umax );
    GCPnts_AbscissaPoint algo(adaptorCurve, len*alpha, umin);
    if (algo.IsDone()) {
        double par = algo.Parameter();
        adaptorCurve.D1( par, point, tangent );
        // normalize tangent to length of the curve
        tangent = len*tangent/tangent.Magnitude();
    }
    else {
        throw tigl::CTiglError("EdgeGetPointTangent: Cannot compute point on curve.", TIGL_MATH_ERROR);
    }
}

Standard_Real ProjectPointOnWire(const TopoDS_Wire& wire, gp_Pnt p)
{
    double smallestDist = DBL_MAX;
    double alpha  = 0.;
    int edgeIndex = 0;

    // find edge with closest dist to point p
    BRepTools_WireExplorer wireExplorer;
    int iwire = 0;
    for (wireExplorer.Init(wire); wireExplorer.More(); wireExplorer.Next(), iwire++) {
        Standard_Real firstParam, lastParam;
        TopoDS_Edge edge = wireExplorer.Current();
        Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, firstParam, lastParam);

        GeomAPI_ProjectPointOnCurve proj(p, curve, firstParam, lastParam);
        if (proj.NbPoints() > 0 && proj.LowerDistance() < smallestDist) {
            smallestDist = proj.LowerDistance();
            edgeIndex = iwire;
            alpha = proj.LowerDistanceParameter();
        }
    }

    // compute partial length of wire until projection point is reached
    wireExplorer.Init(wire);
    double partLength = 0.;
    for (int iwire = 0; iwire <= edgeIndex; ++iwire) {
        Standard_Real firstParam;
        Standard_Real lastParam;
        TopoDS_Edge edge = wireExplorer.Current();
        Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, firstParam, lastParam);
        GeomAdaptor_Curve adaptorCurve(curve, firstParam, lastParam);
        if (iwire == edgeIndex) {
            lastParam = alpha;
        }

        partLength += GCPnts_AbscissaPoint::Length(adaptorCurve, firstParam, lastParam);
        wireExplorer.Next();
    }

    // return relative coordinate
    double normalizedLength = partLength/GetWireLength(wire);
    if (normalizedLength > 1.0) {
        normalizedLength = 1.0;
    }
    else if (normalizedLength < 0.0) {
        normalizedLength = 0.0;
    }
    return normalizedLength;
}

gp_Pnt GetCentralFacePoint(const TopoDS_Face& face)
{
    // compute point on face
    Standard_Real umin, umax, vmin, vmax;

    gp_Pnt p;

    Handle_Geom_Surface surface = BRep_Tool::Surface(face);
    BRepTools::UVBounds(face, umin, umax, vmin, vmax);
    Standard_Real umean = 0.5*(umin+umax);
    Standard_Real vmean = 0.5*(vmin+vmax);


    // compute intersection of u-iso line with face boundaries
    Handle_Geom2d_Curve uiso = new Geom2d_Line(
                gp_Pnt2d(umean,0.),
                gp_Dir2d(0., 1.)
                );

    TopExp_Explorer exp (face,TopAbs_EDGE);
    std::list<double> intersections;
    for (; exp.More(); exp.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(exp.Current());
        Standard_Real first, last;

        // Get geomteric curve from edge
        Handle_Geom2d_Curve hcurve = BRep_Tool::CurveOnSurface(edge, face, first, last);
        hcurve = new Geom2d_TrimmedCurve(hcurve, first, last);

        Geom2dAPI_InterCurveCurve intersector(uiso, hcurve);
        for (int ipoint = 0; ipoint < intersector.NbPoints(); ++ipoint) {
            gp_Pnt2d p = intersector.Point(ipoint+1);
            intersections.push_back(p.Y());
        }
    }

    // remove duplicate solutions defined by tolerance
    double tolerance = 1e-5;
    intersections.sort();
    intersections.unique(IsSame((vmax-vmin)*tolerance));

    // normally we should have at least two intersections
    // also the number of sections should be even - else something is really strange
    //assert(intersections.size() % 2 == 0);
    if (intersections.size() >= 2) {
        std::list<double>::iterator it = intersections.begin();
        double int1 = *it++;
        double int2 = *it;
        vmean = (int1 + int2)/2.;
    }

    surface->D0(umean, vmean, p);

    return p;
}

ListPNamedShape GroupFaces(const PNamedShape shape, tigl::ShapeGroupMode groupType)
{
    ListPNamedShape shapeList;
    if (!shape) {
        return shapeList;
    }

    if (groupType == tigl::NAMED_COMPOUNDS) {
        BRep_Builder b;
        TopTools_IndexedMapOfShape faceMap;
        std::map<PNamedShape, TopoDS_Shape> map;
        TopExp::MapShapes(shape->Shape(),   TopAbs_FACE, faceMap);
        if (faceMap.Extent() == 0) {
            // return the shape as is
            shapeList.push_back(shape);
            return shapeList;
        }
        
        for (int iface = 1; iface <= faceMap.Extent(); ++iface) {
            TopoDS_Face face = TopoDS::Face(faceMap(iface));
            PNamedShape origin = shape->GetFaceTraits(iface-1).Origin();

            std::map<PNamedShape, TopoDS_Shape>::iterator it = map.find(origin);
            if (it == map.end()) {
                TopoDS_Compound c;
                b.MakeCompound(c);
                b.Add(c, face);
                map[origin] = c;
            }
            else {
                TopoDS_Shape& c = it->second;
                b.Add(c, face);
            }
        }

        // create Named Shapes
        std::map<PNamedShape, TopoDS_Shape>::iterator it;
        for (it = map.begin(); it != map.end(); ++it) {
            PNamedShape  origin     = it->first;
            TopoDS_Shape toposhape  = it->second;
            PNamedShape curshape;
            if (origin) {
                curshape = PNamedShape(new CNamedShape(toposhape, origin->Name()));
                curshape->SetShortName(origin->ShortName());
            }
            else {
                curshape = PNamedShape(new CNamedShape(toposhape, shape->Name()));
                curshape->SetShortName(shape->ShortName());
            }
            // set the original face traits
            CBooleanOperTools::AppendNamesToShape(shape, curshape);

            // make shells
            curshape = CBooleanOperTools::Shellify(curshape);
            shapeList.push_back(curshape);
        }
    }
    else if (groupType == tigl::WHOLE_SHAPE) {
        shapeList.push_back(shape);
    }
    else if (groupType == tigl::FACES) {
        // store each face as an own shape
        TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(shape->Shape(), TopAbs_FACE, faceMap);
        if (faceMap.Extent() == 0) {
            // return the shape as is
            shapeList.push_back(shape);
            return shapeList;
        }
        
        for (int iface = 1; iface <= faceMap.Extent(); ++iface) {
            TopoDS_Face face = TopoDS::Face(faceMap(iface));
            const CFaceTraits& traits = shape->GetFaceTraits(iface-1);

            PNamedShape faceShape;
            if (traits.Origin()) {
                faceShape = PNamedShape(new CNamedShape(face, traits.Origin()->Name()));
                faceShape->SetShortName(traits.Origin()->ShortName());
            }
            else {
                faceShape = PNamedShape(new CNamedShape(face, shape->Name()));
                faceShape->SetShortName(shape->ShortName());
            }
            faceShape->SetFaceTraits(0, shape->GetFaceTraits(iface-1));
            shapeList.push_back(faceShape);
        }
        
    }
    return shapeList;
}

// projects a point onto the line (lineStart<->lineStop) and returns the projection parameter
Standard_Real ProjectPointOnLine(gp_Pnt p, gp_Pnt lineStart, gp_Pnt lineStop)
{
    return gp_Vec(lineStart, p) * gp_Vec(lineStart, lineStop) / gp_Vec(lineStart, lineStop).SquareMagnitude();
}

// Returns the coordinates of the bounding box of the shape
void GetShapeExtension(const TopoDS_Shape& shape,
                       double& minx, double& maxx,
                       double& miny, double& maxy,
                       double& minz, double& maxz)
{
    Bnd_Box boundingBox;
    BRepBndLib::Add(shape, boundingBox);
    boundingBox.Get(minx, miny, minz, maxx, maxy, maxz);
}

// Returns a unique Hashcode for a specific geometric component
int GetComponentHashCode(tigl::ITiglGeometricComponent& component)
{
    const TopoDS_Shape& loft = (*component.GetLoft()).Shape();
    if (!loft.IsNull()) {
        return loft.HashCode(2294967295);
    }
    else {
        return 0;
    }
}


TopoDS_Edge GetEdge(const TopoDS_Shape &shape, int iEdge)
{
    TopTools_IndexedMapOfShape edgeMap;
    TopExp::MapShapes(shape, TopAbs_EDGE, edgeMap);
    
    if (iEdge < 0 || iEdge >= edgeMap.Extent()) {
        return TopoDS_Edge();
    }
    else {
        return TopoDS::Edge(edgeMap(iEdge+1));
    }
}

Handle_Geom_BSplineCurve GetBSplineCurve(const TopoDS_Edge& e)
{
    double u1, u2;
    Handle_Geom_Curve curve = BRep_Tool::Curve(e, u1, u2);
    curve = new Geom_TrimmedCurve(curve, u1, u2);
    
    // convert to bspline
    Handle_Geom_BSplineCurve bspl =  GeomConvert::CurveToBSplineCurve(curve);
    return bspl;
}

TopoDS_Face BuildFace(const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3, const gp_Pnt& p4)
{
    BRepBuilderAPI_MakeEdge me1(p1, p2);
    TopoDS_Edge e1 = me1.Edge();
    BRepBuilderAPI_MakeEdge me2(p3, p4);
    TopoDS_Edge e2 = me2.Edge();

    TopoDS_Face face = BRepFill::Face(e1, e2);
    return face;
}

TopoDS_Face BuildFace(const TopoDS_Wire& wire1, const TopoDS_Wire& wire2)
{
    TopoDS_Wire closedWire = CloseWires(wire1, wire2);
    TopoDS_Face face = BRepBuilderAPI_MakeFace(closedWire);
    return face;
}

// Method for finding the intersection point of a face and an edge
bool GetIntersectionPoint(const TopoDS_Face& face, const TopoDS_Edge& edge, gp_Pnt& dst)
{
    BRepExtrema_ExtCF edgeFaceIntersect(edge, face);
    if (!edgeFaceIntersect.IsDone()) {
        throw tigl::CTiglError("Error intersecting edge and face in CTiglCommon::getIntersectionPoint!");
    }
    double minDistance = std::numeric_limits<double>::max();
    int minIndex = 0;
    for (int i=1; i <= edgeFaceIntersect.NbExt(); i++) {
        double distance = edgeFaceIntersect.SquareDistance(i);
        if (distance < minDistance) {
            minDistance = distance;
            minIndex = i;
        }
    }

    if (minDistance <= Precision::Confusion()) {
        dst = edgeFaceIntersect.PointOnEdge(minIndex);
        return true;
    }
    return false;
}

// Method for finding the intersection point of a face and a wire (containing edges)
bool GetIntersectionPoint(const TopoDS_Face& face, const TopoDS_Wire& wire, gp_Pnt& dst)
{
    BRepTools_WireExplorer wireExp;
    for (wireExp.Init(wire); wireExp.More(); wireExp.Next()) {
        const TopoDS_Edge& edge = wireExp.Current();
        if (GetIntersectionPoint(face, edge, dst)) {
            return true;
        }
    }
    return false;
}

TopoDS_Shape CutShapes(const TopoDS_Shape& shape1, const TopoDS_Shape& shape2)
{
    BRepAlgoAPI_Section splitter(shape1, shape2, Standard_False);
    splitter.Approximation(Standard_True);
    splitter.Build();
    if (!splitter.IsDone()) {
        throw tigl::CTiglError("Error cutting shapes in CTiglCommon::cutShapes!");
    }
    TopoDS_Shape result = splitter.Shape();
    return result;
}

// Method for getting all wires from the passed shape (list of edges)
void MakeWiresFromConnectedEdges(const TopoDS_Shape& shape, TopTools_ListOfShape& wireList)
{
    TopExp_Explorer myEdgeExplorer (shape, TopAbs_EDGE);
    Handle(TopTools_HSequenceOfShape) Edges = new TopTools_HSequenceOfShape();

    while (myEdgeExplorer.More()) {
        Edges->Append(TopoDS::Edge(myEdgeExplorer.Current()));
        myEdgeExplorer.Next();
    }

    // connect all connected edges to wires and save them in container Edges again
    double tolerance = 1.0e-7;
    ShapeAnalysis_FreeBounds::ConnectEdgesToWires(Edges, tolerance, false, Edges);
    int numWires = Edges->Length();

    std::vector<TopoDS_Wire> Wires;

    // filter duplicated wires
    for (int wireID=1; wireID <= numWires; wireID++) {
        bool found = false;
        TopoDS_Wire wire = TopoDS::Wire(Edges->Value(wireID));
        for (std::vector<TopoDS_Wire>::size_type i = 0; i < Wires.size(); i++) {
            if (Wires[i].HashCode(200000) == wire.HashCode(200000)) {
                    found = true;
            }
        }

        if (!found) {
            Wires.push_back(wire);
            wireList.Append(wire);
        }
    }
}

void GetListOfShape(const TopoDS_Shape& shape, TopAbs_ShapeEnum type, TopTools_ListOfShape& result)
{
    TopTools_IndexedMapOfShape typeMap;
    TopExp::MapShapes(shape, type, typeMap);
    for (int i = 1; i <= typeMap.Extent(); i++) {
        result.Append(typeMap.FindKey(i));
    }
}

// Method for finding all directly and indirectly connected edges
// The method loops over the passed edgeList and checks for each element if it
// is connected to the passed edge. When an edge is found it is removed from
// the edgeList and added to the targetList. Additionally for this edge all
// connected edges are also added to the targetList by recursively calling this
// method. Finally all directly or indirectly connected edges to the passed
// edge are moved from the edgeList to the targetList
void FindAllConnectedEdges(const TopoDS_Edge& edge, TopTools_ListOfShape& edgeList, TopTools_ListOfShape& targetList)
{
    TopTools_ListIteratorOfListOfShape edgeIt;
    TopoDS_Vertex tempVertex;
    TopoDS_Edge connectedEdge;

    // finished gets true as soon as no more connected edge is found in the edgeList
    bool finished = false;
    // loop until all connected e
    while (!finished) {
        connectedEdge.Nullify();
        // iterate over all edges in edgeList, and break when a connected edge was found
        for (edgeIt.Initialize(edgeList); edgeIt.More(); edgeIt.Next()) {
            TopoDS_Edge testEdge = TopoDS::Edge(edgeIt.Value());

            // check if edges are connected
            if (CheckCommonVertex(edge, testEdge)) {
                connectedEdge = testEdge;
                break;
            }
        }
        // check if a connected edge was found
        if (!connectedEdge.IsNull()) {
            // remove connected edge from edgeList
            edgeList.Remove(edgeIt);
            // append connected edge to connectedEdge list
            targetList.Append(connectedEdge);
            // append all edges which are connected to the new edge to the list
            FindAllConnectedEdges(connectedEdge, edgeList, targetList);
        } else {
            finished = true;
        }
    }
}

// Method for checking if two edges have a common vertex (same position)
bool CheckCommonVertex(const TopoDS_Edge& e1, const TopoDS_Edge& e2)
{
    TopoDS_Vertex v1First, v1Last, v2First, v2Last;
    TopExp::Vertices(e1, v1First, v1Last);
    TopExp::Vertices(e2, v2First, v2Last);

    gp_Pnt p1First = BRep_Tool::Pnt(v1First);
    gp_Pnt p1Last = BRep_Tool::Pnt(v1Last);
    gp_Pnt p2First = BRep_Tool::Pnt(v2First);
    gp_Pnt p2Last = BRep_Tool::Pnt(v2Last);

    if (p1First.Distance(p2First) < Precision::Confusion() ||
        p1First.Distance(p2Last) < Precision::Confusion() ||
        p1Last.Distance(p2First) < Precision::Confusion() ||
        p1Last.Distance(p2Last) < Precision::Confusion()) {
        return true;
    } else {
        return false;
    }
}

TopoDS_Wire SortWireEdges(const TopoDS_Wire& wire, bool closed)
{
    // Inspired by ShapeAnalysis_Wire::CheckOrder
    ShapeExtend_WireData checkWire(wire);

    // Determine correct order of edges in wire
    ShapeAnalysis_WireOrder wireOrder(Standard_True /* 3D */, Precision::Confusion());
    // Use ShapeAnalysis_Edge here, because this returns the correct vertex because it interprets the orientation!!!
    ShapeAnalysis_Edge checkEdge;

    // iterate over the number of edges
    int numEdges = checkWire.NbEdges();
    for (int i  = 1; i <= numEdges; i++) {
        TopoDS_Edge edge = checkWire.Edge(i);
        TopoDS_Vertex v1 = checkEdge.FirstVertex(edge);
        TopoDS_Vertex v2 = checkEdge.LastVertex(edge);
        gp_Pnt p1 = BRep_Tool::Pnt(v1);
        gp_Pnt p2 = BRep_Tool::Pnt(v2);
        // add edges to wireOrder class
        wireOrder.Add(p1.XYZ(), p2.XYZ());
    }
    wireOrder.Perform(closed);

    if (!wireOrder.IsDone()) {
        throw tigl::CTiglError("Error: wireOrder could not be determined in CTiglCommon::SortWireEdges!");
    }

    // get status from wireOrder
    int status = wireOrder.Status();
    TopoDS_Wire fixedWire;

    if (status == 0) {
        // wire is ok
        fixedWire = wire;
    } else if (status == 1 || status == -1) {
        // wire is out of order and edges are probably inversed
        BRepBuilderAPI_MakeWire makeWire;
        for (int i = 1; i <= numEdges; i++) {
            int edgeIndex = wireOrder.Ordered(i);
            bool reverse = false;
            if (edgeIndex < 0) {
                reverse = true;
                edgeIndex = -edgeIndex;
            }
            TopoDS_Edge edge = TopoDS::Edge(checkWire.Edge(edgeIndex));
            if (reverse) {
                edge.Reverse();
            }
            makeWire.Add(edge);
        }
        fixedWire = makeWire.Wire();
    } else {
        throw tigl::CTiglError("Error: failure determining wire order in CTiglCommon::SortWireEdges!");
    }
    return fixedWire;
}

TopoDS_Wire CloseWire(const TopoDS_Wire& wire)
{
    // get the list of end vertices
    TopTools_ListOfShape endVertices;
    GetEndVertices(wire, endVertices);

    // determine number of end vertices
    int numEndVertices = endVertices.Extent();

    // check if wire is already closed
    if (numEndVertices == 0) {
        return wire;
    }

    // check if we have exatcly two end vertices
    if (numEndVertices != 2) {
        throw tigl::CTiglError("Error: invalid number of end vertices found in CTiglCommon::closeWire!");
    }

    // next generate an edge between the end vertices
    TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(TopoDS::Vertex(endVertices.First()), TopoDS::Vertex(endVertices.Last()));

    // build new wire
    BRepBuilderAPI_MakeWire makeWire;
    makeWire.Add(wire);
    makeWire.Add(edge);
    if (!makeWire.IsDone()) {
        throw tigl::CTiglError("Error: unable to build closed wire in CTiglCommon::closeWire!");
    }
    TopoDS_Wire result = SortWireEdges(makeWire.Wire(), true);
    return result;
}

// determine direction from wire1 by using end vertices
TopoDS_Wire CloseWires(const TopoDS_Wire& wire1, const TopoDS_Wire& wire2)
{
    // get a map of vertices with connected edges
    TopTools_IndexedDataMapOfShapeListOfShape vertexMap;
    TopTools_ListOfShape endVertices1, endVertices2;

    TopExp::MapShapesAndAncestors(wire1, TopAbs_VERTEX, TopAbs_EDGE, vertexMap);
    // filter out all vertices which have only a single edge connected (end vertices)
    for (int i = 1; i <= vertexMap.Extent(); i++) {
        const TopTools_ListOfShape& edgeList = vertexMap.FindFromIndex(i);
        if (edgeList.Extent() == 1) {
            endVertices1.Append(vertexMap.FindKey(i));
        }
    }

    // check for correct number of end vertices
    if (endVertices1.Extent() != 2) {
        throw tigl::CTiglError("Error: Unable to close wires because invalid number of end-vertices found!");
    }

    TopoDS_Vertex& v1 = TopoDS::Vertex(endVertices1.First());
    TopoDS_Vertex& v2 = TopoDS::Vertex(endVertices1.Last());
    gp_Pnt p1 = BRep_Tool::Pnt(v1);
    gp_Pnt p2 = BRep_Tool::Pnt(v2);
    gp_Vec dir(p1, p2);

    return CloseWires(wire1, wire2, dir);
}

TopoDS_Wire CloseWires(const TopoDS_Wire& wire1, const TopoDS_Wire& wire2, const gp_Vec& dir)
{
    // get a map of vertices with connected edges
    TopTools_IndexedDataMapOfShapeListOfShape vertexMap;
    TopTools_ListOfShape endVertices1, endVertices2;

    TopExp::MapShapesAndAncestors(wire1, TopAbs_VERTEX, TopAbs_EDGE, vertexMap);
    // filter out all vertices which have only a single edge connected (end vertices)
    for (int i = 1; i <= vertexMap.Extent(); i++) {
        const TopTools_ListOfShape& edgeList = vertexMap.FindFromIndex(i);
        if (edgeList.Extent() == 1) {
            endVertices1.Append(vertexMap.FindKey(i));
        }
    }
    vertexMap.Clear();
    TopExp::MapShapesAndAncestors(wire2, TopAbs_VERTEX, TopAbs_EDGE, vertexMap);
    // filter out all vertices which have only a single edge connected (end vertices)
    for (int i = 1; i <= vertexMap.Extent(); i++) {
        const TopTools_ListOfShape& edgeList = vertexMap.FindFromIndex(i);
        if (edgeList.Extent() == 1) {
            endVertices2.Append(vertexMap.FindKey(i));
        }
    }

    // check if number of end vertices is correct
    int numEndVertices1 = endVertices1.Extent();
    int numEndVertices2 = endVertices2.Extent();
    if (numEndVertices1 != 2 || numEndVertices2 != 2) {
        throw tigl::CTiglError("Error: Unable to close wires because invalid number of end-vertices found!");
    }

    // sort end vertices according to direction vector
    TopoDS_Vertex vUpper1 = TopoDS::Vertex(endVertices1.First());
    TopoDS_Vertex vLower1 = TopoDS::Vertex(endVertices1.Last());
    gp_Pnt pUpper = BRep_Tool::Pnt(vUpper1);
    gp_Pnt pLower = BRep_Tool::Pnt(vLower1);
    gp_Vec vDiff(pLower, pUpper);
    if (vDiff.Normalized().Dot(dir.Normalized()) < 0) {
        vUpper1 = TopoDS::Vertex(endVertices1.Last());
        vLower1 = TopoDS::Vertex(endVertices1.First());
    }
    TopoDS_Vertex vUpper2 = TopoDS::Vertex(endVertices2.First());
    TopoDS_Vertex vLower2 = TopoDS::Vertex(endVertices2.Last());
    pUpper = BRep_Tool::Pnt(vUpper2);
    pLower = BRep_Tool::Pnt(vLower2);
    vDiff = gp_Vec(pLower, pUpper);
    if (vDiff.Normalized().Dot(dir.Normalized()) < 0) {
        vUpper2 = TopoDS::Vertex(endVertices2.Last());
        vLower2 = TopoDS::Vertex(endVertices2.First());
    }

    // next build the edges and combine the wires to a single wire
    TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(vUpper1, vUpper2);
    TopoDS_Edge edge2 = BRepBuilderAPI_MakeEdge(vLower1, vLower2);

    // put all edges into a large list
    TopTools_ListOfShape edgeList;
    edgeList.Append(edge1);
    edgeList.Append(edge2);
    TopoDS_Iterator edgeIt;
    for (edgeIt.Initialize(wire1); edgeIt.More(); edgeIt.Next()) {
        edgeList.Append(edgeIt.Value());
    }
    for (edgeIt.Initialize(wire2); edgeIt.More(); edgeIt.Next()) {
        edgeList.Append(edgeIt.Value());
    }

    BRepBuilderAPI_MakeWire makeWire;
    makeWire.Add(edgeList);
    if (!makeWire.IsDone()) {
        throw tigl::CTiglError("Error: error during creation of closed wire in CTiglCommon::closeWires!");
    }

    TopoDS_Wire result = SortWireEdges(makeWire.Wire(), true /* closed */);

    return result;
}

void GetEndVertices(const TopoDS_Shape& shape, TopTools_ListOfShape& endVertices)
{
    // get a map of vertices with connected edges
    TopTools_IndexedDataMapOfShapeListOfShape vertexMap;
    TopExp::MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, vertexMap);

    // filter out all vertices which have only a single edge connected (end vertices)
    for (int i = 1; i <= vertexMap.Extent(); i++) {
        const TopTools_ListOfShape& edgeList = vertexMap.FindFromIndex(i);
        if (edgeList.Extent() == 1) {
            endVertices.Append(vertexMap.FindKey(i));
        }
    }
}
