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

#ifndef TIGLCOMMONFUNCTIONS_H
#define TIGLCOMMONFUNCTIONS_H

#include "tigl_internal.h"
#include "CCPACSImportExport.h"
#include "Standard.hxx"
#include "gp_Pnt.hxx"
#include "gp_Vec.hxx"
#include "TopoDS_Shape.hxx"
#include "PNamedShape.h"
#include "ListPNamedShape.h"
#include "CCPACSConfiguration.h"
#include "CTiglAbstractPhysicalComponent.h"
#include <TopoDS_Edge.hxx>
#include "TopoDS_Face.hxx"
#include <Geom_BSplineCurve.hxx>
#include "TopTools_ListOfShape.hxx"

#include <map>
#include <string>

typedef std::map<std::string, PNamedShape> ShapeMap;

// calculates a wire's circumfence
TIGL_EXPORT Standard_Real GetWireLength(const class TopoDS_Wire& wire);

TIGL_EXPORT Standard_Real GetEdgeLength(const class TopoDS_Edge& edge);

// returns a point on the wire (0 <= alpha <= 1)
TIGL_EXPORT gp_Pnt WireGetPoint(const TopoDS_Wire& wire, double alpha);
TIGL_EXPORT void WireGetPointTangent(const TopoDS_Wire& wire, double alpha, gp_Pnt& point, gp_Vec& normal);

TIGL_EXPORT gp_Pnt EdgeGetPoint(const TopoDS_Edge& edge, double alpha);
TIGL_EXPORT void EdgeGetPointTangent(const TopoDS_Edge& edge, double alpha, gp_Pnt& point, gp_Vec& normal);

// calculates the alpha value for a given point on a wire
TIGL_EXPORT Standard_Real ProjectPointOnWire(const TopoDS_Wire& wire, gp_Pnt p);

// projects a point onto the line (lineStart<->lineStop) and returns the projection parameter
TIGL_EXPORT Standard_Real ProjectPointOnLine(gp_Pnt p, gp_Pnt lineStart, gp_Pnt lineStop);

// returns the number of edges of the current shape
TIGL_EXPORT unsigned int GetNumberOfEdges(const TopoDS_Shape& shape);

// returns the number of faces of the current shape
TIGL_EXPORT unsigned int GetNumberOfFaces(const TopoDS_Shape& shape);

TIGL_EXPORT TopoDS_Edge GetEdge(const TopoDS_Shape& shape, int iEdge);

TIGL_EXPORT Handle_Geom_BSplineCurve GetBSplineCurve(const TopoDS_Edge& e);

// Returns the number of subshapes, if the shape is a compound
TIGL_EXPORT unsigned int GetNumberOfSubshapes(const TopoDS_Shape& shape);

// returns the central point of the face
TIGL_EXPORT gp_Pnt GetCentralFacePoint(const class TopoDS_Face& face);

// puts all faces with the same origin to one TopoDS_Compound
// Maps all compounds with its name in the map
TIGL_EXPORT ListPNamedShape GroupFaces(const PNamedShape shape, tigl::ShapeGroupMode groupType);

// Returns the coordinates of the bounding box of the shape
TIGL_EXPORT void GetShapeExtension(const TopoDS_Shape& shape,
                                   double& minx, double& maxx,
                                   double& miny, double& maxy,
                                   double& minz, double& maxz);

// Returns a unique Hashcode for a specific geometric component based on its loft
TIGL_EXPORT int GetComponentHashCode(tigl::ITiglGeometricComponent&);

// Method for building a face out of 4 points
TIGL_EXPORT TopoDS_Face BuildFace(const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3, const gp_Pnt& p4);

// Method for building a face out of two wires
TIGL_EXPORT TopoDS_Face BuildFace(const TopoDS_Wire& wire1, const TopoDS_Wire& wire2);

// Method for finding the intersection point of a face and an edge
TIGL_EXPORT bool GetIntersectionPoint(const TopoDS_Face& face, const TopoDS_Edge& edge, gp_Pnt& dst);

// Method for finding the intersection point of a face and a wire (containing edges)
TIGL_EXPORT bool GetIntersectionPoint(const TopoDS_Face& face, const TopoDS_Wire& wire, gp_Pnt& dst);

// Method for cutting two shapes, resulting in the common geometry (e.g. intersection edges)
TIGL_EXPORT TopoDS_Shape CutShapes(const TopoDS_Shape& shape1, const TopoDS_Shape& shape2);

// Method for getting all wires from the passed shape (list of edges)
TIGL_EXPORT void MakeWiresFromConnectedEdges(const TopoDS_Shape& shape, TopTools_ListOfShape& wireList);

// Method for getting a list of subshapes of a passed geometry
TIGL_EXPORT void GetListOfShape(const TopoDS_Shape& shape, TopAbs_ShapeEnum type, TopTools_ListOfShape& result);

// Method for finding all directly and indirectly connected edges
TIGL_EXPORT void FindAllConnectedEdges(const TopoDS_Edge& edge, TopTools_ListOfShape& edgeList, TopTools_ListOfShape& targetList);

// Method for checking whether two edges have a common vertex (regarding the position)
TIGL_EXPORT bool CheckCommonVertex(const TopoDS_Edge& e1, const TopoDS_Edge& e2);

// Method for sorting the edges of a wire
TIGL_EXPORT TopoDS_Wire SortWireEdges(const TopoDS_Wire& wire, bool closed);

// Method for creating a face from an opened wire
TIGL_EXPORT TopoDS_Wire CloseWire(const TopoDS_Wire& wire);

// Method for closing two wires to a single one,
// The method determines a direction vector based on the end vertices of wire1
// and calls the second closeWires method
TIGL_EXPORT TopoDS_Wire CloseWires(const TopoDS_Wire& wire1, const TopoDS_Wire& wire2);

// Method for closing two wires to a single one,
// the passed vector is used to define the upper and lower end vertices of the wires
TIGL_EXPORT TopoDS_Wire CloseWires(const TopoDS_Wire& wire1, const TopoDS_Wire& wire2, const gp_Vec& dir);

// Method for searching all vertices which are only connected to a single edge
TIGL_EXPORT void GetEndVertices(const TopoDS_Shape& shape, TopTools_ListOfShape& endVertices);

#endif // TIGLCOMMONFUNCTIONS_H
