/* 
* Copyright (C) 2007-2015 German Aerospace Center (DLR/SC)
*
* Created: 2015-05-11 Jonas Jepsen <Jonas.Jepsen@dlr.de>
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

#include "CCPACSWingSparSegment.h"

#include "CTiglError.h"
#include "CTiglLogging.h"
#include "tiglcommonfunctions.h"
#include "CCPACSWingComponentSegment.h"

#include "BRep_Builder.hxx"
#include "BRepBndLib.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "Bnd_Box.hxx"
#include "TopExp_Explorer.hxx"

#include "BRepBuilderAPI_MakeEdge.hxx"

namespace tigl 
{

CCPACSWingSparSegment::CCPACSWingSparSegment(CCPACSWingSpars* spars)
    : sparsNode(spars),
      isvalid(false),
      m_auxGeomValid(false)
{
    Reset();
}

void CCPACSWingSparSegment::Reset()
{
    //children.Reset();
    Invalidate();
}

// Update internal segment data
void CCPACSWingSparSegment::Update(void)
{
    if (isvalid) {
        return;
    }
    BuildGeometry();
    isvalid = true;
}

void CCPACSWingSparSegment::ReadCPACS(TixiDocumentHandle tixiHandle, const std::string &sparSegmentXPath)
{
    Reset();
    
    // check path
    if ( tixiCheckElement(tixiHandle, sparSegmentXPath.c_str()) != SUCCESS) {
        LOG(ERROR) << "SparSegment " << sparSegmentXPath << " not found in CPACS file!" << std::endl;
        return;
    }

    // Get UID
    char * uidStr = NULL;
    if ( tixiGetTextAttribute(tixiHandle, sparSegmentXPath.c_str(), "uID", &uidStr) != SUCCESS ) {
        throw tigl::CTiglError("No UID given for wing spar segment " + sparSegmentXPath + "!", TIGL_UID_ERROR);
    }
    
    // Get name
    char * nameStr = NULL;
    std::string tmpPath = sparSegmentXPath + "/name";
    if ( tixiGetTextElement(tixiHandle, tmpPath.c_str(), &nameStr) != SUCCESS ) {
        throw tigl::CTiglError("No name given for wing spar segment " + tmpPath + "!", TIGL_UID_ERROR);
    }

    // Get description
    char * descrStr = NULL;
    tmpPath = sparSegmentXPath + "/description";
    if ( tixiGetTextElement(tixiHandle, tmpPath.c_str(), &descrStr) != SUCCESS ) {
        throw tigl::CTiglError("No description given for wing spar segment " + tmpPath + "!", TIGL_UID_ERROR);
    }

    // read spar cross section
    std::string crossSecPath = sparSegmentXPath + "/sparCrossSection";
    sparCrossSection.ReadCPACS(tixiHandle, crossSecPath.c_str());

    // read spar position uids
    std::string sparPosUIDsPath = sparSegmentXPath + "/sparPositionUIDs";
    sparPositionUIDs.ReadCPACS(tixiHandle, sparPosUIDsPath.c_str());

    uid = uidStr;
    name = nameStr;
    description = descrStr;
    isvalid = false;
}


void CCPACSWingSparSegment::Invalidate()
{
    isvalid = false;
    m_auxGeomValid = false;
}

bool CCPACSWingSparSegment::IsValid() const
{
    return isvalid;
}

// Gets the spar segment uid
const std::string& CCPACSWingSparSegment::GetUID(void) const
{
    return uid;
}

CCPACSWingSparCrossSection& CCPACSWingSparSegment::GetSparCrossSection()
{
    return sparCrossSection;
}

CCPACSWingSparPositionUIDs& CCPACSWingSparSegment::GetSparPositionUIDs()
{
    return sparPositionUIDs;
}

// Method for returning midplane line of spar
TopoDS_Wire CCPACSWingSparSegment::GetSparMidplaneLine()
{
    if (!m_auxGeomValid) {
        BuildAuxiliaryGeometry();
        m_auxGeomValid = true;
    }
    return m_sparMidplaneLine;
}

// Builds the cutting geometry for the spar as well as the midplane line
void CCPACSWingSparSegment::BuildAuxiliaryGeometry()
{
    // get assigned componentsegment
    CCPACSWingComponentSegment& componentSegment = sparsNode->GetWingComponentSegment();

    // build compound for cut geometry
    BRep_Builder cutCompBuilder;
    cutCompBuilder.MakeCompound(m_sparCutGeometry);

    // get bounding box of loft
    TopoDS_Shape loft = componentSegment.GetLoft()->Shape();
    Bnd_Box bbox;
    BRepBndLib::Add(loft, bbox);
    double bboxSize = sqrt(bbox.SquareExtent());

    // get rotation from sparCrossSection
    double rotation = sparCrossSection.GetRotation() * M_PI/180.0;
    if (rotation == 0) {
        throw CTiglError("Invalid Spar definition: angle of 0 degrees to midplane line!");
    }

    // container for all midplane points of the spar segment
    BRepBuilderAPI_MakeWire sparMidplaneLineBuilder;

    // up-vector of spar, initialized at first spar segment
    gp_Vec upVec;

    // iterate over all spar position points
    for (int i=1; i < sparPositionUIDs.GetSparPositionUIDCount(); i++) {
        std::string innerPositionUID = sparPositionUIDs.GetSparPositionUID(i);
        std::string outerPositionUID = sparPositionUIDs.GetSparPositionUID(i+1);

        // get spar position objects from spars node
        CCPACSWingSparPosition& innerSparPosition = sparsNode->GetSparPosition(innerPositionUID);
        CCPACSWingSparPosition& outerSparPosition = sparsNode->GetSparPosition(outerPositionUID);

        // get the inner and outer midplane points
        double innerEta = innerSparPosition.GetEta();
        double innerXsi = innerSparPosition.GetXsi();
        gp_Pnt innerPoint = componentSegment.GetMidplanePoint(innerEta, innerXsi);
        gp_Pnt outerPoint = componentSegment.GetMidplanePoint(outerSparPosition.GetEta(), outerSparPosition.GetXsi());

        // add midplane line to spar line
        TopoDS_Wire sparMidplaneLinePart = componentSegment.GetMidplaneLine(innerSparPosition.GetEta(), innerSparPosition.GetXsi(),
                                                                            outerSparPosition.GetEta(), outerSparPosition.GetXsi());
        sparMidplaneLineBuilder.Add(sparMidplaneLinePart);



        // handle special cases for eta == 0 or eta == 1:
        // - extend spar to bounding box in order to handle x-rotation of sections correctly
        if (innerSparPosition.GetEta() <= Precision::Confusion()) {
            gp_Vec sparDir(outerPoint, innerPoint);
            innerPoint.Translate(bboxSize * sparDir.Normalized());
        }
        if (outerSparPosition.GetEta() >= (1-Precision::Confusion())) {
            gp_Vec sparDir(innerPoint, outerPoint);
            outerPoint.Translate(bboxSize * sparDir.Normalized());
        }

        // BUG #149 and #152
        // because of issues with the spar up vectors in adjacent component
        // segments the up vector is temporarily set to the z direction
        upVec = gp_Vec(0,0,1);
        /*
        // determine up-vector based on midplane line of inner spar point
        if (i == 1) {
            double eta = innerSparPosition.GetEta();
            gp_Pnt pl = componentSegment.GetMidplanePoint(eta, 0);
            gp_Pnt pt = componentSegment.GetMidplanePoint(eta, 1);
            gp_Vec midplaneLine(pl, pt);
            // determine default segment, in case of inner/outer eta value
            // (required for extended eta line)
            std::string defaultSegmentUID;
            if (eta < 0.5) {
                defaultSegmentUID = componentSegment.GetInnerSegmentUID();
            } else {
                defaultSegmentUID = componentSegment.GetOuterSegmentUID();
            }
            gp_Vec leDir = componentSegment.GetLeadingEdgeDirectionYZ(pl, defaultSegmentUID);
            // determine up-vector by rotating the midplaneLine by the defined rotation angle,
            // invert the result because after the rotation the vector shows downwards
            upVec = -1 * midplaneLine.Rotated(gp_Ax1(gp_Pnt(0,0,0), leDir), rotation);
            upVec.Normalize();
        }
        */

        // Compute points for spar face used for cutting with loft
        gp_Pnt p1 = innerPoint.Translated(bboxSize * upVec);
        gp_Pnt p2 = innerPoint.Translated(-bboxSize * upVec);
        gp_Pnt p3 = outerPoint.Translated(bboxSize * upVec);
        gp_Pnt p4 = outerPoint.Translated(-bboxSize * upVec);

        // build face for cutting with loft
        TopoDS_Shape sparCutFace = BuildFace(p1, p2, p3, p4);

        // add face to split geometry compound
        cutCompBuilder.Add(m_sparCutGeometry, sparCutFace);
    }
    // store spar midplane line (required for rib position computation)
    m_sparMidplaneLine = sparMidplaneLineBuilder.Wire();
}

TopoDS_Shape CCPACSWingSparSegment::GetSparGeometry(bool relativeToWing)
{
    Update();

    if (relativeToWing) {
        return geometry;
    } else {
        return sparsNode->GetWingComponentSegment().GetWing().GetWingTransformation().Transform(geometry);
    }
}

TopoDS_Shape CCPACSWingSparSegment::GetSparCutGeometry()
{
    if (!m_auxGeomValid)
    {
        BuildAuxiliaryGeometry();
        m_auxGeomValid = true;
    }
    return m_sparCutGeometry;
}

// Builds the geometry
void CCPACSWingSparSegment::BuildGeometry(void) {
    CCPACSWingComponentSegment& componentSegment = sparsNode->GetWingComponentSegment();

    // build compound for spar geometry
    TopoDS_Compound compound;
    BRep_Builder builder;
    builder.MakeCompound(compound);

    // get bounding box of loft
    TopoDS_Shape loft = componentSegment.GetLoft()->Shape();
    Bnd_Box bbox;
    BRepBndLib::Add(loft, bbox);
    double bboxSize = sqrt(bbox.SquareExtent());

    // iterate over all spar cut faces
    TopoDS_Shape sparCutGeometry = GetSparCutGeometry();
    TopExp_Explorer exp;
    for (exp.Init(sparCutGeometry, TopAbs_FACE); exp.More(); exp.Next())
    {
        TopoDS_Face sparCutFace = TopoDS::Face(exp.Current());

        // intersect spar cut face with loft
        TopoDS_Shape sparCutEdges = CutShapes(loft, sparCutFace);

        // build wires out of connected edges
        TopTools_ListOfShape wireList;
        MakeWiresFromConnectedEdges(sparCutEdges, wireList);

/* ALTERNATIVE IMPLEMENTATION ALLOWING SPAR TO EXIT AND ENTER LOFT
        // make faces out of all closed wires
        TopTools_ListOfShape openedWires;
        TopTools_ListIteratorOfListOfShape wireListIt;
        for (wireListIt.Initialize(wireList); wireListIt.More(); wireListIt.Next())
        {
            const TopoDS_Wire& wire = TopoDS::Wire(wireListIt.Value());
            if (wire.Closed())
            {
                TopoDS_Face sparFace = BRepBuilderAPI_MakeFace(wire);
                builder.Add(compound, sparFace);
            }
            else
            {
                openedWires.Append(wire);
            }
        }

        // close opened wires and build spar faces
        if (openedWires.Extent() == 1)
        {
            TopoDS_Wire sparWire = CTiglCommon::closeWire(TopoDS::Wire(openedWires.First()));
            TopoDS_Shape sparFaces = BRepBuilderAPI_MakeFace(sparWire);
            builder.Add(compound, sparFaces);
        }
        else if (openedWires.Extent() == 2)
        {
            TopoDS_Shape sparFaces = CTiglCommon::buildFace(TopoDS::Wire(openedWires.First()), TopoDS::Wire(openedWires.Last()));
            builder.Add(compound, sparFaces);
        }
        else if (openedWires.Extent() > 2)
        {
            throw CTiglError("Error: unable to build spar geometry: more than two cut edges found on loft!");
        }
*/
        // build face(s) for spar
        TopoDS_Shape sparFaces;
        if (wireList.Extent() == 1) {
            TopoDS_Wire sparWire = CloseWire(TopoDS::Wire(wireList.First()));
            sparFaces = BRepBuilderAPI_MakeFace(sparWire);
        } else if (wireList.Extent() == 2) {
            sparFaces = BuildFace(TopoDS::Wire(wireList.First()), TopoDS::Wire(wireList.Last()));
        } else {
            throw CTiglError("Error: no geometry for spar definition found!");
        }

        // add spar face to compound
        builder.Add(compound, sparFaces);
    }

    // return spar geometry
    geometry = compound;
}

} // namespace tigl
