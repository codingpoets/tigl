/* 
* Copyright (C) 2007-2013 German Aerospace Center (DLR/SC)
*
* Created: 2013-05-28 Martin Siggel <Martin.Siggel@dlr.de>
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

#include "CCPACSWingSpars.h"
#include "CCPACSWingCSStructure.h"
#include "CCPACSWingComponentSegment.h"
#include "CCPACSWingSparPositionUIDs.h"
#include "CCPACSWingSparPosition.h"

#include "CTiglError.h"
#include "CTiglLogging.h"

#include "BRepPrimAPI_MakePrism.hxx"
#include "Geom_TrimmedCurve.hxx"
#include "GC_MakeSegment.hxx"
#include "TopoDS_Edge.hxx"
#include "TopoDS_Wire.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"

namespace tigl 
{

CCPACSWingSpars::CCPACSWingSpars(CCPACSWingCSStructure* parentStructure)
    : structure(parentStructure),
      sparPositions(),
      sparSegments(this)
{
    Reset();
}

void CCPACSWingSpars::Reset()
{
    //cells.Reset();
    Invalidate();
}

void CCPACSWingSpars::ReadCPACS(TixiDocumentHandle tixiHandle, const std::string &sparsXPath)
{
    Reset();
    
    // check path
    if ( tixiCheckElement(tixiHandle, sparsXPath.c_str()) != SUCCESS) {
        LOG(ERROR) << "Wing Spars " << sparsXPath << " not found in CPACS file!" << std::endl;
        return;
    }

    // read sparPositions
    std::string tempString;
    tempString = sparsXPath + "/sparPositions";
    if (tixiCheckElement(tixiHandle, tempString.c_str()) == SUCCESS){
        sparPositions.ReadCPACS(tixiHandle, tempString);
    }

    // read sparSegments
    tempString = sparsXPath + "/sparSegments";
    if (tixiCheckElement(tixiHandle, tempString.c_str()) == SUCCESS){
        sparSegments.ReadCPACS(tixiHandle, tempString);
    }
    
    isvalid = true;
}


void CCPACSWingSpars::Invalidate()
{
    isvalid = false;
}

bool CCPACSWingSpars::IsValid() const
{
    return isvalid;
}

CCPACSWingComponentSegment& CCPACSWingSpars::GetWingComponentSegment(void) const
{
    return structure->GetComponentSegment();
}

CCPACSWingSparPosition& CCPACSWingSpars::GetSparPosition(const std::string& uid)
{
    return sparPositions.GetSparPosition(uid);
}

CCPACSWingSparSegments& CCPACSWingSpars::GetSparSegments()
{
    return sparSegments;
}

TopoDS_Shape CCPACSWingSpars::GetLoft()
{
    /*// loop over all spar segments
    for (int i = 0; i < sparSegments.GetSparSegmentCount(); i++) {
        CCPACSWingSparSegment sparSegment = sparSegments.GetSparSegment(i);
        // get an ordered list of all the spar points
        std::vector<double> etas;
        std::vector<double> xsis;
        std::vector<gp_Pnt> sparPoints;
        CCPACSWingSparPositionUIDs positionUIDs = sparSegment.GetSparPositionUIDs();
        for (int j = 0; j < positionUIDs.GetSparPositionUIDCount(); j++) {
            // get the spar position from the uid
            std::string sparPosUID = positionUIDs.GetSparPositionUID(j);
            CCPACSWingSparPosition sparPosition = sparPositions.GetSparPosition(sparPosUID);
            // get eta and xsi coordinates of the spar point
            double eta = sparPosition.GetEta();
            double xsi = sparPosition.GetXsi();
            etas.push_back(eta);
            xsis.push_back(xsi);
            // get xyz coordinates for the spar point
            gp_Pnt p = structure->GetComponentSegment().GetPoint(eta, xsi);
            sparPoints.push_back(p);
        }
        // create spar segment wire
        BRepBuilderAPI_MakeWire wireBuilder;
        for (int j = 1; j < positionUIDs.GetSparPositionUIDCount(); j++) {
            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(sparPoints[j-1], sparPoints[j]);
            wireBuilder.Add(edge);
        }
        //return wireBuilder.Wire();
    }*/

    // TODO: get two spar points
    double eta1, eta2, eta3, eta4, xsi1, xsi2, xsi3, xsi4;
    std::string sparPosUID1 = "sparPos1";
    std::string sparPosUID2 = "sparPos2";
    std::string sparPosUID3 = "sparPos4";
    std::string sparPosUID4 = "sparPos3";
    CCPACSWingSparPosition sparPos1 = sparPositions.GetSparPosition(sparPosUID1.c_str());
    CCPACSWingSparPosition sparPos2 = sparPositions.GetSparPosition(sparPosUID2.c_str());
    CCPACSWingSparPosition sparPos3 = sparPositions.GetSparPosition(sparPosUID3.c_str());
    CCPACSWingSparPosition sparPos4 = sparPositions.GetSparPosition(sparPosUID4.c_str());
    eta1 = sparPos1.GetEta();
    eta2 = sparPos2.GetEta();
    eta3 = sparPos3.GetEta();
    eta4 = sparPos4.GetEta();
    xsi1 = sparPos1.GetXsi();
    xsi2 = sparPos2.GetXsi();
    xsi3 = sparPos3.GetXsi();
    xsi4 = sparPos4.GetXsi();
    gp_Pnt aPnt1 = structure->GetComponentSegment().GetPoint(eta1, xsi1);
    gp_Pnt aPnt2 = structure->GetComponentSegment().GetPoint(eta2, xsi2);
    gp_Pnt aPnt3 = structure->GetComponentSegment().GetPoint(eta3, xsi3);
    gp_Pnt aPnt4 = structure->GetComponentSegment().GetPoint(eta4, xsi4);
    std::cout << aPnt1.X() << aPnt2.X() << aPnt3.X() << aPnt4.X() << endl;
    std::cout << aPnt1.Y() << aPnt2.Y() << aPnt3.Y() << aPnt4.Y() << endl;
    std::cout << aPnt1.Z() << aPnt2.Z() << aPnt3.Z() << aPnt4.Z() << endl;
    // TODO: get the outer points corresponding to the spar points
    /*gp_Pnt aPnt1(10.0, 0.0, 0.0);
    gp_Pnt aPnt2(10.0, 0.0, 10.0);
    gp_Pnt aPnt3(0.0, 0.0, 10.0);
    gp_Pnt aPnt4(0.0, 0.0, 0.0);
    */

    // TODO: build a surface from those points
    Handle(Geom_TrimmedCurve) aSegment1 = GC_MakeSegment(aPnt1, aPnt2);
    Handle(Geom_TrimmedCurve) aSegment2 = GC_MakeSegment(aPnt2, aPnt3);
    Handle(Geom_TrimmedCurve) aSegment3 = GC_MakeSegment(aPnt3, aPnt4);
    Handle(Geom_TrimmedCurve) aSegment4 = GC_MakeSegment(aPnt4, aPnt1);

    TopoDS_Edge aEdge1 = BRepBuilderAPI_MakeEdge(aSegment1);
    TopoDS_Edge aEdge2 = BRepBuilderAPI_MakeEdge(aSegment2);
    TopoDS_Edge aEdge3 = BRepBuilderAPI_MakeEdge(aSegment3);
    TopoDS_Edge aEdge4 = BRepBuilderAPI_MakeEdge(aSegment4);

    TopoDS_Wire aWire = BRepBuilderAPI_MakeWire(aEdge1, aEdge2, aEdge3, aEdge4);

    TopoDS_Face myFaceProfile = BRepBuilderAPI_MakeFace(aWire);
    gp_Vec aPrismVec(1, 1, 1);
    TopoDS_Shape myBody = BRepPrimAPI_MakePrism(myFaceProfile, aPrismVec);
    return myBody;
}

} // namespace tigl
