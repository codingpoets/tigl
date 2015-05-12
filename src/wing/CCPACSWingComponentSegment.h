/* 
* Copyright (C) 2007-2013 German Aerospace Center (DLR/SC)
*
* Created: 2010-08-13 Markus Litz <Markus.Litz@dlr.de>
* Changed: $Id$ 
*
* Version: $Revision$
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
* @file
* @brief  Implementation of CPACS wing ComponentSegment handling routines.
*/

#ifndef CCPACSWINGCOMPONENTSEGMENT_H
#define CCPACSWINGCOMPONENTSEGMENT_H

#include <string>

#include "tigl_config.h"
#include "tigl_internal.h"
#include "tixi.h"
#include "CCPACSWingConnection.h"
#include "CCPACSWingCSStructure.h"
#include "CTiglPoint.h"
#include "CTiglAbstractSegment.h"

#include "TopoDS_Shape.hxx"
#include "TopoDS_Wire.hxx"
#include "Geom_BSplineSurface.hxx"
#include "CCPACSMaterial.h"

namespace tigl
{

class CCPACSWingSegment;

typedef std::vector<const CCPACSMaterial*>    MaterialList;
typedef std::vector<CCPACSWingSegment*>       SegmentList;

class CCPACSWing;

class CCPACSWingComponentSegment : public CTiglAbstractSegment
{

public:
    // Constructor
    TIGL_EXPORT CCPACSWingComponentSegment(CCPACSWing* aWing, int aSegmentIndex);

    // Virtual Destructor
    TIGL_EXPORT virtual ~CCPACSWingComponentSegment(void);

    // Invalidates internal state
    TIGL_EXPORT void Invalidate(void);

    // Read CPACS segment elements
    TIGL_EXPORT void ReadCPACS(TixiDocumentHandle tixiHandle, const std::string & segmentXPath);

    // Returns the wing this segment belongs to
    TIGL_EXPORT CCPACSWing& GetWing(void) const;

    // Gets a point in relative wing coordinates for a given eta and xsi
    TIGL_EXPORT gp_Pnt GetPoint(double eta, double xsi);

    // Get the eta xsi coordinate from a segment point (given by seta, sxsi)
    TIGL_EXPORT void GetEtaXsiFromSegmentEtaXsi(const std::string &segmentUID, double seta, double sxsi, double &eta, double &xsi);

    // Gets the volume of this segment
    TIGL_EXPORT double GetVolume();

    // Gets the surface area of this segment
    TIGL_EXPORT double GetSurfaceArea();

    // Gets the fromElementUID of this segment
    TIGL_EXPORT const std::string & GetFromElementUID(void) const;

    // Gets the toElementUID of this segment
    TIGL_EXPORT const std::string & GetToElementUID(void) const;

    // Returns the segment to a given point on the componentSegment and the nearest point projected onto the loft.
    // Returns null if the point is not an that wing, i.e. deviates more than 1 cm from the wing
    TIGL_EXPORT const CTiglAbstractSegment* findSegment(double x, double y, double z, gp_Pnt& nearestPoint);

    TIGL_EXPORT TiglGeometricComponentType GetComponentType(){ return TIGL_COMPONENT_WINGCOMPSEGMENT | TIGL_COMPONENT_SEGMENT | TIGL_COMPONENT_LOGICAL; }

    TIGL_EXPORT MaterialList GetMaterials(double eta, double xsi, TiglStructureType);

    TIGL_EXPORT CCPACSWingCSStructure& GetStructure();

    // returns a list of segments that belong to this component segment
    TIGL_EXPORT SegmentList& GetSegmentList();
        
    // creates an (iso) component segment line 
    TIGL_EXPORT TopoDS_Wire GetCSLine(double eta1, double xsi1, double eta2, double xsi2, int NSTEPS=101);

    // creates a straight line between two eta xsi coordinates
    TIGL_EXPORT TopoDS_Edge GetEtaXsiLine(double eta1, double xsi1, double eta2, double xsi2);

    // calculates the intersection of a segment iso eta line with a component segment line (defined by its start and end point)
    // returns the xsi coordinate of the intersection
    TIGL_EXPORT void GetSegmentIntersection(const std::string& segmentUID, double csEta1, double csXsi1, double csEta2, double csXsi2, double eta, double& xsi);

    // Helper method for getting the start segment index
    TIGL_EXPORT int GetStartSegmentIndex();

    // Return the shape for the componentsegment midplane
    TIGL_EXPORT TopoDS_Shape GetMidplaneShape();

    // [[CAS_AES]] added getter for midplane point
    TIGL_EXPORT gp_Pnt GetMidplanePoint(double eta, double xsi);

    // [[CAS_AES]] added getter for eta direction of midplane (no
    //             X-component)
    TIGL_EXPORT gp_Vec GetMidplaneEtaDir(double eta) const;

    // [[CAS_AES]] added getter for midplane normal vector
    TIGL_EXPORT gp_Vec GetMidplaneNormal(double eta);

    // [[CAS_AES]] added getter for the normalized leading edge direction
    TIGL_EXPORT gp_Vec GetLeadingEdgeDirection(const std::string& segmentUID) const;
    TIGL_EXPORT gp_Vec GetLeadingEdgeDirection(const gp_Pnt& point, const std::string& defaultSegmentUID = "") const;

    // [[CAS_AES]] added getter for the normalized trailing edge direction
    TIGL_EXPORT gp_Vec GetTrailingEdgeDirection(const std::string& segmentUID) const;
    TIGL_EXPORT gp_Vec GetTrailingEdgeDirection(const gp_Pnt& point, const std::string& defaultSegmentUID = "") const;

    // [[CAS_AES]] added getter for the normalized leading edge direction in the YZ plane
    TIGL_EXPORT gp_Vec GetLeadingEdgeDirectionYZ(const std::string& segmentUID) const;
    TIGL_EXPORT gp_Vec GetLeadingEdgeDirectionYZ(const gp_Pnt& point, const std::string& defaultSegmentUID = "") const;

    // [[CAS_AES]] added getter for leading edge point
    TIGL_EXPORT gp_Pnt GetLeadingEdgePoint(double eta) const;

    // [[CAS_AES]] added getter for trailing edge point
    TIGL_EXPORT gp_Pnt GetTrailingEdgePoint(double eta) const;

    // [[CAS_AES]] added getter for inner chordline point
    TIGL_EXPORT gp_Pnt GetInnerChordlinePoint(double xsi) const;

    // [[CAS_AES]] added getter for outer chordline point
    TIGL_EXPORT gp_Pnt GetOuterChordlinePoint(double xsi) const;

    // [[CAS_AES]] added getter for the midplane line between two eta-xsi points
    TIGL_EXPORT TopoDS_Wire GetMidplaneLine(double etaStart, double xsiStart, double etaEnd, double xsiEnd);

    // Returns the segment to a given point on the componentSegment.
    // Returns null if the point is not an that wing!
    // [[CAS_AES]] added const
    TIGL_EXPORT const std::string findSegment(double x, double y, double z) const;

    // [[CAS_AES]] added getter for inner segment UID
    TIGL_EXPORT std::string GetInnerSegmentUID() const;

    // [[CAS_AES]] added getter for outer segment UID
    TIGL_EXPORT std::string GetOuterSegmentUID() const;

    // [[CAS_AES]] added method for checking whether segment is contained in componentSegment
    TIGL_EXPORT bool IsSegmentContained(const CCPACSWingSegment& segment) const;

    // Signals if a structure is defined in the component segment
    TIGL_EXPORT bool HasStructure() const;

protected:
    // Cleanup routine
    void Cleanup(void);

    // Update internal segment data
    void Update(void);

    // Builds the loft between the two segment sections
    PNamedShape BuildLoft(void);

    // [[CAS_AES]] added method for building wires for eta-, leading edge-, trailing edge-lines
    void BuildLines(void);

    // Returns an upper or lower point on the segment surface in
    // dependence of parameters eta and xsi, which range from 0.0 to 1.0.
    // For eta = 0.0, xsi = 0.0 point is equal to leading edge on the
    // inner wing profile. For eta = 1.0, xsi = 1.0 point is equal to the trailing
    // edge on the outer wing profile. If fromUpper is true, a point
    // on the upper surface is returned, otherwise from the lower.
//  gp_Pnt GetPoint(double eta, double xsi, bool fromUpper);

private:
    // get short name for loft
    std::string GetShortShapeName(void);

    // Copy constructor
    CCPACSWingComponentSegment(const CCPACSWingComponentSegment& );

    // Assignment operator
    void operator=(const CCPACSWingComponentSegment& );

    std::vector<int> findPath(const std::string& fromUid, const::std::string& toUID, const std::vector<int>& curPath, bool forward) const;

    void UpdateProjectedLeadingEdge();

    // create short name
    std::string MakeShortName();

private:
    std::string          name;                 /**< Segment name                            */
    std::string          fromElementUID;       /**< Inner segment uid (root                 */
    std::string          toElementUID;         /**< Outer segment uid (tip)                 */
    CCPACSWing*          wing;                 /**< Parent wing                             */
    TopoDS_Shape         loft;                 /**< The loft between two sections           */
    double               myVolume;             /**< Volume of this segment                  */
    double               mySurfaceArea;        /**< Surface area of this segment            */
    TopoDS_Shape         upperShape;           /**< Upper shape of this componentSegment    */
    TopoDS_Shape         lowerShape;           /**< Lower shape of this componentSegment    */
    TopoDS_Wire          projLeadingEdge;      /**< (Extended) Leading edge projected into y-z plane */
    SegmentList          wingSegments;         /**< List of segments belonging to the component segment */
    Handle(Geom_Surface) upperSurface;
    Handle(Geom_Surface) lowerSurface;
    // [[CAS_AES]] added wires used for eta/xsi computations
    TopoDS_Wire          etaLine;                  // 2d version (in YZ plane) of leadingEdgeLine
    TopoDS_Wire          extendedEtaLine;          // 2d version (in YZ plane) of extendedLeadingEdgeLine
    TopoDS_Wire          leadingEdgeLine;          // leading edge as wire
    TopoDS_Wire          extendedLeadingEdgeLine;  // leading edge extended along y-axis, see CPACS spar geometry definition
    TopoDS_Wire          trailingEdgeLine;         // trailing edge as wire
    TopoDS_Wire          extendedTrailingEdgeLine; // trailing edge extended along y-axis, see CPACS spar geometry definition

    bool                 surfacesAreValid;
    CCPACSWingCSStructure structure;
    bool                 hasStructure;

};

} // end namespace tigl

#endif // CCPACSWINGCOMPONENTSEGMENT_H
