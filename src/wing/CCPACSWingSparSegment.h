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

#ifndef CCPACSWINGSPARSEGMENT_H
#define CCPACSWINGSPARSEGMENT_H

#include <string>
#include "tigl_internal.h"
#include "tixi.h"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Wire.hxx"
#include "TopoDS_Compound.hxx"

#include "CCPACSWingSparPositionUIDs.h"
#include "CCPACSWingSparCrossSection.h"

namespace tigl 
{
class CCPACSWingSpars;
class CCPACSWingSparSegment
{
public:
    TIGL_EXPORT CCPACSWingSparSegment(CCPACSWingSpars* spars);
    
    TIGL_EXPORT void Reset();
    
    TIGL_EXPORT void ReadCPACS(TixiDocumentHandle tixiHandle, const std::string& sparSegmentXPath);
    TIGL_EXPORT void Invalidate();
    TIGL_EXPORT bool IsValid() const;

    // Gets the spar segment uid
    TIGL_EXPORT virtual const std::string& GetUID(void) const;

    // Get spar cross section
    TIGL_EXPORT CCPACSWingSparCrossSection& GetSparCrossSection(void);

    // Get spar position uids
    TIGL_EXPORT CCPACSWingSparPositionUIDs& GetSparPositionUIDs(void);

    // Method for returning midplane line of spar
    TIGL_EXPORT TopoDS_Wire GetSparMidplaneLine();

    // Builds the cutting geometry for the spar as well as the midplane line
    TIGL_EXPORT void BuildAuxiliaryGeometry();

private:
    std::string uid;
    std::string name;
    std::string description;
    CCPACSWingSparCrossSection sparCrossSection;
    CCPACSWingSparPositionUIDs sparPositionUIDs;
    TopoDS_Shape geometry;
    TopoDS_Wire m_sparMidplaneLine;

    bool m_auxGeomValid;
    TopoDS_Compound m_sparCutGeometry;

    CCPACSWingSpars* sparsNode;
    bool isvalid;
};

} // namespace tigl

#endif // CCPACSWINGSPARSEGMENT_H
