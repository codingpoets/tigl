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

#include "CCPACSWingSparSegments.h"

#include "CTiglError.h"
#include "CTiglLogging.h"


namespace tigl 
{

CCPACSWingSparSegments::CCPACSWingSparSegments(CCPACSWingSpars* parent)
    : parent(parent)
{
    Reset();
}

void CCPACSWingSparSegments::Reset()
{
    //children.Reset();
    Invalidate();
}

void CCPACSWingSparSegments::ReadCPACS(TixiDocumentHandle tixiHandle, const std::string &sparSegmentsXPath)
{
    Reset();
    
    // check path
    if ( tixiCheckElement(tixiHandle, sparSegmentsXPath.c_str()) != SUCCESS) {
        LOG(ERROR) << "Wing Spar Segments " << sparSegmentsXPath << " not found in CPACS file!" << std::endl;
        return;
    }
    
    int nsparSegments = 0;
    if (tixiGetNamedChildrenCount(tixiHandle, sparSegmentsXPath.c_str(), "sparSegment", &nsparSegments) != SUCCESS) {
        // no spar segment found
        return;
    }

    for (int isparSegment = 1; isparSegment <= nsparSegments; ++isparSegment) {
        std::stringstream stream;
        stream << sparSegmentsXPath << "/" << "sparSegment[" << isparSegment << "]";

        // check path
        if ( tixiCheckElement(tixiHandle, stream.str().c_str()) == SUCCESS) {
            CCPACSWingSparSegment * sparSegment = new CCPACSWingSparSegment(parent);
            sparSegment->ReadCPACS(tixiHandle, stream.str().c_str());
            sparSegments.push_back(sparSegment);
        }
    }
    
    isvalid = true;
}


void CCPACSWingSparSegments::Invalidate()
{
    isvalid = false;
}

bool CCPACSWingSparSegments::IsValid() const
{
    return isvalid;
}

int CCPACSWingSparSegments::GetSparSegmentCount() const
{
    return sparSegments.size();
}

CCPACSWingSparSegment& CCPACSWingSparSegments::GetSparSegment(int index) const
{
    if (index < 1 || index > GetSparSegmentCount()) {
        throw CTiglError("Illegal index in CCPACSWingSparSegments::GetSparSegment", TIGL_INDEX_ERROR);
    }

    return *sparSegments.at(index-1);
}

// Gets a spar segment by uid.
CCPACSWingSparSegment& CCPACSWingSparSegments::GetSparSegment(const std::string& sparSegmentUID)
{
    for (CCPACSWingSparSegmentContainer::size_type i = 0; i < sparSegments.size(); i++) {
        if (sparSegments[i]->GetUID() == sparSegmentUID) {
            return (CCPACSWingSparSegment &) (*(sparSegments[i]));
        }
    }
    throw CTiglError("Error: Invalid uid in CCPACSWingSparSegments::GetSparSegment", TIGL_UID_ERROR);
}

} // namespace tigl
