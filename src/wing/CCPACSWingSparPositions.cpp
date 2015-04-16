/* 
* Copyright (C) 2007-2015 German Aerospace Center (DLR/SC)
*
* Created: 2015-04-16 Jonas Jepsen <Jonas.Jepsen@dlr.de>
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

#include "CCPACSWingSparPositions.h"

#include "CTiglError.h"
#include "CTiglLogging.h"


namespace tigl 
{

CCPACSWingSparPositions::CCPACSWingSparPositions()
{
    Reset();
}

void CCPACSWingSparPositions::Reset()
{
    //children.Reset();
    Invalidate();
}

void CCPACSWingSparPositions::ReadCPACS(TixiDocumentHandle tixiHandle, const std::string &sparPositionsXPath)
{
    Reset();
    
    // check path
    if ( tixiCheckElement(tixiHandle, sparPositionsXPath.c_str()) != SUCCESS) {
        LOG(ERROR) << "Wing Spar Positions " << sparPositionsXPath << " not found in CPACS file!" << std::endl;
        return;
    }
    
    int nsparPositions = 0;
    if (tixiGetNamedChildrenCount(tixiHandle, sparPositionsXPath.c_str(), "sparPosition", &nsparPositions) != SUCCESS) {
        // no spar positions found
        return;
    }

    for (int isparPosition = 1; isparPosition <= nsparPositions; ++isparPosition) {
        std::stringstream stream;
        stream << sparPositionsXPath << "/" << "sparPosition[" << isparPosition << "]";

        // check path
        if ( tixiCheckElement(tixiHandle, stream.str().c_str()) == SUCCESS) {
            CCPACSWingSparPosition * sparPosition = new CCPACSWingSparPosition();
            sparPosition->ReadCPACS(tixiHandle, stream.str().c_str());
            sparPositions.push_back(sparPosition);
        }
    }
    
    isvalid = true;
}


void CCPACSWingSparPositions::Invalidate()
{
    isvalid = false;
}

bool CCPACSWingSparPositions::IsValid() const
{
    return isvalid;
}

int CCPACSWingSparPositions::GetSparPositionCount() const
{
    return sparPositions.size();
}

CCPACSWingSparPosition& CCPACSWingSparPositions::GetSparPosition(int index) const
{
    if (index < 1 || index > GetSparPositionCount()) {
        throw CTiglError("Illegal index in CCPACSWingSparPositions::GetSparPosition", TIGL_INDEX_ERROR);
    }

    return *sparPositions.at(index-1);
}

// Gets a spar position by uid.
CCPACSWingSparPosition& CCPACSWingSparPositions::GetSparPosition(const std::string& sparPositionUID)
{
    for (CCPACSWingSparPositionContainer::size_type i = 0; i < sparPositions.size(); i++) {
        if (sparPositions[i]->GetUID() == sparPositionUID) {
            return (CCPACSWingSparPosition &) (*(sparPositions[i]));
        }
    }
    throw CTiglError("Error: Invalid uid in CCPACSWingSparPositions::GetSparPosition", TIGL_UID_ERROR);
}

} // namespace tigl
