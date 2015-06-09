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

#include "CCPACSWingSparPositionUIDs.h"

#include "CTiglError.h"
#include "CTiglLogging.h"

namespace tigl 
{

CCPACSWingSparPositionUIDs::CCPACSWingSparPositionUIDs()
{
    Reset();
}

void CCPACSWingSparPositionUIDs::Reset()
{
    //children.Reset();
    Invalidate();
}

void CCPACSWingSparPositionUIDs::ReadCPACS(TixiDocumentHandle tixiHandle, const std::string &sparPositionUIDsXPath)
{
    Reset();
    
    // check path
    if ( tixiCheckElement(tixiHandle, sparPositionUIDsXPath.c_str()) != SUCCESS) {
        LOG(ERROR) << "Wing Spar Position UIDs " << sparPositionUIDsXPath << " not found in CPACS file!" << std::endl;
        return;
    }
    
    int nsparPositions = 0;
    if (tixiGetNamedChildrenCount(tixiHandle, sparPositionUIDsXPath.c_str(), "sparPositionUID", &nsparPositions) != SUCCESS) {
        // no spar positions found
        return;
    }

    for (int isparPosition = 1; isparPosition <= nsparPositions; ++isparPosition) {
        std::stringstream stream;
        stream << sparPositionUIDsXPath << "/" << "sparPositionUID[" << isparPosition << "]";

        // check path
        if ( tixiCheckElement(tixiHandle, stream.str().c_str()) == SUCCESS) {
            char* sparPositionUID = NULL;
            tixiGetTextElement(tixiHandle, stream.str().c_str(), &sparPositionUID);

            sparPositionUIDs.push_back((std::string)sparPositionUID);
        }
    }

    isvalid = true;
}


void CCPACSWingSparPositionUIDs::Invalidate()
{
    isvalid = false;
}

bool CCPACSWingSparPositionUIDs::IsValid() const
{
    return isvalid;
}

int CCPACSWingSparPositionUIDs::GetSparPositionUIDCount() const
{
    return sparPositionUIDs.size();
}

std::string CCPACSWingSparPositionUIDs::GetSparPositionUID(int index) const
{
    if (index < 1 || index > GetSparPositionUIDCount()) {
        throw CTiglError("Illegal index in CCPACSWingSparPositionUIDs::GetSparPositionUID", TIGL_INDEX_ERROR);
    }

    return sparPositionUIDs.at(index-1);
}

} // namespace tigl
