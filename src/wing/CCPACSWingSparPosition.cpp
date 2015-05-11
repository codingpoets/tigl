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

#include "CCPACSWingSparPosition.h"

#include "CTiglError.h"
#include "CTiglLogging.h"


namespace tigl 
{

CCPACSWingSparPosition::CCPACSWingSparPosition()
{
    Reset();
}

void CCPACSWingSparPosition::Reset()
{
    //children.Reset();
    Invalidate();
}

void CCPACSWingSparPosition::ReadCPACS(TixiDocumentHandle tixiHandle, const std::string &sparPositionXPath)
{
    Reset();
    
    // check path
    if ( tixiCheckElement(tixiHandle, sparPositionXPath.c_str()) != SUCCESS) {
        LOG(ERROR) << "SparPosition " << sparPositionXPath << " not found in CPACS file!" << std::endl;
        return;
    }

    // Get UID
    char * nameStr = NULL;
    if ( tixiGetTextAttribute(tixiHandle, sparPositionXPath.c_str(), "uID", &nameStr) != SUCCESS ) {
        throw tigl::CTiglError("No UID given for wing spar position " + sparPositionXPath + "!", TIGL_UID_ERROR);
    }

    double etaRead, xsiRead;

    // get parameters of spar position
    std::string tmpString;
    tmpString = sparPositionXPath + "/elementUID";
    if ( tixiCheckElement(tixiHandle, tmpString.c_str()) == SUCCESS) {
        LOG(WARNING) << "In " << sparPositionXPath << ": Spar positionings via elements is currently not supported by TiGL. Please use eta/xsi definitions.";
        etaRead = 0.;
    }
    else {
        tmpString = sparPositionXPath + "/eta";
        if ( tixiGetDoubleElement(tixiHandle, tmpString.c_str(), &etaRead) != SUCCESS) {
            throw tigl::CTiglError("No eta positioning given for wing spar position " + sparPositionXPath + "!", TIGL_ERROR);
        }
    }

    tmpString = sparPositionXPath + "/xsi";
    if ( tixiGetDoubleElement(tixiHandle, tmpString.c_str(), &xsiRead) != SUCCESS) {
        throw tigl::CTiglError("No xsi positioning given for wing spar position " + sparPositionXPath + "!", TIGL_ERROR);
    }
    
    uid = nameStr;
    eta = etaRead;
    xsi = xsiRead;
    isvalid = true;
}


void CCPACSWingSparPosition::Invalidate()
{
    isvalid = false;
}

bool CCPACSWingSparPosition::IsValid() const
{
    return isvalid;
}

// Gets the spar position uid
const std::string& CCPACSWingSparPosition::GetUID(void) const
{
    return uid;
}

double CCPACSWingSparPosition::GetEta()
{
    return eta;
}

double CCPACSWingSparPosition::GetXsi()
{
    return xsi;
}

} // namespace tigl
