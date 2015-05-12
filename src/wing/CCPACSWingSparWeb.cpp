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

#include "CCPACSWingSparWeb.h"

#include "CTiglError.h"
#include "CTiglLogging.h"

namespace tigl 
{

CCPACSWingSparWeb::CCPACSWingSparWeb()
{
    Reset();
}

void CCPACSWingSparWeb::Reset()
{
    Invalidate();
}

void CCPACSWingSparWeb::ReadCPACS(TixiDocumentHandle tixiHandle, const std::string &sparWebXPath)
{
    Reset();
    
    // check path
    if ( tixiCheckElement(tixiHandle, sparWebXPath.c_str()) != SUCCESS) {
        LOG(ERROR) << "Wing Spar Web " << sparWebXPath << " not found in CPACS file!" << std::endl;
        return;
    }
    
    // read relative position
    double pos;
    std::string relPosPath = sparWebXPath + "/relPos";
    if (tixiCheckElement(tixiHandle, relPosPath.c_str()) == SUCCESS) {
        tixiGetDoubleElement(tixiHandle, relPosPath.c_str(), &pos);
    }
    
    relPos = pos;
    isvalid = true;
}


void CCPACSWingSparWeb::Invalidate()
{
    isvalid = false;
}

bool CCPACSWingSparWeb::IsValid() const
{
    return isvalid;
}

} // namespace tigl
