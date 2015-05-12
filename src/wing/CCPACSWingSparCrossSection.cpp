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

#include "CCPACSWingSparCrossSection.h"

#include "CTiglError.h"
#include "CTiglLogging.h"

namespace tigl 
{

CCPACSWingSparCrossSection::CCPACSWingSparCrossSection()
{
    Reset();
}

void CCPACSWingSparCrossSection::Reset()
{
    web1.Reset();
    Invalidate();
}

void CCPACSWingSparCrossSection::ReadCPACS(TixiDocumentHandle tixiHandle, const std::string &sparCrossSectionXPath)
{
    Reset();
    
    // check path
    if ( tixiCheckElement(tixiHandle, sparCrossSectionXPath.c_str()) != SUCCESS) {
        LOG(ERROR) << "Wing Spar Cross Section " << sparCrossSectionXPath << " not found in CPACS file!" << std::endl;
        return;
    }
    
    // read rotation
    double rot;
    std::string rotationpath = sparCrossSectionXPath + "/rotation";
    if (tixiCheckElement(tixiHandle, rotationpath.c_str()) == SUCCESS) {
        tixiGetDoubleElement(tixiHandle, rotationpath.c_str(), &rot);
    }
    
    // check web path
    std::string web1path = sparCrossSectionXPath + "/web1";
    web1.ReadCPACS(tixiHandle, web1path.c_str());

    rotation = rot;
    isvalid = true;
}


void CCPACSWingSparCrossSection::Invalidate()
{
    isvalid = false;
}

bool CCPACSWingSparCrossSection::IsValid() const
{
    return isvalid;
}

double CCPACSWingSparCrossSection::GetRotation()
{
    return rotation;
}

} // namespace tigl
