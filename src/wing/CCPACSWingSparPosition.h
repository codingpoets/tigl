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

#ifndef CCPACSWINGSPARPOSITION_H
#define CCPACSWINGSPARPOSITION_H

#include <string>
#include "tigl_internal.h"
#include "tixi.h"

namespace tigl 
{

class CCPACSWingSparPosition
{
public:
    TIGL_EXPORT CCPACSWingSparPosition();
    
    TIGL_EXPORT void Reset();
    
    TIGL_EXPORT void ReadCPACS(TixiDocumentHandle tixiHandle, const std::string& sparPositionXPath);
    TIGL_EXPORT void Invalidate();
    TIGL_EXPORT bool IsValid() const;

    // Gets the spar position uid
    TIGL_EXPORT virtual const std::string& GetUID(void) const;

    TIGL_EXPORT double GetEta();
    TIGL_EXPORT double GetXsi();

private:
    std::string uid;
    double eta;
    double xsi;
    bool isvalid;
};

} // namespace tigl

#endif // CCPACSWINGSPARPOSITION_H
