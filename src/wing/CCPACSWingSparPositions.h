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

#ifndef CCPACSWINGSPARPOSITIONS_H
#define CCPACSWINGSPARPOSITIONS_H

#include <vector>
#include "tigl_internal.h"
#include "tixi.h"

#include "CCPACSWingSparPosition.h"

namespace tigl 
{

class CCPACSWingSparPositions
{
private:
    // Typedef for a CCPACSWingSparPosition container to store the spar positions of a wing component segment.
    typedef std::vector<CCPACSWingSparPosition*> CCPACSWingSparPositionContainer;

public:
    TIGL_EXPORT CCPACSWingSparPositions();
    
    TIGL_EXPORT void Reset();
    
    TIGL_EXPORT void ReadCPACS(TixiDocumentHandle tixiHandle, const std::string& sparPositionsXPath);
    TIGL_EXPORT void Invalidate();
    TIGL_EXPORT bool IsValid() const;

    // Returns the total count of spar positions for the wing component segment
    TIGL_EXPORT int GetSparPositionCount(void) const;

    // Returns the spar position for a given index.
    TIGL_EXPORT CCPACSWingSparPosition& GetSparPosition(int index) const;

    // Returns the spar position for a given UID.
    TIGL_EXPORT CCPACSWingSparPosition& GetSparPosition(const std::string& UID);


private:
    CCPACSWingSparPositionContainer sparPositions;

    bool isvalid;
};

} // namespace tigl

#endif // CCPACSWINGSPARPOSITIONS_H
