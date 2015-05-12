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

#ifndef CCPACSWINGSPARSEGMENTS_H
#define CCPACSWINGSPARSEGMENTS_H

#include <vector>
#include "tigl_internal.h"
#include "tixi.h"

#include "CCPACSWingSparSegment.h"

namespace tigl 
{
class CCPACSWingSpars;
class CCPACSWingSparSegments
{
private:
    // Typedef for a CCPACSWingSparSegment container to store the spar segments of a wing component segment.
    typedef std::vector<CCPACSWingSparSegment*> CCPACSWingSparSegmentContainer;

public:
    TIGL_EXPORT CCPACSWingSparSegments(CCPACSWingSpars* parent);
    
    TIGL_EXPORT void Reset();
    
    TIGL_EXPORT void ReadCPACS(TixiDocumentHandle tixiHandle, const std::string& sparSegmentsXPath);
    TIGL_EXPORT void Invalidate();
    TIGL_EXPORT bool IsValid() const;

    // Returns the total count of spar segments for the wing component segment
    TIGL_EXPORT int GetSparSegmentCount(void) const;

    // Returns the spar position for a given index.
    TIGL_EXPORT CCPACSWingSparSegment& GetSparSegment(int index) const;

    // Returns the spar position for a given UID.
    TIGL_EXPORT CCPACSWingSparSegment& GetSparSegment(const std::string& UID);


private:
    CCPACSWingSparSegmentContainer sparSegments;
    CCPACSWingSpars* parent;

    bool isvalid;
};

} // namespace tigl

#endif // CCPACSWINGSPARSEGMENTS_H
