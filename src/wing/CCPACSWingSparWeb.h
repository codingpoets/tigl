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

#ifndef CCPACSWINGSPARWEB_H
#define CCPACSWINGSPARWEB_H

#include <string>
#include "tigl_internal.h"
#include "tixi.h"

namespace tigl
{

class CCPACSWingSparWeb
{
public:
    TIGL_EXPORT CCPACSWingSparWeb();

    TIGL_EXPORT void Reset();

    TIGL_EXPORT void ReadCPACS(TixiDocumentHandle tixiHandle, const std::string& sparWebXPath);
    TIGL_EXPORT void Invalidate();
    TIGL_EXPORT bool IsValid() const;

    TIGL_EXPORT double GetRelPos(void);

private:
    double relPos;

    bool isvalid;
};

} // namespace tigl

#endif // CCPACSWINGSPARWEB_H
