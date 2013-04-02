/*
* Copyright (C) 2007-2013 German Aerospace Center (DLR/SC)
*
* Created: 2012-11-08 Martin Siggel <Martin.Siggel@dlr.de>
* Changed: $Id$
*
* Version: $Revision$
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
/**
* @file
* @brief  Implementation of an abstract class for phyical components (fuselage, wing, wheels etc...).
*/

#include "CTiglAbstractPhysicalComponent.h"
#include "CTiglError.h"

#include "TDocStd_Document.hxx"
#include "TCollection_ExtendedString.hxx"
#include "TCollection_HAsciiString.hxx"
#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFApp_Application.hxx"
#include "XCAFDoc_DocumentTool.hxx"
#include "TDataStd_Name.hxx"
#include "TDataXtd_Shape.hxx"

namespace tigl {

    CTiglAbstractPhysicalComponent::CTiglAbstractPhysicalComponent()
        : CTiglAbstractGeometricComponent()
    {
        childContainer.clear();
    }

    // Adds a child to this geometric component.
    void CTiglAbstractPhysicalComponent::AddChild(CTiglAbstractPhysicalComponent* componentPtr)
    {
        if (componentPtr == 0) {
            throw CTiglError("Error: Null pointer for component in CTiglAbstractGeometricComponent::AddChild", TIGL_NULL_POINTER);
        }

        childContainer.push_back(componentPtr);
    }

    // Returns a pointer to the list of children of a component.
    CTiglAbstractPhysicalComponent::ChildContainerType& CTiglAbstractPhysicalComponent::GetChildren(void)
    {
        return childContainer;
    }

    // Resets the geometric component.
    void CTiglAbstractPhysicalComponent::Reset(void)
    {
        CTiglAbstractGeometricComponent::Reset();
        childContainer.clear();
        SetParentUID("");
    }

    // Returns the parent unique id
    std::string& CTiglAbstractPhysicalComponent::GetParentUID(void)
    {
        return parentUID;
    }

    // Sets the parent uid.
    void CTiglAbstractPhysicalComponent::SetParentUID(const std::string& parentUID)
    {
        this->parentUID = parentUID;
    }

    void CTiglAbstractPhysicalComponent::SetSymmetryAxis(const std::string& axis)
    {
        CTiglAbstractGeometricComponent::SetSymmetryAxis(axis);

        ChildContainerType::iterator it = childContainer.begin();
        for(; it != childContainer.end(); ++it){
            CTiglAbstractPhysicalComponent * pChild = *it;
            if(pChild) pChild->SetSymmetryAxis(axis);
        }
    }


    TDF_Label CTiglAbstractPhysicalComponent::ExportDataStructure(Handle_XCAFDoc_ShapeTool &myAssembly, TDF_Label& label)
    {
        // Define layer for the wing
//        TDF_Label wingLabel= label.NewChild();
//        Handle_XCAFDoc_ShapeTool hShapeTool = XCAFDoc_DocumentTool::ShapeTool(wingLabel);

        // This component
        TDF_Label aLabel = myAssembly->AddShape(GetLoft(), false);
        TDataStd_Name::Set (aLabel, GetUID().c_str());
        Handle_XCAFDoc_ShapeTool hShapeTool = XCAFDoc_DocumentTool::ShapeTool(aLabel);

        // Other (sub)-components
        ChildContainerType::iterator it = childContainer.begin();
        for(; it != childContainer.end(); ++it){
            CTiglAbstractPhysicalComponent * pChild = *it;
            if(pChild) TDF_Label newLabel = pChild->ExportDataStructure(hShapeTool, aLabel);
        }

        return aLabel;
    }

} // namespace tigl

