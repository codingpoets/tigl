/* 
* Copyright (C) 2007-2015 German Aerospace Center (DLR/SC)
*
* Created: 2015-05-04 Jonas Jepsen <Jonas.Jepsen@dlr.de>
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

#include "test.h"
#include "tiglcommonfunctions.h"
#include "TopoDS_Edge.hxx"
#include "TopoDS_Face.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "gp_Pnt.hxx"

TEST(TiglCommon, GetEdgeLength_on_line)
{
    BRepBuilderAPI_MakeEdge edgeBuilder(gp_Pnt(0,0,1), gp_Pnt(1,0,1));
    ASSERT_EQ(1, GetEdgeLength(edgeBuilder.Edge()));
    
    BRepBuilderAPI_MakeEdge edgeBuilder2(gp_Pnt(2,0,1), gp_Pnt(1,2,1));
    ASSERT_NEAR(2.23607, GetEdgeLength(edgeBuilder2.Edge()), 0.00001);
}


TEST(TiglCommon, BuildFace_from_4_gp_Pnt)
{
    gp_Pnt p1(0,0,0);
    gp_Pnt p2(0,0,1);
    gp_Pnt p3(0,1,0);
    gp_Pnt p4(0,1,1);
    TopoDS_Face face;
    face = BuildFace(p1, p2, p3, p4);
    ASSERT_FALSE(face.IsNull());
}


