/*
* Copyright (C) 2007-2013 German Aerospace Center (DLR/SC)
*
* Created: 2013-10-29 Martin Siggel <Martin.Siggel@dlr.de>
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

#include "TIGLViewerLoggerHTMLDecorator.h"
#include <string>

TIGLViewerLoggerHTMLDecorator::TIGLViewerLoggerHTMLDecorator(tigl::ITiglLogger* logger) 
    : _mylogger(logger)
{
}

TIGLViewerLoggerHTMLDecorator::~TIGLViewerLoggerHTMLDecorator() {
    if(_mylogger) {
        delete _mylogger;
    }
}

void TIGLViewerLoggerHTMLDecorator::LogMessage(TiglLogLevel level, const char * message){
    if(!_mylogger){
        return;
    }
    
    std::string newmsg;
    switch(level) {
    case TILOG_ERROR:
        newmsg = "<font color=\"red\">" + std::string(message) + "</font>";
        break;
    case TILOG_WARNING:
        newmsg = "<font color=\"orange\">" + std::string(message)  + "</font>";
        break;
    default:
        newmsg = message;
        break;
    }
    newmsg = "<i>" + newmsg + "</i>";
    _mylogger->LogMessage(level, newmsg.c_str());
}
void TIGLViewerLoggerHTMLDecorator::SetVerbosity(TiglLogLevel vlevel)
{
}
