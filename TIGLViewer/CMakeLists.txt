# TIGLViewer cmake project
# author: Martin Siggel

# set up QT

OPTION(TIGL_VIEWER "Builds the TiGL Viewer program (requires OpenGL and Qt4)" ON)

set(PROGNAME "TIGLViewer")


if(TIGL_VIEWER)

include(SearchQt)
# disable automatic linking of qtmain
if (POLICY CMP0020)
  cmake_policy(SET CMP0020 OLD)
endif (POLICY CMP0020)

#build only tiglviewer if qt was found
if(QT4_FOUND OR Qt5Widgets_FOUND)

include(tiglmacros)
AddToCheckstyle()

# required opencascade libraries
set (TV_OCC_LIBS 
    TKOpenGl 
    ${OCC_VISU_LIBS} ${OCC_LIBS}
)

find_package(OCE 0.15 COMPONENTS ${TV_OCC_LIBS} QUIET)
if(OCE_FOUND)
  set(OpenCASCADE_LIBRARIES ${TV_OCC_LIBS})
  set(OpenCASCADE_INCLUDE_DIR ${OCE_INCLUDE_DIRS})
else(OCE_FOUND)
  set(OpenCASCADE_LIBRARIES "")
  find_package(OpenCASCADE REQUIRED COMPONENTS ${TV_OCC_LIBS})
endif(OCE_FOUND)

# we need opengl in order to build tiglviewer
find_package(OpenGL REQUIRED)

if(CMAKE_SYSTEM_NAME STREQUAL Linux)
  # set library search path to local lib directory during installation
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
  set(RELDIR "../${CMAKE_INSTALL_LIBDIR}")

  if (IS_ABSOLUTE ${CMAKE_INSTALL_LIBDIR})
    # get relative path to tiglviewer binary
    file(RELATIVE_PATH RELDIR ${CMAKE_INSTALL_PREFIX}/bin "${CMAKE_INSTALL_LIBDIR}")
  endif()
  
  set(CMAKE_INSTALL_RPATH \$ORIGIN/${RELDIR})

  if (NOT RELDIR STREQUAL "../lib")
    message(STATUS ${RELDIR})
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:\$ORIGIN/../lib")
  endif()
endif() # linux‚

# add all tigl directories to include path
set(TIGL_INCLUDES ${PROJECT_SOURCE_DIR}/src/)
subdirlist(SUBDIRS ${PROJECT_SOURCE_DIR}/src/)
foreach(subdir ${SUBDIRS})
  set(TIGL_INCLUDES ${TIGL_INCLUDES} ${PROJECT_SOURCE_DIR}/src/${subdir})
endforeach()

include_directories(${TIGL_INCLUDES})

set(tv_src
    CommandLineParameters.cpp
    ISession_Point.cpp
    ISession_Direction.cpp
    ISession_Text.cpp
    main.cpp
    TIGLViewerContext.cpp
    TIGLViewerDocument.cpp
    TIGLViewerInputoutput.cpp
    TIGLViewerInternal.h.cpp
    TIGLViewerWidget.cpp
    TIGLViewerWindow.cpp
    TIGLViewerSettingsDialog.cpp
    TIGLScriptProxy.cpp
    TIXIScriptProxy.cpp
    TIGLScriptEngine.cpp
    TiglViewerConsole.cpp
    TIGLViewerSettings.cpp
    TIGLViewerControlFile.cpp
    TIGLViewerEtaXsiDialog.cpp
    CHotsoseMeshReader.cpp
    TIGLAISTriangulation.cpp
    TIGLViewerErrorDialog.cpp
    TIGLViewerLogHistory.cpp
    TIGLViewerLogRedirection.cpp
    TIGLViewerLoggerHTMLDecorator.cpp
    TIGLViewerFuseDialog.cpp
    TIGLViewerShapeIntersectionDialog.cpp
    TIGLViewerScreenshotDialog.cpp
    TIGLViewerDrawVectorDialog.cpp
    TIGLViewerScopedCommand.cpp
)
	
# normal header files
set(tv_hdr
    CommandLineParameters.h
    ISession_Point.h
    ISession_Direction.h
    ISession_Text.h
    TIGLViewer.h
    TIGLViewerSettings.h
    TIGLViewerControlFile.h
    TIGLAISTriangulation.h
    CHotsoseMeshReader.h
    TIGLViewerLogHistory.h
    TIGLViewerLoggerHTMLDecorator.h
    TIGLViewerScopedCommand.h
)

# qobject header files that have to be moced
set(tv_hdr_qobj
	TIGLDebugStream.h
	TIGLViewerContext.h
	TIGLViewerDocument.h
	TIGLViewerInputoutput.h
	TIGLViewerInputoutput.h
	TIGLViewerWidget.h
	TIGLViewerWindow.h
	TIGLViewerSettingsDialog.h
	TIGLScriptProxy.h
	TIXIScriptProxy.h
	TIGLScriptEngine.h
	TiglViewerConsole.h
	TIGLViewerEtaXsiDialog.h
	TIGLViewerErrorDialog.h
	TIGLViewerDrawVectorDialog.h
	TIGLViewerLogRedirection.h
	TIGLViewerFuseDialog.h
	TIGLViewerShapeIntersectionDialog.h
	TIGLViewerScreenshotDialog.h
)

set(tv_ui_comp
	src/TIGLViewerWindow.ui
	src/TIGLViewerSettingsDialog.ui
	src/TIGLViewerEtaXsiDialog.ui
	src/TIGLViewerErrorDialog.ui
	src/TIGLViewerFuseDialog.ui
	src/TIGLViewerShapeIntersectionDialog.ui
	src/TIGLViewerScreenshotDialog.ui
	src/TIGLViewerDrawVectorDialog.ui
)
	
foreach ( _fname ${tv_src})
   set(tv_src_comp ${tv_src_comp} src/${_fname})
endforeach()
foreach ( _fname ${tv_hdr_qobj})
   set(tv_hdr_qobj_comp ${tv_hdr_qobj_comp} src/${_fname})
endforeach()
foreach ( _fname ${tv_hdr})
   set(tv_hdr_comp ${tv_hdr_comp} src/${_fname})
endforeach()

# moc header files
QT4_WRAP_CPP(tv_hdr_MOC ${tv_hdr_qobj_comp})
QT4_WRAP_UI(tv_ui_HEADERS ${tv_ui_comp})
# compile icons
QT4_ADD_RESOURCES( tv_rcrs src/TIGLViewer.qrc )

add_definitions(-DCSFDB -DQOCC_STATIC)

include_directories( ${OpenCASCADE_INCLUDE_DIR} ${TIXI_INCLUDE_DIRS})

# opengl
include_directories(${OPENGL_INCLUDE_DIR})

# boost
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/boost_1_55_0)

include_directories( ${PROJECT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)

# include path for ui files
include_directories( ${CMAKE_CURRENT_BINARY_DIR} )

source_group(MOCS REGULAR_EXPRESSION "moc_")
source_group(RESOURCES REGULAR_EXPRESSION "qrc_")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if(WIN32)
  # comment this line if you want to see an additional console window
  set(TYPE WIN32)

  # add program icon
  set(tv_icon_rc src/TIGLViewer.rc)
endif(WIN32)

if(APPLE)
   set(TYPE MACOSX_BUNDLE)

   set( MACOSX_BUNDLE_SHORT_VERSION_STRING ${TIGL_VERSION} )
   set( MACOSX_BUNDLE_VERSION ${TIGL_VERSION} )
   set( MACOSX_BUNDLE_LONG_VERSION_STRING Version ${TIGL_VERSION} )

   # add program icon
   set(MACOSX_BUNDLE_ICON_FILE ${PROGNAME}.icns )
endif()

add_executable(${PROGNAME} ${TYPE}
		${tv_src_comp} 
		${tv_hdr_MOC} 
 		${tv_hdr_qobj_comp}        
		${tv_hdr_comp} 
		${tv_ui_comp}
		${tv_rcrs}
		${tv_icon_rc}
)

target_link_libraries(${PROGNAME} 
	TIGL_static
	${TIXI_LIBRARIES} 
	${CMake_QT_LIBRARIES} 
	${QT_QTMAIN_LIBRARY} 
	${OPENGL_LIBRARIES} 
	${OpenCASCADE_LIBRARIES} 
	${EXTRA_LIBRARIES} 
)

if(GLOG_FOUND)
  include_directories(${GLOG_INCLUDE_DIR})
  target_link_libraries (${PROGNAME} ${GLOG_LIBRARIES})
endif(GLOG_FOUND)

#--------------------------------------------------------------------------------
# Now the installation stuff below
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Install the TIGLViewer application, on Apple, the bundle is at the root of the
# install tree, and on other platforms it'll go into the bin directory.
INSTALL(TARGETS ${PROGNAME} 
    BUNDLE DESTINATION . COMPONENT viewer
    RUNTIME DESTINATION bin COMPONENT viewer
    )


set(APPS "bin/${PROGNAME}")
IF(APPLE)
  set(APPS "${PROGNAME}.app")

  # TKOpenGL is loaded as a plugin, hence we must manually install it
  INSTALL(FILES ${TKOpenGl_OCCLIB} DESTINATION ${PROGNAME}.app/Contents/MacOS/ COMPONENT viewer)
  # copy icon
  INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/gfx/${PROGNAME}.icns DESTINATION ${PROGNAME}.app/Contents/Resources/ COMPONENT viewer)
ENDIF(APPLE)
IF(WIN32)
  set(APPS "bin/${PROGNAME}.exe")

  # install runtime libs
  FILE(TO_CMAKE_PATH "$ENV{SYSTEMROOT}" SYSTEMROOT)

ENDIF(WIN32)

SET(CMAKE_INSTALL_OPENMP_LIBRARIES "ON")
INCLUDE(InstallRequiredSystemLibraries)

	
# set directories to look for dependencies
set(DIRS ${TIXI_LIB_PATH} ${OpenCASCADE_DLL_DIRECTORY}  ${OpenCASCADE_LINK_DIRECTORY} ${QT_BINARY_DIR} ${TIGL_ADD_LIB_PATH})

if (WIN32 OR APPLE)
  if(Qt5Widgets_FOUND)
    if(APPLE)
      set(QT_PLUGINS  "qcocoa;qtaccessiblewidgets")
    endif()
    if(WIN32)
      set(QT_PLUGINS  "qwindows")
    endif()

    include(DeployQt5)
    INSTALL_QT5_EXECUTABLE(${APPS} "${QT_PLUGINS}" "" "${DIRS}" "" "" viewer)
  else()
    include(DeployQt4)
    INSTALL_QT4_EXECUTABLE(${APPS} "" "" "${DIRS}" "" "" viewer)
  endif()
endif()

else (QT4_FOUND OR Qt5Widgets_FOUND)
message(STATUS "No qt was found. TIGLViewer will not be build!")
endif(QT4_FOUND OR Qt5Widgets_FOUND)

endif(TIGL_VIEWER)


