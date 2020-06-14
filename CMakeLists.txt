# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(arich)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

set(VGM_INCLUDE_DIR /home/burmist/VGM/vgm-install/include/)
set(VGM_LIBRARY_DIR /home/burmist/VGM/vgm-install/lib/)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories("/home/leonid/home2/KEK/arich_standalone_G4/xmlarichdata/include/")
include_directories($ENV{XMLARICHDATA_BASF2_INCLUDE})
include_directories($ENV{XMLARICHDATAHOME}/include/)
include_directories(${VGM_INCLUDE_DIR})


EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")



#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
#----------------------------------------------------------------------------
# Find ROOT (optional package)
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /etc/root/cmake/)
#include(FindROOT.cmake)
find_package(ROOT)
if(ROOT_FOUND)
  add_definitions(-DG4ANALYSIS_USE_ROOT)
else()
  #set(ROOTSYS /home/leonid/root_v5.34.34/root_v5.34.34-install)
  #set(ROOTSYS /home/leonid/home2/root_v6.08.06/root_v6.08.06-install/)
  #set(ROOT_VERSION 5.34/34)
  #set(ROOT_VERSION 6.08/06)
  set(ROOT_INCLUDE_DIR ${ROOTSYS}/include/root/)
  set(ROOT_LIBRARY_DIR ${ROOTSYS}/lib/root/)
  set(ROOT_LIBRARIES -LIBPATH:${ROOT_LIBRARY_DIR} -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGeom -pthread -lm -ldl -rdynamic)
endif()

#----------------------------------------------------------------------------
# Add shared library
#
##Not working
#add_library(testlib SHARED IMPORTED)
#set_property(TARGET testlib PROPERTY IMPORTED_LOCATION "/home/leonid/home2/KEK/arich_standalone_G4/xmlarichdata/obj/libxmlarichdata.so")
#set(MYANA_INCLUDE_DIRS "/home/leonid/home2/KEK/arich_standalone_G4/xmlarichdata/include/")
#set(MYANA_LIBRARY      "/home/leonid/home2/KEK/arich_standalone_G4/xmlarichdata/lib/libxmlarichdata.so")
#include_directories("${MYANA_INCLUDE_DIRS}")
#target_link_libraries(myapplication "${MYANA_LIBRARY}")



#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
MESSAGE( STATUS "Geant4_INCLUDE_DIR:         " ${Geant4_INCLUDE_DIR} )
MESSAGE( STATUS "ROOT_INCLUDE_DIR:           " ${ROOT_INCLUDE_DIR} )

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR}
		    ${MYANA_INCLUDE_DIRS})
                    
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(arich arich.cc ${sources} ${headers})

EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})

target_link_libraries(arich ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} "$ENV{XMLARICHDATAHOME}/lib/libxmlarichdata.so" "/home/burmist/VGM/vgm-install/lib/libBaseVGM.so" "/home/burmist/VGM/vgm-install/lib/libClhepVGM.so" "/home/burmist/VGM/vgm-install/lib/libGeant4GM.so" "/home/burmist/VGM/vgm-install/lib/libRootGM.so" "/home/burmist/VGM/vgm-install/lib/libXmlVGM.so")

MESSAGE( STATUS " --> $ENV{XMLARICHDATAHOME}/lib/libxmlarichdata.so" )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  run.mac
  vis.mac
  vis2.mac
  runarich.bash
  tessellatedSolidVertices.dat
  geomID7_VGM_vis.C
  geometry.root
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#file(COPY ${PROJECT_SOURCE_DIR}/anaarich DESTINATION ${PROJECT_BINARY_DIR}/../)
#file(COPY ${PROJECT_SOURCE_DIR}/anaarich DESTINATION ${PROJECT_BINARY_DIR}/)

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(xxx DEPENDS arich)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS arich DESTINATION bin)
