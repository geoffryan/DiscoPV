cmake_minimum_required(VERSION 3.3)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

if(NOT DEFINED CMAKE_MACOSX_RPATH)
    set(CMAKE_MACOSX_RPATH 0)
endif()

ADD_PARAVIEW_PLUGIN(DiscoReader "1.0"
    SERVER_MANAGER_XML DiscoReader.xml
    SERVER_MANAGER_SOURCES vtkDiscoReader.cxx
    REQUIRED_ON_SERVER)
