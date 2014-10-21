# This plugin can be loaded directly as the XML or compiled into a binary and
# then the binary can be loaded in ParaView.
cmake_minimum_required(VERSION 2.8)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

INCLUDE_DIRECTORIES(
  ${VTK_INCLUDE_DIRS}
  )

add_paraview_plugin(
  zCFDReader "1.0"
  SERVER_MANAGER_XML zCFDReader.xml
  SERVER_MANAGER_SOURCES vtkzCFDReader.cxx
  REQUIRED_ON_SERVER)
