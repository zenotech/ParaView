# This plugin can be loaded directly as the XML or compiled into a binary and
# then the binary can be loaded in ParaView.
cmake_minimum_required(VERSION 3.12.4)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

SET(Boost_NO_SYSTEM_PATHS ON)
FIND_PACKAGE(Boost 1.50 REQUIRED COMPONENTS python)

vtk_module_load(vtkhdf5)

INCLUDE_DIRECTORIES(
  ${VTK_INCLUDE_DIRS}
  ${vtkhdf5_INCLUDE_DIRS}
  ${ParaView_SOURCE_DIR}/ThirdParty
  )

add_paraview_plugin(
  zCFDReader "1.0"
  SERVER_MANAGER_XML zCFDReader.xml
  SERVER_MANAGER_SOURCES vtkzCFDReader.cxx
  REQUIRED_ON_SERVER
  )

target_link_libraries(zCFDReader LINK_PRIVATE
  ${vtkhdf5_LIBRARIES}
  vtksys vtkjsoncpp
  Boost::python
)

if(PARAVIEW_USE_MPI)
  include(vtkMPI)
  vtk_mpi_link(zCFDReader)
endif()