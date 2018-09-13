# This plugin can be loaded directly as the XML or compiled into a binary and
# then the binary can be loaded in ParaView.
cmake_minimum_required(VERSION 2.8)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

SET(Boost_ADDITIONAL_VERSIONS "1.50" ) 
SET(Boost_NO_SYSTEM_PATHS ON)
FIND_PACKAGE(Boost 1.50 REQUIRED COMPONENTS python)
LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS} )
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )


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
  ${Boost_PYTHON_LIBRARY} ${PYTHON_LIBRARIES}
  )

if(PARAVIEW_USE_MPI)
  include(vtkMPI)
  vtk_mpi_link(zCFDReader)
endif()