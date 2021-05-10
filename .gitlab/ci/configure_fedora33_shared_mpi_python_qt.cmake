set(PARAVIEW_BUILD_EXAMPLES ON CACHE BOOL "")

set(PARAVIEW_ENABLE_ADIOS2 ON CACHE BOOL "")
set(PARAVIEW_ENABLE_FFMPEG ON CACHE BOOL "")
set(PARAVIEW_ENABLE_FIDES ON CACHE BOOL "")
set(PARAVIEW_ENABLE_GDAL ON CACHE BOOL "")
# set(PARAVIEW_ENABLE_OPENTURNS ON CACHE BOOL "")
set(PARAVIEW_ENABLE_PDAL ON CACHE BOOL "")
set(PARAVIEW_ENABLE_VISITBRIDGE ON CACHE BOOL "")
set(PARAVIEW_ENABLE_XDMF3 ON CACHE BOOL "")

set(PARAVIEW_PLUGINS_DEFAULT ON CACHE BOOL "")

# External package settings.
# set(PARAVIEW_BUILD_WITH_EXTERNAL ON CACHE BOOL "")
# set(VTK_MODULE_USE_EXTERNAL_VTK_libharu OFF CACHE BOOL "")
# set(VTK_MODULE_USE_EXTERNAL_VTK_cgns OFF CACHE BOOL "")
# set(VTK_MODULE_USE_EXTERNAL_VTK_gl2ps OFF CACHE BOOL "")

include("${CMAKE_CURRENT_LIST_DIR}/configure_fedora33.cmake")
