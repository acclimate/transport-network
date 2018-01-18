cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(cmake/helpers.cmake)

set_default_build_type(Release)

project(acclimate_transport)

add_executable(acclimate_transport src/main.cpp)

target_include_directories(acclimate_transport PRIVATE include)
target_compile_options(acclimate_transport PRIVATE -std=c++11)

set_advanced_cpp_warnings(acclimate_transport)
set_build_type_specifics(acclimate_transport)
add_git_version(acclimate_transport
  FALLBACK_VERSION 1.0.0)

option(PARALLELIZATION "" ON)
if(PARALLELIZATION)
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

option(WITH_TQDM "" ON)
if(WITH_TQDM)
  target_compile_definitions(acclimate_transport PRIVATE ACCLIMATE_TRANSPORT_WITH_TQDM)
  target_include_directories(acclimate_transport PRIVATE lib/tqdm.cpp/include)
endif()

#include NetCDF C library
find_package(NETCDF REQUIRED)
message(STATUS "NetCDF include directory: ${NETCDF_INCLUDE_DIR}")
message(STATUS "NetCDF library: ${NETCDF_LIBRARY}")
target_link_libraries(acclimate_transport netcdf)

#include NetCDF C++ library
find_package(NETCDF_CPP4 REQUIRED)
message(STATUS "NetCDF_c++4 include directory: ${NETCDF_CPP4_INCLUDE_DIR}")
message(STATUS "NetCDF_c++4 library: ${NETCDF_CPP4_LIBRARY}")
target_link_libraries(acclimate_transport netcdf_c++4)

#include GDAL library
find_package(GDAL REQUIRED)
message(STATUS "GDAL include directory: ${GDAL_INCLUDE_DIR}")
message(STATUS "GDAL library: ${GDAL_LIBRARY}")
mark_as_advanced(GDAL_CONFIG GDAL_INCLUDE_DIR GDAL_LIBRARY)
target_link_libraries(acclimate_transport gdal)
target_include_directories(acclimate_transport PRIVATE ${GDAL_INCLUDE_DIR})

add_cpp_tools(acclimate_transport)
