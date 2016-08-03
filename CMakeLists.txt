# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( example )


cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_CXX_FLAGS "-std=c++1y -stdlib=libc++")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

find_package(Eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})


# include helper file
include( ${CGAL_USE_FILE} )

add_definitions(-DCGAL_EIGEN3_ENABLED)



# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package
include_directories( BEFORE ../include )


# Creating entries for target: example
# ############################

add_executable( iter  iteration.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS iter )

# Link the executable to CGAL and third-party libraries
target_link_libraries(iter   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

