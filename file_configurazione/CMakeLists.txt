CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

PROJECT(PACS)

SET(CMAKE_CXX_FLAGS "-std=c++0x -Wall ${CMAKE_CXX_FLAGS}")
set(CMAKE_AUTOMOC ON)

INCLUDE_DIRECTORIES(/usr/include)
INCLUDE_DIRECTORIES(/usr/local/lib)
INCLUDE_DIRECTORIES(/opt/local/include)
INCLUDE_DIRECTORIES(/sw/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Include Eigen3

FIND_PACKAGE(PkgConfig)
PKG_CHECK_MODULES(EIGEN3 REQUIRED eigen3)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})

#FIND_PACKAGE(Eigen3)
#IF (EIGEN3_FOUND)
#INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
#ENDIF (EIGEN3_FOUND)


# Include Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)


# Include Getfem
if (GETFEM_LIBRARIES AND GETFEM_INCLUDE_DIRS)
# in cache already
	set(GETFEM_FOUND TRUE)
else (GETFEM_LIBRARIES AND GETFEM_INCLUDE_DIRS)
	find_path(GETFEM_INCLUDE_DIR
	NAMES
		getfem/getfem_mesh.h
	PATHS
		"${BASE_DIR}/include"
		/usr/include
		/usr/local/include
		/opt/local/include
		/sw/include
	)

find_library(GETFEM_LIBRARY
	NAMES
		getfem
	PATHS
		"${BASE_DIR}/lib"
		/usr/lib
		/usr/local/lib
		/opt/local/lib
		/sw/lib

	)
	
set(GETFEM_INCLUDE_DIRS
	${GETFEM_INCLUDE_DIR}
	)
set(GETFEM_LIBRARIES
	${GETFEM_LIBRARY}
	)

if (GETFEM_INCLUDE_DIRS AND GETFEM_LIBRARIES)
	set(GETFEM_FOUND TRUE)
endif (GETFEM_INCLUDE_DIRS AND GETFEM_LIBRARIES)

if (GETFEM_FOUND)
	if (NOT GetFEM_FIND_QUIETLY)
	message(STATUS "Found GetFEM: ${GETFEM_LIBRARIES}")
	endif (NOT GetFEM_FIND_QUIETLY)
else (GETFEM_FOUND)
	if (GetFEM_FIND_REQUIRED)
		message(FATAL_ERROR "Could not find GetFEM")
	endif (GetFEM_FIND_REQUIRED)
endif (GETFEM_FOUND)

# show the GETFEM_INCLUDE_DIRS and GETFEM_LIBRARIES variables only in the advanced view
mark_as_advanced(GETFEM_INCLUDE_DIRS GETFEM_LIBRARIES)
	
endif (GETFEM_LIBRARIES AND GETFEM_INCLUDE_DIRS)



#Include for BLAS library
find_package ( BLAS )
if (BLAS_FOUND)
  if (BLA_VENDOR MATCHES "^ACML")
    set (HAVE_ACML 1)#TODO: for Windows only
    set (HAVE_BLAS_ACML 1)
    set (HAVE_BLAS 1)
  elseif (BLA_VENDOR MATCHES "^Intel")
    set (HAVE_MKL 1)#TODO: for Windows only
    set (HAVE_BLAS_MKL 1)
    set (HAVE_BLAS 1)
  elseif (BLA_VENDOR MATCHES "^ATLAS")
    set (HAVE_BLAS_ATLAS 1)
    set (HAVE_BLAS 1)
  elseif( (BLA_VENDOR STREQUAL "Generic") OR (BLA_VENDOR STREQUAL "All") )
    set (HAVE_BLAS 1)
  else()
    message (WARNING "Unknown BLAS vendor: ${BLA_VENDOR}")
    set(BLAS_FOUND FALSE)
  endif()
endif(BLAS_FOUND)

#check for LAPACK library
find_package ( LAPACK )
if (LAPACK_FOUND)
  if ((BLA_VENDOR MATCHES "^ACML") OR (BLA_VENDOR MATCHES "^Intel") OR (BLA_VENDOR MATCHES "^ATLAS") OR ( BLA_VENDOR STREQUAL "Generic" ) OR (BLA_VENDOR STREQUAL "All"))
    set (HAVE_LAPACK 1)
  else()
    message (WARNING "Unknown LAPACK vendor: ${BLA_VENDOR}")
    set(LAPACK_FOUND FALSE)
  endif()
endif(LAPACK_FOUND)


#Include for QHULL library
set(QHULL_MAJOR_VERSION 6)

if(QHULL_USE_STATIC)
	set(QHULL_RELEASE_NAME qhullstatic)
	set(QHULL_DEBUG_NAME qhullstatic_d)
else(QHULL_USE_STATIC)
	set(QHULL_RELEASE_NAME qhull qhull${QHULL_MAJOR_VERSION})
	set(QHULL_DEBUG_NAME qhull_d qhull${QHULL_MAJOR_VERSION}_d qhull_d${QHULL_MAJOR_VERSION})
endif(QHULL_USE_STATIC)

find_file(QHULL_HEADER
	NAMES libqhull/libqhull.h qhull.h
	HINTS "${QHULL_ROOT}" "$ENV{QHULL_ROOT}" "${QHULL_INCLUDE_DIR}"
	PATHS "$ENV{PROGRAMFILES}/QHull" "$ENV{PROGRAMW6432}/QHull"
	PATH_SUFFIXES qhull src/libqhull libqhull include)

set(QHULL_HEADER "${QHULL_HEADER}" CACHE INTERNAL "QHull header" FORCE )

if(QHULL_HEADER)
	get_filename_component(qhull_header ${QHULL_HEADER} NAME_WE)
	if("${qhull_header}" STREQUAL "qhull")
		set(HAVE_QHULL_2011 OFF)
		get_filename_component(QHULL_INCLUDE_DIR ${QHULL_HEADER} PATH)
	elseif("${qhull_header}" STREQUAL "libqhull")
		set(HAVE_QHULL_2011 ON)
		get_filename_component(QHULL_INCLUDE_DIR ${QHULL_HEADER} PATH)
		get_filename_component(QHULL_INCLUDE_DIR ${QHULL_INCLUDE_DIR} PATH)
	endif()
else(QHULL_HEADER)
	set(QHULL_INCLUDE_DIR "QHULL_INCLUDE_DIR-NOTFOUND")
endif(QHULL_HEADER)

set(QHULL_INCLUDE_DIR "${QHULL_INCLUDE_DIR}" CACHE PATH "QHull include dir." FORCE)

find_library(QHULL_LIBRARY
	NAMES ${QHULL_RELEASE_NAME}
	HINTS "${QHULL_ROOT}" "$ENV{QHULL_ROOT}"
	PATHS "$ENV{PROGRAMFILES}/QHull" "$ENV{PROGRAMW6432}/QHull"
	PATH_SUFFIXES project build bin lib)

find_library(QHULL_LIBRARY_DEBUG
	NAMES ${QHULL_DEBUG_NAME} ${QHULL_RELEASE_NAME}
	HINTS "${QHULL_ROOT}" "$ENV{QHULL_ROOT}"
	PATHS "$ENV{PROGRAMFILES}/QHull" "$ENV{PROGRAMW6432}/QHull"
	PATH_SUFFIXES project build bin lib)

if(NOT QHULL_LIBRARY_DEBUG)
	set(QHULL_LIBRARY_DEBUG ${QHULL_LIBRARY})
endif(NOT QHULL_LIBRARY_DEBUG)

set(QHULL_INCLUDE_DIRS ${QHULL_INCLUDE_DIR})

set(QHULL_LIBRARIES optimized ${QHULL_LIBRARY} debug ${QHULL_LIBRARY_DEBUG})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qhull DEFAULT_MSG QHULL_LIBRARY QHULL_INCLUDE_DIR)

mark_as_advanced(QHULL_LIBRARY QHULL_LIBRARY_DEBUG QHULL_INCLUDE_DIR)

if(QHULL_FOUND)
	set(HAVE_QHULL ON)
	if(NOT QHULL_USE_STATIC)
		add_definitions("-Dqh_QHpointer")
		if(MSVC)
			add_definitions("-Dqh_QHpointer_dllimport")
		endif(MSVC)
	endif(NOT QHULL_USE_STATIC)

	message(STATUS "QHULL found (include: ${QHULL_INCLUDE_DIRS}, lib: ${QHULL_LIBRARIES})")
endif(QHULL_FOUND)

SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install-dir CACHE PATH "" FORCE)


ADD_SUBDIRECTORY(src)

ADD_SUBDIRECTORY(test)

