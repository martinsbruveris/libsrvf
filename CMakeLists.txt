cmake_minimum_required(VERSION 2.6)
project(LIBSRVF)

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

find_library(GSL gsl)

set(CMAKE_CXX_FLAGS "-fPIC")

add_subdirectory(src)
add_subdirectory(tests EXCLUDE_FROM_ALL)

install (DIRECTORY matlab DESTINATION share/libsrvf )

set (CPACK_GENERATOR "TGZ")
if (WIN32)
  list(APPEND CPACK_GENERATOR "NSIS")
endif (WIN32)
set(CPACK_SOURCE_GENERATOR "ZIP;TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR 0)
set(CPACK_PACKAGE_VERSION_MINOR 1)
set(CPACK_PACKAGE_VERSION_PATCH 0)
include(CPack)
