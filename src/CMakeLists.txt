include_directories(${LIBSRVF_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR})

add_library(srvf STATIC c_api.cc dpnbhd.cc fileio.cc functions.cc interp.cc matchselect.cc matrix.cc opencurves.cc partialmatch.cc plf.cc qmap.cc reparam.cc rotate.cc rotselect.cc srvf.cc util.cc)

if(GSL)
  target_link_libraries(srvf gsl gslcblas)
  set_target_properties(srvf PROPERTIES COMPILE_FLAGS "-DUSE_GSL")
endif()

include (GenerateExportHeader)
GENERATE_EXPORT_HEADER( srvf
  BASE_NAME srvf
  EXPORT_MACRO_NAME SRVF_EXPORT
  EXPORT_FILE_NAME srvf_export.h
  STATIC_DEFINE SRVF_BUILT_AS_STATIC
)

install(TARGETS srvf 
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static
  RUNTIME DESTINATION bin
)
