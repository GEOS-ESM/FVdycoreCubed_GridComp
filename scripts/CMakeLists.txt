if(HYDROSTATIC)
   set(CFG_NONHYDROSTATIC FALSE)
else()
   set(CFG_NONHYDROSTATIC TRUE)
endif()

# Did we build for AMD Rome hardware (aka EPYC)?
cmake_host_system_information(RESULT proc_decription QUERY PROCESSOR_DESCRIPTION)
if (${proc_decription} MATCHES "EPYC")
  set(CFG_BUILT_ON_ROME TRUE)
else ()
  set(CFG_BUILT_ON_ROME FALSE)
endif ()

configure_file(fv3_setup fv3_setup @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/fv3_setup DESTINATION bin)

configure_file(fv3.j fv3.j COPYONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/fv3.j DESTINATION bin)

configure_file(.FV3_VERSION .FV3_VERSION @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/.FV3_VERSION DESTINATION etc)
