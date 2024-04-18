esma_set_this ()

option(BUILD_PYFV3_INTERFACE "Build pyFV3 interface" OFF)

set (srcs
  sw.f90 jw.f90 testcases_3_4_5_6_stand_alone.f90
  GetWeightsC2C.F90
  GetWeights.F90
  CubeHalo.F90
  Cube2LatLon.F90 LatLon2Cube.F90 AppGridCreate.F90 FV_StateMod.F90
  AdvCore_GridCompMod.F90
  DynCore_GridCompMod.F90 CreateInterpWeights_GridCompMod.F90
  StandAlone_DynAdvCore_GridCompMod.F90
  CubeToLatLonRegridder.F90
  LatLonToCubeRegridder.F90
  CubeToCubeRegridder.F90
  CubeToLatLon.F90
  CubeGridPrototype.F90
  GEOS_FV3_Utilities.F90
  fv_regrid_c2c.F90
  fv_regrid_c2c_bin.F90
  fv_regridding_utils.F90
  rs_scaleMod.F90
  )
if (BUILD_PYFV3_INTERFACE)
  list (APPEND srcs
  python/interface/cffi_lib/interface.f90
  python/interface/cffi_lib/interface.c)
endif ()

if (BUILD_PYFV3_INTERFACE)

  message(STATUS "Building pyFV3 interface")

  add_definitions(-DRUN_PYFV3)
  
  # The Python library creation requires mpiexec/mpirun to run on a
  # compute node. Probably a weird SLURM thing?
  find_package(MPI REQUIRED)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  # Set up some variables in case names change
  set(PYFV3_INTERFACE_LIBRARY     ${CMAKE_CURRENT_BINARY_DIR}/libpyFV3_interface_py.so)
  set(PYFV3_INTERFACE_HEADER_FILE ${CMAKE_CURRENT_BINARY_DIR}/pyfv3_interface_py.h)
  set(PYFV3_INTERFACE_FLAG_HEADER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/python/interface/cffi_lib/fv_flags.h)
  set(PYFV3_INTERFACE_SRCS        ${CMAKE_CURRENT_SOURCE_DIR}/python/interface/cffi_lib/interface.py)

  # This command creates the shared object library from Python
  add_custom_command(
    OUTPUT ${PYFV3_INTERFACE_LIBRARY}
    # Note below is essentially:
    #  mpirun -np 1 python file
    # but we use the CMake options as much as we can for flexibility
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PYFV3_INTERFACE_FLAG_HEADER_FILE} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${Python3_EXECUTABLE} ${PYFV3_INTERFACE_SRCS}
    BYPRODUCTS ${PYFV3_INTERFACE_HEADER_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    MAIN_DEPENDENCY ${PYFV3_INTERFACE_SRCS}
    COMMENT "Building pyFV3 interface library with Python"
    VERBATIM
    )

  # This creates a target we can use for dependencies and post build
  add_custom_target(generate_python_interface_library DEPENDS ${PYFV3_INTERFACE_LIBRARY})

  # Because of the weird hacking of INTERFACE libraries below, we cannot
  # use the "usual" CMake calls to install() the .so. I think it's because
  # INTERFACE libraries don't actually produce any artifacts as far as
  # CMake is concerned. So we add a POST_BUILD custom command to "install"
  # the library into install/lib
  add_custom_command(TARGET generate_python_interface_library
    POST_BUILD
    # We first need to make a lib dir if it doesn't exist. If not, then
    # the next command can copy the script into a *file* called lib because
    # of a race condition (if install/lib/ isn't mkdir'd first)
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
    # Now we copy the file (if different...though not sure if this is useful)
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PYFV3_INTERFACE_LIBRARY}" ${CMAKE_INSTALL_PREFIX}/lib
    )

  # We use INTERFACE libraries to create a sort of "fake" target library we can use
  # to make libFVdycoreCubed_GridComp.a depend on. It seems to work!
  add_library(pyfv3_interface_py INTERFACE)

  # The target_include_directories bits were essentially stolen from the esma_add_library
  # code...
  target_include_directories(pyfv3_interface_py INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # stubs
    # modules and copied *.h, *.inc
    $<BUILD_INTERFACE:${esma_include}/${this}>
    $<INSTALL_INTERFACE:include/${this}>
    )
  target_link_libraries(pyfv3_interface_py INTERFACE ${PYFV3_INTERFACE_LIBRARY})

  # This makes sure the library is built first
  add_dependencies(pyfv3_interface_py generate_python_interface_library)

  # This bit is to resolve an issue and Google told me to do this. I'm not
  # sure that the LIBRARY DESTINATION bit actually does anything since
  # this is using INTERFACE
  install(TARGETS pyfv3_interface_py
    EXPORT ${PROJECT_NAME}-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )

endif ()

set(dependencies MAPL GFTL_SHARED::gftl-shared GMAO_hermes GEOS_Shared esmf OpenMP::OpenMP_Fortran)

if (BUILD_PYFV3_INTERFACE)
  esma_add_library (${this}
    SRCS ${srcs}
    SUBCOMPONENTS fvdycore
    DEPENDENCIES ${dependencies}
    DEPENDENCIES ${GFDL}
    DEPENDENCIES pyfv3_interface_py) # Make the main library depend on the Python library
else ()
  esma_add_library (${this}
    SRCS ${srcs}
    SUBCOMPONENTS fvdycore
    DEPENDENCIES ${dependencies}
    DEPENDENCIES ${GFDL})
endif ()

if (FV_PRECISION STREQUAL R4)
   set (GFDL fms_r4)
elseif (FV_PRECISION STREQUAL R4R8) # FV is R4 but FMS is R8
   get_target_property (extra_incs fms_r4 INTERFACE_INCLUDE_DIRECTORIES)
   target_include_directories(${this} PRIVATE
   $<BUILD_INTERFACE:${extra_incs}>
   )
   set (GFDL fms_r8)
elseif (FV_PRECISION STREQUAL R8)
   set (GFDL fms_r8)
   string(REPLACE " " ";" tmp ${FREAL8})
   foreach(flag ${tmp})
     target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
   endforeach()
endif ()

if (FV_PRECISION MATCHES R4)
  target_compile_definitions (${this} PRIVATE -DSINGLE_FV -DOVERLOAD_R4)
elseif (FV_PRECISION MATCHES R4R8) # FV is R4 but FMS is R8
  target_compile_definitions (${this} PRIVATE -DSINGLE_FV -DOVERLOAD_R4)
endif ()

message(STATUS "Building FV as ${FV_PRECISION}")

#set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

if (CRAY_POINTER)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${CRAY_POINTER})
endif()

add_definitions (-DSPMD -DMAPL_MODE -DFVREGRID_MAPL_MODE)

foreach(flag ${tmp})
   target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()

ecbuild_add_executable (
  TARGET StandAlone_FV3_Dycore.x
  SOURCES StandAlone_FV3_Dycore.F90
  LIBS ${this} OpenMP::OpenMP_Fortran)

ecbuild_add_executable (
  TARGET rs_scale.x
  SOURCES rs_scale.F90
  LIBS ${this} OpenMP::OpenMP_Fortran)

ecbuild_add_executable (
   TARGET StandAlone_AdvCore.x
   SOURCES StandAlone_AdvCore.F90
   LIBS ${this} OpenMP::OpenMP_Fortran)

ecbuild_add_executable (
   TARGET StandAlone_DynAdvCore.x
   SOURCES StandAlone_DynAdvCore.F90
   LIBS ${this} OpenMP::OpenMP_Fortran)

ecbuild_add_executable (
   TARGET c2c.x
   SOURCES c2c.F90
   LIBS ${this} OpenMP::OpenMP_Fortran)

ecbuild_add_executable (
   TARGET interp_restarts.x
   SOURCES interp_restarts.F90
   LIBS ${this} OpenMP::OpenMP_Fortran)

ecbuild_add_executable (
   TARGET interp_restarts_bin.x
   SOURCES interp_restarts_bin.F90
   LIBS ${this} OpenMP::OpenMP_Fortran)


add_subdirectory(scripts)
