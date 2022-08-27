esma_set_this ()

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

set(dependencies MAPL GFTL_SHARED::gftl-shared GMAO_hermes GEOS_Shared esmf OpenMP::OpenMP_Fortran)

esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS fvdycore
  DEPENDENCIES ${dependencies}
  DEPENDENCIES ${GFDL})

if (FV_PRECISION STREQUAL R4)
   set (GFDL fms_r4)
elseif (FV_PRECISION STREQUAL R4R8) # FV is R4 but FMS is R8
   get_target_property (extra_incs fms_r4 INCLUDE_DIRECTORIES)
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


set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

if (CRAY_POINTER)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${CRAY_POINTER})
endif()

add_definitions (-DSPMD -DMAPL_MODE)

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
