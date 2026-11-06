include(add_pfunit_ctest)

set(state_source "${CMAKE_CURRENT_SOURCE_DIR}/../FV_StateMod.F90")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${state_source})
file(READ "${state_source}" source)

# Compile the private production routines unchanged, with only external services stubbed.
function(extract_block source_var first last output)
  string(FIND "${${source_var}}" "${first}" begin)
  if (begin LESS 0)
    message(FATAL_ERROR "Cannot find ${first} in ${source_var}")
  endif ()
  string(SUBSTRING "${${source_var}}" ${begin} -1 block)
  string(FIND "${block}" "${last}" end)
  if (end LESS 0)
    message(FATAL_ERROR "Cannot find ${last} after ${first}")
  endif ()
  string(LENGTH "${last}" length)
  math(EXPR length "${end} + ${length}")
  string(SUBSTRING "${block}" 0 ${length} block)
  set(${output} "${block}" PARENT_SCOPE)
endfunction()

extract_block(source "\n  type T_FVDYCORE_STATE" "\n  end type T_FVDYCORE_STATE" split_state)
set(split_routines "")
foreach(routine FV_InitSplits FV_SaveSplits adjust_fv3_splits)
  extract_block(source "\nsubroutine ${routine}(" "\nend subroutine ${routine}" block)
  string(APPEND split_routines "${block}\n")
endforeach()
configure_file(adaptive_split_restart_support.F90.in adaptive_split_restart_support.F90 @ONLY)

separate_arguments(real8_flags NATIVE_COMMAND "${FREAL8}")
foreach(precision r4 r8)
  set(test FVdycoreCubed.adaptive_split_restart.${precision})
  # Separate generated sources and modules allow both precisions to build in parallel.
  set(test_source "${CMAKE_CURRENT_BINARY_DIR}/${precision}/test_AdaptiveSplitRestart.pf")
  configure_file(test_AdaptiveSplitRestart.pf "${test_source}" COPYONLY)
  add_pfunit_ctest(${test}
    TEST_SOURCES "${test_source}"
    OTHER_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/adaptive_split_restart_support.F90"
    LABELS ESSENTIAL)
  if (precision STREQUAL r8)
    target_compile_options(${test} PRIVATE ${real8_flags})
  endif ()
  if (TARGET build-tests)
    add_dependencies(build-tests ${test})
  endif ()
endforeach()
