cmake_minimum_required(VERSION 3.17)
cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0054 NEW)

project(
  fv-dynamics-driver
  Fortran C
  )

find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_package(MPI REQUIRED)

set(exec dyn-driver)

add_custom_target(
  fv_dynamics_interface_py_h
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../fv_dynamics_interface.py)

add_executable(
  ${exec}
  input/grid_bounds.f90
  input/domain_dim.f90
  input/input_scalars.f90
  input/input_arrays.f90
  ../fv_dynamics_interface.c
  ../fv_dynamics_interface.f90
  main.f90)
add_dependencies(${exec} fv_dynamics_interface_py_h)
target_link_libraries(${exec} MPI::MPI_Fortran MPI::MPI_C)
target_include_directories(${exec} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${exec} -L${CMAKE_CURRENT_BINARY_DIR} fv_dynamics_interface_py)

set_property(
  TARGET ${exec}
  APPEND
  PROPERTY ADDITIONAL_CLEAN_FILES fv_dynamics_interface_py.c fv_dynamics_interface_py.h fv_dynamics_interface_py.o libfv_dynamics_interface_py.so
)
