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
  geos_gtfv3_interface_py_h
  # mpirun -np 1
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
  ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../geos_gtfv3_interface.py)

add_executable(
  ${exec}
  input/grid_bounds.f90
  input/domain_dim.f90
  input/input_scalars.f90
  input/input_arrays.f90
  ../geos_gtfv3_interface.c
  ../geos_gtfv3_interface.f90
  main.f90)
add_dependencies(${exec} geos_gtfv3_interface_py_h)
target_link_libraries(${exec} MPI::MPI_Fortran MPI::MPI_C)
target_include_directories(${exec} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${exec} -L${CMAKE_CURRENT_BINARY_DIR} geos_gtfv3_interface_py)

set_property(
  TARGET ${exec}
  APPEND
  PROPERTY ADDITIONAL_CLEAN_FILES geos_gtfv3_interface_py.c geos_gtfv3_interface_py.h geos_gtfv3_interface_py.o libgeos_gtfv3_interface_py.so
)
