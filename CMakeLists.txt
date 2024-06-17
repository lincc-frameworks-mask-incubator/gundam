cmake_minimum_required(VERSION 3.16)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
find_package(OpenMP)

set(module_name cflibfor)
set(PYF_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/gundam/${module_name}.pyf")
set(FORTRAN_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/gundam/${module_name}.f90")

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

# Generate the C and Fortran wrapper files
add_custom_command(
  OUTPUT "${module_name}module.c" "${module_name}-f2pywrappers2.f90"
  COMMAND "${PYTHON_EXECUTABLE}" -m numpy.f2py "${PYF_FILE}"
)

# Create the Python extension module
python_add_library(${module_name} MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/${module_name}module.c"
  "${CMAKE_CURRENT_BINARY_DIR}/${module_name}-f2pywrappers2.f90"
  "${FORTRAN_SRC}" WITH_SOABI
)

# Link the Fortran object library
target_link_libraries(${module_name} PRIVATE fortranobject)

# Link OpenMP if available
if (OpenMP_Fortran_FOUND)
  target_link_libraries(${module_name} PRIVATE OpenMP::OpenMP_Fortran)
endif()

# Install the module
install(TARGETS ${module_name} DESTINATION ${SKBUILD_PROJECT_NAME})
