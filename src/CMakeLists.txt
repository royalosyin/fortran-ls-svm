enable_language(Fortran)

set(sources flssvm_main.f90 flssvm_predict.f90 flssvm_train.f90 flssvm_utilities.f90)
set(Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
add_executable(flssvm ${sources})
