
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

FILE(GLOB F90SOURCES ${CMAKE_CURRENT_LIST_DIR}/*.f90)
add_library(p1galerkin ${F90SOURCES})

target_include_directories(p1galerkin INTERFACE
   $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
   $<INSTALL_INTERFACE:mod>
)

target_link_libraries(p1galerkin PRIVATE
   globals linalg vtk geometry ${LAPACK_BLAS_LIBRARIES})

include(GNUInstallDirs)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
   DESTINATION mod
   COMPONENT library
)
install(TARGETS p1galerkin
   EXPORT p1galerkin-targets
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   COMPONENT library
)
