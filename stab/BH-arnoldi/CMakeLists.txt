set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/build")

add_executable(arnoldi-DB "${CMAKE_SOURCE_DIR}/stab/BH-arnoldi/iohb.c" "${CMAKE_SOURCE_DIR}/stab/BH-arnoldi/dspmvc.f" "${CMAKE_SOURCE_DIR}/stab/BH-arnoldi/arnoldi.cpp")

if(WITH_ARPACK)
  target_compile_definitions(arnoldi-DB PRIVATE ARPACK=ON)
endif(WITH_ARPACK)

target_link_libraries(arnoldi-DB stab fem vec ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${ARPACK_LIBRARY})

target_include_directories(arnoldi-DB PUBLIC veclib femlib stab include "${CMAKE_SOURCE_DIR}/stab/BH-arnoldi" ${MPI_INCLUDE_PATH})
