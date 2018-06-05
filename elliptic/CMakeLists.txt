# need this to ensure that the executables below are built in the correct directory
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/build")

add_executable(elliptic "${CMAKE_SOURCE_DIR}/elliptic/helmholtz.cpp" "${CMAKE_SOURCE_DIR}/elliptic/drive.cpp")

target_link_libraries(elliptic src fem vec ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES})

target_include_directories(elliptic PUBLIC veclib femlib src include ${MPI_INCLUDE_PATH})
