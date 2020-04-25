# -- Definitions used to compile dns solver.

set (dns_src
  ${CMAKE_CURRENT_LIST_DIR}/dnsanalysis.cpp
  ${CMAKE_CURRENT_LIST_DIR}/fieldforce.cpp
  ${CMAKE_CURRENT_LIST_DIR}/integrate.cpp
  ${CMAKE_CURRENT_LIST_DIR}/integrates.cpp
  ${CMAKE_CURRENT_LIST_DIR}/nonlinear.cpp
  ${CMAKE_CURRENT_LIST_DIR}/drive.cpp
)

include_directories (veclib femlib ${MPI_C_INCLUDE_PATH})

set (import_serial_libs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set (import_par_libs    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
                        ${MPI_C_LIBRARIES})

# -- Default/serial executable (no MPI).

add_executable             (dns ${dns_src})
target_include_directories (dns PUBLIC src dns)
target_link_libraries      (dns src fem vec ${import_serial_libs})

# -- Add multi-process (_mp) executable if possible.

if (MPI_FOUND)
  add_executable             (dns_mp ${dns_src})
  target_compile_definitions (dns_mp PRIVATE -DMPI)
  target_include_directories (dns_mp PUBLIC src dns)
  target_link_libraries      (dns_mp src fem_mp vec ${import_par_libs})
endif (MPI_FOUND)
