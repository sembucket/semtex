# -- Definitions used to compile elliptic solver.

set (ell_src
  ${CMAKE_CURRENT_LIST_DIR}/helmholtz.cpp
  ${CMAKE_CURRENT_LIST_DIR}/drive.cpp
)

include_directories (veclib femlib ${MPI_C_INCLUDE_PATH})

set (import_serial_libs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set (import_par_libs    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
    			${MPI_C_LIBRARIES})

# -- Default/serial executable (no MPI).

add_executable             (elliptic ${ell_src})
target_include_directories (elliptic PUBLIC src elliptic)
target_link_libraries      (elliptic src fem vec ${import_serial_libs})

# -- Add multi-process (_mp) executable if possible.

if (MPI_FOUND)
  add_executable             (elliptic_mp ${ell_src})
  target_compile_definitions (elliptic_mp PRIVATE -DMPI)
  target_include_directories (elliptic_mp PUBLIC src)
  target_link_libraries      (elliptic_mp src fem_mp vec ${import_par_libs})
endif (MPI_FOUND)
