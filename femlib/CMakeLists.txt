BISON_TARGET(Parser initial.y ${CMAKE_CURRENT_BINARY_DIR}/initial.c)

set (fem_lib_src
${BISON_Parser_OUTPUTS}
polyops.c  operators.c  matops.F  polylib.c  filter.c
temfftd.F  canfft.f  fourier.c
mapping.c  family.c  message.c
netlib.f
)
   
add_library (fem STATIC ${fem_lib_src})

target_include_directories(fem PUBLIC "${CMAKE_SOURCE_DIR}/veclib" "${CMAKE_SOURCE_DIR}/femlib" "${CMAKE_SOURCE_DIR}/include" ${MPI_INCLUDE_PATH})
