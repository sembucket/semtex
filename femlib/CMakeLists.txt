# -- NB:
#
#    The only semtex routine to directly require bison/yacc   is initial.y.
#    The only semtex routine to directly require MPI (if set) is message.c.
#
#    The Fortran source files below are written in Fortran77.  Those
#    with extension .F will be run through the C preprocessor prior to
#    compilation; those with extension .f will be compiled without
#    preprocessing.

BISON_TARGET (Parser initial.y ${CMAKE_CURRENT_BINARY_DIR}/initial.c)

set (fem_lib_src
  ${BISON_Parser_OUTPUTS}
  polyops.c  operators.c   polylib.c  filter.c
  fourier.c  mapping.c     family.c   message.c
  temfftd.F  matops.F      canfft.f   netlib.f  
)

add_library (fem STATIC ${fem_lib_src})

if (MPI_FOUND)
  add_library (fem_mp STATIC ${fem_lib_src})
  target_compile_definitions (fem_mp PRIVATE -DMPI)	
endif()
