add_executable (enumerate ${CMAKE_SOURCE_DIR}/utility/enumerate.cpp)
add_executable (compare   ${CMAKE_SOURCE_DIR}/utility/compare.cpp  )
add_executable (slit      ${CMAKE_SOURCE_DIR}/utility/slit.c       )
add_executable (chop      ${CMAKE_SOURCE_DIR}/utility/chop.c       )
add_executable (preplot   ${CMAKE_SOURCE_DIR}/utility/preplot.c    )
add_executable (avgdump   ${CMAKE_SOURCE_DIR}/utility/avgdump.c    )
add_executable (moden     ${CMAKE_SOURCE_DIR}/utility/moden.c      )
add_executable (noiz      ${CMAKE_SOURCE_DIR}/utility/noiz.c       )
add_executable (rstress   ${CMAKE_SOURCE_DIR}/utility/rstress.c    )
add_executable (xplane    ${CMAKE_SOURCE_DIR}/utility/xplane.c     )
add_executable (sem2tec   ${CMAKE_SOURCE_DIR}/utility/sem2tec.c    )
add_executable (sem2vtk   ${CMAKE_SOURCE_DIR}/utility/sem2vtk.c    )
add_executable (convert   ${CMAKE_SOURCE_DIR}/utility/convert.c    )
add_executable (repeatz   ${CMAKE_SOURCE_DIR}/utility/repeatz.c    )
add_executable (repeatxy  ${CMAKE_SOURCE_DIR}/utility/repeatxy.c   )
add_executable (meshpr    ${CMAKE_SOURCE_DIR}/utility/meshpr.cpp   )
add_executable (wallmesh  ${CMAKE_SOURCE_DIR}/utility/wallmesh.cpp )
add_executable (traction  ${CMAKE_SOURCE_DIR}/utility/traction.cpp )
add_executable (calc      ${CMAKE_SOURCE_DIR}/utility/calc.cpp     )
add_executable (addfield
	       ${CMAKE_SOURCE_DIR}/utility/tensorcalcs.c
	       ${CMAKE_SOURCE_DIR}/utility/addfield.cpp )
add_executable (sem2nek   ${CMAKE_SOURCE_DIR}/utility/sem2nek.cpp  )
add_executable (project   ${CMAKE_SOURCE_DIR}/utility/project.cpp  )
add_executable (transform ${CMAKE_SOURCE_DIR}/utility/transform.cpp)
add_executable (phase     ${CMAKE_SOURCE_DIR}/utility/phase.cpp    )
add_executable (lowpass   ${CMAKE_SOURCE_DIR}/utility/lowpass.cpp  )
add_executable (interp    ${CMAKE_SOURCE_DIR}/utility/interp.cpp   )
add_executable (modep     ${CMAKE_SOURCE_DIR}/utility/modep.cpp    )
add_executable (integral  ${CMAKE_SOURCE_DIR}/utility/integral.cpp )
add_executable (massmat   ${CMAKE_SOURCE_DIR}/utility/massmat.cpp  )
add_executable (probe     ${CMAKE_SOURCE_DIR}/utility/probe.cpp    )
add_executable (eneq      ${CMAKE_SOURCE_DIR}/utility/eneq.cpp     )
add_executable (stressdiv ${CMAKE_SOURCE_DIR}/utility/stressdiv.cpp)
add_executable (rectmesh  ${CMAKE_SOURCE_DIR}/utility/rectmesh.cpp )
add_executable (mapmesh   ${CMAKE_SOURCE_DIR}/utility/mapmesh.cpp  )
add_executable (repmesh   ${CMAKE_SOURCE_DIR}/utility/repmesh.cpp  )
add_executable (nonlin    ${CMAKE_SOURCE_DIR}/utility/nonlin.cpp   )
add_executable (meshplot  ${CMAKE_SOURCE_DIR}/utility/meshplot.cpp )

add_custom_command (TARGET probe POST_BUILD
			  COMMAND bash -c "ln -sf probe probeline"
			  COMMAND bash -c "ln -sf probe probeplane" )

target_compile_definitions (preplot PRIVATE IRIS=ON)
target_compile_options     (preplot PRIVATE "-g")

target_link_libraries (enumerate  src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (compare    src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (slit       src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (chop       src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (preplot    src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (avgdump        fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (moden          fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (noiz           fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (rstress        fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (xplane         fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (sem2tec        fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (sem2vtk        fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (convert        fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (repeatz        fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (repeatxy       fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (meshpr     src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (wallmesh   src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (traction   src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (calc           fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (addfield   src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (sem2nek    src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (project        fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (transform  src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (phase      src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (lowpass    src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (interp     src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (modep      src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (integral   src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (massmat    src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (probe      src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (eneq       src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (stressdiv  src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (rectmesh       fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (mapmesh    src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (repmesh    src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (nonlin     src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries (meshplot   src fem vec
		      ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})		      

target_include_directories (enumerate PUBLIC veclib femlib src)
target_include_directories (compare   PUBLIC veclib femlib src)
target_include_directories (slit      PUBLIC veclib femlib src)
target_include_directories (chop      PUBLIC veclib femlib src)
target_include_directories (avgdump   PUBLIC veclib femlib)
target_include_directories (moden     PUBLIC veclib femlib)
target_include_directories (noiz      PUBLIC veclib femlib)
target_include_directories (rstress   PUBLIC veclib femlib)
target_include_directories (xplane    PUBLIC veclib femlib)
target_include_directories (sem2tec   PUBLIC veclib femlib)
target_include_directories (sem2vtk   PUBLIC veclib femlib)
target_include_directories (convert   PUBLIC veclib femlib)
target_include_directories (repeatz   PUBLIC veclib femlib)
target_include_directories (repeatxy  PUBLIC veclib femlib)
target_include_directories (meshpr    PUBLIC veclib femlib src)
target_include_directories (wallmesh  PUBLIC veclib femlib src)
target_include_directories (traction  PUBLIC veclib femlib src)
target_include_directories (calc      PUBLIC veclib femlib)
target_include_directories (addfield  PUBLIC veclib femlib src utility)
target_include_directories (sem2nek   PUBLIC veclib femlib src)
target_include_directories (project   PUBLIC veclib femlib src)
target_include_directories (transform PUBLIC veclib femlib src)
target_include_directories (phase     PUBLIC veclib femlib src)
target_include_directories (lowpass   PUBLIC veclib femlib src)
target_include_directories (interp    PUBLIC veclib femlib src)
target_include_directories (modep     PUBLIC veclib femlib src)
target_include_directories (integral  PUBLIC veclib femlib src)
target_include_directories (massmat   PUBLIC veclib femlib src)
target_include_directories (probe     PUBLIC veclib femlib src)
target_include_directories (eneq      PUBLIC veclib femlib src)
target_include_directories (stressdiv PUBLIC veclib femlib src)
target_include_directories (rectmesh  PUBLIC veclib femlib src)
target_include_directories (mapmesh   PUBLIC veclib femlib src)
target_include_directories (repmesh   PUBLIC veclib femlib src)
target_include_directories (nonlin    PUBLIC veclib femlib src)
target_include_directories (meshplot  PUBLIC veclib femlib)

# -- Copy (non-compiled) shell scripts.

configure_file (${CMAKE_SOURCE_DIR}/utility/addquick
	       addquick.sh        COPYONLY)
configure_file (${CMAKE_SOURCE_DIR}/utility/index
	       	index.sh           COPYONLY)
configure_file (${CMAKE_SOURCE_DIR}/utility/index00
	       	index00.sh         COPYONLY)
configure_file	(${CMAKE_SOURCE_DIR}/utility/linedist
		linedist.sh        COPYONLY)
configure_file	(${CMAKE_SOURCE_DIR}/utility/pline
		pline.sh           COPYONLY)
configure_file (${CMAKE_SOURCE_DIR}/utility/rayavg
	       	rayavg.sh          COPYONLY)
configure_file (${CMAKE_SOURCE_DIR}/utility/resubmit
	       	resubmit.sh        COPYONLY)
configure_file (${CMAKE_SOURCE_DIR}/utility/run_job
	       	run_job.sh         COPYONLY)
configure_file (${CMAKE_SOURCE_DIR}/utility/run_job_example
	       	run_job_example.sh COPYONLY)
configure_file (${CMAKE_SOURCE_DIR}/utility/save
	       	save.sh            COPYONLY)
configure_file (${CMAKE_SOURCE_DIR}/utility/xfield
	       	xfield.sh          COPYONLY)

