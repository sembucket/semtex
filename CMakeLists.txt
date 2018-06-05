cmake_minimum_required(VERSION 3.2)

project(semtex C CXX Fortran)

# pull in the pre-processor flags
option(WITH_MPI    "Build with MPI"                 OFF)
option(WITH_FLIP   "Build dog with half-period map" OFF)
option(WITH_ARPACK "Build dog with ARPACK"          OFF)

# dependencies
find_package(BISON)
find_package(BLAS)
find_package(LAPACK)
if(WITH_MPI)
  find_package(MPI)
endif(WITH_MPI)
if(EXISTS "${CMAKE_SOURCE_DIR}/stab")
  if(WITH_ARPACK)
    find_library(ARPACK_LIBRARY NAMES "arpack" HINTS ENV LD_LIBRARY_PATH PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
    if(ARPACK_LIBRARY)
      message("...found ARPACK: ${ARPACK_LIBRARY}")
    else()
      message(FATAL_ERROR "...have not found ARPACK, terminating build.")
    endif()
  endif(WITH_ARPACK)
endif()

add_definitions(-Wno-parentheses -Wno-unused-value)

# build the semtex libraris
add_subdirectory(veclib)
add_subdirectory(femlib)
add_subdirectory(src)
if(EXISTS "${CMAKE_SOURCE_DIR}/stab")
  add_subdirectory(stab)
endif()

# build the executables
add_executable(dns  "")
if(EXISTS "${CMAKE_SOURCE_DIR}/scat")
  add_executable(scat "")
endif()
if(EXISTS "${CMAKE_SOURCE_DIR}/stab")
  add_executable(dog  "${CMAKE_SOURCE_DIR}/stab/integrate.cpp" "${CMAKE_SOURCE_DIR}/stab/stabanalysis.cpp" "${CMAKE_SOURCE_DIR}/stab/drive.cpp")
endif()

# setting pre-processor compiler directives based on cmake flags (with defaults)
if(WITH_MPI)
  target_compile_definitions(fem  PRIVATE MPI=ON   )
endif(WITH_MPI)

if(EXISTS "${CMAKE_SOURCE_DIR}/stab")
  if(WITH_FLIP)
    target_compile_definitions(dog PRIVATE FLIP=ON  )
  endif(WITH_FLIP)
  if(WITH_ARPACK)
    target_compile_definitions(dog PRIVATE ARPACK=ON)
  endif(WITH_ARPACK)
endif()

# call subdirectory cmake files to build bring in additional sources files and build additional binaries
include(elliptic/CMakeLists.txt)
include(dns/CMakeLists.txt)
include(utility/CMakeLists.txt)
if(EXISTS "${CMAKE_SOURCE_DIR}/scat")
  include(scat/CMakeLists.txt)
endif()
if(EXISTS "${CMAKE_SOURCE_DIR}/stab")
  include(stab/BH-arnoldi/CMakeLists.txt)
endif()

# link libraries to the binaries
target_link_libraries(dns  src fem vec ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES})
if(EXISTS "${CMAKE_SOURCE_DIR}/scat")
  target_link_libraries(scat src fem vec ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES})
endif()
if(EXISTS "${CMAKE_SOURCE_DIR}/stab")
  target_link_libraries(dog stab fem vec ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${ARPACK_LIBRARY})
endif()

# add header files for the binaries
target_include_directories(dns  PUBLIC veclib femlib src dns  include ${MPI_INCLUDE_PATH})
if(EXISTS "${CMAKE_SOURCE_DIR}/scat")
  target_include_directories(scat PUBLIC veclib femlib src scat include ${MPI_INCLUDE_PATH})
endif()
if(EXISTS "${CMAKE_SOURCE_DIR}/stab")
  target_include_directories(dog  PUBLIC veclib femlib     stab include ${MPI_INCLUDE_PATH})
endif()

# include the tests
enable_testing()
include(test/CMakeLists.txt)
