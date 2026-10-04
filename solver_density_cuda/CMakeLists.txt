cmake_minimum_required(VERSION 3.22)
##project(forge LANGUAGES CUDA)
# Enable C language because FindHDF5 runs C tests
project(forge LANGUAGES C CXX CUDA)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES "75;86;89;90")
endif()
##
##set(CMAKE_CXX_COMPILER /usr/bin/g++ )
##add_definitions( -std=c++11 -Wall )
###add_definitions( -std=c++2a )
####set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

#cuda
set(dyn_libs_cuda cuda rt dl cudart_static cublas)

# HighFive はヘッダオンリー。まずはローカル(同梱)を探し、無ければシステムを探す

find_path(HIGHFIVE_ROOT
  NAMES highfive/highfive.hpp highfive/H5Easy.hpp HighFive/H5Easy.hpp
  PATHS
    ${CMAKE_SOURCE_DIR}/third_party/HighFive/include
    ${CMAKE_SOURCE_DIR}/../third_party/HighFive/include
    /usr/local/include
    /usr/include
  NO_CACHE)

if(NOT HIGHFIVE_ROOT)
  # Manual fallback
  if(EXISTS ${CMAKE_SOURCE_DIR}/third_party/HighFive/include/highfive/highfive.hpp)
    set(HIGHFIVE_ROOT ${CMAKE_SOURCE_DIR}/third_party/HighFive/include)
  endif()
endif()

if(NOT HIGHFIVE_ROOT)
  message(FATAL_ERROR "HighFive headers not found. Did you init submodule (git submodule update --init --recursive)? Looked in ${CMAKE_SOURCE_DIR}/third_party/HighFive/include")
endif()
include_directories(${HIGHFIVE_ROOT})

# ---- Core dependencies (placed early so include dirs reach CUDA compilation) ----
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
find_package(YAML-CPP REQUIRED)
find_package(METIS QUIET CONFIG)
if(NOT METIS_FOUND)
  find_path(METIS_INCLUDE_DIR NAMES metis.h PATHS /usr/include /usr/local/include)
  find_library(METIS_LIBRARY NAMES metis PATHS /usr/lib /usr/local/lib)
  if(METIS_INCLUDE_DIR AND METIS_LIBRARY)
    add_library(METIS::metis UNKNOWN IMPORTED)
    set_target_properties(METIS::metis PROPERTIES IMPORTED_LOCATION ${METIS_LIBRARY} INTERFACE_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIR})
    set(METIS_FOUND TRUE)
    message(STATUS "Found METIS (manual): ${METIS_LIBRARY}")
  else()
    message(FATAL_ERROR "METIS not found. Install libmetis-dev or set METIS_DIR.")
  endif()
endif()
include_directories(${HDF5_INCLUDE_DIRS})


#hdf

#kdtree c++ dir


# --------------
# *** forge ***
# --------------
file( GLOB forge_cppfiles 
    main.cpp 
    gradient.cpp 
    mesh/mesh.cpp 
    output/output.cpp 
    update.cpp 
    variables.cpp 
    boundaryCond.cpp 
    conjugateWall.cpp 
    conjugate/solidFem2d.cpp 
    convectiveFlux.cpp 
    input/solverConfig.cpp 
    input/calcWallDistance_kdtree.cpp 
    setStructualVariables.cpp 
    dependentVariables.cpp 
)

add_executable(forge ${forge_cppfiles})
target_include_directories(forge PRIVATE ${CMAKE_SOURCE_DIR} ${HIGHFIVE_ROOT} ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})


add_subdirectory(cuda_forge)
add_subdirectory(probe )


target_compile_options(forge PUBLIC
             $<$<COMPILE_LANGUAGE:CXX>:-O3 -fopenmp>)

target_link_libraries(forge PRIVATE cuda_forge probe)
target_link_libraries(forge PRIVATE ${HDF5_LIBRARIES} yaml-cpp ${dyn_libs_cuda})

# --------------------------------------
# *** CHT 固体 FE の同値試験用 CLI ***
# 参照実装 (tools/solid_fem2d.py) と突き合わせるためだけの小さな実行ファイル。
# 本体 (forge) には conjugate/solidFem2d.cpp が入っているので、ここは main だけ足す。
# --------------------------------------
add_executable(solid_fem2d_tool conjugate/solid_fem2d_tool.cpp conjugate/solidFem2d.cpp)
target_include_directories(solid_fem2d_tool PRIVATE ${CMAKE_SOURCE_DIR} ${HIGHFIVE_ROOT})
target_link_libraries(solid_fem2d_tool ${HDF5_LIBRARIES})
target_compile_options(solid_fem2d_tool PRIVATE -O3)

# ----------------------
# *** mesh conversion***
# ----------------------
file( GLOB mesh_cppfiles 
    mesh/convertGmshToForge.cpp 
    mesh/mesh.cpp 
    variables.cpp 
    boundaryCond.cpp 
    input/solverConfig.cpp 
    input/calcWallDistance_kdtree.cpp 
)

add_executable(convertGmshToForge  ${mesh_cppfiles})
target_include_directories(convertGmshToForge PRIVATE ${CMAKE_SOURCE_DIR} ${HIGHFIVE_ROOT} ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(convertGmshToForge ${HDF5_LIBRARIES} yaml-cpp)
target_link_libraries(convertGmshToForge cuda_forge )

# ---------------------
# *** mesh partition***
# ---------------------
file( GLOB mesh_part_cppfiles 
    mesh/mesh_partition.cpp 
)


add_executable(mesh_part ${mesh_part_cppfiles})
set_target_properties(mesh_part PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(mesh_part PRIVATE ${CMAKE_SOURCE_DIR} ${HIGHFIVE_ROOT})
if(TARGET METIS::metis)
  target_link_libraries(mesh_part ${HDF5_LIBRARIES} METIS::metis yaml-cpp)
else()
  target_link_libraries(mesh_part ${HDF5_LIBRARIES} metis yaml-cpp)
endif()


# (core dependency discovery moved earlier)

# ---- Optional kdtree ----
find_path(KDTREE_INCLUDE_DIR NAMES kd_tree.hpp kdtree.h
          PATHS /usr/include /usr/local/include ${CMAKE_SOURCE_DIR}/third_party/kdtree/include)
find_library(KDTREE_LIBRARY   NAMES kdtree
             PATHS /usr/lib /usr/local/lib ${CMAKE_SOURCE_DIR}/third_party/kdtree/lib)
if(KDTREE_INCLUDE_DIR AND KDTREE_LIBRARY)
  message(STATUS "Found kdtree: ${KDTREE_LIBRARY}")
  add_library(kdtree::kdtree UNKNOWN IMPORTED)
  set_target_properties(kdtree::kdtree PROPERTIES
    IMPORTED_LOCATION "${KDTREE_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${KDTREE_INCLUDE_DIR}")
  set(HAVE_KDTREE TRUE)
else()
  message(WARNING "kdtree not found. Features depending on it will be disabled.")
  set(HAVE_KDTREE FALSE)
endif()

if(HAVE_KDTREE)
  target_link_libraries(probe PUBLIC kdtree::kdtree)
  target_compile_definitions(probe PUBLIC HAVE_KDTREE=1)
  target_link_libraries(forge PRIVATE kdtree::kdtree)
  target_compile_definitions(forge PRIVATE HAVE_KDTREE=1)
  target_link_libraries(convertGmshToForge PRIVATE kdtree::kdtree)
  target_compile_definitions(convertGmshToForge PRIVATE HAVE_KDTREE=1)
endif()
