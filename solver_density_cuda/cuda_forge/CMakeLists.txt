project(cuda_forge LANGUAGES CUDA)

# 共通 species データ (data/species/forge_species_v1.yaml) をヘッダへ埋め込む。
# input/speciesDB.cpp が起動時にパースする (実行時にファイルを読まない; plans/active/thermophysics-solver-owned-species-db.md #4)。
set(FORGE_SPECIES_DATA_YAML ${CMAKE_CURRENT_SOURCE_DIR}/../data/species/forge_species_v1.yaml)
# 輸送データ (CEA trans.inp 由来; tools/cea_trans_to_forge_transport.py) も同じヘッダへ埋め込む (plan #5t2)。
set(FORGE_TRANSPORT_DATA_YAML ${CMAKE_CURRENT_SOURCE_DIR}/../data/species/forge_transport_v1.yaml)
set(FORGE_SPECIES_DATA_GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
add_custom_command(
    OUTPUT ${FORGE_SPECIES_DATA_GEN_DIR}/forge_species_data.hpp
    COMMAND ${CMAKE_COMMAND} -DIN=${FORGE_SPECIES_DATA_YAML} -DOUT=${FORGE_SPECIES_DATA_GEN_DIR}/forge_species_data.hpp
            -DTRANS=${FORGE_TRANSPORT_DATA_YAML}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/embed_species_data.cmake
    DEPENDS ${FORGE_SPECIES_DATA_YAML} ${FORGE_TRANSPORT_DATA_YAML} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/embed_species_data.cmake
    COMMENT "Embedding species data ${FORGE_SPECIES_DATA_YAML}"
    VERBATIM)

add_library(cuda_forge STATIC
    ${FORGE_SPECIES_DATA_GEN_DIR}/forge_species_data.hpp
    cudaConfig.cuh
    cudaWrapper.cu
    cudaWrapper.cuh
    calcStructualVariables_d.cu
    calcStructualVariables_d.cuh
    dependentVariables_d.cu
    dependentVariables_d.cuh
    update_d.cu
    update_d.cuh
    qAccumulator_d.cu
    qAccumulator_d.cuh
    qAccumulator.hpp
    boundaryCond_d.cu
    boundaryCond_d.cuh
    ransBoundary_d.cu
    ransBoundary_d.cuh
    ransWallFunction_d.cu
    ransWallFunction_d.cuh
    wmlesWallModel_d.cu
    wmlesWallModel_d.cuh
    nodeWallDirichlet_d.cu
    weakIsothermalWall_d.cu

    nodeWallDirichlet_d.cuh
    bodyForce_d.cu
    bodyForce_d.cuh
    wallLaw_d.cuh
    calcGradient_d.cu
    calcGradient_d.cuh
    limiter_d.cu
    limiter_d.cuh
    convection/convectiveFlux_d.cu
    convection/convectiveFlux_d.cuh
    scalarTransport_d.cu
    scalarTransport_d.cuh
    ransTransport_d.cu
    ransTransport_d.cuh
    speciesTransport_d.cu
    speciesTransport_d.cuh
    passiveTransport_d.cuh
    passiveKernels_d.cuh
    passiveLimiter_d.cuh
    limiterFunctions_d.cuh
    chemistry_d.cu
    chemistry_d.cuh
    chemistry_mech_io.hpp
    chemistrySource_d.cuh
    condensationTransport_d.cu
    condensationTransport_d.cuh
    condensationUpdateLimiter_d.cuh
    condensationRealizability_d.cuh
    condensationProperties_d.cuh
    condensationEOS_d.cuh
    condensationSource_d.cu
    condensationSource_d.cuh
    ransSource_d.cu
    ransSource_d.cuh
    axisymmetricSource_d.cu
    axisymmetricSource_d.cuh
    ducrosSensor_d.cu
    ducrosSensor_d.cuh
    turbulent_viscosity_d.cu
    turbulent_viscosity_d.cuh
 
    #convection/AUSM_d.cu
    #convection/AUSM_d.cuh
    viscousFlux_d.cu
    viscousFlux_d.cuh
    # V0 の評価の経路 (plans/active/architecture-float-state-double-geometry.md §4.2c、既定 off。ホスト側の段取りと h5 の書き出し)
    geomAbDiag.cpp
    geomAbDiag.hpp

    implicitCorrection_d.cu
    implicitCorrection_d.cuh
    timeIntegration_d.cu
    timeIntegration_d.cuh
    setDT_d.cu
    setDT_d.cuh
    periodicNode_d.cu
    periodicNode_d.cuh
    #derived_atomic_functions.cu
    #derived_atomic_functions.cuh
    fluct_variables_d.cu
    fluct_variables_d.cuh
 
    gasProperties_d.cu
    gasProperties_d.cuh
    thermo_d.cu
    thermo_d.cuh
    ../input/speciesDB.cpp
    ../input/speciesDB.hpp
    ../input/speciesTransportDB.cpp
    ../input/speciesTransportDB.hpp
    transition_d.cu
    transition_d.cuh
    tracerTransport_d.cu
    tracerTransport_d.cuh
    residualMonitor_d.cu
    residualMonitor_d.cuh
 

)


# speciesDB.cpp (C++) は flowFormat.hpp 経由で CUDA のヘッダ (vector_types.h) を読むので、forge・変換器と同じく CUDA の include を足す
target_include_directories(cuda_forge PRIVATE .. ${CMAKE_CURRENT_SOURCE_DIR} ${FORGE_SPECIES_DATA_GEN_DIR} ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# Inherit CUDA_ARCHITECTURES from the top-level CMakeLists.txt.
if(CMAKE_CUDA_ARCHITECTURES)
    set_property(TARGET cuda_forge PROPERTY CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
endif()

target_compile_options(cuda_forge PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3 -lineinfo>)
