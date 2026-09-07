project(cuda_forge LANGUAGES CUDA)

add_library(cuda_forge STATIC
    cudaConfig.cuh
    cudaWrapper.cu
    cudaWrapper.cuh
    calcStructualVariables_d.cu
    calcStructualVariables_d.cuh
    dependentVariables_d.cu
    dependentVariables_d.cuh
    update_d.cu
    update_d.cuh
    boundaryCond_d.cu
    boundaryCond_d.cuh
    ransBoundary_d.cu
    ransBoundary_d.cuh
    ransWallFunction_d.cu
    ransWallFunction_d.cuh
    wmlesWallModel_d.cu
    wmlesWallModel_d.cuh
    nodeWallDirichlet_d.cu
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
    chemistry_d.cu
    chemistry_d.cuh
    chemistry_mech_io.hpp
    chemistrySource_d.cuh
    condensationTransport_d.cu
    condensationTransport_d.cuh
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
    residualMonitor_d.cu
    residualMonitor_d.cuh
 

)


target_include_directories(cuda_forge PRIVATE .. ${CMAKE_CURRENT_SOURCE_DIR} )

# Inherit CUDA_ARCHITECTURES from the top-level CMakeLists.txt.
if(CMAKE_CUDA_ARCHITECTURES)
    set_property(TARGET cuda_forge PROPERTY CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
endif()

target_compile_options(cuda_forge PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3 -lineinfo>)
