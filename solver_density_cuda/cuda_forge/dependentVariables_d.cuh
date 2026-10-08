#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"
#include "thermo_d.cuh"

__global__ void dependentVariables_d
(
 // gas properties
 int thermalMethod ,
 flow_float gamma , flow_float cp ,

 // thermally-perfect (thermalMethod==2) 用化学種データ
 const SpeciesThermo* sp , int nSpecies , flow_float** roY ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* k   ,
 flow_float* omega,
 flow_float* T   ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 flow_float* gam_array ,
 flow_float* cp_array

);

void dependentVariables_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// TP / 凝縮の温度反転の下限 DEPVAR_TMIN [K] (床事象のカウンタが記録・判定に使う; floorEvents_d.cu)
double dependentVariablesTminTP();
