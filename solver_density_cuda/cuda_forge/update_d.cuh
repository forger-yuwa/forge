#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void updateVariablesOuter_d
( 
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

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,
 flow_float* roKN ,
 flow_float* roOmegaN ,
 
 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM ,
 double* qacc_ro , double* qacc_roUx , double* qacc_roUy , double* qacc_roUz , double* qacc_roe ,
 int* qaccAdopt
);

void updateVariablesOuter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

__global__ void updateVariablesInner_d
( 
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

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 
);

void updateVariablesInner_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

__global__ void applyScalarImplicitCorrection_d
(
 geom_int nCells_all , geom_int nCells,

 flow_float* ro,
 flow_float* roUx,
 flow_float* roUy,
 flow_float* roUz,
 flow_float* roe,

 flow_float* roN,
 flow_float* roUxN,
 flow_float* roUyN,
 flow_float* roUzN,
 flow_float* roeN,

 flow_float* dq_ro,
 flow_float* dq_roUx,
 flow_float* dq_roUy,
 flow_float* dq_roUz,
 flow_float* dq_roe,
 double* qacc_ro , double* qacc_roUx , double* qacc_roUy , double* qacc_roUz , double* qacc_roe
);

void applyScalarImplicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// block DPLUR の補正反映: Q = Q_baseline(roN) + dq_block_old を 1 回 commit する。
// 古典 DPLUR では sweep 中に Q を更新せず、全 sweep 後にこの kernel でまとめて反映する。
__global__ void applyBlockImplicitCorrection_d
(
 geom_int nCells_all , geom_int nCells,

 flow_float* ro,
 flow_float* roUx,
 flow_float* roUy,
 flow_float* roUz,
 flow_float* roe,

 flow_float* roN,
 flow_float* roUxN,
 flow_float* roUyN,
 flow_float* roUzN,
 flow_float* roeN,

 flow_float* dq_block_0,
 flow_float* dq_block_1,
 flow_float* dq_block_2,
 flow_float* dq_block_3,
 flow_float* dq_block_4,
 geom_int* axis_flag,
 flow_float updateGuardAlpha ,
 double* qacc_ro , double* qacc_roUx , double* qacc_roUy , double* qacc_roUz , double* qacc_roe
);

void applyBlockImplicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// SST (k-ω) segregated point-implicit 更新（消散項の対角を陰化）。
__global__ void applySSTPointImplicit_d
(
 geom_int nCells,
 geom_float* vol,
 flow_float* dt_local,
 flow_float implicit_relax,
 flow_float* roK,
 flow_float* roOmega,
 flow_float* roKN,
 flow_float* roOmegaN,
 flow_float* res_roK,
 flow_float* res_roOmega,
 flow_float* src_jac_k,
 flow_float* src_jac_omega,
 flow_float* transport_diag_k,
 flow_float* transport_diag_omega,
 flow_float unsteady_diag
);

void applySSTPointImplicit_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// dual-time: 物理時間 BDF 項を残差に加える / in-place commit / 時間レベルシフト。
void addUnsteadyTimeTerm_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var,
                                   flow_float a, flow_float b, flow_float c, int include_scalar);
void applyBlockImplicitCorrectionInPlace_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
void shiftDualTimeLevels_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
