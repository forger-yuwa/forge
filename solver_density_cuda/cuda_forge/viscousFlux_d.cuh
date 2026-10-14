#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void viscousFlux_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 flow_float mu , flow_float thermCond,

 // variables
//flow_float* convx , flow_float* convy , flow_float* convz,
// flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* ro   ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* Ps  ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* Ts,
 
 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,

 flow_float* dTdx  , flow_float* dTdy , flow_float* dTdz
);

void viscousFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns);

// V0 の評価の経路 (plans/active/architecture-float-state-double-geometry.md §4.2c) 用: 内部面の粘性流束カーネルを、
// 本番 (viscousFlux_d_wrapper) と同じ引数の組み立てで 1 回起動する。残差の書き先を res[5] = {ro, roUx, roUy, roUz, roe}
// に差し替え、ge_x/y/z が非 nullptr なら面ごとの差 e を cc1 − cc0 の引き算の代わりに読む。faceFlux (非 nullptr、
// [6*nPlanes]) には atomicAdd の直前の面の流束を書く (成分の並びは viscousFlux_d の引数の説明)。W-I 実力診断は渡さない。起動後に同期する。
void viscousFluxInternalArm_d(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var ,
                              flow_float* const res[5],
                              const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,
                              flow_float* faceFlux);
