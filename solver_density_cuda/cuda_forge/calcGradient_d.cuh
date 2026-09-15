#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void calcGradient_1_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* roe  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx , flow_float* droUxdy, flow_float* droUxdz,
 //flow_float* droUydx , flow_float* droUydy, flow_float* droUydz,
 //flow_float* droUzdx , flow_float* droUzdy, flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,

 flow_float* divU
);


// calcGradient_b_d は calcGradient_d.cu 内でのみ定義・呼出される (単一 TU)。
// 宣言はここに置かない (シグネチャは同ファイル内の定義を正とする)。

__global__ void calcGradient_2_d
( 
 // mesh structure
 geom_int nCells,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* roe  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx  , flow_float* droUxdy , flow_float* droUxdz,
 //flow_float* droUydx  , flow_float* droUydy , flow_float* droUydz,
 //flow_float* droUzdx  , flow_float* droUzdy , flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,

 flow_float* divU
);

void calcGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 原始量の AoS パック (stride 8: ro,Ux,Uy,Uz,P,T,0,0)。calcGradient_d_wrapper (gradLSQ==2) が applyBconds 後に組み、
// LSQ 勾配とリミッタの近傍 gather が 6 配列 (6 セクタ) の代わりに 1 セクタで読む。未構築なら nullptr。
const flow_float* prim_pack_device_ptr();
