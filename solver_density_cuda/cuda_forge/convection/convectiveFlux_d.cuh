#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"
#include "cuda_forge/thermo_d.cuh"

__device__ flow_float calcRratio(flow_float phiC, flow_float phiD, 
                                 flow_float dphidx, flow_float dphidy, flow_float dphidz,
                                 flow_float cdx, flow_float cdy, flow_float cdz
                                 );

__device__ flow_float MINMOD(flow_float rf);

__device__ flow_float TVD_MUSCL(flow_float rf) ;

__device__ flow_float firstUp(flow_float rf) ;

__device__ flow_float MUSCL(flow_float phiC, flow_float phiD, 
                            flow_float dphidx, flow_float dphidy, flow_float dphidz,
                            flow_float cdx, flow_float cdy, flow_float cdz
                           );

__device__ flow_float betaPls_slau(flow_float M);

__device__ flow_float betaMns_slau(flow_float M);

__global__ void HLLE_d
(
 int conv_scheme, int limit_scheme,

 flow_float ga,

 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,
 geom_int nNormal_halo_Planes, geom_int* normal_halo_planes_d,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
flow_float* massflux,

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

 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   ,

 flow_float* limiter_ro,
 flow_float* limiter_Ux,
 flow_float* limiter_Uy,
 flow_float* limiter_Uz,
 flow_float* limiter_P,

 flow_float* ducros,

 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz

);


__global__ void SLAU_d
( 
 flow_float ga,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_int* normal_ghst_planes_d,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 flow_float* massflux,

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
 
 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz

);

__device__ flow_float sign_sano(flow_float x);

__device__ flow_float betaPls(flow_float M, flow_float alpha);

__device__ flow_float betaMns(flow_float M, flow_float alpha);

__device__ flow_float fPls(flow_float M);

__device__ flow_float fMns(flow_float M);

__global__ void AUSMp_UP_d
( 
 flow_float ga,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

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
 
 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,
 flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz

);

void convectiveFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns);

// 局所帳簿ダンプ (FORGE_DUMP_LEDGER、既定 off・出力専用。plan tooling-nozzle-sern-3d §5.1 R5h)
void ledgerBeginAssemble(mesh& msh);
void ledgerCapture(mesh& msh, variables& var, const char* tag, bool residual);
void ledgerFlushFaces();

// 厚さ 0 の板の自由端の近傍の速度再構成の重み w [nCells] を device へ上げる (space.zeroThicknessEdgeVelocity、起動時に 1 回)。
// 空なら何もしない (無効 = 旧経路とビット同一)。plan convection-zero-thickness-edge-reconstruction §4
void zteSetVelocityWeightDevice(const std::vector<flow_float>& w);

#include "cuda_forge/convection/farfieldFace.hpp"   // 遠方境界 farfield の面の値配列
