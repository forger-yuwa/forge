#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

//#include <cuda_runtime.h>
//#include <cublas_v2.h>

#include <thrust/extrema.h>
#include <thrust/device_ptr.h>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
//#include "cuda_forge/derived_atomic_functions.cuh"
//#include "cuda_forge/atomic.cuh"

#include "input/solverConfig.hpp"
#include "variables.hpp"

//#include "atomic.cuh"

__global__ void setCFL_pln_d
( 
 flow_float dt,
 //flow_float dt_pseudo,
 flow_float visc,
 flow_float* vis_turb  ,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 //flow_float* cfl_pseudo ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 //plane variables
 flow_float* cfl_pln
 //flow_float* cfl_pseudo_pln
);

__global__ void setCFL_cell_d
( 
 int dtControl, flow_float cfl_target, flow_float cfl_pseudo_target,
 int dualTime, int unsteady,
 flow_float dt,
 //flow_float dt_pseudo,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes,geom_int* cell_planes_index, geom_int* cell_planes,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 flow_float* dt_local  ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 flow_float* cfl_pln
 //flow_float* cfl_pseudo_pln 
);

// dt 適応と console 出力は独立に制御する (両者を束ねない)。
//   adaptDt : dtControl==1 のとき cfg.dt を max cfl から適応する (時間精度に効く)。explicit unsteady では毎ステップ要。
//             定常 implicit では dt_local=cfl_pseudo·dx/λ で cfg.dt が打ち消され不影響なので間引いてよい。
//   printCfl: host 読みした max cfl を cfg.monitorCflMax に格納する (console モニタ行が表示。unsteady==1 の monitor step のみ true にする)。
// host 読み出し (thrust::max_element の D2H+同期) は adaptDt||printCfl のときだけ発生する。両 false なら同期ゼロ。
// 既定は両 true (従来挙動)。
void setDT_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var ,
                     bool adaptDt = true , bool printCfl = true);
