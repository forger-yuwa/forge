#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

struct deltas {
    flow_float del, delp, delm;
    flow_float rik;
};

__device__ flow_float venkata_limiter(flow_float delp, flow_float delm , flow_float volume);
//__device__ flow_float barth_Jespersen_limiter(flow_float delp, flow_float delm , flow_float volume);
__device__ flow_float barth_Jespersen_limiter(deltas delta);
__device__ flow_float nishikawa_r1_limiter(deltas delta);

__device__ deltas calcDeltaIJ(geom_float pcx , geom_float pcy, geom_float pcz, 
                              geom_float pc2x , geom_float pc2y, geom_float pc2z, 
                              flow_float dudx, flow_float dudy,flow_float dudz,
                              flow_float delu_max, flow_float delu_min);

__device__ deltas calcDeltaIJ_dash(geom_float pcx , geom_float pcy, geom_float pcz, 
                              geom_float pc2x , geom_float pc2y, geom_float pc2z, 
                              flow_float dudx, flow_float dudy,flow_float dudz,
                              flow_float delu_max, flow_float delu_min);


void limiter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 診断専用 (plan limiter-inlet-column-oscillation §5.1 #5b'、既定 off): ψ の保存と入口集合への差し替え
void limiterPsiSave(solverConfig& cfg, mesh& msh, variables& var, flow_float xmax);
void limiterPsiOverride(bool on);
const flow_float* limiterPsiSaved(int k);

// run 終了時に 1 度呼び、診断の末尾取りこぼしを回収する (plan §4.35)
void limiterDiag_finalize(solverConfig& cfg);
