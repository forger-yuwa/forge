#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

void gasProperties_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 試験用 (FORGE_TRANSPORT_PROBE; physProp.transport が無ければ false を返して何もしない)。
//   セル経路と同じ組成・評価関数で double の μ・λ・展開後 X (nCells_all×nReal) を書く。
bool gasPropertiesTransportProbe_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                                           double* mu_d, double* lam_d, double* Xreal_d);
//   与えた double の輸送種 Y (K×nSpecies) と T で同じ評価関数を呼ぶ。
bool gasPropertiesTransportProbeStates_d_wrapper(solverConfig& cfg, int K, const double* Y_d, const double* T_d,
                                                 double* mu_d, double* lam_d, double* Xreal_d);
// 試験用 (#5t2-3): 同じ Y・T を float に丸めて表引きの経路 (transport_mix_Y_tab) を呼ぶ (表が無ければ false)。
bool gasPropertiesTransportProbeStatesTab_d_wrapper(solverConfig& cfg, int K, const double* Y_d, const double* T_d,
                                                    float* mu_d, float* lam_d);
// 試験用 (#5t2-3, FORGE_TRANSPORT_TABLE_PROBE): 表の単体値 (kind 0: 実種の μ・λ、kind 1: 組の η)。
bool transportTableSingles_d_wrapper(int N, const int* kind_d, const int* idx_d, const float* T_d, float* v0_d, float* v1_d);
