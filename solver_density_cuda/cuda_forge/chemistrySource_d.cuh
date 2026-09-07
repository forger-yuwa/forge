#pragma once
// 有限速度化学ソース項の host 側インタフェース (chemistry_d.cu)。理論・設計: methods/chemistry.md。
#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"
#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"
#include "cuda_forge/chemistry_d.cuh"

// chemistry.enabled==1 かつ TP 多成分のとき true。
bool chemistryEnabled(const solverConfig& cfg);

// 反応機構ファイルを読み ReactionTable を device へアップロードする。thermo_init_db の後に 1 度呼ぶ。
// chemistry.enabled==0 では何もしない。
void chemistry_init(solverConfig& cfg);

// host 側の反応表 (未初期化なら nullptr)。
const ReactionTable* chemistry_table_host();

// assembleResidual 内 (speciesTransport の直後) で呼ぶ: res_roY{s} += Vω_s, res_roe += VQ̇,
// src_jac_Y{s} += max(0,−∂ω_s/∂(ρY_s))。chemistry off / 単成分では no-op。
void chemistrySource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
