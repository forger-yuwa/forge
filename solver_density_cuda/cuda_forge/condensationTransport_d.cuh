#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "cuda_forge/condensationProperties_d.cuh"   // CondPropOpts / condProps_make
#include "variables.hpp"

// 非平衡凝縮 (Phase 1): 凝縮種ごとの 4 モーメント (ρg,ρQ2,ρQ1,ρQ0) を、汎用スカラ輸送コア
// scalarTransport_d (ScalarTransportDesc) を再利用して受動スカラーとして移流する。
// Phase 1 は移流のみ (拡散なし・ソース=0)。核生成/成長ソースは Phase 2。methods/condensation/ 参照。
// 凝縮無効 (var.nCondSpeciesRegistered < 1) のときは全 wrapper が no-op で従来経路を保つ。

// device の rog (液相質量分率の保存量 ρg_sp) ポインタ配列 (flow_float*[nCondSpecies]) を 1 度だけ構築。
// 二相 EOS (dependentVariables) が液相質量分率 g_sp=ρg_sp/ρ を読むために使う。allocVariables 後に呼ぶ。
void condensationInit_d(solverConfig& cfg, variables& var);

// 二相 EOS へ渡す device rog 配列ポインタ。凝縮無効時は nullptr。
flow_float** cond_rog_device_ptr();

// config → kernel 値渡しの凝縮物性オプション (plans/active/condensation-air.md, condensation-kantrowitz-carrier.md)
inline CondPropOpts cond_prop_opts(const solverConfig& cfg)
{
    CondPropOpts o;
    o.latentLowT = cfg.condN2LatentLowT; o.psatLowT = cfg.condN2PsatLowT; o.liquidCp = cfg.condN2LiquidCp;
    o.gasKgasModel = (cfg.condVaporMassFraction > 0.0) ? 1 : 0;   // CPG carrier (空気) は空気 Sutherland
    o.sigmaScale = cfg.condSigmaScale; o.Yw = cfg.condVaporMassFraction;
    return o;
}
int          cond_num_species();

// 原始量 φ = ρφ/ρ を全セル (ghost 含む) について更新する。スカラ移流の上流値に使う。
void condensationPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 凝縮モーメント ghost を埋める (入口 inlet_* は dry=0 の Dirichlet、他は Neumann zero-gradient)。
void condensationBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var);

// 全 bcond をループして凝縮モーメント ghost を埋める (applySpeciesBoundaries と同形)。
void applyCondensationBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 凝縮モーメント移流残差を組み立てる (res_/transport_diag/src_jac をゼロ初期化してから集計)。
void condensationTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 凝縮モーメントの時間積分 (scalarTimeIntegration_d をモーメントごとに呼ぶ)。
void condensationTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// RK ステップ/ステージ始点の保存 (ro*_N / ro*_M)。NS の updateVariablesOuter/Inner に対応。
void condensationUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void condensationUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 相変化ソース (核生成+成長) を res_ro<φ> に加え、point-implicit 線形化を src_jac へ書く (Phase 2)。
// assembleResidual の condensationTransport の直後に呼ぶ (res/src_jac は transport がゼロ初期化済)。
// condensationSource_d.cu が実装。
void condensationSource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
