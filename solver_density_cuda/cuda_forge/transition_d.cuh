#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// Langtry–Menter 2009 γ–Re_θt 遷移モデル (methods/turbulence/theory.md §11, plan turbulence-transition-lm2009)。
// 式の順序・下限は SU2 8.x (trans_sources.hpp / trans_correlations.hpp MENTER_LANGTRY / CTransLMSolver.cpp) に合わせる。
// 保存量 roGamma=ργ, roReth=ρ·Re_θt。SST 側は gammaEff (=max(γ, γ_sep)) だけを読む。
// 全 wrapper は turbulence.transition: none (var.transitionRegistered==0) で no-op。

// 受付条件の検査 (node / 低 Re SST / 非軸対称 / DES なし / sstEnergyIncludesK 0 / 定常陰解法)。違反は例外。
void transitionValidateConfig(const solverConfig& cfg);

// 原始量 γ=ργ/ρ∈[1e-4,1], Re_θt=ρRe_θt/ρ≥20。入力に遷移量が無かった初回だけ γ=1, Re_θt=自由流相関(局所 Tu) で初期化する。
void transitionPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// node 入口 (scalarDirichletPin==1) ノードを γ=1, Re_θt=自由流相関(局所 Tu) にピンする。applyRansScalarBoundaries の後に呼ぶ。
void applyTransitionBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 移流 + 拡散 (γ: μ+μt, Re_θt: 2(μ+μt)) の残差と輸送対角。residual/diag のゼロ化込み。
void transitionTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// ソース (P_γ−E_γ, P_θt)・陰的対角・gammaEff。ransSource_d_wrapper の**前**に呼ぶ (SST が同じ反復の gammaEff を読む)。
void transitionSource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// segregated point-implicit 更新 (SST 更新の直後)。周期 node は root→member ミラー込み。
void applyTransitionPointImplicit_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
