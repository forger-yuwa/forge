#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// 受動トレーサ (排気率 ξ; physProp.tracer: exhaust)。保存量 roXi を汎用スカラ輸送コア
// scalarTransport_d (ScalarTransportDesc) で移流し、入口 (inlet_*) では bcond floats の Xi を Dirichlet、
// 他境界は Neumann (zero-gradient)。凝縮モーメント (condensationTransport_d) と同形の受動スカラ。
// 拡散は 0 (汎用拡散は μ 係数で Sc を持たず、化学種 Fick 拡散は多成分専用カーネルのため; 混合層の
// ξ は移流のみ。必要なら plan §4.5 の混合平均 Sc を別途実装)。
// 無効 (var.tracerRegistered == 0) のときは全 wrapper が no-op で従来経路を保つ (ゼロコスト)。
// plans/active/thermophysics-cea-mole-fraction-species.md §4.5。

// 原始量 Xi = roXi/ρ を全セル (ghost 含む) について更新する。roXi は [0, ρ] にクランプ (Xi ∈ [0,1])。
void tracerPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// ghost 充填 (inlet_* は bvar Xi の Dirichlet + node 境界ノードピン、他は Neumann)。
void tracerBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var);
void applyTracerBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// node 入口ピン (scalarDirichletPin==1) ノードの res_roXi / src_jac_Xi を 0 化する (cell では no-op)。
void tracerPinResidual_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 移流残差 (res_roXi / transport_diag_Xi / src_jac_Xi をゼロ初期化してから集計)。
void tracerTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 時間積分 (RK / point-implicit)。
void tracerTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// RK ステップ/ステージ始点の保存 (roXiN / roXiM)。
void tracerUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void tracerUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
