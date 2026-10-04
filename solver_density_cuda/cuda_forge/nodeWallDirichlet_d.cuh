#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// =============================================================================
// node (median-dual) の壁強境界 (Dirichlet) 一式 — 「状態層」(methods/boundary.md
// 「node 等温壁の壁ノード温度ピン」/ discretization §7.2.1)。
//
// 3 層構造の分担:
//   [状態層]   本ファイル: 壁ノードの状態ピン (u=0 / T=Tw) + 壁ノード残差射影。
//              壁関数・壁モデルの有無に依存しない (wallTreatmentSST / wallModelLES を参照しない)。
//   [流束層]   viscousFlux_d: W-I 双対面流束の評価切替 (Tau_Wall / Qw_Wall マーカ,
//              -1=解像 / それ以外=モデル置換)。
//   [モデル層] ransWallFunction_d (SST) / wmlesWallModel_d (WMLES): マーカと bvar を書く。
//
// 陰解法整合: 状態ピンと対で block-DPLUR の対角ブロック行 decouple が必須
// (timeIntegration_d.cu, wall_flag → 運動量 3 行 / iso_wall_flag → エネルギー行)。
// decouple 無しの状態ピンは implicit で数 step 発散する (2026-07-20 実測)。
// =============================================================================

// 壁ノード状態を u=0 に確定する (IC 確定後に一度 + 各 step)。roe から運動エネルギーを除去。
void enforceWallNoSlip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 壁ノードの Dirichlet 行残差射影 (assembleResidual 末尾・軸射影の後)。
// 運動量 3 行は常時。ω (SST) / roe (等温ピン・WMLES 等温) は wrapper が cfg から判断して渡す。
// 旧名 zeroWallMomentumResidual_d_wrapper (運動量以外もゼロ化するため改名, 2026-07-20)。
void zeroWallDirichletResiduals_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 素の node 等温壁 (非 WMLES) の壁ノード温度ピンが有効か (node && nodeWallDirichlet && wall_isothermal あり)
bool nodeIsothermalPinActive(const solverConfig& cfg, const mesh& msh);

// 1 bcond 分の壁ノード温度状態ピン (T/roe/P/sonic を bvar Ts に整合)。WMLES 等温壁も共用。
void pinWallNodeTemperature_bcond(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , variables& var);

// 素の node 等温壁 (非 WMLES) 全 bcond の壁ノード温度ピン (applyBconds 位相で毎 step)
void applyNodeIsothermalWallPin(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 素の node 等温壁の壁ノード res_roe ゼロ化 (zeroWallDirichletResiduals から呼ばれる)
void zeroNodeIsothermalEnergyResidual(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
// 残差内訳の診断 (interfaceDiag!=0 のみ。読むだけ): 等温壁ノードの res_roe を bvar `name` に写す。
void captureNodeIsothermalEnergyResidual(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var ,
                                         const std::string& name , const std::string& srcField = "res_roe");

// SST 壁関数の熱的閉包 (methods/turbulence §6.5(f), 出力層): 断熱壁の壁面出力温度 (bvar Ts) を
// Crocco 型 T_aw に置換する。状態・保存量は不変 (applyRansScalarBoundaries の後に呼ぶ —
// Taw_diag が当該 step 値である必要)。
bool sstThermalWallFunctionActive(const solverConfig& cfg);
void applySstThermalWallFunction(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
