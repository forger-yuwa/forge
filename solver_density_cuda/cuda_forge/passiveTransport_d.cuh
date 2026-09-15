#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

#include <string>
#include <vector>

// 受動種 (排気トレーサ roXi・凝縮モーメント rog/roQ2/roQ1/roQ0) を化学種輸送経路に乗せる基盤
// (plans/active/species-passive-scalar-unification.md §4.1–4.3, passiveScalarScheme 1)。
// 実装は speciesTransport_d.cu (化学種カーネルを共用するため同一 TU)。受動種の順序は
//   q = 0             : トレーサ (tracer 有効時。上限 ρ [ξ<=1] を持つ唯一の受動種)
//   q = t .. t+4n−1   : 凝縮モーメント (registerCondensation の順 {g, Q2, Q1, Q0} × 種)
// 受動種は熱力学・ΣY 再正規化・ΣJ=0 補正・speciesImplicitCoupling の予測/commit には入らない。
// passiveScalarScheme 0 では passiveInit_d 以外は呼ばれない (旧経路はビット不変)。

// 受動種の device ポインタ配列を構築する (allocVariables + speciesInit_d/condensationInit_d の後に 1 度)。
void passiveInit_d(solverConfig& cfg, variables& var);
int  passive_count();          // nPassive (= tracer + 4·nCondSpecies)。0 なら受動種なし
int  passive_tracer_index();   // トレーサの q (無ければ -1)
int  passive_moment_index0();  // 最初のモーメントの q (無ければ -1)
const std::vector<std::string>& passive_cons_names();   // "roXi", "rog_0", ...
const std::vector<std::string>& passive_prim_names();   // "Xi",   "g_0",   ...

// convectiveFlux (SLAU S3) 用アクセサ: 原始量・勾配・リミッタのポインタ配列と面値バッファ。
flow_float** passive_P_device_ptr();
flow_float** passive_dPdx_device_ptr();
flow_float** passive_dPdy_device_ptr();
flow_float** passive_dPdz_device_ptr();
flow_float** passive_limiter_device_ptr();
flow_float*  passive_Pface_alloc(int nPlanes);   // scheme 1 かつ speciesFaceReconstruction>=2 のときだけ確保
flow_float*  passive_Pface_device_ptr();

// 有効判定: passiveScalarScheme==1 かつ nPassive>0。
bool passiveSchemeEnabled(const solverConfig& cfg);
// 化学種・受動種の更新に使う dt_local の倍率 (scalarCflMax; 既定 1.0)。timeIntegration 11 のみ <1 になり得る。
flow_float scalarDtScale(const solverConfig& cfg);

// 原始量 φ = ρφ/ρ (全セル, ghost 含む; クランプなし)。受動種 [q0, q0+nq)。
void passivePrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq);

// 勾配 ∇φ (Green-Gauss; node 周期半割面は除外し後段 periodicGradientGather で合併)。speciesFaceReconstruction>=1 のとき。
void passiveGradient_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// 受動種ごとの無次元化 Venkat リミッタ ψ_P (limiter_d.cu が実装; limiter_d_wrapper から呼ぶ)。
void passiveLimiter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 境界 (bcond ごと): inlet_* は Dirichlet (bvar が nullptr ならゼロ; node は境界ノードをピンし scalarDirichletPin を立てる)、
// 他は Neumann。化学種の species_dirichlet_boundary_d / species_neumann_boundary_d を受動種ポインタで呼ぶ。
void passiveBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var,
                               int q, flow_float* bvar_dirichlet);

// 移流残差 [q0, q0+nq): res_<cons>/transport_diag_<prim>/src_jac_<prim> をゼロ初期化してから、
// S3 (SLAU, speciesFaceReconstruction>=2) は convectiveFlux が書いた Pface で、それ以外は 1 次風上 (scalarTransportResidualMulti_d) で集計する。
void passiveAdvection_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq);
// トレーサの Fick 拡散 (粘性 run のみ; D = μ/(ρSc) + μ_t/(ρSc_t))。
void passiveDiffusion_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q);
// node 入口ピンノードの res/src_jac を 0 化 (全受動種; 移流・ソース集計の後に呼ぶ)。cell では no-op。
void passivePinResidual_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 陰解法 (timeIntegration 11) の受動種増分: passiveImplicitCoupling 1 のとき scalar-DPLUR sweep で dq_<cons>_old に
// δ(ρφ) を作る (緩和 passiveImplicitRelax・ピン行・周期 dq ミラー)。呼び出し前に N ベースラインを取ること。
// 戻り値: 増分が dq_<cons>_old に入っていれば true (coupling 0 では false = point-implicit を使う)。
bool passiveDPLURIncrement_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq);
// DPLUR 増分 dq_<cons>_old を commit する: ρφ = ρφ_N + δ (floor なし; 後段 passiveBounds)。
void passiveCommitIncrement_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q);
// トレーサの更新 (scheme 1): point-implicit (relax) または DPLUR 増分の commit → 上下限 0<=ρξ<=ρ (更新済み ρ) と補正収支。
void passiveTracerUpdate_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// 更新確定時の上下限と補正収支 [q0, q0+nq) (トレーサは上限 ρ、モーメントは下限 0 のみ)。
void passiveBounds_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq);
// 補正収支の要約を stdout に出す (monitorInterval ごと): 前回出力からの平均/step と全期間積算、相対量 (積算|Δ|/∫ρφ dV)。
void passiveFloorCorrLog_d_wrapper(solverConfig& cfg, int iStep);
// 単体試験・後処理用: 積算 (lo, hi, abs, total) を host へ (nPassive×4)。
std::vector<double> passiveFloorCorrTotals();

// node 周期: 受動種の状態を root→member でミラー (更新の直後に呼ぶ)。
void passiveMirrorPeriodic_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// ===== dual-time (§4.4): 受動種の物理時間レベル <cons>P/PP と BDF 項 (passiveScalarScheme 1 のみ; 旧経路 0 は不変) =====
void passiveInitDualTimeLevels_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void passiveShiftDualTimeLevels_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void passiveAddUnsteadyTimeTerm_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                                          flow_float a, flow_float b, flow_float c);
