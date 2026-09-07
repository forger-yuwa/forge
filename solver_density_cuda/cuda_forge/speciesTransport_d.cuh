#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// 多成分化学種輸送 (M2)。汎用スカラ輸送コア scalarTransport_d (ScalarTransportDesc) を
// 化学種ごとに再利用し、保存質量分率 ρY_s を流れと共に移流する (M2 は移流のみ。拡散は M4)。
// 単成分 (var.nSpeciesRegistered < 2) のときは全て no-op で、M1 と同一経路を保つ。

// device の roY ポインタ配列 (flow_float*[nSpecies]) を 1 度だけ構築する。
// allocVariables 後 (c_d["roY{s}"] が確定後) に呼ぶこと。
void speciesInit_d(solverConfig& cfg, variables& var);

// dependentVariables_d_wrapper へ渡す device roY 配列ポインタ。
// 単成分時は nullptr (混合則 thermo は Y={1} に縮退)。
flow_float** species_roY_device_ptr();
// 化学反応ソース (chemistry_d.cu) 用: res_roY{s} / src_jac_Y{s} の device ポインタ配列 (単成分は nullptr)。
flow_float** species_resroY_device_ptr();
flow_float** species_srcjac_device_ptr();

// face 整合再構成 (speciesFaceReconstruction==1) 用: Y{s}/∇Y{s} の device ポインタ配列と勾配計算。
flow_float** species_Y_device_ptr();
flow_float** species_dYdx_device_ptr();
flow_float** species_dYdy_device_ptr();
flow_float** species_dYdz_device_ptr();
flow_float** species_limiterY_device_ptr();
void speciesGradient_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// S3: convectiveFlux 用 face 組成バッファ確保 + 同一面組成での species 移流。
flow_float* species_Yface_alloc(int nPlanes);
void speciesAdvectionFaceY_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 原始質量分率 Y_s = ρY_s/ρ を全セル (ghost 含む) について更新する。
void speciesPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種 ghost を Neumann (zero-gradient) で埋める。bcond ごとに呼ぶ applySpeciesBoundaries から使用。
void speciesBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var);

// 全 bcond をループして化学種 ghost を埋める (applyRansScalarBoundaries と同形)。
void applySpeciesBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種移流残差を組み立てる (res_roY{s} / transport_diag をゼロ初期化してから集計)。
// 粘性 (viscMethod!=0) かつ nSpecies>=2 のとき M4 の Fick 拡散 + ΣJ=0 補正 + エンタルピー拡散も加える。
void speciesTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// TP 多成分気体の組成-エネルギー整合補正 (speciesTimeIntegration 直後に呼ぶ)。
// roe[ic] += Σ_s (roY_s[ic] - roYN_s[ic]) * h_s(T[ic]) により組成変化に伴う roe ずれを補正し、
// speciesPrimitive (Newton 反転) が発散しないようにする。thermalMethod!=2 / 単成分では no-op。
void speciesEnergyCorrection_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種の時間積分 (scalarTimeIntegration_d を化学種ごとに呼ぶ)。
void speciesTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 緩和整合 scalar-DPLUR (speciesImplicitCoupling==1 かつ nSpecies>=2) が有効か。
bool speciesImplicitCoupled(solverConfig& cfg, variables& var);

// 緩和整合 scalar-DPLUR ソルバ。凍結残差に対し δ(ρY_s) を nStepInner 回 Jacobi sweep で緩和し
// (流れ block と同一 implicitRelax/sweep)、ρY_s=ρY_s^N+δ(ρY_s) を commit する。
// 呼び出し前に speciesUpdateOuter で ρY_s^N=ρY_s を取り、呼び出し後に renormalize/primitive すること。
void speciesImplicitDPLURSolve_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 案C: block-triangular roe↔roY coupling (speciesImplicitCoupling==2 かつ nSpecies>=2) が有効か。
bool speciesEOSCoupled(solverConfig& cfg, variables& var);

// 案C 予測+移項: species scalar-DPLUR で δ(ρY_s)* を予測 → 組成接空間 z_s=ρδY_s へ射影
// (Σz_s=0) → 解析 EOS-JVP δp_Y を評価 → flow エネルギー残差 res_roe へクロス作用 A_QY δY を
// 移項する。flow block 解 (blockDPLURSolve) の直前に呼ぶ。commit はしない (z_s を dq_roY{s}_old に残す)。
// 呼び出し前に speciesUpdateOuter で ρY_s^N=ρY_s を取ること。
void speciesEOSCrossPredictInject_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 案C 最終 commit: flow 密度更新 δρ=ρ-ρ^N も含め ρY_s = ρY_s^N + z_s + Y_s^N δρ (Σδ(ρY_s)=δρ)。
// flow commit (applyBlockImplicitCorrection) の直後に呼ぶ。呼び出し後 renormalize/primitive すること。
void speciesEOSFinalCommit_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種の実現可能性・再正規化: ρY_s>=0 にクランプし Σ_s ρY_s = ρ となるよう再スケール (ΣY_s=1)。
void speciesRenormalize_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// RK ステップ/ステージ始点の保存 (roY{s}N / roY{s}M)。NS の updateVariablesOuter/Inner に対応。
void speciesUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void speciesUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
