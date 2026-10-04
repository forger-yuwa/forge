#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"
#include <string>
#include <vector>

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

// node 入口ピン (scalarDirichletPin==1) ノードの res_roY{s}/src_jac_Y{s} を 0 化する (化学ソース集計の後に呼ぶ)。
// cell / 単成分では no-op。入口 Dirichlet の境界ノード値は speciesBoundary_d_wrapper (node 分岐) が毎 step ピンする。
void speciesPinResidual_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 化学種移流残差を組み立てる (res_roY{s} / transport_diag をゼロ初期化してから集計)。
// 粘性 (viscMethod!=0) かつ nSpecies>=2 のとき M4 の Fick 拡散 + ΣJ=0 補正 + エンタルピー拡散も加える。
void speciesTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 試験用 (FORGE_TRANSPORT_PROBE): Fick 拡散と同じ組成 (凝縮 carrier は気相組成) と thermo_Dmix_species_f で、セルごとの
// 分子拡散係数 D_s を D_d[s*nCells_all + ic] に書く (化学種なしは false)。
bool speciesDmixProbe_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, float* D_d);
// lump を含む化学種拡散の縮約が有効か (起動時に speciesInit_d が決める; plan thermophysics-solver-owned-species-db #7b)。
bool speciesLumpDiffusionActive();

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
// 二相拡散 (plan condensation-two-phase-transport §4.2, #4e; condTwoPhaseDiffusionActive の構成だけ):
//   twoPhaseDiffusion_d_wrapper: 面流束を 1 回組み、化学種・液 ρg・Q・エネルギーの残差と点対角に足す (凝縮モーメントの残差ゼロ化の後に呼ぶ)。
//   speciesRenormalizeTwoPhase_d_wrapper: speciesRenormalize と同じ係数 ρ/ΣρY を液 ρg・Q にも掛ける (非分割更新の commit の後に呼ぶ)。
void twoPhaseDiffusion_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void speciesRenormalizeTwoPhase_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// 収束受入の独立残差監査 (#4f (4)): 呼び出し側が assembleResidual を回した直後に呼ぶ。二相系の各成分を double で組み直して判定し
// [twophase-audit] 行に出す (final=false で r0 を記録、true で VERDICT)。res_rog/Q はソースだけの値に置き換わる (次の組立てで戻る)。
void twoPhaseAudit_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int iStep, bool final);
// 診断 D1 (FORGE_DIAG_TP_FACES; plan condensation-two-phase-default §5.1 #4): 組立の前処理・後処理を 1 回ずつ通した状態から、全通常面の
// OFF (species_diffusion_d の診断用の写し) と ON (本番の tp_build_face_in + tp_face_flux<float> と double 参照) の面作用素を評価し、
// OFF の写しを本番カーネルの拡散寄与 (退避用の 0 初期化配列に流したもの) と節点で照合する。状態・本番の残差配列は書かない。
// face: 面ごとのブロック [nFaces][width]、node: 節点ごとのブロック [nNodes][width] (名前は h5 のパス)。kind 0 float / 1 double / 2 int。
struct TpFaceDiagData {
    struct Block { std::string name; int width = 1; int kind = 0; std::vector<float> f; std::vector<double> d; std::vector<int> i; };
    long nFaces = 0, nNodes = 0; int nSpecies = 0, iw = -1;
    std::vector<Block> face, node;
    std::vector<std::string> summary;   // [tp-faces] 要約行 (ログと h5 属性)
};
bool twoPhaseFaceDiag_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, TpFaceDiagData& out, std::string& why);
// 再正規化の受入ゲート (#1b-pre): 更新ごとの成分別相対補正 C_q,n と局所係数偏差の max を、区間 (final=false; 前回ログからの全更新) と
// 末尾窓 (final=true; 実更新数 N の最後の ceil(0.1N) 更新、κ = 2 n_s ε₃₂ で VERDICT) で [renorm-gate] 行に出す。TP carrier 凝縮のみ。
void renormGateLog(const solverConfig& cfg, int iStep, bool final);
// #1b-r2 (condTwoPhaseDiag 3): 組立の直後に double の組立 B を作る (状態・組立 A は不変)。twoPhaseDiagBGet は device の B 残差 [c*NQ+q]・尺度 A・
// 液の内訳 [c*3+{移流,拡散,ソース}] (NQ = 化学種数 + 5; 成分 q: 化学種, 蒸気 n, 液 n+1, Q2 n+2, Q1 n+3, Q0 n+4) を返す。twoPhaseAuditR0: 開始時監査の max|r|。
void twoPhaseDiagB_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
bool twoPhaseDiagBGet(const double** r, const double** A, const double** brk, int* NQ);
const std::vector<double>& twoPhaseAuditR0();

// RK ステップ/ステージ始点の保存 (roY{s}N / roY{s}M)。NS の updateVariablesOuter/Inner に対応。
void speciesUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
void speciesUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// ===== dual-time (plan species-passive-scalar-unification §4.4; chem e296f0d0 の移植) =====
// 物理時間レベルの初期化 P = PP = 現在値 (起動時・旧形式 restart)。
void speciesInitDualTimeLevels_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// 物理 step 冒頭のレベルシフト roY_PP ← roY_P ← roY (流れの shiftDualTimeLevels と同時)。
void speciesShiftDualTimeLevels_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// 化学種残差に BDF 項 −(V/Δt)(a ρY − b ρY^P + c ρY^PP) を加え、輸送対角に V a/Δt を足す。
// 周期 gather の後に合併体積で一度だけ呼ぶ (係数 a,b,c は呼び出し側 = 全系共有の履歴契約から決める)。
void speciesAddUnsteadyTimeTerm_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                                          flow_float a, flow_float b, flow_float c);
