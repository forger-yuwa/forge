#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include <vector>
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "cuda_forge/condensationProperties_d.cuh"   // CondPropOpts / condProps_make
#include "cuda_forge/condensationTables_d.cuh"       // CondTablesF (float 経路の物性表)
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
const CondTablesF& cond_tables_device();   // condensationInit_d が構築 (valid=0 なら float 経路は使わない)

// H2O の潜熱の気液ペア (plans/active/thermophysics-solver-owned-species-db.md §4.8, #10)。凝縮 ON・condModel 1 のときだけ
// 種 DB の解決結果 (speciesDB_current().condensed: ペアの気相種と共通データの H2O(L)) と datum (thermoHrefTemp, TP のみ) から
// 1 回だけ作り、host と device に 1 つずつ置いて参照を返す (キャッシュ)。それ以外は空の参照 (N2 と凝縮 OFF は使わない)。
CondLatentRef cond_latent_pair_for(const solverConfig& cfg);

// config → kernel 値渡しの凝縮物性オプション (plans/accepted/condensation-air.md, condensation-kantrowitz-carrier.md)
inline CondPropOpts cond_prop_opts(const solverConfig& cfg)
{
    CondPropOpts o;
    o.latentLowT = cfg.condN2LatentLowT; o.psatLowT = cfg.condN2PsatLowT; o.liquidCp = cfg.condN2LiquidCp;
    o.gasKgasModel = (cfg.condVaporMassFraction > 0.0) ? 1 : 0;   // CPG carrier (空気) は空気 Sutherland
    o.sigmaScale = cfg.condSigmaScale; o.Yw = cfg.condVaporMassFraction;
    o.h2oLatent = cond_latent_pair_for(cfg);
    return o;
}
int          cond_num_species();

// 二相拡散 (condTwoPhaseDiffusion; plan condensation-two-phase-transport §4.2, #4e) が実際に働く構成か (cfg だけで決まる)。
// TP carrier の凝縮 (condGasSpecies ≥ 0、多成分 TP) かつ粘性あり。CPG carrier・pure 凝縮・Euler・凝縮 OFF は false (現行経路)。
// 併用不可の設定 (dual-time 等) は main の起動時検査が拒否する。
inline bool condTwoPhaseDiffusionActive(const solverConfig& cfg)
{
    return cfg.condTwoPhaseDiffusion == 1 && cfg.condensation == 1 && cfg.nCondSpecies >= 1 && cfg.thermalMethod == 2
        && cfg.condGasSpecies >= 0 && cfg.nSpecies >= 2 && cfg.viscMethod != 0;
}
// 収束受入の独立残差監査を行う構成か (#4f (4) / #1b-pre (1))。二相拡散が働く run は新作用素、condAuditResidual 1 の
// TP carrier 凝縮 run (二相拡散 OFF) は旧作用素を監査する。
inline bool condResidualAuditActive(const solverConfig& cfg)
{
    if (condTwoPhaseDiffusionActive(cfg)) return true;
    return cfg.condAuditResidual == 1 && cfg.condensation == 1 && cfg.nCondSpecies == 1 && cfg.thermalMethod == 2
        && cfg.condGasSpecies >= 0 && cfg.nSpecies >= 2;
}
// 二相拡散の起動時検査とログ (main が bcond 読込後に 1 回呼ぶ)。併用不可の設定は理由を出して終了する。
void condTwoPhaseDiffusionValidate(const solverConfig& cfg);
// 蒸気の残差 res_roYv = res_roY_w − res_rog (監視・residual_history の rms_roYv 列; 周期集約の後に呼ぶ)。
void twoPhaseVaporResidual_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// 非分割更新 (定常; 設計メモ §6.1 の vl_limit_commit)。化学種の更新 (水は commit しない) と流れの更新の後、
// 凝縮モーメントの N 退避 (condensationUpdateOuter) の後に呼ぶ。液 ρg・Q と総水分 ρY_w を commit し、再正規化 (係数を液・Q にも) と
// 受動種の最後の砦 (floor の記録)・周期ミラーまで行う。化学種の周期ミラーと Y の同期は呼び出し側。
void twoPhaseUpdate_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// 化学種の更新が commit した水 ρY_w を更新前 (roY_w^N) に戻す (水は twoPhaseUpdate が蒸気 + 液から commit する)。
void twoPhaseHoldWater_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// 二相更新の監視 (θ・保留量・状態補正) を 1 行出す (condCorrectionLog から monitorInterval ごと)。
void twoPhaseUpdateLog(solverConfig& cfg, int iStep);
// 二相更新の診断 (condTwoPhaseDiag, #1b-r1; 読むだけ)。twoPhaseDiagRenormPtr: 再正規化がセルごとの f−1・Δ(ρY_w)・Δ(ρg) を書く先 (無効なら nullptr)。
// twoPhaseDiagWrite: 終了時に末尾 200 更新の θ = 0 (または θ < 1) のセルを CSV に書く (run ディレクトリ)。
double* twoPhaseDiagRenormPtr();
bool    twoPhaseDiagInWindow(const solverConfig& cfg);   // 診断の窓 (末尾 200 更新) の更新か (次に行う更新の番号で)
void    twoPhaseDiagWrite(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// θ の全更新を覆う集計 (#1b-pre (3)): kind 0 = θ_src (condensationSource の直後), 1 = 更新の θ (モーメント更新の直後)。計上だけ。
void condThetaScan_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int kind);

// 原始量 φ = ρφ/ρ を全セル (ghost 含む) について更新する。スカラ移流の上流値に使う。
void condensationPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 凝縮モーメント ghost を埋める (入口 inlet_* は dry=0 の Dirichlet、他は Neumann zero-gradient)。
void condensationBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var);

// 全 bcond をループして凝縮モーメント ghost を埋める (applySpeciesBoundaries と同形)。
void applyCondensationBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// モーメント実現可能性射影 (Q1²≤Q0Q2, Q2²≤Q1Q3) の作動数: device カウンタと、読み出して 0 に戻す host 関数 (monitor ログ用)。
int* condRealizViolCounter();
int  condRealizViolReadReset(int* degenerate = nullptr);   // 戻り値: 最近点射影の作動数, *degenerate: 退化の単分散再初期化数
void condensationRealizabilityProject_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);   // Q1/Q2 だけの射影 (EOS 後)
double* condClampBudget(int s);                    // 種 s の成分別収支スロット (device, 8 doubles)
std::vector<double> condClampBudgetTotals(int nSpecies);   // [s*8 + 2k(+1)]: g,Q0,Q1,Q2 の符号付き/絶対 ∫Δq dV (全期間)
// 理由別の補正量監視 (plan condensation-two-phase-transport §4.3, #2): 種 s の device スロット (COND_REASON_N doubles;
// 配置は condensationRealizability_d.cuh の COND_REASON_*)。凝縮なし・s 範囲外は nullptr。
double* condCorrReasons(int s);
// 化学種再正規化 (species_renormalize_d) に渡す計上先: TP carrier の凝縮 run なら種 0 のスロットと凝縮種 index (*iw)、それ以外は nullptr・-1。
double* condCorrReasonsForRenormalize(const solverConfig& cfg, variables& var, int* iw);
// 理由別の補正量を monitorInterval ごとに 1 行/種で出す (区間値 + 累積 + 総液量比; 数値補正が総液量比 1e-6 超なら WARN)。凝縮なしは no-op。
void condCorrectionLog_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int iStep);

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
