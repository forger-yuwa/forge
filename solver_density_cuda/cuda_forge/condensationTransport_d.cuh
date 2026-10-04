#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include <vector>
#include <string>
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

// 二相拡散の実効状態 (plans/active/condensation-two-phase-default.md §4-2)。ソルバ・tools/twophase_state.py・res_*.h5 属性・
// stage_manifest で同じ判定を使う (Python 側を変えたらここも変える)。指定値 (省略 / 0 / 1) は見ない — 構成だけで決まる。
//   Active       : 包絡内 (TP carrier・非平衡・NS・定常陰解法・nCondSpecies 1・非軸対称・非周期 …)。指定 ON で作動する。
//   InactiveA    : 物理が同一 (凝縮 OFF・viscMethod 0 で拡散自体が無い)。どの指定でもビット一致。
//   InactiveB    : モデルが構造的に適用できない (CPG carrier・pure 凝縮・平衡凝縮 EOS 拘束形 condEquilibrium 2)。不活性 + 毎回 WARNING。
//   UnsupportedC : 実装が未対応 (dual-time・RK/陽解法・speciesImplicitCoupling 2・passiveScalarScheme 0・condEquilibrium 1・
//                  condLimiterMode 0・nCondSpecies ≥ 2・軸対称・周期)。指定 ON ならエラー終了。
// 判定順は (a) → (b) → (c) (dual-time + Euler は (a) で通す)。hasPeriodic は bcond に periodic があるか。reason は理由 (英語)。
enum class TwoPhaseDiffusionState { Active, InactiveA, InactiveB, UnsupportedC };
inline const char* twoPhaseDiffusionStateName(TwoPhaseDiffusionState st)
{
    switch (st) {
        case TwoPhaseDiffusionState::Active:       return "active";
        case TwoPhaseDiffusionState::InactiveA:    return "inactive-a";
        case TwoPhaseDiffusionState::InactiveB:    return "inactive-b";
        case TwoPhaseDiffusionState::UnsupportedC: return "unsupported-c";
    }
    return "unknown";
}
inline TwoPhaseDiffusionState condTwoPhaseDiffusionClassify(const solverConfig& cfg, bool hasPeriodic, std::string* reason)
{
    using S = TwoPhaseDiffusionState;
    auto ret = [reason](S st, const char* why) { if (reason != nullptr) *reason = why; return st; };
    // (a) 物理が同一
    if (cfg.condensation != 1 || cfg.nCondSpecies < 1) return ret(S::InactiveA, "condensation is off");
    if (cfg.viscMethod == 0) return ret(S::InactiveA, "inviscid (viscMethod 0): no diffusion at all");
    // (b) モデルが構造的に適用できない (液は拡散しない現行のまま; 既定経路のエネルギー流束は液を蒸気として数える近似が残る)
    if (cfg.condGasSpecies < 0 || cfg.thermalMethod != 2 || cfg.nSpecies < 2)
        return ret(S::InactiveB, (cfg.condVaporMassFraction > 0.0) ? "CPG carrier" : "not a TP carrier (pure condensible)");
    if (cfg.condEquilibrium == 2)
        return ret(S::InactiveB, "EOS-constrained equilibrium condensation (condEquilibrium 2; the liquid is an EOS state, rog is not transported)");
    // (c) 実装が未対応
    if (cfg.unsteady == 1 && cfg.dualTime == 1)
        return ret(S::UnsupportedC, "cannot be combined with dual-time (unsteady 1, dualTime 1): the coupled FCT of vapour and liquid is not designed yet");
    if (cfg.unsteady != 0 || cfg.timeIntegration != 11)
        return ret(S::UnsupportedC, "requires steady implicit pseudo-time (unsteady 0, timeIntegration 11)");
    if (cfg.speciesImplicitCoupling == 2)
        return ret(S::UnsupportedC, "cannot be combined with speciesImplicitCoupling 2 (EOS cross coupling commits water before the liquid)");
    if (cfg.passiveScalarScheme != 1) return ret(S::UnsupportedC, "requires passiveScalarScheme 1 (moments on the passive-scalar path)");
    if (cfg.condEquilibrium != 0)
        return ret(S::UnsupportedC, "requires non-equilibrium condensation (condEquilibrium 0; the relaxation form condEquilibrium 1 is not supported)");
    if (cfg.condLimiterMode != 1) return ret(S::UnsupportedC, "requires condLimiterMode 1 (the source enters the residual without theta)");
    if (cfg.nCondSpecies != 1) return ret(S::UnsupportedC, "supports one condensing species (nCondSpecies 1)");
    if (cfg.isAxisymmetric != 0)
        return ret(S::UnsupportedC, "is not verified on axisymmetric meshes (no ON Navier-Stokes test yet; plan condensation-two-phase-default #7)");
    if (hasPeriodic)
        return ret(S::UnsupportedC, "is not verified with periodic boundaries (no ON Navier-Stokes test yet; plan condensation-two-phase-default #7)");
    return ret(S::Active, "inside the verified envelope");
}
// 二相拡散 (condTwoPhaseDiffusion; plan condensation-two-phase-transport §4.2, #4e) が実際に働く構成か (cfg だけで決まる)。
// 指定 ON かつ実効状態が Active。周期境界は cfg から分からないが、指定 ON + 周期は condTwoPhaseDiffusionValidate が起動時に止める。
inline bool condTwoPhaseDiffusionActive(const solverConfig& cfg)
{
    return cfg.condTwoPhaseDiffusion == 1 && condTwoPhaseDiffusionClassify(cfg, false, nullptr) == TwoPhaseDiffusionState::Active;
}
// 収束受入の独立残差監査を行う構成か (#4f (4) / #1b-pre (1))。二相拡散が働く run は新作用素、condAuditResidual 1 の
// TP carrier 凝縮 run (二相拡散 OFF) は旧作用素を監査する。
inline bool condResidualAuditActive(const solverConfig& cfg)
{
    if (condTwoPhaseDiffusionActive(cfg)) return true;
    return cfg.condAuditResidual == 1 && cfg.condensation == 1 && cfg.nCondSpecies == 1 && cfg.thermalMethod == 2
        && cfg.condGasSpecies >= 0 && cfg.nSpecies >= 2;
}
// 二相拡散の起動時検査とログ (main が bcond 読込後に 1 回呼ぶ)。実効状態を判定して cfg.condTwoPhaseDiffusionState /
// condTwoPhaseDiffusionEffective に書き、[twophase] 行を出す。指定 ON + 未対応 (c) は理由を出して終了する。hasPeriodic: bcond に periodic があるか。
void condTwoPhaseDiffusionValidate(solverConfig& cfg, bool hasPeriodic);
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
// #4h: 二相の更新ごと・成分ごとの補正計測。End は実現可能性クランプの後 (更新を閉じて履歴へ)、Log は区間 (false) / 末尾窓 ceil(0.1N) と κ の VERDICT (true)。
void    twoPhaseCorrGateEnd();
void    twoPhaseCorrGateLog(const solverConfig& cfg, int iStep, bool final);
void    twoPhaseDiagWrite(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);
// θ の全更新を覆う集計 (#1b-pre (3)): kind 0 = θ_src (condensationSource の直後), 1 = 更新の θ (モーメント更新の直後)。計上だけ。
void condThetaScan_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int kind);

// ---- 診断 G3-b (更新写像の収支; FORGE_DIAG_TP_UPDATE=<出力 h5>、既定 off; plan condensation-two-phase-default §5.1 #4g3) ----
// main の runTpUpdateDiag だけが tpuBegin/tpuArm を呼ぶ。それ以外の経路では tpuSlots() は常に nullptr、tpuSnap は何もしない
// (本番カーネルの記録分岐は nullptr で無効、演算は変えない)。スロット配置は twoPhaseUpdateDiag_d.cuh。
//   tpuBegin: 対応構成を検査し (why に拒否理由)、前処理用・更新用の 2 本のスロット (NaN 初期化) を確保する。
//   tpuArm(phase): 0 = 記録しない、1 = 前処理のバッファへ、2 = 更新のバッファへ。
//   tpuSnap(label): 記録中なら 5 成分 (ρY_w, ρg, ρQ2, ρQ1, ρQ0) の格納値 (実節点) を写して label 付きで残す (操作の境界の照合用)。
//   tpuCollect: スロットと写しを host へ移す (終了時に 1 回)。
struct TpuDiagData {
    long n = 0;                                      // 実節点数 (msh.nCells)
    int iw = -1;                                     // 総水分の化学種 index
    std::vector<double> pre, upd;                    // [TPU_NSLOT][n] (前処理 / 更新)
    std::vector<std::string> snapLabel;              // 写しの名前 (起きた順)
    std::vector<int> snapPhase;                      // 写しを取ったときの phase (1 前処理 / 2 更新)
    std::vector<std::vector<float>> snap;            // [nsnap][5·n] (成分 c の節点 i は c·n + i)
};
bool    tpuBegin(solverConfig& cfg, mesh& msh, variables& var, std::string& why);
void    tpuArm(int phase);
double* tpuSlots();
void    tpuSnap(const char* label);
bool    tpuCollect(TpuDiagData& out);

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
