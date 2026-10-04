#pragma once
// =============================================================================
// twoPhaseUpdateDiag_d.cuh — 診断 G3-b (更新写像の収支; FORGE_DIAG_TP_UPDATE=<出力 h5>、既定 off) の記録スロットの配置
//   plans/active/condensation-two-phase-default.md §5.1 #4g3、設計の判断 notes/reviews/2026-10-04-twophase-g3-design-diagnose.md。
//
// 本番カーネルの中で、診断ポインタ (tpuSlots(); 既定 nullptr) が非 null のときだけ、そのカーネルが**実際に読んだ値・書いた値**を
// 節点ごとの double スロットへ写す (積算は後処理で double)。本番の演算式・評価順は変えない (記録は既存の値を読むだけ)。
//   成分 c: 0 = ρY_w (総水分), 1 = ρg (液), 2 = ρQ2, 3 = ρQ1, 4 = ρQ0 (ρv = ρY_w − ρg は後処理で double で作る)。
//   操作 op ごと・成分ごとに before / after (操作が読んだ格納値と書いた格納値; float の格納値を double にしたもの)。
//   commit 操作 (SPC・TPC・CMC) は加えて δq_limited (double; 制限後の増分を float の入力から double で組んだ値)、
//   δq_limited の非融合 float 評価 (参考; 本番は FMA で融合しうるので「本番の float の増分」は存在しないことがある)、
//   格納前の float 候補 cand_f (下限・丸め補正の前)、double 候補 cand_d = before + δ、格納で失われた増分 (cand_f == before かつ δ ≠ 0 のとき δ)。
//   恒等式: after − before = δ + C_round + C_floor、C_round = cand_f − cand_d、C_floor = after − cand_f (commit 内の床・vround)。
// スロット配置: buf[slot·n + ic] (n = msh.nCells、実節点のみ)。書かれないスロットは NaN のまま (その操作は経路に無い)。
// main.cpp (g++) からも include する (定数と名前表だけ; __device__ 関数は __CUDACC__ のときだけ)。
// =============================================================================
#include "flowFormat.hpp"

#define TPU_NC 5   // 成分数 (w, g, Q2, Q1, Q0)

// 操作 (更新の中で起きる順; 前処理の実現可能性クランプも同じ番号で前処理用バッファへ書く)
enum TpuOp {
    TPU_OP_SPC = 0,   // 化学種の commit ρY = max(ρY_N + δ, 0) (水だけ記録; speciesTransport_d.cu species_commit_correction_d)
    TPU_OP_TPC,       // 二相の非分割 commit (twophase_vl_update_d / tp_vl_update; 総水分・液・Q)
    TPU_OP_CMC,       // OFF のモーメント commit (cond_moment_update_limited_passive_d; 液・Q)
    TPU_OP_RNF,       // 再正規化の中の負値の 0 化 (水だけ記録)
    TPU_OP_RNS,       // 再正規化の係数の乗算 (OFF: 水; ON: 水・液・Q)
    TPU_OP_PFL,       // 受動種の床 (passive_bounds_d; 液・Q)
    TPU_OP_ARH,       // 受動種の φ_N δρ 項 (passive_add_rho_term_d; OFF の液・Q。補正でなく更新の一部)
    TPU_OP_LIM,       // 受動種の増分スケーリング θ_b (passive_limit_increment_d; OFF の液・Q)
    TPU_OP_RZL,       // 実現可能性クランプの下限 (ρg < 0 → 0、ρQn < 0 → 0)
    TPU_OP_RZU,       // 実現可能性クランプの上限 (ρg > ρY_w → ρY_w)
    TPU_OP_RZP,       // モーメント射影 (Q1・Q2; 表現できない塵の Q1 = Q2 = 0 を含む)
    TPU_OP_RZR,       // 液滴消滅 (物理; 数値補正とは別に集計)
    TPU_NOP
};
enum TpuX { TPU_X_DELTA_D = 0, TPU_X_DELTA_F, TPU_X_CAND_F, TPU_X_CAND_D, TPU_X_LOST, TPU_NX };
#define TPU_NBA (TPU_NOP*TPU_NC*2)
#define TPU_BA(op, c, k) ((((op)*TPU_NC + (c))*2) + (k))          // k: 0 before, 1 after
#define TPU_XS(ci, c, x) (TPU_NBA + (((ci)*TPU_NC + (c))*TPU_NX) + (x))   // ci: 0 SPC, 1 TPC, 2 CMC
#define TPU_RN_FACTOR (TPU_NBA + 3*TPU_NC*TPU_NX)                // 再正規化の係数 ρ/ΣρY (double)
// 診断 #4pj (射影の打ち消しの判別 A/B; plan condensation-two-phase-default §5.1 #4pj):
//   二相 DPLUR の分母 (tp_dplur_prep_d)。行 q: 0 = 蒸気, 1 = 液 ρg, 2 = ρQ2, 3 = ρQ1, 4 = ρQ0。k: 0 = 本番の分母 (tp_denoms)、
//   1 = sweep が実際に使った分母 (FORGE_DIAG_TP_COMMON_DIAG=1 のときだけ液・Q の行が節点ごとの max(D_g, D_Q2, D_Q1, D_Q0) に置き換わる)。
#define TPU_DIAG_BASE (TPU_RN_FACTOR + 1)
#define TPU_DIAG(q, k) (TPU_DIAG_BASE + (q)*2 + (k))
//   実現可能性クランプ (float 実体) の射影が使う量: T、ρ_l(T) (射影と同じ cond_rho_cond)、Q3 = ρg/((4/3)π ρ_l) (ρg は下限・上限の後)、
//   射影の分岐 kind: −1 = 分岐外 (ρQ0 ≤ 0 または ρg ≤ 0; g = 0 の塵は後段の消滅処理)、0 = 変更なし、1 = 最近点射影、
//   2 = 退化の単分散再初期化、3 = Q3 が表現できず (x, y を作れず) Q1 = Q2 = 0 にする枝。
#define TPU_RZ_T      (TPU_DIAG_BASE + 2*(2 + 3))
#define TPU_RZ_RHOL   (TPU_RZ_T + 1)
#define TPU_RZ_Q3     (TPU_RZ_T + 2)
#define TPU_RZ_KIND   (TPU_RZ_T + 3)
#define TPU_NSLOT     (TPU_RZ_T + 4)

static inline const char* tpuOpName(int op)
{
    static const char* nm[TPU_NOP] = {"SPC_species_commit", "TPC_twophase_commit", "CMC_moment_commit", "RNF_renorm_negfloor",
        "RNS_renorm_scale", "PFL_passive_floor", "ARH_add_rho_term", "LIM_limit_increment", "RZL_realiz_lower",
        "RZU_realiz_upper", "RZP_realiz_projection", "RZR_droplet_removal"};
    return (op >= 0 && op < TPU_NOP) ? nm[op] : "?";
}
static inline const char* tpuCompName(int c)
{
    static const char* nm[TPU_NC] = {"w", "g", "Q2", "Q1", "Q0"};
    return (c >= 0 && c < TPU_NC) ? nm[c] : "?";
}
static inline const char* tpuXName(int x)
{
    static const char* nm[TPU_NX] = {"delta_d", "delta_f_unfused", "cand_f", "cand_d", "lost"};
    return (x >= 0 && x < TPU_NX) ? nm[x] : "?";
}

#ifdef __CUDACC__
__device__ inline void tpu_put(double* s, int slot, geom_int n, geom_int ic, double v) { s[(size_t)slot*(size_t)n + (size_t)ic] = v; }
__device__ inline void tpu_ba(double* s, int op, int c, geom_int n, geom_int ic, double before, double after)
{
    tpu_put(s, TPU_BA(op, c, 0), n, ic, before);
    tpu_put(s, TPU_BA(op, c, 1), n, ic, after);
}
// commit の記録 (ci: 0 SPC, 1 TPC, 2 CMC)。before/cand_f/after は格納 (または格納前) の float、delta_d は double の制限後増分。
__device__ inline void tpu_commit(double* s, int op, int ci, int c, geom_int n, geom_int ic,
                                  float before, double delta_d, float delta_f, float cand_f, float after)
{
    tpu_ba(s, op, c, n, ic, (double)before, (double)after);
    const double cand_d = (double)before + delta_d;
    tpu_put(s, TPU_XS(ci, c, TPU_X_DELTA_D), n, ic, delta_d);
    tpu_put(s, TPU_XS(ci, c, TPU_X_DELTA_F), n, ic, (double)delta_f);
    tpu_put(s, TPU_XS(ci, c, TPU_X_CAND_F), n, ic, (double)cand_f);
    tpu_put(s, TPU_XS(ci, c, TPU_X_CAND_D), n, ic, cand_d);
    tpu_put(s, TPU_XS(ci, c, TPU_X_LOST), n, ic, (cand_f == before && delta_d != 0.0) ? delta_d : 0.0);
}
#endif
