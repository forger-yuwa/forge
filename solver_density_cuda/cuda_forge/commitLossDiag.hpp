#pragma once

// 診断 FORGE_DIAG_COMMIT_LOSS=<N> (plans/active/architecture-float-state-double-geometry.md §4.5・§6 V4 の記録、既定 off)。
//
// 定常の陰解法の commit (Q_after = Q_before + dq_req) で、要求した更新 dq_req のうち丸めで消えた分を数える。
// 外側の step の番号 step = iStep + 1 (res_<step>.h5 と同じ数え方。residual_history.csv の step 列は iStep = step − 1) が
// N の倍数の step だけ、commit のカーネルの直後に Q_before・dq_req・Q_after を読んで量ごとに double で集計する。
// 状態には書かない。commit のカーネル (SST・化学種) には増分の写しを書く省略できる引数を足したが、off では nullptr で、
// 演算の式は変えていない。off では main の分岐 1 つだけで、カーネルの起動も確保もしない。
//
// 比べる値 (どれも commit のカーネルの直後。その後の等温壁のピン・SST の E_t 補正・化学種の再正規化・次の step の EOS の床の前):
//   - ρ・ρu_x・ρu_y・ρu_z・ρE (block DPLUR / スカラー DPLUR の commit、update_d.cu):
//       Q_before = roN (commit の基準。ro = roN + dq)、dq_req = dq_block_old_k / dq_*_old (implicitRelax を掛けた最終の補正)、
//       Q_after = commit 直後の ro..roe。軸の射影は本番で無効 (axis_flag = nullptr)。
//       qAccumulatorFP64 と updateGuardAlpha > 0 の経路は計らない (起動時に書く)。
//   - ρk・ρω (SST の point-implicit、update_d.cu): Q_before = カーネル直前の roK・roOmega (その場で足すので基準はこれ)、
//       dq_req = カーネルが計算した dk・dw (implicitRelax 込み。ピンで上書きする節点は 0 = 増分なし)、Q_after = カーネル直後。
//   - ρY_s (化学種): Q_before = roY{s}N、Q_after = commit 直後 (再正規化の前)。dq_req は
//       speciesImplicitCoupling 0 (point-implicit) で speciesImplicitRelax·δ (積は double で厳密)、1 (scalar DPLUR) で dq_roY{s}_old。
//       2 (EOS 結合) は計らない (起動時に書く)。
// 床に当たった節点 (Q_before + dq_req < 床。ρk は 0、ρω は 1e-20、ρY は 0) は丸めでなく床なので n_clip に数え、
// n_req 以下の集計から外す。
//
// 集計 (量ごと、領域ごと。領域は interior = どの境界にも属さない節点 (bnode_flag 0)、boundary = 境界の節点):
//   n_req  : dq_req ≠ 0 の節点数         n_lost : dq_req ≠ 0 なのに Q_after と Q_before がビット単位で同じ節点数
//   S_req  : Σ|dq_req|                  S_act  : Σ|Q_after − Q_before| (double で引く)
//   S_lostfrac = Σ|dq_req − (Q_after − Q_before)| / S_req
// 出力: ログに 1 step 1 行 ([commitLoss] step ... 量 n_req/n_lost S_act/S_req ...、interior だけ) と、
//       実行ディレクトリの commit_loss.csv (列: step,qty,n_req,n_lost,S_req,S_act,S_lostfrac,n_clip,region)。

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

namespace commitLossDiag {

namespace detail {
extern bool g_enabled;   // 環境変数で有効にした run
extern bool g_active;    // この step で計る (beginStep 〜 endStep の間)
}

// 有効か (main の step の頭で見る分岐)。
inline bool enabled() { return detail::g_enabled; }
// この step で計るか (commit の呼び出し側が見る)。
inline bool active() { return detail::g_active; }

// 起動時に 1 回 (main、時間更新の前)。環境変数を読み、有効なら集計の配列と増分の写しの置き場を確保する。
// 定常の陰解法 (isImplicit 1・unsteady 0・gpu 1) 以外で有効にされたら理由を出して止める。
void init(solverConfig& cfg, mesh& msh, variables& var);

// 定常の外側の step の頭と終わり (main の advanceImplicitSteady)。iStep は 0 始まり。
void beginStep(int iStep);
void endStep();

// 終了時に CSV を閉じて確保を返す。
void finalize();

// ---- commit の直後に呼ぶ (active() のときだけ) ----
// ρ..ρE。dq は 5 本 (ρ, ρu_x, ρu_y, ρu_z, ρE の順)。guardAlpha は block の updateGuardAlpha (スカラー DPLUR は 0)。
void flowCommit(solverConfig& cfg, mesh& msh, variables& var, flow_float* const* dq, flow_float guardAlpha);

// SST: カーネルの直前に Q_before を写し、カーネルが dk・dw を書く置き場を返す。直後に sstAfter。
void sstBefore(mesh& msh, variables& var, flow_float** dk, flow_float** dw);
void sstAfter(mesh& msh, variables& var);

// 化学種 coupling 0: point-implicit のカーネルが relax·δ を書く置き場 (double、全種で使い回す)。
double* speciesDqBuffer();
// 化学種 s の commit の直後 (coupling 0: dq は speciesDqBuffer、coupling 1: dq_roY{s}_old)。floor はカーネルの下限。
void speciesCommitPointImplicit(mesh& msh, int s, const flow_float* roYN, const flow_float* roY, const double* dq, double floor);
void speciesCommitDPLUR(mesh& msh, int s, const flow_float* roYN, const flow_float* roY, const flow_float* dq, double floor);

}  // namespace commitLossDiag
