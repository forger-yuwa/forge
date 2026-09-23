#pragma once

// ============================================================================
// 共役熱伝達 (CHT) の界面契約 — **符号・定義の正本**
//
// 仕様は methods/boundary.md「共役熱伝達 (CHT)」、設計判断は
// plans/active/boundary-conjugate-heat-transfer.md (codex plan レビュー 3 巡)。
// 散文で符号を守らないこと: ここのコメントと methods/boundary.md が唯一の正本で、
// 他所は必ずこれを引く。
//
// 【符号規約】
//   * 面流束 F^E は **流体 CV から外向きを正**。
//   * 拘束反力 C は **流体への供給を正** (node 等温壁は温度ピン後に res_roe を 0 化するので、
//     壁 CV が実際に授受した熱はこの C に入る)。
//   * 界面熱量 (**連成に渡す正本**) は
//         Q_f = Σ_{壁面} F^E − C        [W]  (平面 2D は W/m)
//     で、**固体へ入る向きを正**とする。
//     検算: Σ F^E = 80, C = −20 なら Q_f = 100。
//     定常の Dirichlet 行では C = −R^raw。**過渡では C = D_t(V E) − R^raw** であり、
//     定常式を瞬時入熱に使ってはならない。
//
// 【診断量 (本ヘッダが提供するのはこちら。Q_f そのものではない)】
//   拘束反力 C の採取は plans/active/tooling-energy-balance-diagnostics.md が提供する。
//   それが入るまでの間、および入った後も精度診断として、次を**固体向き正**で出す。
//     q_compact = k_eff (T_1 − T_w) / d_1     … コンパクト差分形。SU2 CHT の界面転送と同じ形で、
//                                               界面抵抗 (D_f = k_eff A / d_1) と同じ離散。
//     q_recon   = − qwall                     … 再構成勾配形。viscousFlux_d.cu が res_roe に入れた値
//                                               (qwall は「壁→流体が正」なので符号を反転する)。
//     q_2nd     = k_eff dT/dn|_0 (3 点片側)   … 第一・第二内部点を使う 2 次片側差分。
//   **3 つは一致しない**。実測 (case/48 run_0014, 壁法線に整列した node メッシュ):
//   q_recon は q_compact と 2.6e-8 相対で一致する (この配置では再構成勾配がコンパクト差分に帰着する) が、
//   q_2nd は数 % ずれる。**どれを見ているかを必ず明示する**。
//
// 【第一内部点 (T_1, d_1) の定義】
//   tools/check_wall_resolution.py と同一規則にする (数字が食い違うと切り分け不能になるため):
//     1. 壁 DOF ごとに、その DOF に属する境界面の面ベクトルを**合算**して単位法線 n̂ を作る。
//     2. 隣接 DOF のうち |d·n̂|/|d| (整列度) が最大のものを第一内部点とする。
//     3. d_1 = |d·n̂|。整列度が alignMin 未満なら**評価不能** (角・斜交で第一内部点が定まらない)。
//        評価不能を 0 距離や合格に変換しない。
//   値の位置は node モードでは**ノード座標** (T の DOF 位置。壁ノードは壁面上に乗る)、
//   cell モードでは centCoords を使う。
// ============================================================================

#include <string>
#include <vector>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

namespace conjugateWall {

// 壁 DOF ごとの第一内部点 (bplane 順)。
struct FirstInterior {
    std::vector<geom_int> jdof;   // 第一内部点の DOF index (評価不能は -1)
    std::vector<geom_int> jdof2;  // 第二内部点 (2 次片側差分用。見つからなければ -1)
    std::vector<double>   d1;     // 法線方向距離 [m] (評価不能は NaN)
    std::vector<double>   d2;     // 第二内部点までの法線方向距離 [m]
    std::vector<double>   align;  // 整列度 |d·n̂|/|d|
    std::vector<double>   nx, ny, nz; // 単位法線
    std::vector<char>     ok;     // 1 = 評価可
};

// 第一内部点マップを作る (physID ごとにキャッシュ。メッシュは実行中に変わらない)。
const FirstInterior& firstInterior(const solverConfig& cfg, const mesh& msh, const bcond& bc);

// 壁面ダンプ用の界面診断を bc.diagVar に詰める (host のみ。デバイス bvar は触らない)。
// 呼び出しは bc.copyVariables_bplane_D2H() の**後**。cfg.interfaceDiag != 1 なら何もしない。
void fillInterfaceDiagnostics(const solverConfig& cfg, const mesh& msh, variables& var, bcond& bc);

// 温度を拘束する壁どうしが CV を共有していて、**そこに与える壁温が食い違う**構成を検出する。
// node の温度ピンは bcond ごとに順に適用されるので (nodeWallDirichlet_d.cu)、角ノードでは
// **最後に適用した bcond が勝つ** = 壁温が設定順で決まってしまう。CHT では隣接壁が別々の
// 壁温を持つので必ず踏む。**壁温を陽に扱っている run (wallProfile か interfaceDiag が有効) では
// 起動時エラーにする**。それ以外の run は挙動を変えないため何もしない。
// 呼び出しは applyWallProfiles の後 (per-face Ts が確定した状態で見る)。
void checkWallTemperatureSharing(const solverConfig& cfg, const mesh& msh);

// ソルバ内 CHT (Phase 2a): `conjugate:` ブロック + bcond `ints: {conjugate: 1}` の壁で、
// interval step ごとに壁温を更新する。
//
//     g_f = k_eff / d_1  [W/m2K]  (流体側の第一内部点までのコンダクタンス)
//     g_s = 1 / R_tot,  R_tot = t/k_s + R_back
//
// **`flux: q_eff` (既定, 保存形)** — plan boundary-conjugate-heat-transfer §4.2 の更新式:
//
//     (g_s + D_f) T_w^{k+1} = g_s T_b + q_eff(T_w^k) + D_f T_w^k,   D_f = g_f (初期推定)
//
// 収束すると g_s (T_w - T_b) = q_eff、すなわち**保存形の界面熱量**と固体の 1 次元法則が釣り合う。
// D_f は収束速度だけを決める (固定点は D_f に依らない)。`output: {interfaceDiag: 1}` が要る。
//
// **`flux: q_compact` (旧実装)** — 抵抗加重平均 T_w^{new} = (g_f T_1 + g_s T_b)/(g_f + g_s)。
// これは q_compact = k_eff (T_1 - T_w)/d_1 の固定点であり、**保存形 q_eff とは一致しない**:
// 差は壁半 CV 内の粘性加熱 tau.u と流動仕事で、第一層厚 d_1 に比例する。
// 実測 (case/48, d_1=3.0 um): 同一状態の G-cons が q_eff で **1.77 % (FAIL)**、q_compact で
// 0.000017 % (更新式の固定点なので恒等)。前縁では節点差が +26 % に達する
// (plan boundary-conjugate-heat-transfer §5.1 #66、run `case/48.flat_plate_cooled_m4/run_0026_cht_qeff_gcons`)。
// したがって **q_compact は A/B のときだけ使う**。
//
// 初版の制限 (いずれも起動時に拒否): node 以外、dual-time、`mode != local1d`、背面断熱、
// `flux: q_eff` で `interfaceDiag != 1`。
// 面内伝導が要る場合は外部ループ (tools/cht_loop.py + solid_shell.py) を使う。
// 界面の収束判定 (G-if) の素材は run 直下の `conjugate_history.csv` に出る
// (step, physID, Tw 統計, max|dTw|, 未緩和の界面残差 res_abs_Wm2 / res_max_W / res_rel, q_total)。
void initConjugateWalls(const solverConfig& cfg, const mesh& msh);   // 起動時の検査 (拒否条件)
void updateConjugateWalls(const solverConfig& cfg, mesh& msh, variables& var, int iStep);
bool conjugateActive(const solverConfig& cfg, const mesh& msh);
// 収束した壁温を run ディレクトリに残す (再開時は wall_profile_<physID>.csv にコピーして使う)。
void writeConjugateState(const solverConfig& cfg, const mesh& msh, int iStep);
// 再開時の累積 step オフセット (出力側が壁ダンプに `step_abs` を書くために使う)。
int stepOffsetForOutput();

} // namespace conjugateWall
