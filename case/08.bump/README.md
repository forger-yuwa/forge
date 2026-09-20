# case/08.bump — 2D 非粘性 bump チャネル

x∈[0,3], y∈[0,1] の 2D チャネル。下壁中央に円弧 bump (peak 0.1 @x=1.5)、上下壁 slip、非粘性。
低マッハ亜音速 (inlet Pt=120193 Pa / Tt=302.557 K、outlet Ps=101325 Pa → Pt/Ps≈1.19、M~0.5)。
スムーズな亜音速場でメッシュ収束しやすく、スキーム/離散化 (cell vs node)・CFL・陽/陰解法の確認や
SU2 とのライン比較の基準ケースとして使う。形状元: [mesh/bump.geo](mesh/bump.geo)。

検証フロー一式: [verify/](verify/) (低/高マッハ × 陽/陰で y≈0.3 の数値解一致と収束非劣化を自動判定)。

## 計算 run 一覧

詳細な考察・設計判断は各 `plans/` 側に置く。本表は「どの検証がどの run か」の一次情報。

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| `run_su2cmp` | **node/cell/SU2 三者ライン比較** (y=0.25, P/T/\|U\|)。SU2 2D Euler ROE/Venkat を基準に forge node(HLLE,収束) / cell(HLLE,plateau) を突き合わせ | node は SU2 と <0.1% 一致、cell は ~0.5% (atomicAdd 床で未収束の刻印)。`bump_line_compare.png` / `COMPARISON_SUMMARY.md` / `bump_euler.cfg` / `mesh/bump2d.geo` | ref |
| `run_0001`〜`run_0018_*` | SLAU 低/高マッハ × 陽/陰 × CFL スイープ (`*_exp_cfl*` / `*_imp_cfl*`) | 陰解法の安定 CFL 上限と陽解法比の速度評価 ([[bump-implicit-vs-explicit-result]]) | ref |
| `run_dual_*` | cell vs node (median-dual) 双対比較 (m2/m3, tri 含む) | 離散化レイアウト検証の入力 | ref |
| `run_0018`〜`run_0024_std8*` | **リミッタの標準ケース回帰**と、その過程で見つかった**退行の切り分け** (保管 `run_dual_hiM_node_m3imp` = node + 陰解法 + `convMethod 2` + `limiter 2` と同一入力。旧キー `LESorRANS: 0 / LESmodel: 0` は `model: "none"` に置換) | **保管 run は PASS (rms_ro 4.0 桁) だが現行バイナリでは全て DIVERGED**: 生産既定 step 217 / ② K=0.05 step 240 / 診断キー無し step 228 / **`build-native` (09-19 10:16、本セッションの変更前) でも step 243** → **本セッションの変更ではない**。step 0 残差は保管とビット一致、有効設定値にも差は無い。`build-passive` (09-18) は削除済み `meshFormat` を要求するので同一 config で遡れない | ref (**plan convection-node-wall-reconstruction §4.31 / 退行は未解決**) |

> 注: 上記 SLAU/dual 系は既存の入力リファレンス群。本表は新規 run 追加・破棄時に同期する。
