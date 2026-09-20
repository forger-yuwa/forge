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
| `run_0025`/`run_0026_std8_nobnd_*` | **禁止フラグ `mesh.bndFirstOrder: 1` を外して再実行** (保管 run から複製する際に落とし忘れていた。codex 2026-09-20 result Major 7) | **step 0 残差も発散も変わらない** (`rms_roUy` 2.01e-4、DIVERGED)。→ 退行の原因は bndFirstOrder ではない。**保管 run との step 0 差は `rms_roUy` 21.70760 vs 0.000200996 = 10⁵ 倍** (`rms_ro` と `rms_roe` だけ一致)。bump は **slip 壁**なので node slip 境界の扱いの変化が疑われる | ref (plan §4.31) |
| `run_0027_std8_pref` / `run_0028`〜`run_0031_std8_*` | **W1i: 「退行」の追跡** (`pRef: 21839.5` 追加 + `convMethod`/`cfl_pseudo` スイープ) | **退行ではなく自由流保存の修正だった**: IC は厳密に一様 (ρ std 2.98e-8, **Uy std 0**) なので `rms_roUy≈0` が正しく、保管 run の 21.708 がスプリアスな壁力。`pRef` で 2.01e-4 → **1.54e-5**。**`convMethod: 0` + `pRef` は `PASS (converged)` 全列 5.1〜5.3 桁**。1 次以外 (cm 1/2, cfl_pseudo 1/3) は全て DIVERGED | ref (**plan §4.33**) |
| `run_0032`〜`run_0035_lim8_*` | **リミッタ修正で 2 次が救えるか** (cm 1/2 × ① / ② K=0.05) | **4 本とも DIVERGED** (cm2 は step 182)。上下が `slip` 壁なので [[node-slip-spurious-flow]] が本命で、リミッタの問題ではない。→ **標準ケース回帰は「node + 2 次で収束するケースが無い」ためブロック継続** | ref (plan §4.33) |
| `run_0036_slip8_noslip` / `run_0037_slip8_slipref` | **slip 壁説の検証** (上下を `slip` → `wall` に変えて `convMethod: 1`) | **両方 DIVERGED** → slip は原因でない。`case/05` Sod は slip 壁 2 枚・`convMethod: 2` で完走するのでも裏づけ。**case/08 の 2 次発散の原因は未特定** | ref (plan §4.33) |
| `run_0038_diag2nd` / `run_0039_diag2nd_dense` | **2 次発散の 1 step 刻み追跡** (`detectNaN: 1`, `outStepStart: 150, outStepInterval: 1`) | **発散は入口の 1 節点から**: NaN は 6076 中 **1 節点**、座標 **(0.0000, 0.5208)** = 入口 `x=0` 上 (壁でも角でもない)。step 195 まで場は収束方向。id 317 の履歴は ρ 0.2604→0.0068→負、Ux 511.8→628.8→856.5 で **ρ↓→Ux↑→Ps↓→ρ↓ の正帰還** | ref (**plan §4.34**) |
| `run_0041_inletuni_cm0` / `run_0042_inletuni_cm1C` | **入口を揃えた 1 次/2 次・①/② 対照** (全て `inlet_uniformVelocity`, 5000 step。前版は 1 次が `inlet_Pressure` で交絡していた) | **1 次 `PASS (converged)`** (`rms_ro` 4.63e-7) / 2 次 ① `NOT CONVERGED` (4.90e-3) / **2 次 ② K=0.05 `NOT CONVERGED` (1.11e-3)**。→ **2 次のプラトーは本物** (交絡でない)、かつ **② は ① より 3〜5 倍良い** (`rms_ro` 4.4 / `roUx` 3.0 / `roe` 4.8 倍) が解決ではない | ref (**plan §4.35 (3)**) |
| `run_0043`/`run_0044_inletchk*`, `run_0045`/`run_0046_inletfix_*` | **`inlet_Pressure` CPG 修正の検証** (境界出力から実 Tt/Pt を逆算、新旧バイナリ) | 修正後 **Tt 293.1500 / Pt 99999.98 = 誤差 0.0000 %**。ただし**修正前も誤差 0.000 %** — この IC は指定全条件と等エントロピーで整合しており内点音速=新音速になるため。**欠陥は過渡・オフデザインでのみ出る**。発散は修正後 step 204 / 前 202 で**変わらない** | ref (plan §4.36) |
| `run_0047_plateau`, `run_0048`〜`run_0062_pl*` | **2 次プラトーの原因特定** (支配節点の特定 → CFL 依存性 → リミッタ A/B → K スイープ) | **支配節点は上壁の衝撃反射点** (x 1.46–2.83 / y 0.81–1.00、下壁 0/50・入口 0/50)。プラトーは **CFL に比例** (cfl 1→6 で 1.18e-3→4.99e-3) = 反復律速のリミットサイクル。**`limiter: 0` なら 2 次でも 5.110e-07** (1 次と同等) = **原因はリミッタ**。K スイープは**収束と有界性を交換するだけ** (K=0.05 は逸脱 ro/Ux/P とも 0 で rms 1.21e-3、K=50 は rms 5.2e-7 だが逸脱 255 万) | ref (**plan §4.37**) |
| **`run_0040_inletuni`** | **入口 BC を `inlet_Pressure` → `inlet_uniformVelocity` に替えて 2 次** (`convMethod: 1`, `limiter: 2` のまま) | **5000 step 完走・NaN なし** (`NOT CONVERGED (stalled/plateau)`, 1.2〜1.4 桁)。→ **真因は case の config**: 入口 M=1.650 (超音速) なのに `inlet_Pressure` (亜音速の全条件入口) を使っていた。forge の欠陥ではない。残る課題は 2 次のプラトー (1 次は 5.1〜5.3 桁) | **active (case/08 node 2 次の正しい入口)** |

> 注: 上記 SLAU/dual 系は既存の入力リファレンス群。本表は新規 run 追加・破棄時に同期する。
