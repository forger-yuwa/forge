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

> 注: 上記 SLAU/dual 系は既存の入力リファレンス群。本表は新規 run 追加・破棄時に同期する。
