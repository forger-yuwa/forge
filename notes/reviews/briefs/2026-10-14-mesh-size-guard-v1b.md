# 諮問ブリーフ: Euler の ni の決め方 — V1b の結果の解釈と plan 段レビューの採否 (plan tooling-nozzle-mesh-size-guard)

AGENTS.md 条件 3 (事前登録の判定が FAIL・判定不能)・条件 5 (codex の Major の採否)・条件 7 (result の解釈)。判断役のスイッチが codex。
前回の諮問 `notes/reviews/2026-10-14-mesh-size-guard-diagnose.md` (事前の判定を実効の節点数で) は採用済み (plan §4.2・§6 V3)。

## 観測事実

### V1b (事前登録は plan §6 V1b、2026-10-14 実施、commit 88675245 + 未 commit の plan/台本、forge バイナリ sha256 6c8cff55…)
case/70 の `problem_m5_air_c2b_euler.yaml` の `mesh_euler.ni` だけを変えた Euler 3 本 (soft + 本段 cfl 1・48,000 step、1 本通し、手元)。領域 x ∈ [−20.02, 48.07] r_t (L = 68.1)。

| run (`case/70.m5_air_wt_demo/`) | ni | Δx 最小 / スロート平均 / 最大 [r_t] | 残差 (`check_convergence --segment`、区間 main) | 出口コア M 末尾 5 枚 | 流量 (無次元) 末尾 5 枚の平均 ± 振幅 |
| --- | --- | --- | --- | --- | --- |
| `run_0015_euler_c2b_ni1800` | 1800 | 0.0117 / 0.0120 / 0.0467 | NOT CONVERGED (stalled/plateau) | 4.999365 (STEADY) | 10537.82 ± 1.4 (OSCILLATING) |
| `run_0016_euler_c2b_ni900` | 900 | 0.0234 / 0.0240 / 0.0935 | NOT CONVERGED (still converging) | 4.999518 (STEADY) | 10538.79 ± 0.29 (OSCILLATING) |
| `run_0017_euler_c2b_ni600` | 600 | 0.0351 / 0.0361 / 0.1403 | NOT CONVERGED (stalled/plateau) | 4.999744 (STEADY) | 10540.34 ± 0.37 (OSCILLATING) |

準定常は `qs_series.py RUN --tail 0.4 --drift 5e-5 --osc 5e-5` (= `check_quasisteady --series-csv`)。NaN: 3 本とも残差の全履歴・最終 res の VALUE に非有限値 0、ρ・P・T > 0。メッシュ品質 PASS。

抽出 (NS `run_0010b_ns_c2b_cont60k` の最終の場を各 Euler に対して `extract_and_merge(.., band_select="edge")`、台本 `case/70.m5_air_wt_demo/euler_ni_interp_ab.py --runs`、出力 `_euler_ni_ab/V1b/`):

| ni | max\|Δδ_E\| 使用域 [0, 22.23] | その x | 流量の比 | Δ(流量の比) | Δ出口コア M |
| --- | --- | --- | --- | --- | --- |
| 1800 | — | — | 1.0076110 | — | — |
| 900 | 6.0e-5 r_t | 10.7 | 1.0075124 | −9.86e-5 | +1.53e-4 |
| 600 | 3.0e-5 r_t | 9.9 | 1.0073679 | −2.43e-4 | +3.79e-4 |

事前の判定 (δ_E ≤ 1e-4 r_t・出口コア M ≤ 2.5e-4・流量の比 ≤ 1e-4、流量と M は STEADY を確かめてから): ni 600 (比 3) は出口 M で **FAIL**。ni 900 (比 2) は δ_E・M は合格、
流量は 3 本とも OSCILLATING で、比 2 の差 9.9e-5 は ni 1800 の振幅 (相対 1.3e-4) の中 → **判定不能**。事前の分岐は「比 2 だけ合格 → Euler の ni = max(min(mesh.ni, 2000), ⌈mesh.ni/2⌉)」。

当方の見立て (仮説): 差は 2 次の格子収束に乗る。出口 M の差 1.53e-4 (900)・3.79e-4 (600) の比 2.48 (2 次なら (9−1)/(4−1) = 2.67)、流量の差 0.97・2.52 の比 2.6。
Richardson で ni 1800 の M の格子誤差 ≈ 5e-5、流量 ≈ 3e-5 (相対)。δ_E の差は 900 > 600 で単調でない (帯の切り替わりの雑音の桁; 窓の 5 枚の平均の雑音は約 1e-4 r_t、accepted plan tooling-nozzle-deltastar-band-root-zone)。
**効くのは NS と Euler の点数の比ではなく、Euler の絶対の Δx (r_t 単位)** と考える。

### 関連の経緯
- case/45 (M6): 粗い Euler 格子 1100 × 65 で出口較正 (`Md_moc_offset`) を決めると、細分格子の NS で出口 M が約 −0.0008 ずれた (`procedures/nozzle-design-workflow.md:161`)。
  出口較正は Euler から決めて設計に入る (`runner_axismach.py:373-377`)。Euler の格子誤差は設計の出口 M に直接入る。
- 生成器は Euler の ni に NS の ni を写す (`intake/wt_card.py:562-563`)。NS の ni は壁際の AR から決まり上限 6000。

### plan 段のレビュー `notes/reviews/2026-10-14-tooling-nozzle-mesh-size-guard-plan.md` (GO-with-changes, C0/M3/m1)
- M1: 「比 3 合格 → 全ケースで上限 2000」は成立しない (配点は領域長・細分で決まる、`mesh2d.py:53`)。case/70 に限定し、ni > 2000 の形状で直接比較を。
- M2: 準定常のゲートが許容差に対応していない (流量の振動が 1e-4 の許容を覆う)。δ_E も複数の Euler 時刻で。覆うなら判定不能。
- M3: 要求カードの生成器の `estimate_mesh` は `_mesh_ar_full` (`wt_card.py:716`) で NS 格子を丸ごと作る → prepare の前に 58 万点を確保。点数の計算も配列を作らずに。
- m4: 明示の `overrides.mesh_euler.ni` を黙って上書きしない (`wt_card.py:1497`)。methods の更新も。

## 問い

**Euler の ni の既定をどう決めるかを 1 つ推奨してほしい。** 候補 (縛られなくてよい):
(A) 事前の分岐どおり max(min(ni_NS, 2000), ⌈ni_NS/2⌉) (比の規則; M1 は一般化できないと言う)、
(B) Euler の Δx の目標 (r_t 単位) で決める: ni_E = min(ni_NS, Δx の最大が case/70 の ni 1800 の値 0.047 r_t (M の格子誤差 約 5e-5) 以下になる最小の ni) — 比でなく絶対の解像度、
(C) 切り離さない (Euler = NS の ni のまま)。節点数の上限 (§4.2) と事前の提示 (§4.3) で「黙って大きな格子」を止める。2D の 58 万点は手元の GPU で 48,000 step 約 10 分。
あわせて、レビューの M1〜M3・m4 の採否 (各 1 行) と、M2 を受けて V1b の流量の判定不能をどう扱うか (窓を長くする / 量を変える / 判定から外す) を。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md`, `_euler_ni_ab/**/*.csv` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,130p' plans/active/tooling-nozzle-mesh-size-guard.md`
- `sed -n '1,80p' notes/reviews/2026-10-14-tooling-nozzle-mesh-size-guard-plan.md`
- `sed -n '40,80p' design/forge_design/meshing/mesh2d.py`
- `sed -n '150,170p' procedures/nozzle-design-workflow.md`
- `sed -n '365,380p' design/forge_design/evaluate/runner_axismach.py`
- `cat case/70.m5_air_wt_demo/_euler_ni_ab/V1b/euler_ni_runs_ab.json | head -80`
