# 諮問ブリーフ: スカラー LSQ S2 物理 A/B の結果解釈 — 収束ゲートが gg・lsq 両方で外れた件と FCT smoke の差 (2026-09-26)

plan: `plans/active/gradient-scalar-lsq-unification.md` (§6 S2 表と共通規則、§5.1 #5)。
結果の run 一覧: `case/{48.flat_plate_cooled_m4,40.nozzle_design_tool,39.periodic_hills,16.nozzle_wys}/README.md` の `run_095x_sglsq_*` 行 (h5 は AWS `~/forge-pgrad-new/case/...` のみ)。

## 決めてほしいこと

1. **収束ゲートの不成立 (4 ケース共通)**: 起点 run の床に対する判定 (共通規則 (ii) の `check_floor_ratio`、0476 の `--from-floor run_0476`、0482 の「0482 と同基準で PASS」) が、**gg・lsq の両方で同程度に**外れた。起点は旧バイナリ・旧既定の run なので、現行バイナリの床は起点の床と違う、と読んでよいか。よいなら S2 の収束ゲートをどう読み替えるか (例: 現行バイナリの gg 双子の初段床を基準にする — 下の表では両双子 PASS)。**結果を見てから規則を変えることになる**ので、その扱いも含めて判断してほしい。
2. **FCT smoke の差が §6 の上限を超えた**: 上限は `run_0486` vs `run_0487` の同一設定反復差の 2 倍 (= 差が出ない前提のノイズ規則) だが、SFR 2 では lsq が Y・ξ の再構成を変えるので離散解そのものが変わる。これを F としてどう判定するか。
3. 上の判断を踏まえ、Phase 1 の S2 を「PASS」「条件付き」「不合格」のどれとして codex result 1 回目 (#5r) に回すか。追加の run が要るなら 1 つに絞って (設定・step 数・合格条件)。

## 読んでよいファイル

- plan 全文、4 case README の `run_095x_sglsq_*` 行
- `case/09.Taylor-Green/_g0_lsq_seam/{s2_setup.py,s2_twin_diff.py,s2_eval_case40.py}`、`solver_density_cuda/tools/check_floor_ratio.py`
- `procedures/recommended-settings.md` §1.0a (`slauWallNormalChi`)
- 巨大ファイル・h5 は読まない (AWS にしか無い)

## 観測事実

共通: バイナリ `f99f236d` (sha256 `98703ca1…`)、`FORGE_CUDA_BLOCKSIZE=128`、双子の差は `mesh.scalarGradient` の 1 行のみ (`s2_setup.py` が起点の入力を複製し `restart_field.py` でビット一致の IC を作る)。起点 res の sha256 は plan §6 表と一致 (case/48 `9a9022d1…`、case/40 `15e47af0…`/nozzle `849660af…`、0476 `c1143289…`、0482 IC `909786ef…`、case/39 `01f90b83…`)。

**起点と現行バイナリの設定差 (両双子に共通)**: `slauWallNormalChi` は 2026-09-26 に node+nodeWallDirichlet+SLAU で既定 auto = 実効 1 (commit `58754955`)。起点 (09-12〜09-17 作成) は実効 0 相当。case/40 は廃止キー `mesh.nodeAxisDirichlet: 1` を削除 (2026-08-16 廃止、書くと起動エラー)。case/39 は起点 config が `slauWallNormalChi: 0` を明示しているのでそのまま。起点作成以降、`solver_density_cuda/cuda_forge` と `input` に 96 commit。

### 物理ゲートと lsq−gg の差 (§6 表の上限)

| ケース | 双子 run | 物理ゲート (既存) | lsq − gg | 上限 | 判定 |
| --- | --- | --- | --- | --- | --- |
| case/48 | `run_0952/0953_*_ext` (24000+24000) | Cf/VD-II 0.969/0.986/0.993 (x 0.3/0.6/0.9、両方)、2St/Cf 1.154–1.159、閉合 1.0157/1.0156、NaN 0、系列 STEADY (drift ≤ 0.09 %) | Cf −0.014/−0.009/−0.010 %、q_w −0.011/−0.006/−0.011 %、δ* −0.012/−0.010/−0.006 %、θ −0.016/−0.010/−0.011 % | 各 ≤ 1 % | 差 PASS・ゲート PASS |
| case/48 F1 (記録) | `run_0954/0955_*_f1probe` (level 2 で 1 step、後処理で F1 を再計算) | 壁ノード 1001 点 k = 0 → F1 = 0 | 第一内層 F1 は両方 1.0 (tanh 飽和)、差 0。全域 L∞ 0.24 (境界層外縁) | 記録 | — |
| case/40 | `run_0950+0952` / `0951+0953` (12000+12000) | η_CF 0.97794/0.97795 (∈ 0.978 ± 0.003)、ṁ 1.29947 (∈ 1.29–1.30)、STEADY | Δη_CF +0.000 %、Δṁ 0.000 %、壁温 L∞ 0.62 K | 0.1 % / 0.2 % / 15 K | 差 PASS・ゲート PASS |
| case/16 0476 | `run_0950/0951_sglsq_s2_0476_*` (+24000) | min Y 1.07e-2、max\|ΣY−1\| 8.1e-8/7.5e-8、ξ ∈ [5.9e-4, 0.99985]、floor 補正累積 rel 8e-28、machmax・pmax STEADY | Y0 1.18e-5、Xi 3.12e-4 | 3.7e-5 / 9.4e-4 | 差 PASS・ゲート PASS (from-floor を除く) |
| case/16 0482 | `run_0952/0953_sglsq_s2_0482_*` (48000) | NaN 0、モーメント最小 0 (負値なし)、cond_series ALL STEADY | onset (中心線 g = 1e-3) 23.23 / 23.23 mm、壁 p/p0 L∞ 0.0029 % | 0.3 mm / 0.3 % | 差 PASS (convergence を除く) |
| FCT smoke | `run_0954/0955_sglsq_s2_fct_*` (0476 入力 + 0486 の dual-time 設定、200 step) | NaN 0、max\|ΣY−1\| 7.5e-8/7.7e-8、ξ ∈ [0, 0.99983]、`[passiveFct] active` | ro 3.4e-5、roY1 2.2e-5、roUy 2.3e-2 (相対 1.2e-3、index 111)、Y0 1.2e-5、Xi 5.2e-6 | ro 3e-6、roY 5e-6、roUy 1e-4 | **差が上限超過** |
| case/39 | `run_0950/0951_sglsq_s2_*` (200000、設定は起点と同一: `slauWallNormalChi: 0` 明示) | 全量 STEADY (`--drift 0.002 --osc 0.005 --tail 0.4`)、r_gradk 0.9844/0.9845・r_gradw 1.0000 (∈ [0.9, 1.1])、NaN なし | Cf_x05/x2/x6 −0.001/−0.01/−0.08 %、x_r/h 5.0339/5.0343、dF1_inf 0.03342/0.03350 (**lsq で STEADY**、gg も今回 STEADY) | 記録のみ | ゲート PASS |

参考 (旧バイナリの起点 vs 現行 gg): 0476 で Y0 5.5e-5・Xi 4.8e-4、0482 で壁 p/p0 L∞ 0.16 %、case/40 で η_CF 0.9754 → 0.9779。いずれも lsq−gg より大きい。

### 収束ゲート (両双子で同程度に不成立)

| ケース | 規則 | gg | lsq | 現行 gg 初段の床を基準にした場合 |
| --- | --- | --- | --- | --- |
| case/48 延長 | 末尾平均 ≤ 1.5× 起点床 | rms_ro 1.85×、roUy 9.1×、roe 1.9× | rms_ro 2.54×、roUy 10.4×、roe 2.7× | gg 1.00–1.03× / lsq 0.99–1.43× (roe 1.425、roUy 1.14) で両方 PASS |
| case/40 延長 | 同上 | rms_roOmega 55.9×、ro 1.59×、roe 1.70× | roOmega 56.7×、ro 1.60× | 両方 PASS (0.98–1.02×) |
| 0476 | `--from-floor run_0476` | 全列 21–64× (peak 170×) | 同じ | (未計算) |
| **case/39** (設定差なし) | 末尾平均 ≤ 1.5× 起点床 | **PASS** (0.99–1.00×) | **PASS** (0.99–1.00×、rms_roK の peak 2.4× @584 → step 2704 で再進入) | — |
| 0482 | `check_convergence` PASS (0482 は 3.8 桁で PASS) | plateau 2.5 桁 (rms_ro 9.9e-8、0482 は 4.7e-9) | plateau 2.5 桁 (8.6e-8) | — |

- `check_convergence`: 全ケース DIVERGED なし、RISING なし (case/48 lsq 延長の rms_roK のみ `falling`)。
- case/48 だけ lsq の床が gg より高い (延長末尾 rms_ro 2.98e-7 vs 2.00e-7、roe 3.7e-1 vs 2.4e-1、≈1.4 倍)。case/40・0476・0482 では床は同程度。

## 期待値と出典

- §6 S2 共通規則 (ii) 末尾平均 ≤ 起点の 1.5 倍、「未達なら同一設定で延長」。0476 は `--from-floor run_0476`、0482 は「0482 と同基準で PASS」。
- FCT smoke の上限は「反復差 0486 vs 0487 の 2 倍」(§6 表)。
- 「差が上限を超えたら自動不合格にも自動合格にもせず格子対診断 → F 判断」(§6 共通規則、codex M5)。

## 仮説

- H1: 収束ゲートの不成立は lsq と無関係で、**設定差の無い case/39 だけが両双子とも起点床に戻った**ことと整合する。起点作成後のバイナリ変更 (特に `slauWallNormalChi` の既定化、case/40 は `nodeAxisDirichlet` の廃止) で床が移動した。根拠は gg 双子も同じだけ外れていること。反証するなら、gg 双子に `slauWallNormalChi: 0` を明示して同じ延長を回し、起点床に戻るかを見る (case/48 は 2 分で済む)。
- H2: case/48 の lsq 床が 1.4 倍高いのは、lsq が第一内層以外の CDkω・F1 (全域 L∞ 0.24 の差) を通じて擬似時間の減衰を変えたため。物理量は 0.016 % 以内で STEADY なので結果には効かない、という見立て (未検証)。
- H3: FCT smoke の差は離散化の変更によるもので、ノイズ規則を当てるのは不適切。200 step の dual-time は起点 (gg の定常解) から lsq の定常解へ向かう過渡の途中なので、差は時間とともに 0476 双子の定常差 (Y0 1.2e-5、Xi 3.1e-4) の方へ向かう、という見立て。代替の規則案: 0476 双子の定常差を上限にする、または FCT の健全性 (有界性・保存・ΣY) だけをゲートにして差は記録にする。
