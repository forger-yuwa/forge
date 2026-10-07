# 諮問: 細分メッシュ生産 pass (#11 ①) の発散 2 回 — 次の一手はメッシュの第 1 セル分布の修正か、起動手順か

日付 2026-10-05。諮問先: codex (diagnose) — `~/.config/forge/diagnose-backend` = codex のため。
エスカレーション条件 2 (手順書 `procedures/divergence-and-startup.md` で 2 回対処して発散) と 4 (plan に無い修正の前)。
plan: `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` §5.1 #11 (事前登録)、§4.9〜4.11、§9 末尾。
作業ツリー: `/home/sano/work/forge-integ-1005` (ブランチ `feature/nozzle-wall-fit-and-pipeline`, commit 754f02d8)。AWS 側も同 commit。
バイナリ: main 統合ビルド (`~/forge-integ/solver_density_cuda/build/forge`, exe_size 40287504、run_0089〜0094 と同一)。`FORGE_CUDA_BLOCKSIZE=128`。

## 目的 (#11 の事前登録)

run_0094 (現メッシュ ni 1250 × nj 65、`wall_first_frac` 4.5e-5) は壁解像 FAIL (y1+ > 1 が 29.4 %、最大 13.39 @ x/r_t −0.32)。
#11 ① は「スロート近傍の第 1 セル細分 (y1+ ≤ 1 目標、壁法線層 AR ≤ 5000) で run_0092 の壁を NS、IC = run_0092 (interp_field)、12000 step」。

## 観測事実

### run_0094 の y1+ の x 分布 (正式ツール `check_wall_resolution.py` と同じ式、step 6000; x は r_t 単位)

| x 帯 | y1/r_t 中央 | y1+ 中央 | y1+ 最大 | y1+>1 |
| --- | --- | --- | --- | --- |
| [−40,−11) | 2.9e-4 | 2.45 | 4.60 | 100 % |
| [−11,−6) | 2.6e-4 | 3.10 | 3.22 | 100 % |
| [−6,−2) | 1.1e-4 | 3.23 | 6.09 | 100 % |
| [−2,0) | 4.9e-5 | 12.06 | 13.39 | 100 % |
| [0,2) | 4.9e-5 | 7.68 | 12.26 | 100 % |
| [2,6) | 8.4e-5 | 3.17 | 5.39 | 100 % |
| [6,12) | 1.5e-4 | 1.45 | 2.01 | 100 % |
| [12,20) | 2.2e-4 | 0.97 | 1.18 | 42 % |
| [20,40) | 3.3e-4 | 0.60 | 0.82 | 0 % |
| [40,96) | 4.1〜4.5e-4 | 0.28〜0.39 | 0.47 | 0 % |

### 作った細分メッシュ (`case/45.isobutane_m6_d155/problem_d155_ns_finemesh_pin.yaml`)

ni 2000 × nj 97、`wall_first_frac` 1.3e-5 (**全域**、旧 4.5e-5 の 1/3.5 — x > 20 の不要な細分を含む)、`wall_first_frac_throat` 4.5e-6、
上流は x ≤ −9 で 1.3e-5 → x ≥ −4 で throat 値、下流は x ≤ 1 で throat 値 → x ≥ 17 で 1.3e-5 (smoothstep、`design/forge_design/meshing/mesh2d.py` L118〜130)。
`throat_refine` 4、`throat_width` 3。品質: AR max 4260・p99 2374・mean 146 (旧 1210・562・51)、skew max 0.44 → `--ar-max 5000` PASS。

### 試行 1: `case/45.isobutane_m6_d155/run_0095_ns_finemesh_pass` (IC 直行)

IC = run_0092 最終場を `interp_field.py` (最近傍) → 本段そのまま (2 次 SLAU `convMethod 1`・`limiter 2`、block-DPLUR、cfl_pseudo 5、implicitRelax 0.7、nStepInner 5、SST 低 Re)。
**step 20 で ro NaN** (detectNaN)。非有限 375 節点、x/r_t 59〜70・r/r_t 9.50〜9.84 (壁際)、T は上限 6000 K に張り付き (IC 1460 K)。
序盤: rms_ro 2.2e-3 → step 15 で 3.0e-4 (下降) → 再上昇、worst は rms_roe (step 16 以降 +)。

### 試行 2: `case/45.isobutane_m6_d155/run_0098_ns_finemesh_pass_staged` (段階起動)

同じ IC に `run_staged_ns(stages="full")`: soft (1 次 cfl 0.5、nStepInner 10、3000 step) → mid (1 次 cfl 1、3000) → 本段 (2 次 cfl 5)。
soft・mid は完走 (各 ~13 s)。**本段 step 35 で ro NaN**。非有限 615 節点、x/r_t 78.6〜93.8・r/r_t 9.97〜10.09 (出口側の壁際)、T → 6000 K (段開始値 1480 K)、
ω/ω_段開始 が 120〜155 倍 (x 86〜88、壁際)。本段の序盤: rms_ro 2.4e-5 から単調増加 (step 10 +0.1 桁、20 +0.2、30 +0.3)、worst は step 20 から rms_roOmega (+1.0 → +2.9 @ 30)。
段の出力ファイル (`_soft_res_*` 等) は残っていない。`stage_manifest.json` は無い。

### 参照 (同じ設定で通った run)

run_0089 (現メッシュ、等エントロピー IC から stages full) → run_0090/0092 (IC = 前 pass の場を interp_field、本段直行 cfl 5) → run_0094 (延長)。いずれも NaN なし。
case/45 の設定は run_0092 と同一 (メッシュ以外)。

## 当方の仮説

H1. **下流の不要に薄い壁セル** (x > 20 で第 1 セルが旧の 1/3.5、y1+ ≈ 0.1) と 2 次・cfl 5 の組合せが壁際の ω (壁値 ∝ ν/y1²) を硬くし、出口側の低密度 BL で T が暴走する。
2 回とも発散が下流 x 55〜94 の壁際 (細分の目的のスロート近傍ではない) であることと整合。
H2. 1 次段 (各 3000 step) が細分メッシュの過渡に短すぎ、本段への切替 (1 次 → 2 次、cfl 1 → 5) が跳びすぎ。

## 提案 (次の一手)

A. メッシュを直す: 下流は旧と同じ 4.5e-5 (y1+ ≤ 0.82 の実績)、上流 (x ≤ −9) は 1.3e-5 (旧 y1+ 2.5〜3.2 → ~0.8)、スロート 4.5e-6、下流ブレンド x 1 → 22。
mesh2d に上流の遠方値 `wall_first_frac_up` を 1 キー追加 (既定 = `wall_first_frac` でビット不変)。起動は試行 2 と同じ stages full。
B. メッシュはそのまま、本段を `stages="ramp"` (cfl ランプ) か cfl_main 2 に下げる。
C. A + B。

## 問い

Q1. H1/H2 のどちらが主因と読むか。データから足りない切り分け (例: 段階起動後の場で x>20 の y1+ と ω 壁値、`res_nan_35.h5` で最初に壊れる変数) は何か。
Q2. A を採るとき、事前登録 (#11 ③ の「細分前後で δ_E(x_F) の差 ≤ 1 %」など) との関係で注意すべきこと (下流を旧メッシュに戻すと δ_E(x_F) の細分感度はスロート近傍の細分分だけになる — それで #11 の目的 [壁解像 PASS と δ_E の格子依存の確認] は満たせるか)。
Q3. 次の試行 (3 回目) の合格/停止条件をどう事前登録するか。

## 読んでよいもの

上記 plan、`case/45.isobutane_m6_d155/README.md`、`case/45.isobutane_m6_d155/{run_finemesh_pin.sh,problem_d155_ns_finemesh_pin.yaml,problem_d155_ns_c2pin_final.yaml,prep_c2pin.py}`、
`design/forge_design/meshing/mesh2d.py`、`design/forge_design/evaluate/runner_axismach.py` (`run_staged_ns`、`prepare_ns`)、`procedures/divergence-and-startup.md`。
run ディレクトリは AWS 上にありローカルには無い (上の数値が全て)。編集は禁止。
