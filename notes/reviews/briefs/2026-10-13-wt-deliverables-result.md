# 諮問ブリーフ: 成果物一式の検証結果の解釈 (plan tooling-nozzle-wt-deliverables、result 段の前)

AGENTS.md「モデル分担とエスカレーション」条件 7 (result 段の解釈を確定する前)。判断役のスイッチが codex。
commit `0785d720` (`feature/nozzle-deliverables`)。plan §6 (事前登録) と §9 末尾 (結果) を読むこと。

## 観測事実 (数値は plan §9 の 2026-10-13「検証の結果」と case/70 README)

- 試験: wt_card 80 ok (+`--slow` 81)・mesh_structure 49・fluent_msh 47・wall_tables 58・deliverables 17・既存 wall_single_bspline 93、FAIL 0。
- V5 の実 run: `case/70.m5_air_wt_demo/run_0006_ns_pass0_iface_60k` (軸対称・等温壁 350 K・SST・`interfaceDiag: 1`・本段 60000 step)。
  `check_convergence --segment` (main): **NOT CONVERGED** (ρ 1.5 桁・ρE 2.0 桁が下降中、他 2.9〜3.1 桁)。壁の総熱量は 12k step (run_0005) 1.51 MW → 60k 0.73 MW。
  同じ 60k の snapshot で、表の q (= 2π `iface_Qf_eff` / A_node) と −`qwall`・`iface_q_compact` の比は r = 25〜150 mm の 9 点で 1.013〜1.026、Σ Q = Σ q·A (相対 0)、
  符号一致、`iface_q_eff` / q = r (9 点とも)。12k の snapshot では比 1.017〜1.049。
  −`qwall` と `iface_q_compact` は 4 桁一致していた。
- 成果物の総合判定は PASS、`thermal_use` は FAIL (未収束)。寸法ゲートは PASS (M5 デモ: 拡大部長 +0.066 mm、他 ≤ 2e-6 mm)。
- Fluent での読み込みは未確認 (ユーザに Fluent が無い、2026-10-13 ユーザ判断で保留し accepted は可)。V8 (skill だけでの再現) は未実施。

## 問い

**この結果から「成果物の書き出し器と熱量の単位は正しい、accepted にしてよい」と解釈してよいか。** 解釈の穴を 1 つに絞って示してほしい。
当方が気にしている候補 (縛られなくてよい): (a) V5 を**未収束の run** の snapshot で判定している — `iface_Qf_eff` は残差と蓄積の補正を含む保存形なので、未収束だと −qwall との比が 1 からずれるはずで、
1.3〜2.6 % が「単位が正しい」根拠として十分か (r 倍のずれの検出には十分だが、5 % の判定が事前登録どおりの意味を持つか)。(b) −qwall と q_compact が 4 桁一致する (独立な 2 形式の照合になっていない可能性)。
(c) 総合 PASS と `thermal_use` FAIL を並べる README の見せ方が、ユーザに「熱量の表を使ってよい」と誤読させないか。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `*.cgns`, `*.stl`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '230,330p' plans/active/tooling-nozzle-wt-deliverables.md` (§6) と `grep -n "検証の結果" -A30 plans/active/tooling-nozzle-wt-deliverables.md`
- `sed -n '1,40p' case/70.m5_air_wt_demo/README.md`
- `grep -n "def thermal_use" -A60 design/forge_design/export/wall_tables.py`
- `grep -n "def cfd_wall_table\|def _wall_block" -A40 design/forge_design/export/wall_tables.py | head -120`
- `sed -n '185,235p' solver_density_cuda/conjugateWall.cpp` と `sed -n '630,660p' solver_density_cuda/cuda_forge/viscousFlux_d.cu`
- `grep -n "_readme" -A40 design/forge_design/export/deliverables.py | head -60`
