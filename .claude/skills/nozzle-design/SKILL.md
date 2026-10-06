---
name: nozzle-design
description: 軸対称風洞ノズル (wind_tunnel_axisym_axismach) の設計・計算の作業手順。(1) 軸長・出口径・入口径・入口条件を指定して非粘性 (MOC) ノズルを新規設計する、(2) 出口・入口・条件を固定して軸長を MOC でパラスタし、有望形状を Euler → NS (境界層補正) → 形状調整 → 凝縮計算する、(3) 作成済みのノズルで Euler/NS を回す — のどれかをユーザが頼んだとき (「ノズルを設計したい」「軸長を振りたい」「このノズルで Euler/NS を回して」「凝縮を見たい」) に使う。
---

# /nozzle-design — 軸対称風洞ノズルの設計・計算

**棲み分け**: 要件の確定 (目的・Pt/Tt/組成・M_design・寸法・凝縮の扱い・コスト承認) と problem YAML・方針メモの作成は
skill [`design-intake`](../design-intake/SKILL.md) が担う (ここでは繰り返さない)。本 skill はその YAML ができた後の
**axismach ノズルの作り方・回し方**だけを扱う。YAML が無い新規の依頼は、まず `design-intake` を通す。

**手順の正本は [`procedures/nozzle-design-workflow.md`](../../../procedures/nozzle-design-workflow.md)**。
ここには AI が判断を誤りやすい点だけを置く。作業前に正本の該当節 (§1 新規設計 / §2 パラスタ〜凝縮 / §3 既存形状) を読む。

## 最初に決めること

1. **どのユースケースか** を確定する (§1 / §2 / §3)。曖昧ならユーザに 1 問だけ聞く。
2. **新規設計 (§1・§2) で problem YAML が未確定なら `/design-intake` に渡す** (要件インタビューはそちらの担当)。
3. **最新の実例を手本にする**: `case/45.isobutane_m6_d155` (M6、2026-10)。生産レシピは
   [`plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md`](../../../plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md) §4・§5.1 #11〜#11f と §9 末尾。
   古い run の YAML をそのまま複製しない (試作の `M_design` 較正値・廃止キーが残っている)。

## 判断の分かれ目

- **既存ノズルで条件 (Pt・Tt・H2O 分率・背圧) だけ変えたい** → 正本 §3a: `solver_density_cuda/tools/rerun_conditions.py REF NEW ...` (形状固定、forge は起動しない)。Pt を変えるなら `--Ps`/`--keep-Ps` 必須、回し方は `RERUN_CONDITIONS.json` の `recommended_stages` (Pt 変更は full + 本段 cfl 1 + `--scale-ic pt` が暫定推奨)。
  **problem YAML の spec/gas を変えて runner/deltastar_loop に通すと壁が作り直される** (それは §1・§2 の設計作業)。
  組成の lump を変えると化学種の定義が変わり restart_field が拒否する → convert_species_field か等エントロピー初期場。

- **探索段階** (軸長・R・M_knot を振る) は Hall 初期線 + 補間壁 (YAML の既定) と粗い格子で十分。
  **形状を詰める段階**で `wall_repr: joint`・`initial_line: cfd`・`Md_moc_offset` を入れる (ユーザ決定 2026-10-05)。
- 壁の評価観点 (ユーザの方針): **壁がなめらか (r″ が連続) / 試験部に圧力波が無い / オーバーシュートが小さい** (凝縮を強めるので)。
  軸 M は形状の生成器にすぎない。Pa 単位の波の追い込みはしない (Euler の数値床と同じ桁)。
- 試験部・膨張部 (MOC の領域) に Euler のずれを返さない。CFD に合わせるのは遷音速のスタート部 (初期線) と出口の 1 係数だけ。
- **出口較正 (`Md_moc_offset`) は生産 NS と同じ格子パラメータの Euler で決める** (粗い Euler 格子だと NS で出口 M が約 −0.0008 ずれた)。
- NS の格子は**壁解像 PASS** (y1+ ≤ 1、`check_wall_resolution.py --over-frac 5`) と**格子ゲート** (δ_E(x_F) の細分前後差 ≤ 1 %) を通したものを使う。
  細分格子は本段 CFL 5 で発散する → 段階起動 + cfl 1・60000 step。

## 落とし穴 (実際に踏んだもの)

- `delta_r_next.csv` の δ は **`delta_E` 列 (未緩和)** を使う。第 2 列 `delta_target` は ω=0.5 緩和後 (取り違えて k_f・r_t を誤って解いた)。
- 出口径は spec に書けない。`solve_rt` で r_t を解いて `spec.r_throat` を手で書き換える。r_t を変えても `r_inlet`・`L_total` は r_t 単位のまま。
- `design_chain` は fold・特性線トポロジ・壁マージンを検査しない → パラスタでは `geometry/moc_diagnostics.py` を併用する。
- Euler の `prepare` と NS の `prepare_ns` は同じ mesh キーを読む (2026-10-06 に `mesh_params()` で統一。以前の Euler は一部のキーを黙って無視していた)。
- 同じ格子の続きは `restart_field.py`、格子が違えば `interp_field.py`。凝縮 ON への引き継ぎは `convert_species_field.py --mode conserve` (restart_field は拒否する)。
- case スクリプトに `/home/sano/work/forge/...` の絶対パスが書かれているものがある (worktree で走らせると別ツリーを読む)。使う前に確認する。
- AWS の共有バイナリは他セッションが作り直すことがある。自分の worktree でビルドし、`RUN_PROVENANCE.txt` の sha256 を記録する。
- 判定は判定ツールの VERDICT で行う。残差 plateau (`NOT CONVERGED`) はこの系の常態 — そのまま書き、**報告する量の時系列**に
  `check_quasisteady.py --series-csv` をかける (`nozzle_report` の末尾幅は VERDICT の代わりにならない)。
- 判定基準は回す前に plan に書く (事前登録)。外れたら閾値を動かさず、延長 1 回 → なお未達なら上位の判断役に諮る (AGENTS.md「モデル分担とエスカレーション」)。

## 未実装 (ユーザに聞かれたら正直に言う)

汎用の MOC パラスタ CLI・壁 CSV 出力/作図 CLI・k_f の YAML キー・凝縮引き継ぎの自動化・外部壁 CSV の取り込みは無い (case スクリプト頼み)。
統一は [`plans/active/tooling-design-problem-campaign-recipe.md`](../../../plans/active/tooling-design-problem-campaign-recipe.md) の残作業。

## 終わったら

`nozzle_report` → pptx ([`procedures/nozzle-design-outputs.md`](../../../procedures/nozzle-design-outputs.md))、case README の run 一覧、plan §9 に記録。
