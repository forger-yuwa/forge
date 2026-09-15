---
name: design-intake
description: 設計キャンペーン (forge_design の MOO/評価) を立ち上げる前の要件インタビュー。目的関数・設計変数と範囲・CFD メニュー (ガス/凝縮/高度/乱流) をユーザと対話で確定し、検証済みの problem YAML と方針メモを残してからキャンペーン起動を提案する。ユーザが「設計したい」「最適化を回したい」「新しいノズル/風洞の検討」を言い出したときに使う。
---

# /design-intake — 設計キャンペーン立ち上げインタビュー

目的: 曖昧な設計要望を、**検証済みの `problem_*.yaml` (契約書) + 方針メモ**に変換する。
対話で決めるのは YAML の中身だけ。決めた後のキャンペーン実行は決定的 (seed 固定・
ledger) であり、ループの中に対話を持ち込まない。

## 原則

- **まず [`design/CAPABILITIES.md`](../../design/CAPABILITIES.md) を読む** (対応状況の正本)。
  記憶で「できる/できない」を答えない。表と食い違う発言をしない。
- 1 ターンに質問は 2–4 個まで。ユーザの答えから決められることは聞き直さない。
- 各決定には**根拠を聞く** (特に dv の min/max)。「なんとなく」の範囲は INFEASIBLE
  だらけの DOE を生む。
- 未対応メニューが選ばれたら**その場で実装を始めない**。§6 の分岐に従う。
- コスト見積もりの承認を得るまでキャンペーンを起動しない。

## インタビュー手順

### 1. 目的の確認

- 何を最大/最小化したいか (η_CF? 推力? 全長? 軸 M 一様性 ε_M? 出口一様性?)
- トレードオフが本当に 2 軸あるか。1 軸なら MOO でなく単目的/直接設計
  (axismach 系) を提案する。
- 制約と目的を区別する (「全長 ≤ X は制約か、短いほど良いのか」)。

### 2. 問題タイプの選定

`probdef.KNOWN_TYPES` から選ぶ。迷ったら:
- スラスタ/ベルノズル形状の MOO → `thruster_bell` (+ `opt/driver.py`)
- 風洞ノズル (目標 M の軸分布を直接設計) → `wind_tunnel_axisym_axismach`
- 既存タイプに合わない → 新タイプ = 実装タスク。§6 へ。

### 3. 仕様固定量 (spec) とガスモデル

- Pt / Tt / 背圧 (p_ambient) / r_throat または M_design / (風洞なら出口径)。
- ガス: 常温空気なら `cpg`。**Tt > 600 K・燃焼ガス・多成分は `semiperfect`** を推奨。
  - semiperfect の罠: **`evaluate.thermo_href_temp: 298.15` 必須** (絶対基準 h だと発散)。
  - 組成はモル分率でもらったら質量分率へ正規化して記録する
    (例: `case/44.vitiated_air_wt/problem_va3_M4.19_Lc8_dry.yaml` のヘッダコメント形式)。
  - 条件 (Pt/Tt/組成) を変えたら **r_inlet を付け替える** (固定値の使い回しは罠)。
- 凝縮の可能性 (H₂O/炭化水素を含み膨張で低温になる) があるかをここで確認。
  あるなら**扱い方をユーザに問う** (既定は a):
  - **a. 後段確認 (既定・推奨)**: キャンペーン本体は dry で回し、パレート勝者に
    対してのみ凝縮評価を行う (dry run の後処理 `Tsat_post` で飽和線接近を診断 →
    必要なら凝縮 ON CFD で再評価)。ループが軽いまま。
  - **b. 制約化**: 「飽和線を越えない (S<1、余裕温度 ΔT を指定)」を評価ゲートに
    追加する。dry CFD + `Tsat_post` 後処理で判定できるため 1 評価のコストは
    ほぼ増えない。凝縮回避が要求仕様のときはこちら。
  - **c. ループ内凝縮 CFD**: 全評価点で凝縮 ON (`tp_species: split_h2o` [= `{mode: lumped, lumps: {MIXDRY: {from: composition, exclude: [H2O]}}, keep: [H2O]}`; `full` で全種独立も可] +
    `condensation`, condEquilibrium: 2 推奨, axismach のみ)。凝縮量・潜熱の影響
    そのものが目的量に効く場合に限る。コスト増を明示して承認を取る。
  - 注意: `condTsat` 診断は凝縮 ON run のみ出力。dry run は後処理 `Tsat_post` を使う。

### 4. 設計変数 (dv) と範囲

- 各 dv の min/max の**物理的根拠**を必ず聞く (製造制約? 既往設計? 理論限界?)。
- 記入後、**範囲の四隅+中心で `Campaign.feasible()` 相当の幾何成立を確認**
  (fold・θ 矛盾を DOE 前に検出)。M6 級では単一多項式壁が fold する既往例あり
  → knot 則を使う。
- dv が 4 個を超えるなら DOE 点数の増加 (≳10×d) とコストを警告する。

### 5. CFD メニュー (evaluate)

CAPABILITIES §2–4 を見ながら選ぶ。主な分岐:
- 乱流: Euler で足りるか (設計ループの大半は Euler+δ* 補正で十分の実績)、
  NS/SST が要るか。
- 凝縮: §3 で確認済みなら `evaluate.condensation` を設定 (axismach のみ対応)。
- 高度/プリューム評価: 現状 🧩 (資産あり・配管未)。必要なら §6 へ。
  「全評価点で必要か、代表点だけ後段で高度評価か」を必ず聞く (コストが桁で違う)。
- 化学反応: **何を知りたいかで段を選ぶ**。
  - 性能の不確かさ幅 → 凍結流 (✅) + CEA 平衡ブラケット (🧩) で挟む。まずこれ。
  - 着火・緩和距離そのもの → 有限速度 (📋・大型)。§6 へ。
- メッシュ: 設計ループは mesh2d (既定)。固定形状の評価だけなら Gmsh/Fluent 経路可
  (Fluent はレシピ 6 項目を solverConfig に必ず織り込む)。

### 6. 未対応メニューの分岐 (CAPABILITIES の状態記号で機械的に)

- ✅ → そのまま使う。
- 🔧 → 条件・罠を明示して使う (該当 case/メモを参照)。
- 🧩 (資産あり配管未) → 移植 plan を `plans/active/` に起こす見積もりを提示し、
  「plan 消化を待つ」か「縮小メニューで先行」かをユーザに選ばせる。
- 📋 (未実装) → 同上。ただし規模感を正直に言う (有限速度反応・CHT は数週間級)。
- どの場合も**インタビュー自体は完走させる** (YAML は書ける範囲で書き、未対応項目は
  TODO コメントで残す)。

### 7. ゲート・許容値

- 収束: `--drop` (既定 3、warm 継続系は 2 の実績)。
- 準定常: 目的量 (軸 M なら tol_M。既定 5e-3 = 0.5%M_d の 1/4)。
- 物理ゲート: η 範囲・Cd 範囲 (既定 0.5<η<1.02, 0.94<Cd<1.005。Cd>1.005 は SUSPECT)。
- 既定から変える場合は理由を YAML コメントに残す。

### 8. コスト見積もりと承認

- 類似 case の実績から 1 評価の壁時計を引く (bell Euler+SST ≈ 34 s/評価 @case/40)。
- 総評価数 = n_doe + n_iter×batch。総時間・ディスク (res_*.h5) を提示し、
  **明示の承認を得てから**起動する。プリューム/NS 系は 1 評価が桁で重くなることを明記。

### 9. 成果物 (この順で作る)

1. `case/NN.<slug>/problem_<tag>.yaml` — 決定事項をヘッダコメントに残す
   (条件の由来・導出値。case/44 の形式に倣う)。
2. **検証**: `cd design && .venv-opt/bin/python -c "from forge_design.probdef import load_problem; load_problem('../case/NN.<slug>/problem_<tag>.yaml')"` が通ること (forge_design は design/ を cwd に import する) + 範囲四隅の feasible 確認。
3. **1 点 smoke 評価** (段階起動込み) → VERDICT を確認してから本キャンペーン提案。
4. 方針メモ: 設計判断 (なぜこの範囲・このメニューか) を `plans/active/design-<slug>.md`
   に記録 (テンプレ `plans/_template.md`)。インタビューの Q&A 要約を含める。
5. case README に run 一覧表を用意 (AGENTS 必須ルール)。
6. 新しい知見 (未対応→実装した等) が出たら **CAPABILITIES.md を同期**する。

## キャンペーン起動 (承認後)

```bash
# thruster_bell MOO の例
design/.venv-opt/bin/python -m forge_design.opt.driver \
    case/NN.<slug>/problem_<tag>.yaml case/NN.<slug>/run_NNNN_moo_<tag> \
    --seed-res <収束済み res_*.h5> --n-doe 24 --n-iter 13 --batch 2
# 長時間になるので tmux / nohup を推奨。終了後は ledger.jsonl の
# status/fail_class 分布を集計して報告する。
```

run 番号は case README の一覧と衝突しない連番を振ること。
