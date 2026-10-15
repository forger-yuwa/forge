# 要求カードの格子を場所ごとの第一層と流れ方向の密度で作る

## メタ

- **area**: `tooling` (`design/forge_design/intake/wt_card.py` の `estimate_mesh`・雛形 `intake/recipes/axismach_*`、`meshing/mesh2d.py` は使うだけ)
- **status**: `draft` (2026-10-15 起票)
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「要求カード」「格子」の節
  - skill `nozzle-design`・`design-intake`
- **related_plans**:
  - [`tooling-generator-defaults-and-launch-gate.md`](tooling-generator-defaults-and-launch-gate.md) — 段取り 2 の「生成器の既定 = 最新の生産レシピ」と `mesh2d` の壁法線の近壁層の既定化 (§6.4)。この plan はカードの格子の作り方だけを先に揃える
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) §5.1 #36r — M6 の H の上流の第一層 (NS の壁量から表を作り直した例)
  - [`tooling-nozzle-mknot-fast.md`](tooling-nozzle-mknot-fast.md) — この枝の元 (M_knot の新しい規則)
- **created**: `2026-10-15`
- **owner**: `M6ノズルの設計` セッション (作業ツリー `forge-cardmesh` / `feature/card-local-mesh`、feature/mknot-fast 7713f84c から)

## 1. 目的

要求カードの雛形 (`axismach_explore_v1`・`axismach_final_v1`) の NS 格子は、第一層を 2 値 (スロート / それ以外) + ブレンドで持ち、流れ方向はスロートの細分だけ、壁に直交する近壁層も無い。
M6 の生産の格子 (case/45 の H: `wall_first_frac_table`・`x_density_table`・`wall_normal_layer`) より古い作り方で、M5.5 (case/71) では最も厳しい場所の薄い第一層が拡大部 8 m 全体に掛かり、
壁際の AR ≤ 5000 のために ni 10,177・約 99 万点になった (ユーザ 2026-10-15「あまりにメッシュ多すぎる、古いやり方でメッシュ作ってない？」「ほかのケースでちゃんと流れ方向に分布つけてたし、壁面に直交するようにメッシュ作ってた」)。
**目標**: カードから作る格子を、場所ごとの第一層の表・流れ方向の密度の表・壁に直交する近壁層で自動で作る。M5.5 は ユーザ決定「a」(2026-10-15) で、見積もりの式の偏りを M6 の実績で補正した約 45 万点の格子で Euler → NS まで回す。

## 2. スコープ

- **やる**: `estimate_mesh` が `wall_first_frac_table`・`x_density_table`・`ni` を作る (作り方 `local_v1`)、`wall_normal_layer` を雛形の既定に、見積もりの補正の係数、記録、M5 デモ (case/70) と M5.5 での確認。
- **やらない**: 見積もりの式そのもの (Eckert の参照温度 + 平板の摩擦係数) の置き換え、`nj` と径方向の配点 (core-grid の範囲)、Euler の格子 (`mesh_euler`、全域 0.005 のまま)。

## 3. 関連 docs と前提 (観測事実)

- **見積もり** (2026-10-15、AWS A `~/m55_mesh_est/`、手元 scratchpad `m55_mesh_est/`; コードは変えずに `estimate_mesh` (7713f84c) の式で計算): M5.5 (r_t 0.0999088 m・M_knot 2.80・nj 97)、y₁⁺ 目標 0.8、AR ≤ 0.97 × 5000 = 4850。

  | 作り方 | ni | 節点 |
  | --- | --- | --- |
  | (A) 今: 第一層 2 値 + ブレンド | 10,177 | 987,169 |
  | (B) 第一層を表、流れ方向は今の分布 | 8,749 | 848,653 |
  | (C) 第一層を表、流れ方向の密度も表 (dx ≤ 4850·y₁·cosθ_w、上限は雛形の ni 2000 相当の間隔) | 3,163 | 306,811 |

  (C) の格子を `mesh2d.generate_axisym_mesh` で実際に作り、AR 最大 4850.0・隣り合う dx の比の最大 1.031 (径方向 1.134)。
  前提を変えると: y₁⁺ 目標 1.0 で 271,503 点、第一層を全域 ×0.5 で 490,335 点。入口の端は y₁,req → 0 (y₁ ∝ s^0.1) なので、断面の 1 点目 (s = 0.011 r_t) より手前はその値で置き換える (生成器の y₁⁺ の検査も入口の節点を除いている)。
- **見積もりの式の偏り** (M6 の H の表は NS の壁量から作ったもの): H の表 ÷ 見積もりの要求 = 配管 2.2〜3.6、縮流部 0.61〜2.2、スロート付近 0.51〜0.61、拡大部 0.50〜0.69。見積もりはスロート・拡大部で第一層を約 2 倍厚く見積もる (壁のせん断を低く見積もる)。配管では逆に薄く見積もる。
  H の表は断熱の NS の y₁⁺ に冷却の倍率 8.5〜9.4 を掛けて y₁⁺_cold ≤ 0.8 にした 300 K 用 (YAML のコメント)。M5.5 は同じ燃焼ガスの組成・壁 350 K。
- M6 の H は今の配点で 46 % の列が AR 4000 超 (既に限界近く)。(C) にしても ni 4,719 → 4,294 (1 割減)。
- `mesh2d` は `wall_first_frac_table` と `wall_first_frac_throat`・`wall_first_up_*` の併用を拒む (`mesh2d.py:100`)。
- `wall_normal_layer` は PR #8 で main に入った (c096d66c)。M6 の H は `[0.02, 0.8]`。カードの雛形には無い。

## 4. 設計方針 (2026-10-15 起案、上位に諮る前)

### 4.1 作り方 `local_v1` (雛形の `mesh.method: local_v1` で選ぶ; 新しい雛形の既定)

1. 設計壁に沿って、今の `estimate_mesh` と同じ式で要求の第一層 y₁,req(s) を出す (y₁⁺ 目標 `Y1PLUS_TARGET` 0.8、入口の端は s = 0.011 r_t の値で頭打ち)。
2. **補正**: y₁,req に区間ごとの係数 c(x) を掛ける。係数は M6 の H の表 ÷ 見積もりから、保守側に取る (表が見積もりより厚い区間は 1.0 = 補正しない): 配管・縮流部 1.0、スロート付近と拡大部 0.5 (H の比の最小 0.50 に合わせる)。区間の境 (縮流部の終わり x = −2 r_t 付近) は滑らかにつなぐ。係数は雛形の `mesh.first_layer_calibration` に表で持ち、出所 (M6 の H の表と見積もりの比、2026-10-15) を記録する。
3. 第一層の表 `wall_first_frac_table`: 0.125 r_t (x < 10) / 1 r_t (x ≥ 10) の点で、前後の区間の c·y₁,req の最小 (下包絡) を局所半径で割る (log 線形補間で要求を超えない)。
4. 流れ方向の密度の表 `x_density_table`: 目標の間隔 h(x) = min(4850·y₁,表·cosθ_w, h_cap(x)) の逆数の上包絡。h_cap は雛形の分布で ni 2000 にしたときの各位置の間隔 (主流の解像の上限)。AR ≤ 4850 と dx ≤ h_cap を満たす最小の ni を二分法で求め、10 の倍数に切り上げる。`throat_refine` は 1。
5. `wall_normal_layer: [0.02, 0.8]` (M6 の H と同じ) を雛形の既定にする。
6. 節点数は `check_node_budget` で判定 (上限を超えたらユーザの了承と `max_nodes`・`max_nodes_reason`; 今のまま)。
7. 記録: card_log の `mesh_estimate` に y₁,req・係数・表・ni・節点数・AR の最大・dx の比の最大。

### 4.2 後で直す仕組み

NS の後に `check_wall_resolution.py` (壁ごと・面積の重み・停止の step) で y₁⁺ を確かめ、超えていれば NS の壁量から表を作り直す (M6 の #36r の台本 `h370_wall_envelope.py` の `design` の手順)。カードの流れでは「NS の後の壁解像の確認 → 必要なら表を直して NS をもう 1 回」を標準の手順として skill に書く。

## 5. 実装ステップ

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~plan 段レビュー~~ | 判断: 2026-10-15 NO-GO (C0/M7/m1) を全件採用 (§6.1)。§4・§6 は書き直して再レビュー | F |
| 3a | **M5.5 の手組みの格子 (先行)** | ユーザ決定「a」(2026-10-15): case/71 の問題 YAML の mesh 節を、正規化した要求 f_req = c(x)·y₁,req/r_w の下包絡の表 (補間の区間の中も検査・物理壁で再検査)・流れ方向の密度の表・`throat_refine 1`・`wall_normal_layer [0.02, 0.8]`・`msh_digits 17`・`max_nodes` (ユーザ了承「a」、約 45 万) で手で組む。係数 c: 配管 1.0、縮流部 0.6、スロート・拡大部 0.5 (M6 の H の表 ÷ 見積もりの比の最小、初期格子の経験の係数)。品質: `check_mesh_quality --ar-max 5000`・`mesh_ar_relax_check`・変換後の第一層厚の照合。生成器の既定は変えない。**実装 (2026-10-15、implementer)**: `case/71.m55_combustion_wt/problem_m55_ns.yaml` (`_mesh/hand_mesh.py build`; 生成は `card_m55.yaml` に `mesh.ni: 3000` を足した `problem_m55_card.yaml`、M_knot 自動 2.80)。ni 4,800 × nj 97 = 465,600 (NS・Euler)。係数の拡大部は比の最小 0.49874 < 0.5 なので 0.49 (0.5 なら ni 4,780)。CFD 前の物理壁の格子 (forge 不起動、A `~/m55_mesh_hand/check/`): 品質 PASS (AR 最大 4849.4)・AR > 1000 の 47,607 セルでスキュー最大 0.0015 (違反 0; 出口面の非平面 0.62 mm は記録だけ)・第一層の実寸 / (c·y₁,req) 最大 0.999999 (超過 0)・見積もりの y₁⁺ 最大 0.797・dx の比の最大 1.033。Euler PASS (AR 最大 6.1)。再現の入力と記録は `case/71.m55_combustion_wt/_mesh/` **判断 (2026-10-15、主セッション)**: (1) 拡大部の係数 0.49 (M6 の比の最小 0.49874 を 2 桁に切り下げ) — 採用。(2) 正規化の半径を max(設計壁, 物理壁 (CONTUR の初期壁)) — 採用 (NS の pass で壁が動けば半径の比だけ第一層がずれる: NS の後の V4 で見る)。(3) FP64 のビルドで回す (第一層/半径が 3.7e−7 まで下がる): AWS A の `~/forge-integ-dflow2-fp64` (2fd36b0c、solver の源は main と同じ) を使う。(4) 探索の雛形 (`axismach_explore_v1`) の物理壁では停止規則が使えない → カードの雛形を `axismach_final_v1` に替えて問題 YAML と `_euler.yaml` を作り直し、手組みの mesh を当て直す (#3b) | O |
| 3b | M5.5 を最終の雛形で作り直し、Euler → NS | カードの `recipe: axismach_final_v1` で問題 YAML・`_euler.yaml` を作り、`hand_mesh.py build` で mesh・mesh_euler を当て直して品質を確かめる。Euler `run_0001_euler` → NS `run_0002_ns` (停止規則、上限 20 万) を AWS A の FP64 のビルドで。**途中 (2026-10-15、implementer、YAML と Euler の準備)**: カード → `axismach_final_v1` で生成 (A `~/forge-cardmesh` 94790833): `case/71.m55_combustion_wt/problem_m55_final_card.yaml`・`_euler.yaml` (M_knot 自動 2.80 で不変、r_t 0.0999084914252905 m、寸法 PASS)。`hand_mesh.py build --plan-item "§5.1 #3b" --axis-cap-frac 0.03` → `problem_m55_final.yaml` (NS)・`problem_m55_final_euler.yaml` (Euler; `wt_card.write_euler_companion` の写し、生成器の写しと mesh 以外が同じことを照合)。係数・正規化・ni 4,800 は #3a と同じ。NS の径方向はユーザ決定 (2026-10-15「上限3%」) で `axis_cap_frac 0.03`、nj 92 (全断面の隣接比の最大 1.1974 ≤ 1.2、nj 91 は 1.2013) → NS 441,600 点、Euler は nj 97・465,600 点のまま。検査 (FP64 の変換器、forge 不起動): NS 品質 PASS (AR 最大 4849.0)・AR 緩和 違反 0・変換後の第一層 実寸/(c·y₁,req) 最大 0.9999997 (超過 0)・見積もりの y₁⁺ 最大 0.797・軸側のセル / r_w 0.030; Euler PASS (AR 最大 6.1)。NS の準備の試行で停止規則 use (上限 200,000)。`run_0001_euler` は準備のみ (A、未起動)。記録 `case/71.m55_combustion_wt/_mesh/hand/final/` **径方向の決定 (ユーザ 2026-10-15「上限3%にしてみますか」)**: 雛形の等比 (壁から軸まで比 1.10〜1.14、nj 97、軸側のセルが局所半径の 11〜12 %) を、`axis_cap_frac 0.03` + 近壁の比 ≤ 1.2 (nj は全断面で比 ≤ 1.2 の最小 = 92) に替えた。軸側のセルは全区間で 3.0 %、節点 441,600。Euler の格子は壁面垂直方向がほぼ一様 (第一層 0.5 %・比 1.014) なので変えない。**投入 (主セッション)**: `run_0001_euler` を AWS A で起動 (FP64 `~/forge-integ-dflow2-fp64`、soft → mid → 本段 cfl 1.0・12,000 step、約 12 ms/step) ユーザが nj 92・上限 3 % の NS 格子を ParaView で確認 (2026-10-15「オッケーそうでした」)。`run_0001_euler` 完走 (227 s、出口の軸の M 5.5002、dM 最大 0.06 % Md、NaN なし; `check_convergence --segment` は NOT CONVERGED = 残差 2.1〜2.5 桁で停滞、M6 の Euler 参照と同じ型 → δ* の参照としては未適格と記録)。`run_0002_ns` を起動 (停止規則 use・上限 20 万、物理壁 single_bspline) **NS の未到達の診断 (2026-10-15、implementer、forge 0 step)**: codex 諮問 `notes/reviews/2026-10-15-m55-ns-unsettled-diagnose.md` (incomplete: 証拠不足) の証拠を A から回収して `case/71.m55_combustion_wt/_runs/` に (h5 なし: 全残差履歴の gz・stage_manifest・solverConfig・bcond・prepare_info・MAIN_STAGE_STARTED・forge_launches・stop_rule.json の全体・帯の抽出・ドライバの log・forge の log の先頭と末尾、Euler の同類)。バイナリは FP64 (`flowFormat.hpp` double)、sha256 3a9e347b…、`~/forge-integ-dflow2-fp64` (HEAD 2fd36b0c、ビルドより新しいソルバの源 0 件)。判定 (原文 `_runs/run_0002_ns/judge/`): (1) 記録の stage_manifest に本段が無い (見張りの `save_record` が OSError で落ち、`run_staged_ns` の `_record(cfg_main, "main")` に届かなかった) ので `--segment` をそのまま回すと S1 → S2 を判定する。本段を足した写しで区間 = 本段 0〜183,998 → `check_convergence --segment` NOT CONVERGED (still converging; ρ 1.8・ρE 2.2・運動量 2.4・k 2.9・ω 2.7 桁)。各列の最小は 157,795〜172,251 step で、直近 2 万 step の log10 の平均はその前の 2 万 step より 0.07〜0.44 桁高い。(2) 停止規則の系列 (1,000〜183,000) の `check_quasisteady --series-csv` DRIFTING (r DRIFTING、θ_r 0.4・0.7 TRANSIENT-UNSETTLED、θ_r 0.95・ṁ_in・ṁ_out・出口コア M STEADY)。(3) res_183000: 非有限 0、ρ・P・T > 0、P > P0 は壁の行 (U = 0) の 298 点で相対 ≤ 5.8e−8、P0 > Pt (1 + 1e−6) が壁から 27 点以上離れた主流〜軸の 50,441 点 (最大 +0.24 %)。(4) 壁解像 (`wall`、area、step 183000) PASS: y1+ > 1 の面積 0.0 %・最大 1.004 (停止の step ではないので V4 の判定ではない)。(5) 実効の設定: 問題 YAML から生成器で作り直した S1・S2・本段の solverConfig と bcond が run の記録と一致 (sha256・cfg_fnv・bcond_fnv)。(6) θ_r の測定帯の A/B (同じ保存場 100k〜183k の 12 枚・同じ Euler res_12000・同じ 3 断面、帯だけ変える; `_runs/run_0002_ns/theta_band_ab.json`): D_θ [%/2 万 step] A 0.97・0.91・0.54、B 2.33・1.96・1.75 → 事前に固定した判定で **「帯の選び方だけで未到達を説明できる」を棄却** (3 断面とも両腕で > 0.05 %)。別に、帯 A (res_1000 から) は後期の縁 y* を使用域の全列で覆っていない (帯 B/A 1.40〜5.52 倍)。容量: A の空き 16 GB、全場 1 枚 100,605,412 B、壁の出力 1 枚 1,170,416 B。記録 `case/71.m55_combustion_wt/_runs/README.md` **判断 (2026-10-15、codex diagnose `notes/reviews/2026-10-15-m55-ns-unsettled-2-diagnose.md`、全件採用)**: 長い延長は保留し、先に既存の保存場で全圧の A/B (A = 出力の原始量・h0 から、B = 保存量 ro・roU・roe・roY からソルバと同じ規約で再構成した状態から、forge 0 step、100k〜183k の 12 枚) を行う (#3c)。manifest の写しに本段を足した判定は採用 (原本を残し「事後復元・本段 0〜183,998」と明記、正常終了の証明ではない)。追加 約 16.8 万 step は到達の予測に使わない。延長に進むときは: 新しい run で同じメッシュの保存量を index restart (前段は回さない)・出力 2,000 step・`evaluate.stop_rule: true`・保持 4 万・**追加 20 万 step を予算の上限** (到達の予測ではない、2 万 step ごとに再評価)・新しい帯と Euler の res_12000 (sha256) を固定して θ_r は別の系列・残差は同じ設定の区間だけ連結 (183k 以後を重ねない)・容量は約 10.5 GB で投入の時に測り直す | O |
| 3c | 全圧の超過の A/B (forge 0 step、**済み: 複合要因、超過の大部分は保存量の状態にある; 3 回目の諮問は thermoFloat の 1 更新の A/B を提案 → ユーザ決定で延長を全部 double にして代える (#3d)**) | 上の判断のとおり。節点ごとの ΔP0/Pt、超過の点数・分位点・位置 (入口・軸・内部)、エネルギーの逆算の残差。判定 (事前に固定): 両経路の差が全時点で ≤ 1e−6 Pt かつ超過が残る → 出力の状態の不整合を主因とする仮説を棄却; B で超過が消え A/B の差が超過を説明 → 保存場の側という仮説を棄却し出力・後処理の経路を直す; 一部だけ → 複合要因として記録。**結果 (2026-10-15、implementer、forge 0 step)**: 台本 `case/71.m55_combustion_wt/_runs/run_0002_ns/p0_ab.py` (A の `~/m55_diag_0002/p0ab/` で実行、ツールは A の `~/forge-cardmesh` b0d528bf = 手元と solver_density_cuda/ の差 0; 原文 `_runs/run_0002_ns/p0_ab.json`・`p0_ab_series.csv`・`p0_ab_series_quasisteady.txt`・`p0_ab.stdout.txt`)。同じ 12 枚 (100k〜180k の 1 万 step ごと + 181k〜183k)・同じ物性の記録 (res の属性 → `resolved_species_fc8ba1a8828147e3_70cc58aa292c9504.yaml`、sha256 70cc58aa…)・datum 298.15 K。B は `dependentVariables_d.cu` の thermoFloat 1・凝縮 off の経路 (ρ の床・float32 の組成の非負化と正規化・e = roe/ρ − \|U\|²/2・T は e_mix の Newton・P = ρ R_mix T・h0 = e_mix + P/ρ + \|U\|²/2、h0_includes_k 0) を写し、全温・全圧は A と同じ式。A は既報 (183k: 50,441 点・最大 +0.23580 %) を再現。(1) \|P0_A − P0_B\|/Pt の全節点の最大は 8.59e−6 (160k)〜1.03e−4 (100k)、1e−6 を超える節点は 7,995〜116,254 (内部が 97〜99.5 %、壁の行 0; x 0〜1 m の割合は時点で 10〜75 %)、中央値 2.5e−8〜2.4e−7。(2) B でも超過が残る: 49,862〜79,964 点 (A 50,161〜80,210、A∩B は A の 97.7〜99.5 %)、B の最大 P0/Pt − 1 は +0.2232〜+0.6166 % (A の最大との差 ≤ 1.4e−5)。位置は A と同じで、入口の列 0・壁の行 0・壁から 27〜91 点・x −1〜8 m (軸の行は B で 100k〜160k に 1,947〜3,782 点、180k〜183k に 51〜107 点; A は 1,947〜3,803 点と 64〜150 点)。(3) A の超過の点で A − B が超過量を説明する割合 (Σ) は 0.11〜0.65 %。(4) エネルギーの逆算の残差 \|e_mix(T_B) − e\| ≤ 2.6e−9 J/kg (T の床・P の床・ρ の床・Y の負値の切り捨て 0 点)。出力の T を保存量の組成と e に当てた差は −35〜+32 J/kg (中央値 0.05〜0.12 J/kg)。T0 > Tt + 1 K は A・B とも入口の列の 9 点 (100k だけ A 59・B 57 点)、B の静温も同じ 9 点で Tt を最大 2.62 K 超える。組成の float32 の丸めの寄与 (B と、double で正規化した B_dbl の差) ≤ 3.4e−7 Pt。**判定 (事前に固定した規則どおり): 複合要因** — 条件 1 (全時点で差 ≤ 1e−6 Pt) も条件 2 (B で超過が消える) も満たさない。差は超過量の 1 % 未満で、超過の大部分は保存量から作り直した状態にも残る (過渡の原因の確定には進まない)。時系列 (最後の時点の A の超過の集合 S、50,441 点を固定、12 枚) の `check_quasisteady --series-csv` (既定の閾値): **DRIFTING** — 超過の点数は A・B とも STEADY、S の平均の ΔP0/Pt は DRIFTING (−5.06e−5 → +1.46e−4、182k で 1.0e−7 下がった 1 回を除いて増加)、S の平均の h0 − h0_in は DRIFTING (−55 → −2,152 J/kg で単調に減少)、軸の近く (S ∩ j ≤ 2) の h0 は OSCILLATING。A と B の列の判定は全て同じ | O |
| 3d | **NS の延長を全部 double で** (ユーザ決定 2026-10-15「どうせやるなら全部doubleかな」) | 3 回目の諮問 (`notes/reviews/2026-10-15-m55-ns-unsettled-3-diagnose.md`) は `thermoFloat` 1/0 の 1 更新の A/B を先に、延長は保留としたが、ユーザの決定で**延長を `physProp.thermoFloat: 0` (FP64 のビルドで状態・熱物性とも double) で回す**。混合精度の経路を延長から外すことで、H2 のうち `thermoFloat` の寄与を延長の結果に混ぜない。条件は 2 回目の諮問のとおり: 新しい run `run_0003_ns_cont` に run_0002_ns の res_183000 を同じメッシュで index restart (`restart_field.py`)、前段は回さない、出力 2,000 step・`evaluate.stop_rule: true`・保持 4 万・追加 20 万 step の予算 (2 万 step ごとに再評価)、帯と Euler の参照 (run_0001_euler/res_12000、sha256) を固定し θ_r は別の系列、判定区間は延長の本段だけ (設定が違うので run_0002_ns とは連結しない)、容量を投入の時に測る。監視: 残差・r・θ_r、S の h0・P0、壁の積分熱流、入口・出口のエネルギー流束、入口の 9 点の静温の超過。入口の面の実効の流入とエネルギー流束は延長の序盤に確かめる (違反なら止めて直す) **投入 (2026-10-15、implementer)**: AWS A `case/71.m55_combustion_wt/run_0003_ns_cont` (台本・記録 `case/71.m55_combustion_wt/_runs/run_0003_ns_cont/`、`cont_0003.py`、A の木 f85e4f1b)。res_183000 (sha256 c29044a8…) を `restart_field.py --keep-src-dtype` で同じ nozzle.h5 に index restart (VERDICT: OK、9 量、SRC とビット一致)、`run_staged_ns(stages="none", stop_rule=True)`。run_0002_ns の本段から変えたのは `physProp.thermoFloat: 0`・`time.deltaT.implicitSolvePrecision: 1` (timeIntegration 11・blockDPLUR 1・lowMachPrecond 0 を確認)・`outStepInterval` 2,000 と、環境変数 `FORGE_TRANSPORT_TABLE=0`・`FORGE_DIAG_FACE_H_DOUBLE=1` (`run_case.sh` の記録の対象に足して `RUN_PROVENANCE.txt` に記録; 4 つとも forge の起動の行で効きを確認)。停止規則は明示・窓 2 万・保持 4 万・上限 200,000 (= 追加分)、Euler の参照 run_0001_euler/res_12000.h5 (sha256 93cb88b0…) に固定、帯は res_2000 から。FP64 のビルド (sha256 3a9e347b…)、26.9 ms/step (run_0002_ns の約 1.9 倍、上限まで約 90 分)、序盤 (step 0〜14,478) の残差・出力に非有限 0。容量: 投入の前の空き 16.36 GB、見込み約 10.4 GB。**残る float**: 化学種の拡散の D_s・h_s・X の float の写し (`speciesTransport_d.cu:240-356`)、リミッタの長さの `sqrtf` (`limiter_d.cu:398`)、XMF の型の宣言。入口の面の確認 (res_183000、forge 0 step、`inlet_face_check.py`): 92 面のうち 91 が流入・逆流 0 (壁との角の 1 点は U = 0 で ṁ = 0)、流入の面の Ht_R − h0_in は相対 1e−13 以下・合計のエネルギー流束 = ṁ·h0_in (相対 −9.5e−14)・mach ≤ 0.066、T0 > Tt + 1 K の 9 点の面も全て流入で境界の状態は Ts_b 1,249.25 K — 事前に置いた違反の目安 (逆流・Ht_R ≠ h0_in・合計 ≠ ṁh0_in・mach ≥ 1) に当たらない。壁から 1〜5 点の面は節点が遅く冷たいので Ps_b が節点の静圧より最大 +0.27 % (面積の割合 1.4e−5) | O |
| 5 | 成果物の対応 (M2) | `export/mesh_structure.py`・`deliverables.py` を壁法線の近壁層の格子でも作れるように (既定化の前提) | O |
| 2 | 実装 | `estimate_mesh` の `local_v1`、雛形 (`axismach_explore_v1`・`axismach_final_v1`) の `mesh.method`・`first_layer_calibration`・`wall_normal_layer`、記録、試験 | O |
| 3 | M5.5 のカードの生成 (AWS A) と Euler → NS | ユーザ了承の節点数で | O |
| 4 | skill の手順 | `nozzle-design`・`design-intake` に格子の作り方と NS の後の壁解像の確認 | O |

## 6. 検証

事前登録 (2026-10-15、上位に諮る前の案):

- **V1 見積もりの再現**: 係数 1.0 で M5.5 の ni が 3,163 (±10)、節点 306,811 前後になる (見積もりの台本と同じ入力)。
- **V2 格子の品質**: M5.5 と M5 デモ (case/70 の c2b のカード) の `local_v1` の格子で `check_mesh_quality.py --bcond --ar-max 5000` が PASS、隣り合う dx の比 ≤ 1.1、壁直交性 (`wall_normal_layer`) の判定 PASS。
- **V3 第一層**: 表の値 ≤ c(x)·y₁,req(x) (全断面)、見積もりの y₁⁺ の最大 ≤ 0.8。
- **V4 NS の後** (M5.5 の NS、停止規則まで): 壁ごとの y₁⁺ > 1 の面積 ≤ 2 % (`--groups <全 no-slip 壁> --weight area --step <停止の step>`)。超えたら §4.2 の手順で表を直す (この plan の不合格ではなく、係数の見直しの材料として記録)。
- 判定不能は合格にしない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-15 | [2026-10-15-tooling-nozzle-card-local-mesh-plan.md](../../notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan.md) | NO-GO (分類 incomplete: 見積もりの台本・出力と M6 の比の元の記録がこの作業ツリーに無い), C0/M7/m1 | **全件採用、進め方を 2 段に分ける**: M5.5 は先にケースだけの手組みの格子 (M6 の H と同じやり方; #3a) で進め、生成器の既定の変更 (#2) は指摘を直して再レビューの後。M1 (実寸の下包絡を半径で割ると補間で要求を超える: 反例 1.0607 倍) → f_req = c·y₁,req/r_w の下包絡を表にし、補間の区間の中も検査、物理壁でも実寸を再検査。M2 (`wall_normal_layer` で成果物 `mesh_structure.json`・STEP の再構築が作られない) → 既定化の前に成果物の対応 (#5)。M3 (補正の係数の根拠: H の表は冷却倍率を掛けた見込み、縮流部の比 0.61〜2.2 に係数 1.0 は保守側でない) → 「初期格子の経験の係数」と位置付け、元の run・step・比の分布を固定して commit、縮流部も比の最小 (0.6) に、補間後の c(x) が比を超えない検査。M4 (V2 の `--bcond` がこの枝の検査器に無い・AR 緩和の条件が未定量) → `runner.mesh_ar_relax_check` (AR > 1000 の全セルで skew ≤ 0.1) と変換後の第一層厚の照合・`msh_digits 17` を V2 に。M5 (V4 を記録だけで閉じられる) → 最終の完了は再調整の後の V4 PASS、全 no-slip 壁の同じ step の出力の存在、`check_convergence --segment`・壁量の系列の `check_quasisteady` を事前登録。M6 (見積もりと実生成の設定の契約・旧い overrides) → 見積もり・節点数・品質・実生成で同じパラメータの組み立て、雛形の旧キーの除去、明示の overrides の衝突は拒否、`mesh.ar_max` に従う。M7 (Euler の ni も変わる) → Euler の ni が NS と結合していることを明記し、M5 で出口 M・流量・δ_r の基準格子との差と準定常を受入に。m8 (参照の解決) → 再現の入力を case/71 に commit、関連 plan は統合の後に参照 |

## 7. 影響範囲

- `design/forge_design/intake/wt_card.py`・雛形 2 本・試験。既存のカード (`mesh.method` を書いていないもの) は新しい既定で格子が変わる (ユーザの方針「標準と決めた方式はコードの既定にする」)。
- case/71 の M5.5 の生成。

## 8. 完了条件

- [ ] §6 V1〜V3 を満たし、V4 を記録
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録
- [ ] `status: done`、§9 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-15` — plan 段レビュー NO-GO (C0/M7/m1) を全件採用。M5.5 は手組みの格子 (#3a) で先行、生成器の既定化は書き直して再レビュー。
- `2026-10-15` — 起票 (ユーザ「古いやり方でメッシュ作ってない？」、決定「a」)。見積もり (3 通りの作り方の ni・節点) と見積もりの式の偏り (M6 の H の表との比) を §3 に。
