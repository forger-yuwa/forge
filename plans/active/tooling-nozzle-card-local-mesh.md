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
| 3a | **M5.5 の手組みの格子 (先行)** | ユーザ決定「a」(2026-10-15): case/71 の問題 YAML の mesh 節を、正規化した要求 f_req = c(x)·y₁,req/r_w の下包絡の表 (補間の区間の中も検査・物理壁で再検査)・流れ方向の密度の表・`throat_refine 1`・`wall_normal_layer [0.02, 0.8]`・`msh_digits 17`・`max_nodes` (ユーザ了承「a」、約 45 万) で手で組む。係数 c: 配管 1.0、縮流部 0.6、スロート・拡大部 0.5 (M6 の H の表 ÷ 見積もりの比の最小、初期格子の経験の係数)。品質: `check_mesh_quality --ar-max 5000`・`mesh_ar_relax_check`・変換後の第一層厚の照合。生成器の既定は変えない。**実装 (2026-10-15、implementer)**: `case/71.m55_combustion_wt/problem_m55_ns.yaml` (`_mesh/hand_mesh.py build`; 生成は `card_m55.yaml` に `mesh.ni: 3000` を足した `problem_m55_card.yaml`、M_knot 自動 2.80)。ni 4,800 × nj 97 = 465,600 (NS・Euler)。係数の拡大部は比の最小 0.49874 < 0.5 なので 0.49 (0.5 なら ni 4,780)。CFD 前の物理壁の格子 (forge 不起動、A `~/m55_mesh_hand/check/`): 品質 PASS (AR 最大 4849.4)・AR > 1000 の 47,607 セルでスキュー最大 0.0015 (違反 0; 出口面の非平面 0.62 mm は記録だけ)・第一層の実寸 / (c·y₁,req) 最大 0.999999 (超過 0)・見積もりの y₁⁺ 最大 0.797・dx の比の最大 1.033。Euler PASS (AR 最大 6.1)。再現の入力と記録は `case/71.m55_combustion_wt/_mesh/` | O |
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
