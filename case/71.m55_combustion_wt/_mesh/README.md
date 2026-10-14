# case/71 の格子の見積もりと手組みの格子の記録

plan [`tooling-nozzle-card-local-mesh`](../../../plans/active/tooling-nozzle-card-local-mesh.md) §3 の見積もりと §5.1 #3a の手組みの格子の再現の入力。
codex plan 段レビュー ([記録](../../../notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan.md)) の「証拠不足」と M3 (比の分布の固定) への対応として commit した。

## 見積もり (2026-10-15、plan §3)

AWS A の `~/m55_mesh_est/` で、`~/forge-mknot` (7713f84c、`estimate_mesh` は変えていない) を import して走らせた使い捨ての台本の写し (パスは書き換えていない)。

| ファイル | 中身 |
| --- | --- |
| `step1_chain.py`・`card_est.yaml`・`step1.log`・`est.json`・`prob_parts.json` | カード (M_knot 2.80 固定) → 問題 → 設計チェーン → `estimate_mesh` |
| `step2_est.py`・`step2_m55*.log`・`m55*_result.json`・`m55_y1req.csv` | 作り方 (A)(B)(C) の最小の ni。`sc*` は入口の端の頭打ち、`t*` は第一層の表の倍率の感度 |
| `plot_est.py` | 見積もりの図 (手元の scratchpad で作った) |
| `step1_c45.py`・`c45_hcheck.py`・`c45/` | M6 (case/45 の H、`c45/problem_h_hall.yaml` = `problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml` の初期線を Hall にした写し) の第一層の表と、同じ式の見積もり (Tw 370 K) の比 |
| `c45_ratio.py`・`c45/h_over_est_ratio.csv`・`c45/h_over_est_ratio_summary.json` | 比の断面ごとの分布と区間ごとの最小 (手組みの格子の係数 c(x) の上限に使う) |

比 (H の表 ÷ 見積もりの要求) の区間ごとの最小: 配管 2.218・縮流部 0.612・スロート付近 0.511・拡大部 0.499 (x = 25.4 r_t)・後半 0.558。
H の表は 300 K 用に断熱 NS の y₁⁺ × 冷却倍率 8.5〜9.4 で作った見込みで、300 K の実測の要求ではない (初期格子の経験の係数として使う)。

## 手組みの格子 (2026-10-15、plan §5.1 #3a)

`hand_mesh.py build` が `../problem_m55_card.yaml` (カードの生成器の出力) の `mesh`・`mesh_euler` を置き換えて `../problem_m55_ns.yaml` を書き、
`hand_mesh.py check` が prepare_ns・prepare と同じ経路で格子を作って変換し、forge を起動せずに検査する (AWS A `~/forge-cardmesh`、
作業ディレクトリ `~/m55_mesh_hand/`)。記録は `hand/`:

| ファイル | 中身 |
| --- | --- |
| `hand/card_gen.log` | カードの生成 (M_knot 自動 2.80、r_t 0.09990877093936801 m) |
| `hand/build.log`・`hand/build_record.json`・`hand/build_record.csv` | 式の照合・係数の検査・第一層の表の検査・ni の決定。CSV は断面 (0.02 r_t おき) の y₁,req・c・f_req・表 |
| `hand/check.log`・`hand/check_record.json` | 品質・AR 緩和・変換後の第一層厚の照合・見積もりの y₁⁺・節点数・dx の比 (NS・Euler) |
| `hand/ns/`・`hand/euler/` | `MESH_QUALITY.txt`・`MESH_AR_RELAX.json`、NS の `first_layer_check.csv` (壁の断面ごとの実寸・c·y₁,req・比) |

係数 c(x): 配管 1.0・縮流部 0.6・スロート付近 0.5・拡大部 0.49 (拡大部は比の最小 0.49874 が 0.5 を下回るので有効 2 桁に切り下げ。0.5 なら ni 4,780)。

## 最終の雛形での作り直し (2026-10-15、plan §5.1 #3b)

カードの雛形を `axismach_final_v1` に替えて (`../card_m55.yaml`; 探索の雛形の物理壁では NS の停止規則が使えない) AWS A `~/forge-cardmesh` (94790833)
で生成 → `../problem_m55_final_card.yaml`・`../problem_m55_final_card_euler.yaml` (M_knot 自動 2.80、r_t 0.0999084914252905 m、寸法 PASS)。
`hand_mesh.py build --plan-item "§5.1 #3b" --axis-cap-frac 0.03` で `mesh`・`mesh_euler` を置き換えて `../problem_m55_final.yaml` (NS) を書き、
Euler 用の写し `../problem_m55_final_euler.yaml` を生成器と同じ `wt_card.write_euler_companion` で書く (生成器の写しと mesh・mesh_euler 以外が同じことを照合)。
係数・正規化・流れ方向の作り方は #3a と同じ (ni 4,800)。**径方向** (ユーザ 2026-10-15「上限3%」): NS は `axis_cap_frac 0.03` (壁の第一セルから等比、
間隔 / r_w が 0.03 で軸まで一様) で、nj は全断面で `mesh2d._radial_fracs_capfixed` の隣接比の最大 ≤ 1.2 の最小 = 92 (比の最大 1.19744、nj 91 では 1.20130;
最も薄い第一セルのスロート付近で決まる)。NS 441,600 点・Euler 465,600 点 (Euler は nj 97 のまま)。検査は FP64 の変換器
(`~/forge-integ-dflow2-fp64`、2fd36b0c) で `hand_mesh.py check` (作業ディレクトリ A `~/m55_final/`)。記録は `hand/final/`:

| ファイル | 中身 |
| --- | --- |
| `hand/final/card_gen.log` | カードの生成 (最終の雛形) |
| `hand/final/build.log`・`build_record.json`・`build_record.csv` | 係数・第一層の表・ni・径方向の nj の決定・Euler の写しの照合 |
| `hand/final/check.log`・`check_record.json` | NS (nj 92): 品質・AR 緩和・変換後の第一層厚・見積もりの y₁⁺・区間ごとの径方向 (軸側のセル / r_w・隣接比)、Euler |
| `hand/final/ns/`・`hand/final/euler/` | `MESH_QUALITY.txt`・`MESH_AR_RELAX.json`、NS の `first_layer_check.csv` |
| `hand/final/check_nj97.log`・`check_record_nj97.json` | 径方向を変える前 (雛形の nj 97 の等比) の同じ検査 (参考) |
| `hand/final/euler_prepare.log` | `run_0001_euler` の準備 (`runner_axismach --prepare-only`) |
| `hand/final/ns_prep_probe.log` | NS の準備の試行 (`deltastar_loop --init-integral --prepare-only`、A の使い捨てのディレクトリ; 停止規則 use・上限 200,000 の確認。消した) |
