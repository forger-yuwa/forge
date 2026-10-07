# codex レビュー: tooling-nozzle-upstream-poly-and-throat-sizing (plan)

- **plan**: [`plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md`](../../plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md)
- **stage**: `plan`
- **date**: 2026-10-07
- **commit**: `5dab1804` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **extra**: `notes/reviews/2026-10-07-upstream-poly-throat-sizing-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

上流の Q 化、生産と共通の δ_r による寸法逆算、B-spline の代数的構成は採用を推奨します。  
ただし、壁の保存・復元、大域スロート探索、寸法逆算の合格条件、U4 の判定条件を実装前に補う必要があります。  
HEAD `9a09c0a5` を確認し、保存データから独立再計算しました。ファイル変更・CFD 実行はありません。

目的は未解決の joint 解析経路を対象としており、既存計画の単なる重複ではありません。非 joint 経路には上流 Hermite 再生成が既にありますが、今回の「x=0 で接続して下流を保持する」構成とは異なります。case/45 の再計算では、最大形状差 **−522.66 µm**、スロート位置差 **−1.457×10⁻⁹ r_t**、半径差 **−3.779×10⁻¹³ r_t**、max｜Q″−H″｜ **1.184×10⁻³**を再現できました。

1. **Major — `poly` を保存・復元する仕様と統合試験が欠けています。**

   **根拠:** 現行の [wall_file_record](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:1268) は `legacy` 表現を必ず「ランプ＋δ_r」として保存し、[復元側](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:1354) も同じ式を使います。Q の係数・区間・`pw_upstream` は保存されません。`SingleBSplinePhysicalWall.REQUIRED_ATTRS` も `_ramp`・`ramp_gate` を必須としています（同ファイル:1036）。

   したがって、壁生成だけを変更すると、`poly`＋`physical_wall_repr: legacy` の保存時に失敗するか、復元後に別の壁になります。`prepare_info.json` に実効キーを書くことだけでは解決しません。

   また、[関連 plan:113](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:113) は代数版で W1・W5 を再実行するとしていますが、旧ランプ・最小二乗版の W3 合格を新しい構成へ引き継ぐことはできません。

   **対案:** §4・§5 に、上流方式と Q の係数・基底・区間の保存、旧形式の読込規則、方式別の必須属性を追加してください。`poly` の区分表現と `single_bspline` の双方で、保存→復元→報告の往復試験を行い、**同じ `poly` 壁を両表現で生成した W3**も再実行してください。U3 の倍精度誤差だけでは、float32 のソルバ入力がビット同一とは保証できません。

2. **Major — Q′ の全実根だけでは、方程式に必要な全域最小を保証できません。**

   **根拠:** [対象 plan:67](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:67) は物理スロートが Q 内にあると仮定しますが、`solve_rt_throat` の定義は **minₓ r_W**です。既存の [delta_r_from_table](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:925) に δ_r′(0)≥0 という制限はありません。

   独立の反例として、case/45 の設計壁に正値の表  
   `δ_r(x)=0.0014−0.000776·x·exp(−(x/0.1)²)`  
   を同じ5次補間経路で渡すと、Q 内の極小根はなく、物理スロートは **x=+0.00155091 r_t** にあります。上流曲率ゲートは **3.106×10⁻⁴ < 5×10⁻³**で通ります。これは CFD 結果ではなく、現行入力形式で作れる幾何の反例です。

   **対案:** Q と下流の S＋δ_r の**全ノット区間**について導関数の実根・端点・継ぎ目を比較し、正半径、一意な大域最小、その前後の単調性を判定してください。δ_r′(0) の正・零・負、下流に追加極値がある負例を U1/U2 に追加し、`SingleBSplinePhysicalWall` も同じ探索器を使う構成を推奨します。

3. **Major — `solve_rt` の変更方針と U0 が両立せず、出口寸法の修正を検証できません。**

   **根拠:** [対象 plan:92](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:92) は CFD 前の `solve_rt` を較正・平滑化済み経路へ変更します。一方、[U0:118](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:118) は変更前との一致を要求しています。

   コードも、旧 [solve_rt:143](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:143) と生産の [integral_delta_r:967](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:967) で異なります。保存設計を用いた再計算では、同じ `r_t=0.0766539 m` に対して：

   | δ_r の経路 | スロート補正 / r_t | 出口補正 / r_t |
   |---|---:|---:|
   | 旧・直接積分 | 0.001013529 | 0.701489173 |
   | 生産・較正＋平滑化 | 0.001389683 | 0.732755147 |

   **同じ寸法での出口半径差は 2.39666 mm**です。これは逆算後の誤差ではありませんが、旧結果への一致を合格条件にできない明確な根拠です。

   **対案:** U0 は明示 `ramp` の壁再現と、変更しない NS 後の出口逆算経路に限定してください。変更する CFD 前の出口逆算には、U2 と同様の**生産壁再生成による出口半径の往復検証**を追加してください。新しいスロート逆算には、NS 後の全分布入力、最終寸法での残差再評価、反復上限時の不合格処理も明記してください。

4. **Major — U4 は参照元の判定条件を欠落させ、比較窓も固定していません。**

   **根拠:** [対象 plan:134](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:134) の dry 条件には、参照元の [monotone plan:261](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-throat-monotone-r2.md:261) にある**壁解像条件と全残差 RISING なし**がありません。

   比較元の選択も結果を変えます。別作業ツリーにある保存 CSV を現行 `check_quasisteady.py` で再判定すると：

   | 保存 run・判定窓 | 出口 M・波・δ_E/δ_C | オーバーシュート |
   |---|---|---|
   | `case/45.isobutane_m6_d155/run_0147_ns_mono_final/`、40000〜60000 | `STEADY` | **`DRIFTING`** |
   | `case/45.isobutane_m6_d155/run_0149_ns_mono_final_ext/` の連結系列、60000〜80000 | `STEADY` | **`STEADY`** |

   両 run の保存済み main 区間の残差判定は **`NOT CONVERGED (stalled/plateau)`**です。残差 CSV は手元になく、残差判定の再実行は `NO residual_history.csv` でした。延長後の準定常判定は再確認できましたが、残差収束とは区別が必要です。run の対応は [case README:105](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:105) にあります。

   **対案:** U4 に、比較元の run・snapshot・判定窓、全残差の扱い、壁解像条件、延長規則を明記してください。メッシュ品質検査は [prepare_ns:1128](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1128) に既にあるので、その VERDICT を投入条件・成果物として指定すれば十分です。cross-mesh 移送後の起動設定と早期 NaN 検査も登録し、旧計画の「番号写像だから段階起動なし」を引き継がないでください。

5. **Minor — 既存の生産 YAML は、記載どおり自動で `poly` にはなりません。**

   **根拠:** [対象 plan:59](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:59) は次回 `poly` になると説明しますが、実際の [生産 YAML:59](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/problem_d155_ns_finemesh_recal_final_mono.yaml:59) には `pw_ramp: [-11, -6]` があります。§4.1 の併記禁止を適用すると、次回は例外になります。

   **対案:** 既定 `poly` というユーザ決定は維持してください。そのうえで、生産採用を U4 後まで保留する対象には `pw_upstream: ramp` を明示し、新規 `poly` テンプレートから `pw_ramp` を除く移行作業を §5 に追加してください。

6. **Minor — §3 の断面積差 −0.270% は再現しません。**

   **根拠:** [対象 plan:47](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:47) の値を、保存壁と δ_r から `100·(Q²/W_ramp²−1)` で再計算すると、x=−7.085 では **−0.246765%**、縮流部全体の最小は **−0.249577%（x=−6.84024）**でした。半径差最大の位置と相対断面積差最大の位置は異なります。

   **対案:** 計算式・比較位置・保存入力を添えて訂正してください。Q 化を採用する判断自体は変わりません。

**推奨は、`poly` を既定とする現方針を維持し、上記 1→2→3→4→5→6 の順で plan を修正してから実装することです。** 実装順は「壁の構築・保存・復元の共通化 → 全域極値判定 → 寸法逆算 → U0〜U3 と関連 W 試験 → U4」としてください。

今回の変更は設計側に閉じており、CUDA・block-DPLUR・周期境界の変更を要求する根拠はありません。検証は現行規則どおり node を主とし、case/45 を使う選択は妥当です。cell 回帰の追加は不要です。

**plan 未反映**です。編集禁止の依頼に従い、指摘のみを返しています。

指摘数: Critical 0 / Major 4 / Minor 2
