# codex レビュー: convection-slau-wall-normal-chi-usage-rule (plan)

- **plan**: [`plans/active/convection-slau-wall-normal-chi-usage-rule.md`](../../plans/active/convection-slau-wall-normal-chi-usage-rule.md)
- **stage**: `plan`
- **date**: 2026-09-25
- **commit**: `e5ed2312` (feature/sern-design)
- **codex**: effort `high`, 6.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m0
- **focus**: 2 回目の plan 段。1 回目 (notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan.md) の Major 7 件への対応が §4・§6・§5.1 で閉じているか、特に §4.3 の解像度キー表・§4.4 の比較式・§6 Q1 の 3 結論が全ケースを覆い測る前に固定されているかに集中してほしい。
- **extra**: `notes/reviews/2026-09-25-convection-slau-wall-normal-chi-usage-rule-plan.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

目的の限定、Roe 比較の保留、Q2 の独立判定は妥当です。  
しかし、格子表は形状不変条件を満たさず、比較式と Q1 の分岐も誤判定を起こします。前回 Major 7 件を「対応済み」と扱える状態ではありません。

1. **Major — §4.3 の表では、固定形状の格子比較が成立しない**

   **根拠:** [`mesh_sern.py:347`](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern.py:347) は、格子上でサンプリングしたランプ高さの最大値から `y_veh` を決めます。同ファイルの `359` 行では、後縁勾配も格子の差分から決めています。したがって、入力の `vehicle_*` を固定しても機体形状は固定されません。

   指定 YAML と表の解像度で、ファイルを保存せず座標を生成して確認しました。

   | 検査量 | 粗 | 中 | 細 |
   |---|---:|---:|---:|
   | 節点数 | 62,067 | 124,158 | 245,163 |
   | `y_veh/H` | 2.746515521 | 2.746525681 | 2.746524518 |
   | `vehicle` 壁から粗格子壁折れ線への最大距離 `/H` | — | `3.97e-4` | `6.56e-4` |

   規定の `1e-6 H` を超えます。ランプ・カウルでも超過しました。

   解像度指定にも抜けがあります。対象 YAML は `t_base` 未指定で既定 `0`、生成結果も `has_wake=False` です。表の `nj_wake` は効かず、`first_wake_frac=t_base/5` の細分化にも意味がありません。また、[`mesh_sern.py:401`](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern.py:401) が使う `first_top_frac` は表に含まれず、機体上面の第一層幅は全水準で `0.02 H` のままでした。

   **対案:** 生産格子の壁折れ線を固定入力として、その線上に節点を追加する格子生成手順を先に確定してください。実際に存在するブロックだけで解像度表を作り、`first_top_frac` と固定点数のフィレット区間も点検すること。§5.1 にこの前提作業を追加し、形状・局所間隔・品質ゲートを通してから Q1 に進むべきです。今回実施したのは座標検査で、CFD・メッシュ品質判定は実施していません。

2. **Major — 生産許容の参照値が旧値になっている**

   **根拠:** [対象 plan:40](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:40) は、sern-3d §8 を出典として `|ΔC_M|≤0.02` としています。しかし、[実際の §8:1918](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1918) は **`0.05`** です。直後に「2026-09-21、ユーザ決定で `0.02→0.05`」と明記されています。

   したがって「差 `0.023` は現行の生産許容を超える」という前提は成立しません。**前回レビューもこの旧値を踏襲しており、訂正が必要です。** ただし、既存 `C_L/C_M` の未定常性が解消されるわけではありません。

   **対案:** 今回の目的である「生産許容との比較」には現行値 `0.05` を採用し、§3・§4.4 の派生閾値・§6 を同期してください。準定常閾値についても、絶対許容とツールの相対閾値を区別して出所を明記すること。

3. **Major — Q1 の「3 結論」は、依然として網羅的でも排他的でもない**

   **根拠:** [対象 plan:135](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:135) の条件を、そのまま人工データで評価しました。振幅はゼロ、許容は記載どおり `0.02` とします。

   | 粗→中→細の値 | 判定結果 |
   |---|---|
   | flag 0: `0, 0, 0`／flag 1: `0.030, 0.018, 0.010` | **どの結論にも入らない** |
   | flag 0: `0, 0.001, 0`／flag 1: `0.005, 0.006, 0.005` | **「許容内」と「判定不能」の両方に入る** |

   前者は粗格子だけ許容超過、細格子間変化は許容内です。後者は全差が許容内ですが、格子間変化の符号が反転します。追加する第4水準の解像度と、追加後にどの水準で判定するかも未指定です。

   **対案:** 優先順位付きの `if/elif/else` として判定表を定義してください。準定常・格子間比較の前提不成立を先に処理し、合格条件にも有意差条件にも入らない場合は「判定不能」とすること。符号反転の扱い、第4水準、複数量の集約順も固定し、上記反例を判定器の試験に含めてください。

4. **Major — `|平均差| + span/2 の和` は、一般には上限にならない**

   **根拠:** [対象 plan:85](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:85) は、平均を中心に半幅 `span/2` を置いています。しかし、非対称な変動では平均と最大・最小の中点は一致しません。

   人工系列で、flag 0 を一定 `−7`、flag 1 を各周期の9点で `−6.987`、1点で `−6.977` とすると、

   - 平均差 `0.014`、半幅 `0.005`、現式の `U=0.019`
   - 実際の最大差 **`0.023`**
   - 末尾窓の前半・後半の平均差 `s=0`

   となり、許容 `0.02` を誤って通します。現行 `check_quasisteady.classify_series` に §4.5 の `C_M` 閾値を与えた結果も **`OSCILLATING`** であり、plan が比較を許す系列です。

   **対案:** 観測された変動帯を比較するなら、各系列の区間 `I_f=[min C_f,max C_f]` を使い、差区間を

   \[
   I_\Delta=[\min C_1-\max C_0,\ \max C_1-\min C_0]
   \]

   としてください。その端点絶対値の最大を上限とすれば、非対称な変動も扱えます。これは観測窓内の上限であり、平均値の統計的不確かさとは区別して記録すること。

5. **Major — 上限が許容を超えただけでは、「flag が許容を超えて解を変える」と結論できない**

   **根拠:** [対象 plan:135](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-usage-rule.md:135) は `U>許容` を有意な flag 差として扱います。

   例えば全格子・両 flag が同じ `−7±0.015` の周期系列なら、平均差も格子間変化もゼロです。それでも現式は `U=0.030>0.02` となり、**全く同じ系列を「flag が解を変える」と判定します**。この系列も指定閾値で `OSCILLATING`、窓間平均差はゼロでした。

   **対案:** 「許容内を保証できない」と「許容超過の差を確認した」を分けてください。上記の差区間について、全体が許容帯内なら「許容内」、全体が許容帯外なら「差が残る」、それ以外は「判定不能」とするのが一貫します。格子間変化にも同じ不確かさの扱いを適用してください。

6. **Major — §4.2 の適用条件を、指定された診断ツールでは測れない**

   **根拠:** [`diag_wall_cv_budget.py:4`](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_wall_cv_budget.py:4) は明記どおり **1次 SLAU** の再計算です。`187–191` 行でも節点値をそのまま使い、2次再構成を行いません。

   さらに [`75–77` 行](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_wall_cv_budget.py:75) は、`--wall-normal-chi` 指定時に元の `chi` を `chi_n` で上書きします。plan の「同オプションで元の `χ=0` を確認する」という使い方は逆です。壁速度ゼロ、隣接接線速度 `1100 m/s`、音速 `700 m/s` の入力では、表示対象の `chi` は指定なしで `0`、指定ありで `1` になりました。

   **対案:** 今回は適用診断を、検証済みの **`convMethod: 0` の起動区間**に限定してください。同じ面状態から元の `χ`、`χ_n`、圧力差、変更前後の質量流束を併記し、対象壁 CV への補充が実際に増えることを条件にするべきです。2次生産場への適用診断は、再構成後の面状態を取得できるまで保留し、この制限と必要なツール作業を §5.1 に反映してください。

7. **Major — 衝撃足ゲートが目的の衝撃を識別しておらず、格子間比較の定義も不足している**

   **根拠:** [`v3sern_foot_perrun.py:42`](/home/sano/work/forge-sern-design/case/46.sern_design/v3sern_foot_perrun.py:42) は、ランプ全域から `max|dp/dx|` を選びます。圧縮と膨張を区別しません。

   保存された [`perrun_flag0.csv:2`](/home/sano/work/forge-sern-design/case/46.sern_design/_r3_m1m3/perrun_flag0.csv:2) と flag 1 の系列では、全72点の `xfoot` が **`−0.000877683 m`**。flag 0 の末尾平均圧力は、その前後で **約99.57→56.60→53.55 kPa** と低下しています。追っているのは上流の圧力低下部であり、これをカウル衝撃足の証拠にはできません。

   両 CSV を指定閾値で再判定した VERDICT は **`OVERALL: ALL STEADY`** でした。これは抽出量の定常性を保証するだけで、対象の識別を保証しません。

   また、[`v3sern_series.py:53`](/home/sano/work/forge-sern-design/case/46.sern_design/v3sern_series.py:53) は単一メッシュの節点番号で両 run を比較します。異なる格子には使えず、「衝撃足 L2」は二つの場の差なので、§6 の各 run のスカラー `C` として引き算する定義も成立していません。

   **対案:** 衝撃を探索する物理区間と圧力上昇の条件を先に固定してください。壁圧は共通の物理座標へ補間し、固定窓・重み・正規化分母で、flag 間と格子間の相対 L2 を直接計算すること。各 run の準定常は絶対位置と窓内壁圧で判定し、「格子1つ」の基準も固定してください。抽出・格子間比較ツールの整備が §5.1 に必要です。

推奨は、**既定0・既知3D構成での限定 opt-in を維持し、上記を plan と判定器へ反映して再レビューしてから Q1・Q2 に着手すること**です。修正順は、許容値と比較・分岐規則 → 固定形状の格子生成 → 衝撃足・適用診断の測定手順です。

なお、既存 `case/46.sern_design/_r3_m1m3/sern2d_conv/flag{0,1}` の本段残差を再判定し、両側とも **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`** を確認しました。この訂正、Roe の保留、Q2 の `D_intra` 先行は維持すべきです。ファイルは変更しておらず、本レビューの提案は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 0
