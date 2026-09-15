# codex レビュー: condensation-source-limiter-steady (result)

- **plan**: [`plans/active/condensation-source-limiter-steady.md`](../../plans/active/condensation-source-limiter-steady.md)
- **stage**: `result` (diff base `f1e40f41`)
- **date**: 2026-09-16
- **commit**: `a0ce5564` (feature/sern-design)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2
- **focus**: 3 回目 result レビュー (NO-GO M1–M2/m3–m5) の採否が §6.1 の記載どおり反映されているか: 小液滴の蒸発継続 (ṙ を max(r30, r_min) で評価) と単体 (k)、残差床ゲート §6 (a)–(c) の妥当性と最終確定バイナリ (d618c0dc) の run_0162–0167/run_0159 の床データ (§9 最上段)、compare_condfix --series-csv による Wys 判定、(h) の Δτ 非依存正規化、文書同期。§9 の主張を run (case/44 run_0162–0167・0159–0161, case/34 lim1d*, case/16 lim1d/lim0d, case/09 run_0063) と単体テストで突き合わせ、accepted にしてよいか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **NO-GO**

主要修正と最終場の比較値は確認できました。  
ただし、蒸発停止条件が一つ残り、定常性判定には **0.2%と20%の取り違え**があります。残差床ゲートも未収束を除外する根拠として不十分です。

## 指摘一覧

### 1. Major — `Q0=0` の液相は、依然として蒸発も消滅もしない

**根拠:** 新しい [condensationSource_d.cuh:333](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSource_d.cuh:333) は `q0 <= 1e-30` で全ソースを0にして戻ります。float版も同じです。一方、[condensationRealizability_d.cuh:53](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:53) は `g > g_rm` なら消滅処理前に戻ります。`g_rm=5e-7` です。

したがって、コード上の反例は **`S=0.5, g=1e-5, Q0=Q1=Q2=0`、輸送なし**です。

- ソースは0。
- 更新増分も0。
- 消滅クランプも無作用。
- 液相が永久に残る。

これは不整合モーメントですが、旧経路は [同ファイル:279](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSource_d.cuh:279) で消滅に向けて減少させていました。新経路の「`Q0=0` はクランプで消滅」という説明は成立しません。なお、単体(a)自身もこの種類の状態を入力していますが、Δτ不変性だけでは全ソース0を検出できません。

**対案:** 非平衡モデルへ渡す前に、この不整合を検出して明示的に拒否してください。平衡場からのrestartでは液滴モーメントの初期化を必須とし、両精度で検出試験を追加するべきです。

### 2. Major — 「全runが0.2%でALL STEADY」は再現しない

**根拠:** 集計用 [final_verify_report_v2.sh:6](/tmp/claude-1000/-home-sano-work-forge/96eabcfb-a7fa-40da-aae8-ef553208f9cd/scratchpad/final_verify_report_v2.sh:6) は `--drift 0.2 --osc 0.2` を指定しています。[判定実装:274](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:274) は比率を直接比較するため、これは **20%** です。

計画の0.2%に対応する `--drift 0.002 --osc 0.002` で、CSVの全報告列を再判定しました。

| run（`case/44.vitiated_air_wt/` 配下） | `x30_j4`の末尾変動 | 再実行VERDICT |
|---|---:|---|
| `run_0164_va3_M4.19_Lc8_noneq_inletTt_cell_lim1d_cfl2/` | 0.5304% | **DRIFTING** |
| `run_0167_va3_M4.19_Lc8_noneq_inletTt_cell_lim1d_cfl2_condFloat0/` | 0.5304% | **DRIFTING** |

両者の総合判定は **`OVERALL: NOT ALL STEADY`**。保存された `ALL STEADY` と、[plan §6](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:138) の条件は一致しません。

**対案:** 閾値を修正し、実行引数・各列の判定を含む全文を保存してください。`x30_j4` は [cond_series_csv.py:36](/home/sano/work/forge/case/44.vitiated_air_wt/cond_series_csv.py:36) の閾値を超えた最初の格子位置なので、今回の変動には格子単位の飛びが含まれます。交差位置を補間して再評価し、実際の移動と抽出の量子化を切り分けるべきです。

### 3. Major — 残差床の比率から「ソースは未収束ではない」とは結論できない

**根拠:** [plan §6(b)](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:135) は、CFL変更に対して流れとモーメントの残差が同方向に動くことを、同一の床機構・未収束の否定に使っています。しかし、比較している状態自体が異なるため、この推論は成立しません。`condFloat` 比較も、流れ・輸送・保存量更新のfloat32丸めを除去する対照ではありません。

末尾20%平均の記載値は概ね再現しました。ただし、その「床」は時間方向に安定していません。

| 対象 | 全4モーメントの末尾20%／直前20%平均 | `check_convergence.py` |
|---|---:|---|
| 上記 `run_0164` | **1.055–1.075** | 全モーメント **RISING** |
| `run_0165_va3_M4.19_Lc8_noneq_inletTt_cell_lim1d_cfl05/` | **0.768–0.819** | 全モーメント **falling** |

`run_0159–0167` は再実行でも全て **`NOT CONVERGED`** でした。また、§6には新しい床ゲートと「全モーメント3桁低下」の旧条件が併存し、どちらを満たせば完了なのか不明確です。

**対案:** CFL間の床比を合否根拠から外し、同一状態で輸送・ソース・合計残差を分離評価してください。輸送速度で正規化した全モーメント残差と、その時間変動に許容値を定め、必要なら流れを固定した対照計算で残差の発生源を確認するべきです。未達の間は§5.1の「回帰完了」を解除してください。

### 4. Minor — 小液滴の単体(k)はdouble実体だけを通している

**根拠:** [test_cond_limiter_steady.cu:332](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_limiter_steady.cu:332) は `condensation_source_d`、337行は `cond_realizability_clamp_d` を呼びます。別実装のfloatソース・float消滅クランプによる一連の試験はありません。

保存ログでは、修正した小液滴が **17ステップで消滅**しています。しかし、3回目レビューで要求された両精度の「ソース→更新→消滅」の確認には届いていません。

**対案:** 同じ反例を `condensation_source_f_d` と `cond_realizability_clamp_f_d` にも通し、液相と全モーメントの消滅を確認してください。

### 5. Minor — 文書同期と検証根拠の更新が未完了

**根拠:**

- [plan §5.1:128](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:128) は4巡目に言及しながら、根拠runを `0151–0161` と記載しています。
- [case/44 README:625](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:625) は旧セットも「検証正本（最終）」「残差床はmode 0と同一」としています。
- [followups F-cf7:59](/home/sano/work/forge/plans/active/condensation-followups.md:59) の本文に「未実装」「平衡緩和形も同時に見直す」が残り、確定スコープと矛盾します。
- `run_0160/0161` の最終場の `max|Δρ|/max|ρ|` は **2.475e-6**。§9の **1.4e-6**、参照する反復差 **1.8e-6**とは一致しません。`g`差 **5.626e-5** は記載に近い値でした。

**対案:** 現行の根拠をrun名・比較式・実測値で統一し、旧セットは履歴として区別してください。今回の未達条件を§5.1へ戻す必要があります。`methods/index.md` は既存文書の更新なので変更不要です。

## 確認できた修正・実測

- `max(r30,r_min)`、単体(h)のΔτ非依存正規化、Wysłouzilの抽出関数共用は実装済みです。
- case/44最終場の差は§9を再現しました。`0163−0162` は **g相対L1=1.638e-3、ΔT最大=0.182 K**、`0165−0164` は **7.579e-4、0.106 K**。最終凝縮域の `condLim=1`、補正G/Q=0も確認しました。
- `case/16.nozzle_wys/run_0470_limiter_regress_wys/lim1d,lim0d` のCSVは報告用抽出値と一致し、正しい0.2%で **`ALL STEADY`**。g場差 **8.776e-6** も再現しました。
- Arthurの `run_0105–0107` 内 `lim1d*` は全対象で **場比較28/28 PASS**。収束判定は別で、**`NOT CONVERGED`** です。
- `case/09.Taylor-Green/run_0063_periodic_seam_cond_source_final2/` のseam比は **0.99999976–1.00000024**。
- 対象135スナップショットの `VALUE/*` にNaN/Infはありませんでした。case/44対象runの保存されたメッシュ品質判定は **PASS**。恒久索引は [case/44 README](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:625) です。

## 推奨

**`in_progress`を維持し、蒸発停止条件、定常性判定、残差床の検証を修正してから再レビューしてください。**

レビュー範囲は `f1e40f41...a0ce5564`。CUDA単体の独立再実行はdriver/runtime不整合で検証不能でした。保存ログの `ALL PASS` は確認しましたが、独立再現済みとは扱っていません。ファイルは変更しておらず、指摘は **plan未反映**です。

指摘数: Critical 0 / Major 3 / Minor 2
