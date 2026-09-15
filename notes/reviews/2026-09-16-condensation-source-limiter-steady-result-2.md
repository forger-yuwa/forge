# codex レビュー: condensation-source-limiter-steady (result)

- **plan**: [`plans/active/condensation-source-limiter-steady.md`](../../plans/active/condensation-source-limiter-steady.md)
- **stage**: `result` (diff base `f1e40f41`)
- **date**: 2026-09-16
- **commit**: `5a3b1369` (feature/sern-design)
- **codex**: effort `high`, 7.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m3
- **focus**: 2 回目 result レビュー (NO-GO M1–M3/m1–m2) の採否が §6.1 の記載どおり反映されているか: 残差ゲートの床基準 (cell mode 1 run_0155 vs mode 0 run_0159 の床)、interp_field の roY 新設と run_0146 の訂正、cell × condFloat 0 (run_0157/0158)、単体 (h) の実現可能性クランプ込み判定と (j) corQ、旧経路の診断リセット、Wys の正式 --series-csv 判定。最終バイナリ (6b500375) での検証セット (case/44 run_0151–0161, case/34 lim1c*, case/16 lim1c/lim0c, case/09 run_0062) が §9 の主張を裏付けるか、accepted にしてよいか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **NO-GO**

残差から Δτ を外す主要実装と、多くの比較値は確認できました。  
ただし、最終修正に蒸発を停止させる欠陥があり、変更後の残差床ゲートも未達です。`accepted` への移動は推奨しません。

## 指摘一覧

### 1. Major — 小液滴が蒸発も消滅もしない停止条件がある

**根拠:** [condensationSource_d.cuh:335](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSource_d.cuh:335) は `r30 < 2*rmin` で全ソースを0にして戻ります。[float版:244](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:244) も同じです。

一方、消滅クランプは [condensationRealizability_d.cuh:53](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:53) で `g > g_rm` なら戻ります。`g_rm` は [condensationTransport_d.cu:206](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:206) の **5e-7** です。

具体的な反例は、次の整合した単分散状態です。

- H2O、`T=250 K`、`ρ=0.1 kg/m³`、`S=0.5`
- `g=1e-5`、`r30=1.5e-9 m`、既定 `rmin=1e-9 m`
- 保存モーメント `q0=7.0965e16`、`q1=1.0645e8`、`q2=0.15967`

コードの物性式から求めた蒸発速度は約 **−5.61e-5 m/s** ですが、新経路はソース0、消滅クランプも無作用になります。輸送のない状態では液相が残り続けます。これは不整合な「数値塵」に限定されません。

**対案:** 消滅クランプの許容質量まで、Δτ 非依存の蒸発速度で減少を継続させてください。更新上限は更新クランプで守り、この反例について両精度の「ソース→更新→消滅」試験を追加すべきです。

### 2. Major — 「旧経路と同じ残差床」は全モーメントで成立していない

**根拠:** [plan §6:133](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:133) は各モーメントについて「3桁低下、または旧経路の床±30 %」を要求しています。

対象は次の2 runです。

- `case/44.vitiated_air_wt/run_0155_va3_M4.19_Lc8_noneq_inletTt_cell_lim1c_cfl2/`
- `case/44.vitiated_air_wt/run_0159_va3_M4.19_Lc8_noneq_inletTt_cell_lim0c_cfl2/`

`residual_history.csv` の**末尾20 %平均**は以下でした。

| 残差 | mode 1：0155 | mode 0：0159 | 比 |
|---|---:|---:|---:|
| `rms_rog_0` | 4.031e-8 | 1.170e-8 | **3.44** |
| `rms_roQ2_0` | 1.471e-4 | 8.285e-5 | **1.78** |
| `rms_roQ1_0` | 1.590e3 | 1.286e3 | 1.24 |
| `rms_roQ0_0` | 2.051e10 | 2.140e10 | 0.959 |

`rog` と `roQ2` が不合格です。`0155` の低下桁数も各々 **0.4／0.3桁**で、代替条件を満たしません。

`check_convergence.py` の再実行は `run_0151–0161` 全て **`NOT CONVERGED`**。特に `0156` は全モーメントが `RISING` 判定でした。`--drift 0.002 --osc 0.002` による `0151–0159` の **`ALL STEADY`** は再現しましたが、残差条件の代替にはなりません。

**対案:** §5.1の回帰完了を解除し、全モーメントの低下桁数・床比較・トレンドを表にしてください。`rog`／`roQ2` の床差を切り分け、ゲートを満たしてから完了扱いにすべきです。

### 3. Minor — Wysłouzil の正式CSVが実際の出口を評価していない

**根拠:** [wys_cond_series.py:11](/home/sano/work/forge/case/16.nozzle_wys/wys_cond_series.py:11) は全領域共通の `|y|` 閾値で中心線候補を選びます。このメッシュでは下流の列が除外され、`gc[-1]` は **x=38.34 mm** の値です。実際の出口は **95 mm** です。

`case/16.nozzle_wys/run_0470_limiter_regress_wys/lim1c/res_48000.h5` では、

| 量 | 保存CSV | 報告用抽出関数による出口 |
|---|---:|---:|
| `g_exit` | 0.00828519 | **0.01085136** |
| `M_exit` | 1.29551 | **1.63251** |

また、CSVには報告対象の壁圧偏差がありません。

報告用の [compare_condfix.py:66](/home/sano/work/forge/case/16.nozzle_wys/compare_condfix.py:66) で全時系列を抽出し直し、正式ツールに渡したところ、`lim1c/lim0c` とも出口量・壁圧偏差を含め **`ALL STEADY`（0.2 %）**でした。したがって、今回は回帰結論を覆す問題ではありません。

**対案:** 報告と検証で抽出関数を共用し、正しい出口・補間onset・壁圧偏差を含むCSVとVERDICTを再保存してください。

### 4. Minor — 単体(h)の「未加工残差」判定が Δτ に依存する

**根拠:** [test_cond_limiter_steady.cu:283](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_limiter_steady.cu:283) が判定しているのは、

`|res| × Δτ / (V × |状態|)`

です。これは相対更新量であり、同じ残差でも Δτ を小さくすれば合格しやすくなります。

実現可能性クランプの呼び出し、補正0の確認、(j)の `corQ=1` 判定は追加されています。ただし、残差判定そのものは本planの問題意識と整合しません。

**対案:** この1セルモデルでは、輸送速度を基準に `|res|/(β V |状態|)` など、Δτ を含まない正規化残差を全4成分で判定してください。

### 5. Minor — 文書と残作業表に訂正漏れが残る

**根拠:**

- [methods/condensation.md:839](/home/sano/work/forge/methods/condensation.md:839) は依然「定常解の一致」と記載し、cellにもnodeの残差床の説明を適用しています。
- [case/44 README:617](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:617) に `0142+0146` の「72000 step」記述が残ります。
- [plan §5.1:128](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:128) と§6.1は `run_0147–0150` を根拠にしていますが、このcheckoutには該当ディレクトリがありません。
- [plans/README.md:27](/home/sano/work/forge/plans/README.md:27) は `draft` のままです。

**対案:** 実在する最終検証セットへ根拠を統一し、未達条件を§5.1へ戻してください。旧runを廃棄・置換したなら台帳に記録すべきです。`methods/index.md` は既存文書の更新なので変更不要です。

## 確認できた修正・実測

- `roY*` 新設と旧経路の診断リセットは実装済み。最終cell runの入力→`res_0.h5` では `roY1`／`rog_0` が保持されています。
- [case/44 最終セットの台帳](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:542) の場差を再現しました。`0153−0151` のg相対L1は **1.7504e-3**、`0156−0155` は **7.7575e-4**。cell両精度の比較値も記載どおりです。これは**未収束の最終保存場の比較**です。
- `case/09.Taylor-Green/run_0062_periodic_seam_cond_source_final/` のseam比は **0.99999976～1.00000024**でした。
- Arthurの場比較は記載どおり。ただし `run_0105_limiter_regress_n2_cell/lim1c_r3` は **27/28 FAIL**、他の対象は28/28 PASSです。これらは収束PASSではありません。
- 対象のcase/44非平衡・Arthur・Wysłouzil保存場にNaN/Infはなく、最終凝縮域の補正診断は0でした。

## 推奨

**`in_progress` を維持してください。** 優先順は、蒸発停止条件の修正、全モーメントの残差ゲート達成、検証CSV・単体判定・文書の同期です。その後に再レビューすることを推奨します。

レビュー範囲は `f1e40f41...5a3b1369`。CUDA単体の再実行はdriver/runtime不整合で検証不能でした。保存ログの `ALL PASS` は確認しましたが、独立再現済みとは扱っていません。ファイルは変更しておらず、指摘は **plan未反映**です。

指摘数: Critical 0 / Major 2 / Minor 3
