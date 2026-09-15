# codex レビュー: thermophysics-cea-mole-fraction-species (result)

- **plan**: [`plans/active/thermophysics-cea-mole-fraction-species.md`](../../plans/active/thermophysics-cea-mole-fraction-species.md)
- **stage**: `result` (diff base `4b357400~1`)
- **date**: 2026-09-16
- **commit**: `a73df90b` (feature/sern-design)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M1/m1
- **focus**: 5 回目 (result-4 の M1/M2/m1 反映後; §6.1 の採否表と §5.1 #28–#29)。前回までの再指摘は不要。実物で確認: M1 変換器の dtype 変換後検査 (tests/unit/test_convert_species_field_fail.py 21 件: 1e39 overflow 拒否・1e30 通過), M2 設計側の範囲外 NASA-9 (design/forge_design/gas/composition.py ResolvedSpeciesDB.species_cp_R/h_RT/s0_R, frozen.py; design/tests/run_gas_tests.py の外部 DB H2O Tlo=300 K でソルバ式独立移植と 6e-15、反転 T 3e-13 K), m1 §5.1 #4b と plans/README。accepted に移してよいか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **NO-GO**

前回の Major 2件は修正・試験結果を確認しました。  
ただし、`interp_field.py` が検査で受理したトレーサを転送せず、成功終了する別の欠陥を再現しました。  
これを修正するまで `accepted` への移動は推奨しません。

## 指摘一覧

### 1. Major — `Xi` で受け取ったトレーサが restart 先へ転送されない

**根拠:** [interp_field.py:92](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:92) は、保存量 `roXi` がなくても原始量 `Xi` があれば必須データセット検査を通します。しかし、[同:136](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:136) の転送処理は `roXi` だけを拾い、`Xi` から再構成しません。

同じ種順序・DB・`tracer: exhaust` を持つ8ノードの合成入力で、実際の書込み分岐までメモリ上で実行しました。

```text
元:   ρ=0.8、Xi=0.3、roXi なし
宛先: roXi/ρ=0.8

species signature OK
終了コード: 0
転送後の宛先 Xi: 0.79999995（期待値 0.3）
転送一覧に roXi なし
```

通常の保存済み run が壊れているという指摘ではありません。**ツール自身が有効と認めた入力で、流入元ラベルを黙って引き継がない**問題です。

**対案:** `roXi` があれば優先し、なければ `ρ×Xi` を転送配列へ追加してください。さらに、検査で要求した保存量が転送配列にも揃っていることを書込み前に確認し、この分岐の回帰試験を追加してください。

### 2. Minor — 完了状態と検証範囲に関する文書の不整合が残る

**根拠:**

- [plans/README.md:27](/home/sano/work/forge/plans/README.md:27) は、修正済みの dtype 検査をまだ「残」としています。plan §5.1 #28 と一致しません。#4b の回帰記録への更新は確認できました。
- [case/44 README:645](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:645) は、凝縮場の差を丸めに帰属させ「固定点の差ではない」と断定しています。対象 `run_0195`–`0197`／`0199` の再判定は **`NOT CONVERGED (stalled/plateau)`**。準定常性と比較許容差だけでは、固定点の同一性や差の原因までは証明できません。

**対案:** 索引を現在の完了状態へ同期し、case/44 は「報告量は比較許容差内。差の原因・固定点の同一性は未確定」としてください。今回の Major は §5.1 に未解決として登録してください。

## 再検証で確認したこと

- **前回 M1:** 変換器 **21/21 PASS**。書込み分岐でも `1e39` は拒否され宛先不変、`1e30` は有限の `float32` として保存されました。
- **前回 M2:** `run_gas_tests.py` は **ALL PASS**。外部 DB の `H2O Tlo=300 K` で独立ソルバ式との差 **6.4e-15**、反転温度差 **3.1e-13 K**。低温試験 **5/5 PASS**、SERN 試験も **ALL PASS**。
- **小型回帰:**  
  `case/16.nozzle_wys/run_0471_tp5_node_euler_va3/`  
  `case/16.nozzle_wys/run_0472_tp2_node_euler_va3_lumped/`  
  両方 **`PASS (converged)`**、`machmax,pmax` は **`OVERALL: ALL STEADY`**、保存済みメッシュ品質は **`VERDICT: PASS`**。[run 一覧](/home/sano/work/forge/case/16.nozzle_wys/README.md:373)
- **未収束回帰:** [case/46 の run 一覧](/home/sano/work/forge/case/46.sern_design/README.md)の `run_0100`–`0103` は **`NOT CONVERGED`**、力系列は指定の0.1%閾値で **`ALL STEADY`**。case/44 の上記4 run も未収束ですが、凝縮系列は0.2%閾値で **`ALL STEADY`**。確認対象の全保存スナップショットに NaN/Inf はありませんでした。

## 推奨

**`active/in_progress` を維持し、指摘1の転送欠落を修正・再検証した後、文書を同期して移動してください。**

対象は指定 diff、`HEAD=a73df90b`。試験の一時入出力はメモリへ置換し、新規ビルド・CFD投入は行っていません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 1 / Minor 1
