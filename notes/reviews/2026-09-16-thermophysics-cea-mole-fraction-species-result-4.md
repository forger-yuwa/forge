# codex レビュー: thermophysics-cea-mole-fraction-species (result)

- **plan**: [`plans/active/thermophysics-cea-mole-fraction-species.md`](../../plans/active/thermophysics-cea-mole-fraction-species.md)
- **stage**: `result` (diff base `4b357400~1`)
- **date**: 2026-09-16
- **commit**: `3fbb9eb5` (feature/sern-design)
- **codex**: effort `high`, 7.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **focus**: 4 回目 (result-3 の M1–M3/m1 反映後; §6.1 の採否表と §5.1 #24–#27、§9 最新項)。前回までの再指摘は不要。実物で確認: M1 変換器の config 駆動 (forge_species.species_signature のトレーサ矛盾拒否、required_conserved; case/46 run_0104 → 自身で roXi 保持、meta enabled false の複製で拒否), M2 total_quantities._TPGas の範囲外処理 (tests/unit/test_tpgas_lowT.py 5 件: .cuh 式の独立移植と 2.7e-15、反例 100/150 K で T 差 4e-13 K), M3 全書込配列の有限性 (tests/unit/test_convert_species_field_fail.py 19 件), m1 docs (procedures/solver-settings.md の ξ, §9 の桁数訂正, plans/README)。accepted に移してよいか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **NO-GO**

重点のトレーサ保持、低温 EOS、既存24試験は修正を確認しました。  
ただし、書込み時の Inf 発生と、設計側の外部 DB 温度域の無視を再現しました。  
この2件を直すまで `accepted` への移動は推奨しません。

## 指摘一覧

### 1. Major — 有限性検査後の型変換で Inf を書き込める

**根拠:** [convert_species_field.py:616](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:616) は float64 配列を検査しますが、[同:649](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:649) の書込先 dtype への変換後は検査しません。

既存の正常試験入力に、有限な float64 の `roOmega=1e39` を追加し、float32 の宛先へ変換しました。ファイル入出力をメモリに置き換え、**実際の書込み分岐まで**実行した結果です。

```text
rc=0
SUMMARY: all checks passed
書込後 VALUE/roOmega: [inf inf inf inf inf inf inf inf]
```

通常 run で発生したという主張ではありません。ただし、異常入力を拒否する契約と、§5.1 #26 の「全書込配列の有限性」は満たしていません。

**対案:** 宛先の dtype を先に読み、全配列を最終 dtype に変換してから有限性・正値性を検査してください。既存データセットの削除は全検査後に行い、`--dry-run` も同じ検査を通してください。試験には型変換による overflow を追加すべきです。

### 2. Major — 設計側は外部 DB の `Tlo`／`Thi` を反映していない

**根拠:** [FrozenGas._per_species:61](/home/sano/work/forge/design/forge_design/gas/frozen.py:61) は、解決済み DB の `Tmid` だけを使い、下限には固定の `T_FLOOR` を使っています。`Thi` の処理もありません。[ResolvedSpeciesDB.cp_mass/h_mass:145](/home/sano/work/forge/design/forge_design/gas/composition.py:145) も同様です。

`full` が受理する配置で、N2/H2O の質量分率を `[0.8,0.2]`、H2O の `Tlo=300 K`、datum を `298.15 K` としました。`FrozenGas.e_sens()` の生成エネルギーを、`thermo_d.cuh` の独立移植式で反転すると次になります。

| 指定温度 | ソルバ式での反転温度 |
|---:|---:|
| 100 K | **100.498255 K** |
| 150 K | **150.341829 K** |
| 250 K | **250.057419 K** |

SERN の [warm_from_run:588](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:588) は、この `e_sens()` で `roe` を生成します。したがって、§4.1 の「同じ解決済み DB を設計・IC・CFD 全経路で使う」は、温度域処理まで含めると未達です。今回直った `_TPGas` 自体への再指摘ではありません。

**対案:** 設計側も種別の `Tlo/Tmid/Thi` に従う cp 固定・h 線形・s° 対数外挿へ統一してください。非標準温度域を持つ外部 DB で、設計側の IC／warm restart が作るエネルギーを独立したソルバ式で検査する試験が必要です。

### 3. Minor — 残作業表と索引が現在の完了記録に揃っていない

**根拠:** [plan §5.1 #4b:180](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:180) は、後段で結果を記録済みの case/46 回帰を依然「残」としています。[plans/README.md:27](/home/sano/work/forge/plans/README.md:27) も、解消を確認した旧トレーサ問題を残作業として列挙しています。

`solver-settings.md` の ξ 説明と、§9 の残差低下桁数の訂正は確認できました。

**対案:** 完了済み回帰を整理し、今回の2件を§5.1へ登録してください。延期済みの `F-sp1` は完了対象外であることと追跡先を維持し、索引を同期してください。

## 検証で確認できたこと

- **試験:** ファイル入出力をメモリへ置き換え、失敗系 **19/19 PASS**、低温 **5/5 PASS**。物性差 `2.66e-15`、低温反例の温度差最大 `4.26e-13 K` を再現。
- **トレーサ:** `case/46.sern_design/run_0104_species_regress_euler_lumped_tracer/` の `res_2000.h5` を同配置へ変換し、`roXi` 差 **0**。メタだけ `enabled: false` にすると書込み前に拒否。
- **小型回帰:**  
  `case/16.nozzle_wys/run_0471_tp5_node_euler_va3/`  
  `case/16.nozzle_wys/run_0472_tp2_node_euler_va3_lumped/`  
  再判定は両方 **`PASS (converged)`**、`machmax,pmax` は **`OVERALL: ALL STEADY`**。保存済みメッシュ品質は **`VERDICT: PASS`**。場差も記載値を再現しました。[run 一覧](/home/sano/work/forge/case/16.nozzle_wys/README.md:373)
- **未収束回帰:** case/46 の `run_0100`–`0103`、case/44 の `run_0195`–`0197`／`0199` は **`NOT CONVERGED (stalled/plateau)`**。指定の力系列0.1%、凝縮系列0.2%は **`OVERALL: ALL STEADY`**。収束解の一致を示すものではありません。対象の全保存スナップショットで非有限値は検出しませんでした。

## 推奨

**`active/in_progress` を維持し、Major 2件を修正・再検証してから移動を判断してください。**

対象は指定 diff、`HEAD=3fbb9eb5`。新規ビルド・CFD投入は行っていません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 1
