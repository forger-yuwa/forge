# codex レビュー: thermophysics-cea-mole-fraction-species (result)

- **plan**: [`plans/active/thermophysics-cea-mole-fraction-species.md`](../../plans/active/thermophysics-cea-mole-fraction-species.md)
- **stage**: `result` (diff base `4b357400~1`)
- **date**: 2026-09-16
- **commit**: `219f3a2f` (feature/sern-design)
- **codex**: effort `high`, 7.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **focus**: 3 回目 (result-2 の M1–M3/m1 反映後; §6.1 の採否表と §5.1 #20–#23、§9 の最新項)。前回までの再指摘は不要。実物で確認してほしい: M1 組成再初期化 (composition.reinit_transport_vector, runner_sern.warm_from_run, convert_species_field.py --mode reinit; case/46 run_0105 [lumped] / run_0106 [full 11 種, tp_species.species で順序固定] の warm start と 2000 step 適応段; 単体 run_sern_frozen_gas_tests の m10_on 目標組成 0.24881767/0.01383296), M2 署名 (runner_sern._species_signature/check_species_compatible と forge_species.species_signature/compare_signatures/required_conserved; 単体の係数摂動・tracer 差・順序矛盾・config 無し), M3 変換器の検査 (tests/unit/test_convert_species_field_fail.py 12 件), m1 docs (methods/thermophysics.md, methods/design/overview.md の ξ, §6 ゲート一本化, 粘性試験の延期 §6/§8)。accepted に移してよいか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **NO-GO**

目標組成への再初期化と、通常 restart の係数照合は改善を確認しました。  
ただし、変換器にはトレーサ消失・低温での温度保存違反・乱流保存量の NaN 見逃しが残っています。  
`accepted` への移動は、この3件の修正後としてください。

## 指摘一覧

### 1. Major — 変換器が config で有効な `roXi` を削除する

**根拠:** [convert_species_field.py:215](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:215) は、トレーサの有無を config よりメタデータ優先で決定します。`species_signature()`／`required_conserved()` を使っておらず、両者の矛盾を拒否しません。

実在する `case/46.sern_design/run_0104_species_regress_euler_lumped_tracer/` は、[solverConfig.yaml:4](/home/sano/work/forge/case/46.sern_design/run_0104_species_regress_euler_lumped_tracer/solverConfig.yaml:4) が `tracer: exhaust`、[species_meta.yaml:9](/home/sano/work/forge/case/46.sern_design/run_0104_species_regress_euler_lumped_tracer/species_meta.yaml:9) が `enabled: false` です。

同 run の `res_2000.h5` → `sern.h5` を `--dry-run` で検査すると、実際に次を返しました。

```text
note: source roXi is dropped (destination has no tracer)
--dry-run: all checks passed
```

通常実行では [同:603](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:603) で既存 `roXi` を削除し、再作成しません。次の solver 起動で排気率がゼロ初期化されます。

**対案:** config から必須保存量を決定し、メタデータとのトレーサ矛盾は書込み前に拒否してください。配置変換で許す種・DBの差と、各 run 内の不整合を分けて検査する必要があります。§5.1 #21 は未完了へ戻してください。

### 2. Major — 温度保存検査がソルバと異なる低温 EOS を使う

**根拠:** 変換器が使う [total_quantities.py:41](/home/sano/work/forge/solver_density_cuda/tools/total_quantities.py:41) の `_TPGas` は、`Tlo`／`Thi` を無視して NASA-9 多項式を評価します。一方、ソルバの [thermo_d.cuh:177](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:177) は、範囲外で cp を固定し、h を線形外挿します。

メモリ上の1セル入力で、同じ内蔵 DB・datum を使い、`reinit` で `[N2,H2O]` の質量分率を `[0.95,0.05]` → `[0.8,0.2]` に変更しました。

| 元のソルバ温度 | 変換器の温度差判定 | 変換後エネルギーをソルバの EOS で反転 |
|---|---:|---:|
| 100 K | 約 `1.3e-13 K`、合格 | **101.38117 K** |
| 150 K | `0 K`、合格 | **150.10185 K** |

いずれも変換器本体は `all checks passed` を返しますが、実際の温度差は既定許容 `0.05 K` を超えます。同じ誤った物性関数で変換と検査を行うため、検査が問題を隠しています。

**対案:** 種ごとの `Tlo`／`Thi`、cp 固定、h 線形外挿をソルバと一致させてください。低温入力の期待エネルギーは独立した参照実装で作り、組成変更後の温度を検証する試験を追加してください。

### 3. Major — 有限性検査から `roK`／`roOmega` が漏れる

**根拠:** [convert_species_field.py:385](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:385) で乱流保存量を読みますが、[同:398](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:398) 以降の有限性検査には含めず、[同:601](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:601) でそのまま書き込みます。

既存の正常試験入力に、メモリ上で `roK=NaN` を追加すると、

```text
rc=0
--dry-run: all checks passed
```

でした。前回問題だった `roe=NaN` は拒否されますが、「全保存量の有限性」を保証する実装にはなっていません。

**対案:** 書き込む配列を一括管理し、全配列の有限性を検査してください。`roK`／`roOmega` の NaN・Inf を失敗系に追加し、§5.1 #22 に残してください。

### 4. Minor — 文書同期に取り残しがある

**根拠:** methods 2文書の ξ 定義、§6 の現行ゲート、粘性試験の延期は修正済みです。ただし、

- [procedures/solver-settings.md:185](/home/sano/work/forge/procedures/solver-settings.md:185) は依然「`lumped` では ξ=Y_EXH」と一般化しています。
- [plan:256](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:256) の変更ログには旧「4.2–4.8桁」が残ります。実測の `rms_roUy` は **3.3桁**です。
- [plans/README.md:27](/home/sano/work/forge/plans/README.md:27) は result 1回目の対応中という記載です。

**対案:** 手順の ξ 定義を methods に揃え、変更ログには訂正を明記し、索引を今回のレビュー段階に更新してください。

## 修正を確認できた点

- **再初期化:**  
  `case/46.sern_design/run_0105_species_regress_euler_lumped_m10_warm/` と  
  `case/46.sern_design/run_0106_species_regress_euler_full_m10_warm/` の `sern.h5` は、再初期化式との最大差がそれぞれ **`1.08e-7`／`7.40e-8`**。full の純排気側は H₂O **0.24881765**、H₂ **0.013832964** でした。入口目標を使う修正は確認できました。[run一覧](/home/sano/work/forge/case/46.sern_design/README.md:144)
- **2000 step 適応段:** 両 run の再判定は **`NOT CONVERGED (still converging)`**。`rms_ro` は4.5桁低下しますが、full の `rms_roXi` は2.6桁です。`res_0`／`res_2000` の非有限値はありません。適応段の検証として扱うことは妥当です。
- **通常 restart の署名:** 設計側・forge ツール側とも、NASA-9 係数 `+1` とトレーサ設定差を検出しました。
- **既存12件:** ファイル入出力をメモリへ置き換え、変換器本体を実行して **12/12 PASS**。上記の追加反例は、その試験範囲外です。
- **小型ケース:**  
  `case/16.nozzle_wys/run_0471_tp5_node_euler_va3/`／`run_0472_tp2_node_euler_va3_lumped/` は再実行した判定で **`PASS (converged)`**、`machmax,pmax` は **`OVERALL: ALL STEADY`**。保存済みメッシュ品質も **`VERDICT: PASS`** です。[run一覧](/home/sano/work/forge/case/16.nozzle_wys/README.md:373)
- **準定常性:** case/46 の `run_0100`–`0103` の力系列は0.1%、case/44 の `run_0195`–`0197`／`0199` の凝縮系列は0.2%で、`check_quasisteady` の **`STEADY`** を再現しました。未収束ケースについて、収束解の一致を証明するものではありません。

## 推奨

**`active/in_progress` を維持し、上記 Major 3件を修正・再検証してから再レビューしてください。**

対象は指定 diff、`HEAD=219f3a2f`。新規ビルド・CFD投入は行っていません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
