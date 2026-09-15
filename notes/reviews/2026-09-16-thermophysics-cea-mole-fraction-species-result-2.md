# codex レビュー: thermophysics-cea-mole-fraction-species (result)

- **plan**: [`plans/active/thermophysics-cea-mole-fraction-species.md`](../../plans/active/thermophysics-cea-mole-fraction-species.md)
- **stage**: `result` (diff base `4b357400~1`)
- **date**: 2026-09-16
- **commit**: `cb3c1267` (feature/sern-design)
- **codex**: effort `high`, 7.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **focus**: 2 回目 (result-1 の M1–M10/m1 反映後; §6.1 の採否表と §5.1 #9–#19、§9 の 2026-09-16 項)。前回指摘の再指摘は不要で、各対応が実際に指摘を解消しているかを実物で確認してほしい: M1 引用符付き species (case/44 run_0193 再生成, run_sern_frozen_gas_tests の safe_load 試験), M2/M4 convert_species_field.py (二相 EOS 反転を run_0196 nozzle.h5 で, stream lump の conserve/reinit を case/46 run_0100/0101 の meta で), M3 restart 照合 (runner_sern.check_species_compatible, forge_species.species_signature / interp_field), M5 periodicNode の roXi (case/09 run_0065/0066), M6/M7 solverConfig (tests/unit/test_solver_config_species.cpp), M8 トレーサ規則と composition.exhaust_fraction, M9 forge_species.vapor_array / axis_csv_va.py, M10 小型 PASS ケース case/16 run_0471/0472・トレーサ差の切り分け case/46 run_0104・作動点変更 run_0105・CEA ゲート改定 (§6 (h) 5e-6)・SERN ξ ゲート改定 (§6)。cell はユーザ指示で対象外。accepted に移してよいか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **NO-GO**

前回の修正の多くは実測で確認でき、小型5種 node ケースの `PASS` も再現しました。  
ただし、作動点変更・restart 照合・変換器の異常値検査に未解消の欠落があります。  
現状のまま `done`／`accepted` へ移すことは推奨しません。

## 指摘一覧

### 1. Major — 作動点変更の「組成再初期化」が実装・検証されていない

**根拠:** [convert_species_field.py:259](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:259) の `transfer_reinit()` は、元の実種質量を行き先へ配分します。行き先が `full` では、各実種の配分先は一つなので、**新作動点の入口組成へ変わりません**。

実際の問題 YAML から m6_on／m10_on の配置を生成し、純排気 `ξ=1` をメモリ上で変換すると、次の結果でした。

| 質量分率 | m6_on 元値 | `reinit` 出力 | m10_on 目標値 |
|---|---:|---:|---:|
| H₂O | 0.24110912 | 0.24110912 | 0.24881767 |
| H₂ | 0.00103232 | 0.00103232 | 0.01383296 |

また、実 run の `run_0102_species_regress_euler_lumped_24k/` → `run_0105_species_regress_euler_lumped_m10_warm/` を変換器の `--mode reinit --dry-run` に通すと、

```text
実種質量の射影損失 1.633e-03 kg/m³ > 許容 2.2e-04
FAILED
```

でした。`run_0105` が検証したのは、変換器を使わず組成分率を持ち越す [warm_from_run:539](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:539) の別経路です。[plan:218](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:218) の完了根拠にはなりません。

**対案:** 保存的な配置変換と作動点変更を分離し、後者では保存済み `ξ` と**目標**入口ベクトルから組成を再構成してください。`full`／`lumped+keep` でも目標組成・ΣY・温度の契約を検証し、§5.1 #12/#18 を未完了へ戻すべきです。

### 2. Major — restart の新しい署名でも、不整合を通してしまう

**根拠:** [forge_species.py:155](/home/sano/work/forge/solver_density_cuda/tools/forge_species.py:155) の署名には NASA-9 係数・温度区切り・トレーサ設定がありません。メタデータのハッシュも、DB 内容の照合を代替しません。

書込みなしの再現試験では、`run_0100` の `EXH.nasa9_low[2]` だけをメモリ上で `+1` しても、`compare_signatures()` は不一致なしの `[]` を返しました。この変更は低温域の cp を約 **340.35 J/kg/K** 変えます。

設計側の [runner_sern.py:436](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:436) も実 config／DB ではなくメタデータを参照します。実際の `run_0100` と、config でトレーサを追加した `run_0104` は、両照合関数を通過しました。さらに、両側にメタデータがなければ、`thermalMethod` を確認せず戻ります。[同:452](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:452)

**対案:** 実際の config と解決済み DB から共通署名を作り、係数・温度区切り・datum・種順序・トレーサを照合してください。メタデータとの矛盾と、TP の照合不能は拒否し、必要な保存量の存在もコピー前に検査してください。§5.1 #11 は完了扱いにできません。

### 3. Major — 変換器の「書込み前 hard fail」は NaN を検出しない

**根拠:** [convert_species_field.py:457](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:457) は温度差を `dT.max() > tolerance` で判定しますが、NaN の比較は偽になります。入力・出力の有限性検査もありません。

`run_0196` の配置を使い、**メモリ上だけに作った1セル入力**で `roe=NaN` を与えると、実際に次の結果になりました。

```text
source T: nan..nan K
check T 保存: max |ΔT| 全セル nan K
--dry-run: 検査 OK
```

通常実行なら、その後の書込みへ進みます。加えて [T_from_e:115](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:115) は、反復上限到達時にも EOS 残差を確認せず温度を返します。

**対案:** `ρ>0`、全保存量・組成・温度の有限性、組成の非負性、EOS 反転残差を明示的に検査してください。NaN、範囲端への張り付き、反転未収束は書込み前に拒否し、失敗系試験を追加してください。

### 4. Minor — M8 修正後の仕様と完了条件が文書に揃っていない

**根拠:** [methods/thermophysics.md:266](/home/sano/work/forge/methods/thermophysics.md:266) は依然「`lumped` では ξ=Y_EXH なので不要」と記述しています。[methods/design/overview.md:783](/home/sano/work/forge/methods/design/overview.md:783) も同様です。実装された `lumped+keep` のトレーサ規則と矛盾します。

また、[plan:219](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:219) には旧 `require_residual_pass=True` と代替ゲートが併記され、未実施の粘性試験も必須検証として残っています。

**対案:** 排気率の説明を「純流入元ラベル種がある場合／ない場合」に統一してください。旧ゲートは履歴として明示し、現行ゲートを一本化すること。`F-sp1` と依存する粘性試験の延期範囲を§6・§8にも反映し、今回の未解消項目を§5.1へ戻してください。

## 再検証で確認できた改善

以下は解消済みの根拠として認めます。

- **M1/M8/M9:** 引用符付き生成 config の再読込、`lumped+keep` のトレーサ有効化、dry `run_0203` の `vapor_array=Y0` を確認。
- **M2/M4の同一作動点変換:** `case/44.vitiated_air_wt/run_0196_va3_M4.19_Lc8_noneq_inletTt_full5/` の湿潤6365セルで、入力反転温度と `res_0` の差は最大 **1.94e-4 K**。`run_0100→0101` の保存的展開・`roXi` 生成、逆方向の再配分も記録値を再現。
- **M5:** `case/09.Taylor-Green/run_0065_periodic_tracer_seam/` と `run_0066_periodic_tracer_interior/` は周期対の ΔXi=0、積分変化 **3.67e-8／4.96e-8**。根拠は[run一覧](/home/sano/work/forge/case/09.Taylor-Green/README.md:35)。
- **M6/M7:** 既存 C++ テストバイナリを既存入力の読込モードで再実行し、15条件すべて期待どおり。
- **M10:** `case/16.nozzle_wys/run_0471_tp5_node_euler_va3/`／`run_0472_tp2_node_euler_va3_lumped/` は **`PASS (converged)`**、準定常判定も **`ALL STEADY`**。最大 |ΣY−1| は **7.92e-8／8.57e-8**。ただし「全列4.2–4.8桁」は不正確で、`rms_roUy` は **3.3桁**です。[run一覧](/home/sano/work/forge/case/16.nozzle_wys/README.md:373)
- **改定ゲート:** `run_0104` の同一 run 内トレーサ差は最大 **1.809e-4**／平均 **4.612e-7**。CEA 物性差は最大 **1.11017e-6** で、改定値 `5e-6` 内でした。

case/44 の凝縮量と case/46 の力は、指定閾値で **`ALL STEADY`** を再現しました。ただし残差判定は引き続き **`NOT CONVERGED`** です。これは報告量の準定常比較の根拠であり、収束解一致の証明ではありません。

## 推奨

**`active/in_progress` を維持し、上の Major 3件を修正・失敗系を含めて検証してから、result レビューを再実施してください。**

対象は `HEAD=cb3c1267`。ファイル変更なし、**plan 未反映**です。Python 単体は書込みブロックを除いて再実行し、新規ビルド・CFD 計算は行っていません。

指摘数: Critical 0 / Major 3 / Minor 1
