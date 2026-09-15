# codex レビュー: thermophysics-cea-mole-fraction-species (result)

- **plan**: [`plans/active/thermophysics-cea-mole-fraction-species.md`](../../plans/active/thermophysics-cea-mole-fraction-species.md)
- **stage**: `result` (diff base `4b357400~1`)
- **date**: 2026-09-16
- **commit**: `c1deeb96` (feature/sern-design)
- **codex**: effort `high`, 9.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M10/m1
- **focus**: 実装 diff (設計チェーン: gas/composition.py の ResolvedSpeciesDB・mole_to_mass・parse_tp_species/resolve_species_layout・species_db_yaml/species_meta、probdef/runner_axismach/runner_sern{,3d} の統一スキーマ化と restart の全種化; forge 本体: input/speciesDB.{hpp,cpp} の host DB 解決、boundaryCond.cpp の X{s}/Y{s} 検証換算、condensationSpecies、rms_roY 列、tracerTransport_d の roXi、tools/convert_species_field.py・forge_species.py・gen_inlet_profile --X・cea --check・interp_field 照合・ParaView 配列解決) が plan §2/§4 と一致するか。検証 (plan §9 の 2026-09-16 回帰 2 項目; case/44 README「species 統一スキーマ」表 run_0190–0199/0202/0203; case/46 README run_0100–0103) が §6 のゲートを満たすか。cell はユーザ指示で対象外 (§6 改定)。既知の制限: トレーサ拡散 0 (F-sp1)、SERN node Euler は残差 plateau で check_convergence PASS 不可 (§6 に読み替えを記載)。accepted に移してよいか。単体: design/tests/run_gas_tests.py, run_sern_frozen_gas_tests.py, solver_density_cuda/tests/unit/test_species_db_host.cpp

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **NO-GO**

力・凝縮量の時系列は再判定でも `ALL STEADY` ですが、受入ゲート全体は満たしていません。  
実際の入力で再現する後処理・restart の不具合と、`roXi` の周期境界・非定常処理の欠落があります。  
現状を `done` として `accepted` に移すことは推奨しません。

## 指摘一覧

Critical はありません。以下は移行前に修正すべき Major です。

### 1. Major — SERN `full` の `NO` が Python 側で `False` になる

**根拠:** [runner_sern.py:223](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:223) は種名を無引用で連結します。実際の [run_0101 の solverConfig.yaml:4](/home/sano/work/forge/case/46.sern_design/run_0101_species_regress_euler_full/solverConfig.yaml:4) を `yaml.safe_load()` すると、species index 6 は文字列 `NO` ではなく `False` になります。

同 run に `species_info()` を実行すると、次の例外を再現しました。

```text
KeyError: species False: MW が species_db.yaml にも内蔵表にも無い
```

`forge_species` を使う後処理・種変換が失敗します。ノズル側も [runner_axismach.py:124](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:124) に同じ生成方法があります。

**対案:** 全生成経路で種名を文字列として安全に YAML 出力し、生成した config を Python 側の実際の読込関数で再読込する試験を追加してください。

### 2. Major — 湿潤 input HDF5 の温度反転が二相 EOS を無視する

**根拠:** [convert_species_field.py:195](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:195) は、`T` のない入力を乾き気相の `e(T)` で反転します。一方、実装の二相 EOS は `e_allvap + g(Rw T−L)` です。[condensationEOS_d.cuh:8](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:8)

実際の `case/44.vitiated_air_wt/run_0196_va3_M4.19_Lc8_noneq_inletTt_full5/nozzle.h5` を同関数で反転すると、node 23724 で次の差が出ます。

| 確認量 | 値 |
|---|---:|
| ソルバの `res_0.h5` の温度 | 257.40765 K |
| 変換器の入力経路で反転した温度 | 203.27890 K |
| 差 | **54.12875 K** |

DB 変更時には、この誤った温度でエネルギー補正を計算します。さらに [同ファイル:266](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:266) は湿潤セルの温度保存を検査せず、乾きセルでも許容差超過を警告だけで通します。§5.1 の「T 保存検査済み」は過大な記述です。

**対案:** ソルバと同じ二相 EOS で入力温度を反転・検査し、保存条件違反は書込み前に失敗させてください。総水量も、`Y_H2O` が既に総水分率であることに合わせ、現在の `ρ(Y_H2O+g)` という二重計上を修正してください。

### 3. Major — restart の名前・DB 照合が仕様を満たしていない

**根拠:** [interp_field.py:50](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:50) は名前列だけを返し、MW・係数・datum を比較しません。読込例外は握り潰して照合を省略します。指摘1の SERN config は、この省略経路に入ります。

同一メッシュの [restart_by_index:438](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:438) も、配列長と `roY` の番号によるコピーだけです。種順序・DB の照合を行うという [methods/thermophysics.md:259](/home/sano/work/forge/methods/thermophysics.md:259) と一致しません。

**対案:** 保存済み species metadata と DB・datum の共通照合を全 restart 経路に適用し、照合不能も通常はエラーにしてください。配列本数が同じという理由で異なる種をコピーしてはいけません。

### 4. Major — 種変換器が SERN の stream lump とトレーサ新設を扱えない

**根拠:** [convert_species_field.py:82](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:82) は、実種の行き先を一意の lump に限定します。実際の `case/46.sern_design/run_0100_species_regress_euler_lumped/species_meta.yaml` を使うと、**同じ配置への変換さえ**次の理由で拒否されます。

```text
destination: 実種 N2 が EXH と AIR の両方に入っている
```

これは正常な SERN の流れ別 lump です。また [同ファイル:304](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:304) は元にある `roXi` をコピーするだけで、旧 `[EXH,AIR]` から `full` へ移す際のトレーサを生成しません。

**対案:** 保存的な種展開と、作動点変更に伴う組成再初期化を別の明示的操作として実装してください。後者では保存済み流れ情報と流入元ラベルを使い、`lumped → full` の `roXi` 初期化まで検査する必要があります。

### 5. Major — `roXi` が周期 node の残差合算・状態同期から漏れている

**根拠:** [periodicNode_d.cu:87](/home/sano/work/forge/solver_density_cuda/cuda_forge/periodicNode_d.cu:87) の追加残差合算対象は SST・species・凝縮だけです。[同ファイル:114](/home/sano/work/forge/solver_density_cuda/cuda_forge/periodicNode_d.cu:114) の保存量同期にも `roXi` がありません。

新しい [registerTracer:96](/home/sano/work/forge/solver_density_cuda/variables.cpp:96) は、これらの既存リストに登録していません。したがって周期 seam では、合併 CV の流れ場と部分 CV のトレーサ残差を組み合わせることになります。

**対案:** `res_roXi` の周期合算と `roXi` の同期を組み込み、非一様トレーサの周期移流で seam の連続性・積分保存を検証してください。

### 6. Major — `roXi` の dual-time「対応」は物理時間積分になっていない

**根拠:** [main.cpp:1616](/home/sano/work/forge/solver_density_cuda/main.cpp:1616) は各擬似時間反復で `roXiN` を更新し、物理時間項なしでトレーサを前進させます。[scalarTransport_d.cu:417](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cu:417) も擬似時間の point-implicit 更新だけです。

このためトレーサの進行量が物理時間刻みではなく、擬似時間刻み・反復回数に依存します。[solver-settings.md:181](/home/sano/work/forge/procedures/solver-settings.md:181) の「dual-time 対応」は誤りです。

**対案:** 当面は `tracer` と dual-time の併用を入力時に拒否してください。対応を掲げるには、物理時間履歴・BDF 残差・対角項と時間精度検証が必要です。

### 7. Major — `condensationSpecies` と `condModel` の不一致を拒否しない

**根拠:** [solverConfig.cpp:711](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:711) は名前・index・種数を検査しますが、物質と `condModel` の対応を検査しません。例えば `[N2,H2O]` に `condensationSpecies: H2O, condModel: 0` を指定しても、この検査を通ります。

物性選択は [condensationProperties_d.cuh:354](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:354) で独立に行われるため、水の保存量に窒素の凝縮物性を適用できます。§4.2 の明示的な拒否条件が未実装です。

**対案:** 解決済み凝縮種と物性モデルの対応を host 入力段階で検査し、不一致・未対応モデルを拒否してください。単体試験は `DB.index()` ではなく config 読込まで通してください。

### 8. Major — `lumped + keep` の排気率契約が成立していない

**根拠:** [composition.py:479](/home/sano/work/forge/design/forge_design/gas/composition.py:479) は `full` のときだけトレーサを有効にします。実際の m6_on、`lumped + keep:[H2O]` の純排気入口は、

```text
Y = [0.7588908814, 0, 0.2411091186]
tracer = False
```

です。現在仕様の `ξ=Y_EXH` では純排気が `ξ=0.75889` となります。また、[methods/design/overview.md:783](/home/sano/work/forge/methods/design/overview.md:783) が実装済みとして記載する `exhaust_fraction(run)` はコードに存在しません。

**対案:** `keep` や任意 lump によって単一種分率が流入元ラベルにならない配置でもトレーサを有効にし、保存メタデータを読む共通アクセサを実装してください。

### 9. Major — 凝縮 OFF の H₂O 後処理が固定組成へ戻り、分圧を約24%誤る

**根拠:** [forge_species.py:78](/home/sano/work/forge/solver_density_cuda/tools/forge_species.py:78) は凝縮 ON の場合だけ `vapor_array` を返します。[axis_csv_va.py:27](/home/sano/work/forge/case/44.vitiated_air_wt/axis_csv_va.py:27) はそれがないと `0.0462711*P` に戻ります。

実際の `case/44.vitiated_air_wt/run_0203_va3_M4.19_Lc8_dry_full_dry_runner/` は H₂O が `Y0` にあるのに `vapor_array=None` です。最終場で比較すると、

```text
固定組成による分圧 / Y0 に基づく分圧 ≈ 0.758764
```

となります。凝縮 OFF の飽和温度・過飽和度診断が誤ります。

**対案:** 凝縮の有効・無効と独立に対象種 `H2O` を名前で解決し、DB の MW を使って分圧を計算してください。

### 10. Major — 検証ゲートの未達が残り、§5.1 に正しく残されていない

**根拠:**

- **小型5種 node TP の PASS がない。** [plan:195](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:195) は独立の必須ゲートです。指定された case/44 の再判定は `NOT CONVERGED (stalled/plateau)` でした。SERN の読み替えは、この小型ケースを免除しません。
- **SERN のトレーサ許容差を超える。** 実際の HDF5 で再計算した結果は次のとおりです。

  | 比較した run | `max |Xi−Y_EXH|` | `1e-6` を超える node 数 |
  |---|---:|---:|
  | `case/46.sern_design/run_0100_species_regress_euler_lumped/` 対 `run_0101_species_regress_euler_full/` | `1.73032e-4` | 4,638 |
  | `case/46.sern_design/run_0102_species_regress_euler_lumped_24k/` 対 `run_0103_species_regress_euler_full_24k/` | `1.73330e-4` | 4,813 |

  平均誤差が小さくても、[plan:208](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:208) の場の許容差を満たしたことにはなりません。
- **CEA 物性差の `rtol 1e-6` も未達。** 保存済み `species_db_cea_va3.yaml` と内蔵表を200–3000 Kで比較すると、H₂O の cp・h 相対差は `1.11017e-6` です。係数検査の結果と物性値ゲートを分けて記録する必要があります。
- **「同じ固定点」は裏付け不足。** [case/44 README:643](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:643) はそう断定しますが、`run_0195/0196` の VERDICT はともに `NOT CONVERGED`。時系列の `ALL STEADY` は報告量の定常性を支持しても、固定点への収束を証明しません。
- §5.1 は restart 照合・種変換を完了扱いにし、小型収束ケースなどを独立した残作業として列挙していません。

**対案:** 未達ゲートと指摘1–9を §5.1 に戻し、小型 node PASS、トレーサ誤差の原因切分け、作動点変更・湿潤 restart を検証してください。`F-sp1` は既知の制限として残せますが、未実施の粘性試験との関係も明示してください。

### 11. Minor — 現在仕様・索引に実装と異なる記述が残る

**根拠:** [methods/thermophysics.md:250](/home/sano/work/forge/methods/thermophysics.md:250) の `loadSpeciesDB()` は実際の関数名と異なり、「内蔵 DB は空パス時のみ」も、外部 DB を内蔵 DB に上書きする実装と異なります。[plans/README.md:27](/home/sano/work/forge/plans/README.md:27) はまだ `draft`、[methods/condensation.md:989](/home/sano/work/forge/methods/condensation.md:989) は追加予定の記述です。

**対案:** 不具合修正後の実装・対応範囲に合わせて同期してください。既定 `pseudo`／SERN `lumped` の維持と、Y入力検査の厳格化も区別して記載してください。

## 推奨

**`active/in_progress` を維持し、実装の欠落と未達ゲートを解消してから result レビューを再実施してください。**

確認した範囲では、SERN の力は `--drift 0.001 --osc 0.001`、case/44 の凝縮量は `0.002` で **`OVERALL: ALL STEADY`** を再現しました。対象最終場の NaN/Inf もありません。これらは有効な回帰証拠ですが、上記の不具合を打ち消しません。

指定 diff を取得して確認しました。Python 単体はファイル生成ブロックを除いて再実行し PASS。書込みを伴う単体部分と C++ 単体の再実行は行っていません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 10 / Minor 1
