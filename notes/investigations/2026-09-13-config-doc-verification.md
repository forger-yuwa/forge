# config_doc.py 検証記録 (2026-09-13)

plan: [tooling-config-self-documentation.md](../../plans/active/tooling-config-self-documentation.md) §6

## 1. 自己試験 `config_doc.py selftest`

```
PASS 節違い (同名キーがトップレベルにあっても見逃さない)
PASS 未知キー (綴り間違い)
PASS 配列キーは誤検知しない
PASS 使用禁止キー (値が空の辞書でも見逃さない)
PASS ドット入りのキー名を入れ子と同一視しない
PASS 未知の節の下でも、正しい節が分かる子は出す
PASS スカラーのキーに辞書
PASS ブロックスカラーの中身をキーと誤認しない (原文復元 + YAML 値一致)
PASS 複数行のフロー形式でも入れ子を取り違えない (原文復元 + YAML 値一致)
PASS キーの後ろのコメントを入れ子開始と誤認しない (原文復元 + YAML 値一致)
PASS YAML 真偽値に見える化学種名 (NO) を壊さない (原文復元 + YAML 値一致)
PASS 空の節を潰さない (原文復元 + YAML 値一致)
PASS 引用符・バックスラッシュを含む文字列 (原文復元 + YAML 値一致)
PASS フロー形式のキーにフルパス付きの注釈が付く
PASS coverage: 解析できなかった読み出し 0 件
PASS 説明なし 0 キー / 抽出 159 キー

VERDICT: PASS
```

期待値と入力はツール内の `SELFTEST` / `SELFTEST_TEXT` にある (過去に実際に壊れた入力をそのまま残してある)。

## 2. 起動ログ (C++ 最小ハーネス)

`solverConfig::read()` だけを呼ぶハーネスを `g++ -std=c++17 -I. harness.cpp input/solverConfig.cpp -lyaml-cpp` で作り、
`case/16.nozzle_wys/run_0456_perf_regress_node2d_cond/solverConfig.yaml` に対して実行した。

| 項目 | 実測 |
| --- | --- |
| 明示キーの行 (`'key' in 'section': value`) | 53 行 (従来書式のまま) |
| 省略キーの行 (`[default] ...`) | 76 行 |
| `FORGE_CONFIG_LOG_DEFAULTS=0` での `[default]` 行 | 0 行 |

## 3. 既存 run の config (無作為 40 run のうち config を持つ 32 本)

対象一覧:

```
case/17.flared_cone/run_keep_slau_v4
case/03.cavity_density/run_keep
case/05.sod_shock_tube/run_0008_diffusion_n2n2
case/16.nozzle_wys/run_0194_user_node_euler_cond
case/44.vitiated_air_wt/run_0039_va2_M4.4_eq
case/44.vitiated_air_wt/run_0055_va2lp_M4.4_dry
case/44.vitiated_air_wt/run_0005_va_R2_LU6_Lc8
case/18.backstep/run_slau
case/15.2D_poiseuille/run_ROE
case/05.sod_shock_tube/run_keep_slau
case/11.Taylor-Green_prism/run_ausm
case/06.mach3_wind_tunnel/run_slau
case/16.nozzle_wys/run_0047_fig3_2d_lam_dry
case/04.laval_nozzle/run_slau
case/44.vitiated_air_wt/run_0089_va2c_M5.0_noneq
case/26.flat_plate_sst/run_0029_node_outletic_subsonic
case/28.cutler_coaxial_jet/run_0001_meshcheck
case/23.axi_nozzle/run_tp_air_axisym
case/43.node_axis_dof/run_0044_eul_van_edgemid
case/17.flared_cone/run_keep_slau_v4.turb0.5_probe
case/44.vitiated_air_wt/run_0033_va2_M4.19_eq
case/05.sod_shock_tube/run_0012_aws_repro_cell
case/23.axi_nozzle/run_cpg_air_axisym
case/44.vitiated_air_wt/run_0036_va2_M4.3_eq
case/44.vitiated_air_wt/run_0072_va2fx_M4.4_eq
case/43.node_axis_dof/run_0048_ns_van_em_lmp2
case/16.nozzle_wys/run_0456_perf_regress_node2d_cond
case/16.nozzle_wys/run_0461_frozen_roe_protect
case/13.nozzle_H/run_KEEP_AUSM
case/10.half_sphere/run_roe
case/44.vitiated_air_wt/run_0103_va2eq2_M4.5_eq2
case/44.vitiated_air_wt/run_0090_va2c_M5.0_eq
```

指摘の内訳 (すべて実在の問題であることを個別に確認済み):

| 種別 | キー | 件数 | 実体 |
| --- | --- | --- | --- |
| 未知 | `time.last.time` | 17 | `last["time"]` を読むコードは無い (コメントアウト済み) |
| 未知 | `time.implicit` (節ごと) | 8 | `config["time"]["implicit"]` を読むコードは無い。`cfl_pseudo` は `time.deltaT` が正 |
| 未知 | `turbulence.kInf` / `omegaInf` | 5 / 5 | 正しくは `kInit` / `omegaInit` (recommended-settings.md の推奨が誤っていた。訂正済み) |
| 未知 | `mesh.nodeWallViscGradFlux` | 1 | ソルバのどこにも無い |
| 廃止 | `mesh.nodeAxisDirichlet` ほか 3 | 各 2 | 2026-08-16 廃止。書くと起動エラー |
| 必須欠落 | `time.unsteady` ほか | 4〜9 | 2019 年前後の古い run。現ソルバでは起動しない |

## 4. 生成物

- `annotate`: 元のテキストには触れず、YAML 構文木のキー位置の前にコメント行を挿すだけ。出力前に
  (1) 挿入行を外すと元テキストに 1 文字違わず戻る (2) YAML として読んだ値が元と一致する、の 2 つを検査する
  (どちらか外れれば出力しない)。(1) だけではブロックスカラー内への `#` 挿入を見逃す。
- `annotate -o <入力自身>`: 拒否して exit 2。
- `template`: YAML として読める (26 トップレベル項目)。廃止・使用禁止キーを含まない。
- `coverage`: 解析できなかった読み出し 0 件 / 抽出 159 キー。
