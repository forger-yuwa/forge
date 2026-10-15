# Forge ソルバー設定リファレンス

> **位置づけ**: 本文書は `solverConfig.yaml` の**運用上の設定リファレンス** (各設定の意味・正しい使い方・注意点) である。
> スキームの理論的背景や実装の詳細は `methods/<area>/` を参照すること (例: 対流スキームは [`../methods/convection/`](../methods/convection/)、
> 時間積分は [`../methods/time_integration/`](../methods/time_integration/))。ここは「どの値を選ぶか」、`methods/` は「なぜそう動くか」を担う。

この文書は `solverConfig.yaml` の主要な数値設定について、正しい使い方と注意点をまとめたものです。解析設定を変更するときは、**必ずこのファイルを参照してから**変更を行うこと。

## convMethod — 対流スキームの次数

| 値 | スキーム | 説明 |
|---|---|---|
| `0` | 1次風上 | 安定性が高い。収束困難なケースの初期化や発散回避に使う |
| `1` | 2次風上 | 通常の高精度計算に使う |
| `2` | 3次 MUSCL | より高精度だが計算コスト増。安定した解が得られてから使う |

## limiter — リミタ

| 値 | リミタ | 説明 |
|---|---|---|
| `0` | なし | **リミタ無しのため最も発散しやすい。** 基本的に単体では使わない |
| `1` | Barth | |
| `2` | Venkatakrishnan | 通常使用するリミタ。`convMethod: 1` または `2` と組み合わせる |

### リミッタの基準値の記録 `limiter_ref.json` (2026-10-16)

`limiterScaled: 1` (既定) の Venkatakrishnan は基準値 ρ・P・a (`space.limiterRoRef`・`limiterPRef`・`limiterARef`) と基準長
L_ref (`limiterRefLength`) で無次元化する。**省略 (0) は起動のたびに開始の場の体積平均** (L_ref は格子の境界箱の対角) で決まるので、
同じ YAML でも起動ごとに作用素が変わりうる (再開・段の切り替えのたび)。forge は起動のたびに、決まった値を作業ディレクトリ (run dir) の
`limiter_ref.json` に**上書き**で書く (一時ファイル → rename; 出力だけで値の決め方は変えない):

```json
{"schema": "forge.limiter_ref/v1", "time": 1792078376, "cfg_fnv": "cdf047de7f18fd9", "limiterScaled": 1, "limiterDiag": 0,
 "computed": true, "limiterRoRef": 0.59298852861870699, "limiterPRef": 192416.41762690683, "limiterARef": 327.13657637494532,
 "limiterRefLength": 9.0296406421583377,
 "auto": {"limiterRoRef": true, "limiterPRef": true, "limiterARef": true, "limiterRefLength": true},
 "venkatK": 0.050000000000000003, "h_i": "sqrt(A_planar)"}
```

- 値は `%.17g` (double を往復できる桁)。log の `[limiter]` の行 (10 桁) と同じ値で、log にも `[limiter] 記録: limiter_ref.json (... %.17g)` を 1 行出す。
- `cfg_fnv` は `solverConfig.yaml` の FNV-1a 64 (`forge_launches.jsonl` と同じ)。`tools/stage_manifest.py` はこれで段と結び付け、値を区間の
  hard キー `limiter.ref` に入れる (値が違う段は別の区間、記録が無く自動の値が分からない段は連結しない)。段の記録 (`StageManifest.add`) は
  **その段の forge の後に**呼ぶ。
- 基準値を使わない起動 (`limiterScaled: 0` かつ `limiterDiag: 0`) も `"computed": false` で書く。書けなかったら log に警告を出して計算は続ける。
- 固定するには値をそのまま `space` に書く (`limiterRoRef: 0.592988528618707` のように Python の `repr` で書けば forge は同じ double を読む;
  case/71 の Euler で自動 → 固定の 2 起動の記録が `%.17g` の文字列まで一致)。設計チェーンは `evaluate.main_line.limiter_ref: {from_run: <run>}` で
  別 run の記録を継承する (`runner_axismach.limiter_ref_from_run`; plan tooling-settle-convergence-and-monitor §4.6)。
- **系列の自動の継承** (2026-10-16、同 plan §5.1 #6b): 設計チェーンの NS (`prepare_ns`) で B0 の本段 (`evaluate.main_line` が mapping) の
  `limiter_ref` を省くと、初期場の元 (`ic_from`; `deltastar_loop` の前の pass・`--ic-from`・`--prev` もここを通る) が NS の run なら、その run の
  本段の起動の `limiter_ref.json` (cfg_fnv が本段の設定と一致・computed) の 4 つを固定値で継承する。記録が無い・別の起動のものなら本段の設定の明示
  ρ・P・a を継承し、`limiterRefLength` は明示があれば写し、無ければ新しい run でも forge が格子から決める (既存の固定値の run は L が自動)。
  初期場の元が Euler・無いときは自動。元の run が NS なのに値が分からない (記録も ρ・P・a の明示も無い)・記録が壊れている・`computed: false`・
  本段の前の段 (S1・S2) で終わった run (`MAIN_STAGE_STARTED` が無く `stage_manifest.json` の最後の段が `main` でない) のときは準備で止まる —
  `limiter_ref: auto` (継承しない自動) か値の明示で通す。段の記録の無い旧い run は `solverConfig.yaml` を本段とみなし、出所に「本段の確認なし」と
  書く。明示・`from_run` は継承より優先、point の本段には効かない。
  出所は `prepare_info.json` の `main_line.limiter_ref_source` と `MAIN_STAGE_STARTED` の `limiter_ref.source` (`runner_axismach.limiter_ref_inherit_from_ic`)。

## 推奨組み合わせ

| 用途 | convMethod | limiter |
|---|---|---|
| 初期化・発散回避 | `0` | `2` |
| 通常の高精度計算 | `1` | `2` |
| 高精度（安定解から） | `2` | `2` |

**`limiter: 0` は原則使用しない。** 意図的に使う場合は理由を明記すること。

## CFL の定義 (`cfl` と `cfl_pseudo`) — 重要

`time.deltaT` には `cfl` と `cfl_pseudo` の 2 つがある。**定常計算 (`unsteady: 0`) と
二重時間刻み (`dualTime: 1`) では、実際の局所時間刻みを決めるのは `cfl` ではなく
`cfl_pseudo` である。** これを取り違えないこと。

根拠 ([`setDT_d.cu`](../solver_density_cuda/cuda_forge/setDT_d.cu) の `setDT_d_wrapper`):

- `unsteady == 0` または `dualTime == 1` のとき、局所刻みは `setDTlocal_pseudo_cell_d` で
  `dt_local = cfl_pseudo * dt / cfl_cell` と設定される。`cfl_cell` は同じ `dt` から計算されるため
  `dt` が相殺し、**実効局所 CFL = `cfl_pseudo`** になる。
- `cfg.cfl` は `cfg.dt` をスケールするだけで、定常計算の積分自体には効かない。したがって定常で `cfl`
  だけを上げても収束は速くならない。このため console モニタ行は定常では `dt`/`maxCFL` を表示しない
  (2026-09-09 以降。[`methods/architecture/overview.md`](../methods/architecture/overview.md) §8.5)。

運用ルール:

- **定常 (`unsteady: 0`) の収束を速めたいときは `cfl_pseudo` を上げる** (例: 陰解法 `blockDPLUR: 1`
  なら `cfl_pseudo` を 20〜50 程度まで)。`cfl` を上げても無意味。
- 非定常 (`unsteady: 1`, `dualTime: 0`) の物理時間刻みは `cfl` (または `dt`) で決まる。
- console モニタ行の `maxCFL` は非定常 (`unsteady: 1`) のみ表示され、そのときは `cfg.dt` 基準の物理 CFL である
  (dual-time では経験的に ≲12 が安定域)。定常には表示されない (以前の `max cfl : ...` 行は `cfg.dt`
  に追従するだけの無意味な値だったため廃止)。

## lineImplicit — 壁法線 line-implicit (block-Thomas) と v2 オプション

```yaml
time:
  deltaT: {..., blockDPLUR: 1, lineImplicit: 1, lineKFreeze: 1, lineDtDirectional: 1}
```

壁法線ライン上の隣接結合 (対流+音響 Jacobian) を point-DPLUR の lag から block 三重対角の
直接解に昇格する。`timeIntegration: 11` + `blockDPLUR: 1` 専用、`lowMachPrecond>=2` 併用不可。
全キー既定 0 = 挙動不変。詳細は [`methods/time_integration/implementation.md`](../methods/time_integration/implementation.md) の
「line-implicit」節と plan ([v1](../plans/accepted/time_integration-line-implicit.md) /
[v2](../plans/accepted/time_integration-line-implicit-viscous-v2.md))。

- **`lineImplicit: 1`**: 本体。v2 (factor/solve 分離) 込みでコスト +32%/subiter (ny160 DDES 実測)。
- **`lineKFreeze: 1`**: dual-time サブ反復間で K/LU を凍結 (収束軌道不変・コスト削減)。dual-time では常用推奨。
- (削除済 2026-09-18) `lineDtWallRelief`: 壁境界半割面の λ も CFL max から除外する診断スイッチ。発散したので削除した (plan config-key-pruning §5.2)。
- **`lineDtDirectional: 1`**: 方向別 dt — 内部 line 面の λ (音響込み) を擬似 dt の CFL max から除外。
- **`lineDtDirectionalCap: R`** (2026-10-09、既定 0 = 上限なし): 方向別 dt の伸びを Δτ ≤ R·Δτ_point に抑える (平均流と SST の両方に効く)。検証中 — plan [time_integration-implicit-thermal-jacobian](../plans/active/time_integration-implicit-thermal-jacobian.md) §4.4。
- **Thomas の配列の並び** (環境変数、2026-10-10): 既定は LAYOUT2 (ラインが隣り合う並び + 前の節点の値をレジスタに持つ。数値は従来とビット一致、case/45 で −7.6 %/step)。`FORGE_LINE_LAYOUT=0` で従来の並び。起動ログの `[line] Thomas の配列の並び: …` で実効の並びを確かめる。詳細は [methods/time_integration/implementation.md](../methods/time_integration/implementation.md) の「Thomas の配列の並び」。
  壁ノードの境界半割面は除外されない (壁 CV 自身の Δτ は境界面律速のまま) 点に注意。
- `lineViscCoupling` / `lineViscousDtRelief`: line 面の粘性結合と粘性 CFL 割引。`lineViscCoupling` は 0 (結合なし、既定) / 1 (スカラー α I を 5 行すべて) /
  2 (薄層の粘性・熱伝導の Jacobian をライン面の K と対角に。連続の行 0、運動量は速度、エネルギーは T と粘性の仕事。等温壁の行は ΔT_w = 0 の拘束)。
  2 は node・`lineImplicit 1`・`blockDPLUR 1`・`timeIntegration 11`・`lowMachPrecond < 2`・`heatCorrSU2 0`・`nodeIsothermalEnergyBC 0`・`wallTreatmentSST 0` 専用 (それ以外は起動時に拒否)。
  検証中 (plan [time_integration-line-viscous-jacobian](../plans/active/time_integration-line-viscous-jacobian.md))。
  圧縮性の壁法線 pseudo-dt 律速は音響 (λ_visc/λ_ac=2ν/(Δn·c)≪1) なので通常は効果僅差 — 既定 off で可。

**使い分け (2026-09-03 実測)**: 定常 M6 ノズル (case/45) では cfl 上限を上げない (streamwise
対流律速) — 使わない。**壁解像 DDES/dual-time (case/39) では point の cfl_pseudo 上限 (2 で発散)
を 8 まで引き上げ、`cfl_pseudo 4 + nSub 13-15` で現行比 ~13% 高速・ωバースト低減** — 採用候補
(長時間 run でのバースト余裕検証は未了)。

## implicitThermalJacobian — エネルギー行の熱伝導 Jacobian (opt-in、2026-10-09)

```yaml
time:
  deltaT: {..., blockDPLUR: 1, implicitThermalJacobian: 1}
```

block DPLUR の粘性の対角は 5 行に同じスカラー (2ν·δ/dcc) を保存量に掛けるので、エネルギー行では ρE にかかり、
ρE がほとんど動かずに ρ と T が入れ替わる等圧のエントロピーのモードを抑えない。ビットマスク (既定 0 = 従来どおりビット同一):

- **ビット 1**: 内部の node 間面でエネルギー行を熱伝導の Jacobian k_face·δ/dcc·(γ/c_p)·∂e/∂Q に置き換える (k_face は残差と同じ式)。
- **ビット 2**: 等温壁の節点のエネルギー行を拘束の行 [−e_w, 0, 0, 0, 1] にする (原因の切り分け用)。
- **ビット 4** (ビット 1 と併用、キー 5): エネルギー行にも従来のスカラーを残し、温度の項はその上に足す (ビット 1 だけは高速・高 μ_t の所で不安定化した)。

node・`timeIntegration 11`・`blockDPLUR 1`・`lowMachPrecond < 2` 専用。`lineViscCoupling 1`・`mesh.nodeIsothermalEnergyBC 1` との併用は起動時に拒否する (`lineViscCoupling 2` とは併用できる)。
LHS だけの変更なので定常解は変わらない。**検証中** (plan [time_integration-implicit-thermal-jacobian](../plans/active/time_integration-implicit-thermal-jacobian.md)) で、既定化はしていない。
式は [`methods/time_integration/implementation.md`](../methods/time_integration/implementation.md) の「エネルギー行の熱伝導 Jacobian」。

## bodyForce — 一様体積力 (周期チャネル駆動)

```yaml
bodyForce: [7.24, 0.0, 0.0]   # [N/m³] top-level キー。既定 [0,0,0] = off (ビット不変)
```

全セルに運動量ソース $f_i V$ とエネルギーソース $(\mathbf{f}\cdot\mathbf{u})V$ を加算する。
主用途は**周期チャネルの駆動** (入口出口が無いため平均圧力勾配の代わりに使う):
平衡で $f_x \delta = \tau_w$、つまり $u_\tau = \sqrt{f_x \delta/\rho}$ を直接指定できる
(δ: チャネル半幅)。発熱の逃げ場が要るため**壁は等温**にすること (断熱では bulk 温度が漸増)。

- 検証: 体積力駆動 Poiseuille 厳密解 (case/24 `run_poiseuille_bf_*`, u_max=fH²/8μ・T 四次分布・q_w=f²H³/24μ)
- 軸対称は未対応 (config 読込で拒否)。陰解法では明示ソース (Jacobian 無し)

## detectNaN — NaN/Inf 検知診断モード

`time.deltaT.detectNaN` (任意, 既定 `0`)。発散の発生ステップと場を特定するための診断オプション。

- `0` (既定): 検査を一切行わない。通常実行はビット不変・性能影響なし。
- `1`: **毎ステップ終端**で保存量 `ro`,`roUx`,`roUy`,`roUz`,`roe` と圧力 `P` (RANS 時は `roK`,`roOmega` も) を内部セル (`nCells`) にわたって検査する。NaN/Inf が 1 つでもあれば、その時点の解を `res_nan_<step>.h5` / `.xmf` に強制ダンプし、最初に該当した変数名を `stderr` に出力して `exit(1)` で停止する。

```yaml
time:
  deltaT:
    detectNaN: 1   # 発散をその場で捕捉してダンプ・停止する
```

毎ステップ GPU リダクション (`thrust::any_of`) が走るため、原因調査時のみ有効化し、本番計算では `0` のままにする。ダンプされた `res_nan_<step>.h5` を ParaView で開けば、どの境界・どのセルで最初に NaN が出たかを直接確認できる。

## lowMachPrecond ≥ 2 — 低マッハ完全前処理の注意 (node・多成分、2026-10-12)

```yaml
time:
  deltaT: {..., blockDPLUR: 1, lowMachPrecond: 2, precondEps: 0.15, speciesPrecondDt: 1}   # キーの置き場所の例 (推奨値は recommended-settings.md)
```

`lowMachPrecond` 2/3 は `timeIntegration 11`・`blockDPLUR 1` 専用 (設計と検証は
[plan](../plans/accepted/time_integration-lowmach-preconditioning.md)、実装は
[`methods/time_integration/implementation.md`](../methods/time_integration/implementation.md) の
`implicit_defect_correction_block_precond_d`)。

- **node × `lowMachPrecond ≥ 2` は、このブランチでは 2026-10-12 の修正より前のバイナリで境界の節点がほぼ凍結していた**。
  前処理の block カーネルに標準の block カーネルの node の行の処理 (境界半割面の粘性対角を課さない・軸の roUy 行・
  壁の運動量 3 行・等温壁の roe 行の拘束) が無く、1 step の試験で入口・出口の節点の更新が P0 の 0.02〜1 % に縮み、
  軸・壁の拘束の成分も 0 にならなかった ([plan architecture-float-state-double-geometry](../plans/active/architecture-float-state-double-geometry.md)
  §6.57・§5.1 #26。chem のブランチは `876ec8b4` で修正済み)。それより前の node の前処理の run の結論は、再計算なしに流用しない。
  等温壁の roe 行の拘束は標準の経路と同じく強いピン (`mesh.nodeIsothermalEnergyBC 0`) のときだけで、弱形式 (`1`) では
  エネルギー行を残す (弱形式の近似対角 `weakIsoDiag` は前処理のカーネルには入っていない)。
- **`time.deltaT.speciesPrecondDt`** (int、既定 `1`): 多成分 (nSpecies ≥ 2) × `lowMachPrecond ≥ 2` の定常擬似時間で、
  化学種の擬似時間刻みを選ぶ。前処理は流れ 5 変数の擬似時間項だけを Γ_c で変え、`dt_local` を前処理スペクトル半径に
  合わせて拡大する。分離更新の ρY_s は前処理されないので、拡大後の Δτ' で進めると組成の前線と密度の前線が擬似時間で
  別の速さで進み、異組成の接触面で P が暴走する (chem のブランチの Cabra H₂/N₂ 浮き上がり火炎、`876ec8b4`。
  `lowMachPrecond 0` と dual-time では起きない)。

  | 値 | 化学種の擬似時間刻み |
  | --- | --- |
  | `0` | 流れと同じ拡大後の Δτ' (旧挙動。A/B 用) |
  | `1` (既定) | 前処理で拡大する前の Δτ (物理スペクトル半径基準) |
  | `2` | Δτ'·β (β = `lowMachBeta`、[1e-6, 1] に収める) |

  効くのは化学種の scalar-DPLUR (`speciesImplicitCoupling` 1/2) だけ。**`speciesImplicitCoupling 0` の point-implicit と
  受動種 (トレーサ・凝縮モーメント) は対象外**で、前処理のときも拡大後の `dt_local` で進む。前処理と多成分を併用するときは
  `speciesImplicitCoupling 1` 以上にする。`lowMachPrecond < 2` では何もしない (ビット不変)。実際に使われた刻みは
  `output.extraFields: [dt_local_sp]` で出せる (多成分のときだけ確保。効いていない run では 0)。
  このブランチでの β < 1・多成分の検証は未了 (§5.1 #26)。

## リスタート手順

前の run の結果から引き継ぐ場合は `valueFileName` に該当の `res_*.h5` を指定する。

```yaml
mesh:
  meshFileName: "axi_nozzle.h5"      # メッシュ (変わらない)
  valueFileName: "res_10000.h5"      # リスタート元の結果ファイル
```

`initial` キーの値はリスタート時には無視される（valueFileName が res_*.h5 を指していれば、そこから値を読む）。

## 段階的な計算戦略

均一初期値から超音速流れを計算する場合、いきなり2次以上で始めると全域超音速流出に発散することがある。次の2段階を踏む：

1. `convMethod: 0, limiter: 2` で安定した解を得る（例：10000ステップ）
2. その結果を `valueFileName` で引き継ぎ、`convMethod: 1, limiter: 2` で計算する

## speciesFaceReconstruction — 多成分 TP の face 組成整合

`time.deltaT.speciesFaceReconstruction` (任意, 既定 `0`)。多成分 thermally-perfect (`thermalMethod: 2`,
`nSpecies > 1`) 専用。face thermo に使う組成の再構成次数を制御する。

| 値 | 動作 | 使いどころ |
|---|---|---|
| `0` | face thermo はセル値 Y (1次, mixed-order) | 既定・ビット不変 |
| `1` | **S2**: Y を ρ と同一リミッタで face へ 2 次再構成し thermo に使用 | **多成分 TP で contact/混合層の圧力振動・残差床を下げたいとき (production 推奨)** |
| `2` | **S3**: S2 に加えて**化学種と受動種 (トレーサ `roXi`・凝縮モーメント) の移流そのもの**を同じ face 値で 2 次に | 2 次風上のスカラ移流が要るとき。**`convMethod ≥ 1` と `speciesImplicitCoupling: 1` が必須** (下記) |

```yaml
time:
  deltaT:
    speciesFaceReconstruction: 1   # 多成分 TP の face thermo 組成を ρ と同次数で整合
```

効果の実測 (case/28 He/空気 coaxial): 圧力振動振幅 cfl2 で −77% / cfl4 で −48%。単成分・CPG では
無効果 (組成を thermo に使わないため)。詳細は [`../methods/convection/theory.md`](../methods/convection/theory.md)

**`2` (S3) の本番化 — 現行 (2026-09-17, plan [species-passive-scalar-unification](../plans/accepted/species-passive-scalar-unification.md))**:
旧記述の「cfl≤2 限定の experimental、高 CFL で発散」は**撤回**した。発散は S3 固有ではなく `speciesImplicitCoupling 0` 固有で、
coupling 1 なら cfl 6 まで安定 (case/28 `run_0064`–`0078`)。使うときの条件は次の 3 つ。

- **`space.convMethod ≥ 1`**: 面値は流れと同じ `interp_dispatch` を通るので、`convMethod 0` では**黙って 1 次に落ちて S3 が不活性**になる。
  `check_solver_config.py` がこの組合せを FAIL にする。
- **`speciesImplicitCoupling: 1`** (受動種の `passiveImplicitCoupling` は自動で 1)。定常は `implicitRelax 0.7` 併用、dual-time は緩和なし。
- **有界性**: 化学種は面正規化 + リミッタ ψ_ρ、受動種はセル局所で無次元化した Venkatakrishnan ψ_P と独立クリップ。
  dual-time ではさらに物理 step 末尾の保存的 FCT (`passiveFct 1`, 既定) が保存と有界性を担う。合否は `check_passive_budget.py` / `check_passive_field.py`。

検証: 時間次数は BDF2/BDF1 × coupling 0/1 の 4 組で次数・収支とも PASS (case/44 `run_0429`–`0456`)、FCT 無効の対照は全モーメント 1.99–2.01。
**既定は `0` のまま** (ユーザ決定 2026-09-17): S3 は凝縮 onset を 0.2 r_t / 0.7 mm 下流に動かし、Wysłouzil では実験からむしろ遠ざかるため、
既存の凝縮回帰と直接比較できない。使うときは run の README にその旨を書く。
の「多成分 TP の face 組成整合」節。

## passiveScalarScheme ほか — 受動スカラ (排気トレーサ・凝縮モーメント) の輸送経路 (2026-09-17)

plan [species-passive-scalar-unification](../plans/accepted/species-passive-scalar-unification.md)。`time.deltaT` 配下。

| キー | 既定 | 意味 |
|---|---|---|
| `passiveScalarScheme` | **1** | 1 = トレーサ `roXi`・凝縮モーメント `rog_s, roQ2_s, roQ1_s, roQ0_s` を化学種の輸送経路 (化学種と同じ勾配・面再構成・移流残差・境界/ピン・周期・拡散 [トレーサのみ Fick]・dual-time BDF) で受動種として解く。熱力学・ΣY 再正規化・ΣJ=0 補正・`speciesImplicitCoupling` 予測/commit には入らない。0 = 旧汎用スカラ経路 (1 次風上・拡散なし・トレーサは primitive 段クランプ; dual-time では物理時間項なし・`condLimiterMode 1` 降格・`tracer` 拒否) — A/B 用でビット不変 |
| `passiveImplicitCoupling` | −1 (自動) | 受動種の陰解法更新: 0 = segregated point-implicit (増分 × `passiveImplicitRelax`), 1 = 化学種と同じ scalar-DPLUR sweep で増分を作り、モーメントは更新クランプ (θ_u) に渡す。自動 = `passiveScalarScheme 1` かつ `speciesFaceReconstruction ≥ 2` で 1、他は 0 |
| `passiveImplicitRelax` | = `implicitRelax` | 受動種 point-implicit 更新の増分緩和 (scheme 1 のみ) |
| `speciesImplicitRelax` | 1.0 | 化学種の segregated (coupling 0) 更新の増分緩和。1.0 で現行と同じ写像 |
| `scalarCflMax` | 無効 | 化学種/受動種の更新だけ擬似 CFL をこの値で頭打ち (物理時間項は変えない)。保険用 |
| `passiveFct` | 1 | **dual-time の受動種 S3 に対する物理 step 末尾の保存的 FCT 補正** (plan §4.7 v5): `passiveScalarScheme 1` + SLAU + `speciesFaceReconstruction ≥ 2` + `timeIntegration 11` + `unsteady 1` + `dualTime 1` のときだけ作動。BDF2 を流束形の BE (前 step の増分を面流束 G/局所 H に分解) として扱い、低次 BE 陰解 (有効質量流束 ṁ^eff) を限界に Zalesak で反拡散流束を面共有の α で制限する (保存; α=1 の面では収束した完全陰的解のまま)。0 = 無効 (A/B; CV ごとの増分制限 θ_b だけになり非保存)。定常・RK では常に無効 |
| `passiveFctPrelimit` | 0 | Zalesak の前制限 (滑らかな極値でも作動するので既定 off; A/B 用) |
| `passiveFctSweeps` / `passiveFctTol` / `passiveFctTolAbs` | 100 / 1e-6 / 1e-30 | 低次 BE 陰解の Jacobi sweep 上限と線形残差の受入 (未達は `[passiveFct] WARNING` と log) |

**S3 (`speciesFaceReconstruction: 2`, 化学種・受動種の 2 次面移流) を使うときの組合せ** (node で検証済, case/28 `run_0064`–`0078`, case/44 `run_0216`–`0223`):
`speciesImplicitCoupling: 1` + `passiveImplicitCoupling` 自動 (=1) + `implicitRelax 0.7` で cfl 6 まで安定。**`speciesImplicitCoupling 0` + S3 は組成せん断層で発散する** (cfl 4 で step ~400)。
**`implicitRelax 0.7` は定常 (擬似時間) 限定** (投入前に `check_solver_config.py` で検査できる): dual-time の非定常計算では sub-iter の残差は下がるのに遅いモードが収束せず、同じ物理時刻の解が sub-iter 数に依存する
(case/44 `run_0399`–`0404`: nSub 40→80 の ρ 差が緩和なしの 1.3e-5 に対し 6.3e-4)。非定常で安定化が要るときは `cfl_pseudo` を下げるか nSub を増やす (2026-09-17, plan species-passive-scalar-unification §5.1 #26)。
S3 は凝縮の固定点を動かす (onset が case/44 で +0.18 r_t、Wysłouzil で +0.72 mm 下流; 前線が鋭くなる) ので、既存の凝縮回帰と直接比較しないこと。既定は `speciesFaceReconstruction 0` のまま。
更新確定時の上下限補正 (0≤ρξ≤ρ, モーメント ≥0) は `passiveFloorCorr_<name>` (level 2) と monitor ログに符号付き/絶対の体積積分 (step 内・累積) で出る。収束時は 0 であること。
**非定常 (dual-time) の合否**は `python3 solver_density_cuda/tools/check_passive_budget.py <run_dir>` (monitor の `[passive]` 積算 [floor / limCorr / FCT の基点逸脱・ピン交換 / 実現可能性クランプの成分別 |Δ|] を総量比 1e-6 で PASS/FAIL) で判定する。
checkpoint には受動種の流束形履歴 (`/CHECKPOINT/<cons>_fctG`, `_fctH`, `passive_fctMeff`) が入り、FCT 有効時の restart はこれが揃わないと全系 BDF1 から再開する。
凝縮モーメントの実現可能性 (許容領域 $x\le1,\ x^2\le y\le\sqrt x$; $x=Q_1/(Q_0r)$, $y=Q_2/(Q_0r^2)$) は更新後に最近点射影 (退化は単分散再初期化) で保証し、作動数と成分別収支を monitor に出す。
凝縮 run (受動種経路の有無を問わず) は monitorInterval ごとに `[cond-corr]` 行で理由別の補正量 (蒸気上限違反・負値 floor・増分制限・受動種 floor・射影・化学種再正規化・液滴消滅) の区間値と累積を総液量比で出し、液滴消滅以外が比 1e-6 を超えると `WARN` を出す。累積は restart で 0 から ([methods/condensation.md](../methods/condensation.md) 実装 §4c)。
TP carrier の凝縮 run では `viscMethod 2` の μ・λ と化学種拡散係数を気相組成 (液を除いた組成) で評価する (同 §7b)。
**二相拡散 `condensation.condTwoPhaseDiffusion`** (既定 0, 2026-10-02, opt-in・CFD 未検証; plan condensation-two-phase-transport #4e, methods/condensation.md 実装 §7c):
`1` で TP carrier 凝縮 (`condGasSpecies` ≥ 0) の NS run に、気相内の分子拡散 (気相基準 z・風上の補正) と全輸送量 (化学種・液 g・Q2/Q1/Q0) 共通の乱流混合 μ_t/Sc_t を
1 回の面流束で足し (エネルギーに Σh_k J_k + h_v J_w − L J_l)、蒸気と液を非分割で更新する (増分は点対角、制限は `condDgMaxStep`/`condDTmaxStep` と蒸気・液の非負)。
**実効状態の判定表** (2026-10-04, plan condensation-two-phase-default §4-2; ソルバ・`check_solver_config.py`・`res_*.h5` 属性・`stage_manifest.json` で共通、Python 側は `tools/twophase_state.py`):
active (包絡内: TP carrier・非平衡・`nCondSpecies 1`・`condLimiterMode 1`・NS・定常陰解法・`passiveScalarScheme 1`・非軸対称・非周期) で指定 ON なら作動。
inactive-a (凝縮 OFF・`viscMethod 0`; dual-time + Euler もここ) は INFO・ビット不変。inactive-b (CPG carrier・pure 凝縮・`condEquilibrium 2`) は毎回 WARNING で実効 0。
unsupported-c (dual-time・陽解法/RK・`speciesImplicitCoupling 2`・`passiveScalarScheme 0`・`condEquilibrium 1`・`condLimiterMode 0`・`nCondSpecies` ≥ 2・**軸対称・周期**) は
指定 ON で起動時にエラー終了 (旧作用素で回すなら `condTwoPhaseDiffusion: 0` を明示)。cell は未検証 (WARNING)。
起動時に `[twophase] condTwoPhaseDiffusion requested <omitted|0|1>, state <…>, effective <0|1>` 行を出し、`res_*.h5` のルート属性
`twophase_diffusion_effective`・`twophase_diffusion_state`・`twophase_diffusion_requested` に記録する (凝縮 run か明示したときだけ)。
`stage_manifest.json` の hard キーに `twophase_diffusion_effective` が入り、OFF の場から ON で継続した段は別区間になる。
**既定はまだ 0**。既定を 1 にした後 (S1-c) は省略も ON と数え、包絡内の明示 0 は `WARNING: legacy operator: liquid not diffused` を出す。
`condTwoPhaseSolver 0` のとき、蒸気・液・Q は `passiveImplicitCoupling`/`speciesImplicitCoupling` の DPLUR ではなく点対角で更新する (非水種は従来どおり)。`condTwoPhaseRelax` (既定 1, 0 < ω ≤ 1) は
蒸気・液・Q の増分の緩和 (前処理の後・制限の前)。#4c の 1D 試験では大きな擬似刻みで核生成の Q0 が周期運動になり ω 0.5 で収束した (高 CFL は保証外)。
計算開始時と終了時に `[twophase-audit]` 行 (格納状態から化学種・蒸気・液・Q の残差を double で組み直した成分ごとの比と、終了時の `VERDICT: PASS | NOT CONVERGED`; エネルギーは対象外) が出る。収束を受け入れる根拠はこの VERDICT と check_convergence の両方。`residual_history.csv` に `rms_roYv` (蒸気の残差 R_w − R_g) が加わり、`check_convergence.py` の検査対象に入る (`rms_roY` 接頭辞)。monitor に `[twophase]` 行 (θ<1 のセル数・最小 θ・保留量・状態補正)。
既定 0 の run の結果・列構成は変わらない。
**`condensation.condTwoPhaseSolver`** (既定: `condTwoPhaseDiffusion 1` のとき **1** (2026-10-04 から、検証済みの組)、0 のとき 0 (使われない); #4g): 二相拡散の蒸気・液・Q の増分の作り方。0 = 点対角 (診断用 opt-in)、1 = 化学種・受動種と同じ緩和整合 scalar-DPLUR
(右辺は全残差、対角は点対角と同じ分母、非対角は流入質量流束、ゼロ開始で `nStepInner` 回、ω = `implicitRelax`; その後 `condTwoPhaseRelax`・θ・commit・再正規化は同じ)。
1 sweep・ω 1 では点対角とビット一致。起動時に `[twophase] condTwoPhaseSolver` 行で実効の implicitRelax・nStepInner・scalarCflMax を出す。
**`condensation.condTwoPhaseNonnegLimit`** (既定: `condTwoPhaseDiffusion 1` のとき **0** (2026-10-04 から、検証済みの組)、0 のとき 1 (使われない); #4h): 0 で二相の非分割更新の共通 θ から蒸気・液の非負制限を外す (θ = dg_max・dT_max だけ)。1 (非負 θ を掛ける) は診断用 opt-in。
2026-10-04 より前に `condTwoPhaseDiffusion: 1` だけを書いた run (Solver 0・NonnegLimit 1 だった) を再現するには `condTwoPhaseSolver: 0`・`condTwoPhaseNonnegLimit: 1` を明示する。
commit は総水分だけ 0 に下限を掛け (`vround` を使わない)、液は再正規化 → 受動種の床 → 実現可能性クランプで固定した総水分に対して 0 ≤ ρg ≤ ρY_w に射影する。
二相拡散 ON の run は更新ごと・成分ごと (ρY_w, ρv, ρg, ρQ2, ρQ1, ρQ0) の補正を `[twophase-corr-gate]` (区間と末尾 ceil(0.1N) 更新の max、κ = 2n_sε₃₂ の VERDICT、段ごとの行) に出す。
**`condensation.condAuditResidual`** (既定 0, 2026-10-02, #1b-pre): `1` で二相拡散 OFF の TP carrier 凝縮 run でも、現行の作用素 (化学種の Fick 拡散、液・Q は移流のみ) を
格納状態から double で組み直す `[twophase-audit]` を開始時・終了時に出す (A/B の A 側の受入用; 二相拡散 ON の run は常に新作用素で監査)。結果は変えない。
凝縮 run の `[cond-corr]` には `renorm components` 行 (再正規化の max|f−1| と ρY_w・ρg・Q2・Q1・Q0 の補正量; 相対は各成分の自分の総量) と
`theta over all updates` 行 (更新の θ<1 と θ_src<1 のセル数を全更新で積算; 区間・累積・直近の更新) が加わる (2026-10-02)。
TP carrier 凝縮 run は再正規化の受入ゲート `[renorm-gate]` を出す: 区間 (前回ログからの全更新) と終了時の末尾窓 (実更新数 N の最後の ceil(0.1N) 更新) について、
更新ごとの成分別相対補正 C_q,n = Σ|q⁺−q⁻|V/Σq⁻V (ρY_w, ρg, ρQ2, ρQ1, ρQ0) と局所係数偏差 max|f−1| の max、終了時は κ = 2 n_s ε₃₂ で `VERDICT`。
`renorm components` 行は記録用 (ゲートではない)。
**`condensation.condTwoPhaseDiag`** (既定 0, 2026-10-02, #1b-r1; 読むだけ・数値は不変): 二相拡散 ON の run で `1` = 末尾 200 更新の θ = 0 のセル
(`2` = θ < 1) について、セル ID・座標・更新前の ρv/ρg・制限前増分・θ を決めた制限・Q 残差・再正規化の f−1 と前後差を記録し、終了時に
`twophase_diag_theta0_cells.csv` (頻度順上位 500) と `twophase_diag_theta0_summary.csv` (固有セル数・更新あたり件数・持続性・制限理由の内訳) を書く。
区間ごとに `[twophase-diag]` 行 (θ 制限のセル×更新数を更新数で割った値と固有セル数) を出す。
`3` (#1b-r2) は 1 に加えて、窓内の各更新で同じ状態・面値・係数・ソース値の double 組立 B (状態は書かない) を作り、乾燥停止セルと成分ごとの残差最大セルについて
R_g の移流/拡散/ソース内訳・前処理分母・組立 A/B の残差・制限理由を `twophase_diag3_cells.csv`、停止集合内外の Q 残差比較を `twophase_diag3_qcompare.csv`、
plan #1b-r2 の判定規則を機械的に当てた要約を `twophase_diag3_summary.csv` に書く。cells CSV の値の列は窓内最大 |値| とその更新番号 (最後の記録は `last_update_*`)。
注意: 受動種/化学種の拡散は `viscMethod != 0` のときだけ加わる (viscMethod 0 は定数粘性ではなく「拡散なし」扱い; 化学種と同じ規約)。

## physProp.viscMethod — 層流の粘性・熱伝導 (2026-09-27 `viscMethod: 2` を置き換え)

| 値 | μ・λ | 備考 |
| --- | --- | --- |
| `0` | 定数 (`visc`, `thermCond` / `thermCondMethod: 1` なら k=μ·cp/`prandtlLam`) | 受動種・化学種の拡散は加わらない (上の注意) |
| `1` | 空気の Sutherland (熱伝導は `thermCondMethod` で選択) | 組成に依らない。CPG・TP どちらでも可 |
| `2` | 種ごとの輸送物性 (`physProp.transport` 必須、CEA 形 frozen 混合則) | `thermalMethod: 2` のみ。[plan §4.3c](../plans/active/thermophysics-solver-owned-species-db.md) |

- **`viscMethod: 2` は `physProp.transport` (実種ごとの出所 `cea` / `kinetic` / `fit` / `custom:<名前>_v<版>`、lump は構成実種ごと) が無いと
  起動時エラー**。以前の `viscMethod: 2` (LJ + Chapman–Enskog の種別値を Wilke の φ で混合し、λ にも同じ φ を使う kinetic 経路) は
  計算から外した。空気の Sutherland で足りるなら `viscMethod: 1` にする。逆に `physProp.transport` を `viscMethod ≠ 2` と併用するのも起動拒否
  (記録と計算が食い違うため)。
- **旧 `viscMethod: 2` の結果を再現するときは旧バイナリで回す** (置き換え前の `e2daaba8` までのコミットでビルドしたもの)。
  新バイナリには旧経路へ戻す設定・環境変数は無い。
- 表引き: `physProp.transport` の μ・λ は既定で ln T の区分 3 次表を float で引く。`FORGE_TRANSPORT_TABLE=0` はその実行だけ同じモデルの
  double 評価に戻す (性能比較・切り分け用; 旧 kinetic 経路には戻らない)。
- `visc` は `viscMethod: 2` でも dt と陰解法対角の剛性見積りに使われるので残す。

## physProp.chemistry — 有限速度化学 (H₂ 燃焼・ノズル化学非平衡)

多成分 TP (`thermalMethod: 2`, `species` ≥2 種) に化学反応ソース項を加える。理論・実装は
[`methods/chemistry.md`](../methods/chemistry.md)、計画は [`plans/active/chemistry-finite-rate-h2.md`](../plans/active/chemistry-finite-rate-h2.md)。

```yaml
physProp: {thermalMethod: 2, species: [H2, O2, H, O, OH, H2O, HO2, H2O2, N2], speciesDBFile: "species_db.yaml",
           thermoHrefTemp: 298.15,
           chemistry: {enabled: 1, mechanismFile: "mech.yaml", jacobianMode: 1, tMaxReaction: 6000.0, freezeBelowT: 0.0}}
```

| キー | 既定 | 意味 |
| --- | --- | --- |
| `enabled` | 0 | 1 で反応ソース $\dot\omega_s$ と反応熱 $\dot Q$ を有効化。0 なら全経路ビット不変 |
| `mechanismFile` | — | 反応機構 (Cantera YAML サブセット: `equation`, `rate-constant {A,b,Ea}`, `type: three-body`, `efficiencies`, `units`)。同梱: `solver_density_cuda/tools/mechanisms/` (Jachimowski 1988 9 種 20 反応 / 13 種 33 反応)。falloff (Troe) は未対応 (Phase 2) |
| `jacobianMode` | 1 | 0: 陽ソースのみ, 1: 対角 point-implicit (`src_jac_Y{s}` に $\max(0,-\partial\dot\omega_s/\partial\rho Y_s)$), 2: 種ブロック (Phase 2, 現状は 1 と同じ) |
| `tMaxReaction` | 6000 | 速度式評価の温度上限 [K] |
| `freezeBelowT` | 0 | この温度未満で反応を凍結 ($\dot\omega=0$) [K]。試験部の低温域で反応評価を省く用 |

- **`mesh.renumber`** (変換時, 既定 `none`, 2026-09-12): `rcm` で `convertGmshToForge` が節点を Reverse Cuthill–McKee で再番号付けする (node では CV 順 = gather の局所性)。A10G 3D 2.37 M 節点で 1 step −3.3 %。順列は `/MESH/RENUMBER_PERM` (new→old) に入り、旧番号の `res_*.h5` は `tools/permute_res_h5.py SRC_res.h5 DST_mesh.h5` で移植する (index コピー不可)。
- **`thermoFloat`** (physProp, 既定 1, 2026-09-12): thermally-perfect の温度反転を float Newton + double 研磨 (収束まで最大 3 段, 通常 1 段) のハイブリッドにする (厳密参照との誤差 ≤1e-8·T, 3D 2.37 M 節点で 1 step −13 %)。`thermoHrefTemp>0` が前提で、datum 無しなら自動で 0 (従来 double Newton) に落ちる。`0` で従来経路。plan performance-3d-node-sst-speedup §4.2-3。
- **`condFloat`** (condensation, 既定 1, 2026-09-13): 凝縮経路の float 実体。相変化ソース (物性は double で作った区分 3 次表を float で評価、範囲外の湿潤セルは旧 double 実体へ退避、核生成は対数空間で ln J ≤ ln 1e35 の上限を本体・摂動の両方に)、面潜熱 (表, g=0 の面はスキップ)、実現可能性 clamp、二相の温度反転 (float Newton + double 研磨、残差 1e-9|e|+0.05 J/kg で成功判定、失敗セルは roe 不変)。平衡形と二温度: ソース kernel は `condEquilibrium`≠0 または `condTwoTemp`=1 で常に double、EOS 反転は `condEquilibrium: 2` (拘束形) だけ double で `condEquilibrium: 1` の湿潤 TP セルは非平衡と同じハイブリッド反転 (面潜熱・clamp・組成構築は常に float 経路)。`0` で従来 double 経路 (ソース診断は変更前とビット一致、二相反転の T は成功判定の統一で ≤1e-3 K 変わり得る)。plan condensation-float-speedup §4.2-6、methods/condensation.md 実装 §9。
- **`condensation.condensationSpecies`** (condensation, 既定なし, 2026-09-16): 凝縮する気相化学種を**名前**で指定する
  (`condensationSpecies: H2O`)。`physProp.species` から大文字小文字無視で index を解決して `condGasSpecies` に入れる。数値
  `condGasSpecies` を併記して食い違えばエラー、数値だけなら範囲検査 (`nSpecies` 超え・単一種で carrier 形はエラー)。起動ログの
  `[species]` 表に `condensing species: H2O (condGasSpecies=1)` と出る。種順序を変えても config を書き直さずに済むので名前を正本にする
  ([plan cea-mole-fraction §2](../plans/accepted/thermophysics-cea-mole-fraction-species.md))。
- **`physProp.tracer: exhaust`** (physProp, 既定 `none`, 2026-09-16; 2026-09-17 改定): 受動トレーサ `roXi` (排気率 ξ∈[0,1])。既定 `passiveScalarScheme 1` では
  **化学種の輸送経路の受動種**として移流 (`speciesFaceReconstruction 2` なら 2 次面再構成)・Fick 拡散 ($D=\mu/(\rho Sc)+\mu_t/(\rho Sc_t)$, 粘性 run のみ)・
  陰解法更新 (`passiveImplicitCoupling` 0 = point-implicit × `passiveImplicitRelax` / 1 = scalar-DPLUR)・**dual-time の BDF 物理時間項** (履歴 `roXiP/PP`) を持つ。
  `passiveScalarScheme 0` (旧経路) では汎用スカラの 1 次風上・拡散なしで、`time.dualTime != 0` との併用は config でエラー (物理時間項なし)。
  SERN で排気/外気の見分けに使う: 輸送種に純粋な流入元ラベル (排気入口 1・外気入口 0 の種, 旧 `[EXH, AIR]` の Y_EXH) があればそれを ξ に使い、無ければ (`full`、`lumped`+`keep`) このトレーサを輸送する (`species_meta.yaml` の `exhaust_fraction` が正本, `methods/thermophysics.md` §5)。
- **bcond `floats: {X0:.., X1:.., ...}`** (多成分 TP の入口, 2026-09-16): 入口組成を**モル分率**で与える。forge が double で検証し
  $Y_k = X_k M_k / \sum_j X_j M_j$ (MW は `speciesDBFile`/内蔵 DB) に換算して従来の `Y{s}` 経路へ流す。**X を 1 つでも書いたら全種必須**
  (既定補完しない)、同じ境界での `X`/`Y` 混在・負値・非有限・総和 0・範囲外 index はエラー終了。`Y{s}` を明示した場合も負値と
  |ΣY−1|>1e-3 はエラー (以前は黙って通した; 未指定種は従来どおり `Y0=1`, 他 0)。起動ログに入口ごとの Y と X (MW から逆算) が出るので
  ここで桁を確認する。`initial` (IC) は文字列のまま; 組成付き IC と `inletProfile` CSV は生成ツール側で Y に換算する
  (`gen_inlet_profile.py --X`, [procedures/inlet-profile.md](inlet-profile.md))。
- **内蔵の化学種** (2026-10-01, [plan #13-2](../plans/active/thermophysics-solver-owned-species-db.md)): `physProp.species` の名前は共通データ
  `solver_density_cuda/data/species/forge_species_v1.yaml` の**全気相種 61 種** (移行前からの N2/O2/CO2/H2O/Ar/He/AIR、CO/H2/OH/H/NO/O、
  CEA `thermo.inp` から生成した 48 種 — 例 `CH4`・`C2H2,acetylene`・`N`・`NH3`・`Kr`・`Xe`) で `speciesDBFile` なしに解決できる。一覧と出典は
  [methods/thermophysics.md](../methods/thermophysics.md)「内蔵 species DB の一覧と出典」。`e-` と CEA の `Air` は内蔵に無い (`Air`/`air` は擬似種 `AIR` の別名)。
  生成 48 種の多くは **LJ を持たない**: `physProp.transport` で `kinetic` を指定する、または LJ の混合平均拡散 (化学種 2 以上・`viscMethod` ≠ 0・
  `speciesDiffusionMethod: 1` = 既定) に使うと起動時エラー (`cea`/`fit` か `speciesDiffusionMethod: 0`、または `speciesDBFile` で LJ を与える)。
  `speciesDBFile` による上書き・追加は従来どおり (同じキーだけ上書き)。既存 config の解決結果・互換性ハッシュは変わらない。
  (LJ の有無は次項の `physProp.ljSource` で決まる。既定でどの集合にも LJ が無いのは D2・D2O・H6F6・N2O4 の 4 種。)
- **`physProp.ljSource`** (2026-10-01, [plan §4.10 #14](../plans/active/thermophysics-solver-owned-species-db.md)): 内蔵種の Lennard-Jones パラメータ
  (σ, ε/k_B) を探す**集合の順序付きリスト**。先頭から探して最初にある集合の値を使う。lump の構成種も同じリストで解決し、
  `speciesDBFile` の種はそのファイルの LJ が優先 (リストに関係なし)。
  - 集合: `gri30` (GRI-Mech 3.0 の transport; Cantera 同梱 `gri30.yaml`)、`svehla1962` (Svehla 1962 NASA TR R-132 Table I(a); 希ガスは粘性フィット行)、
    `legacy_v1` (#14 前の内蔵値の凍結 — GRI と Svehla の混在; 旧 run の再現用)。
  - **既定 (無指定) は `[gri30, svehla1962]`** (2026-10-01 ユーザ決定)。#14 前の値に対して変わるのは H2・OH・H・O・NO・CO (GRI 値へ) と、
    以前 LJ が無かった CEA 生成種 (GRI か Svehla の値が入る)。N2/O2/CO2/H2O/Ar/He/AIR/N/NH3/NO2/N2O は値が同じ。
  - 旧 run を同じ物性で再現・継続するときは `physProp: {..., ljSource: [legacy_v1]}` を明示する (#14 前と同じ互換性ハッシュ)。
  - 効くのは LJ を読む経路だけ: LJ の混合平均拡散 (化学種 2 以上・`viscMethod` ≠ 0・`speciesDiffusionMethod: 1`) の D_ij と、
    `physProp.transport` の `kinetic` の μ・λ。Euler (`viscMethod: 0`)・`cea`/`fit` 輸送・定数 Schmidt 数の run では値は使われない
    (ただし互換性ハッシュには LJ の値が入るので、値が変わる種を含む config はハッシュが変わる)。
  - どの集合にも無い種を LJ を読む使い方に回すと、種名と探した集合を示して起動時エラー。空リスト・未知の集合名・重複もエラー。
  - 解決結果は起動ログ (`[species]   LJ sets searched ...` と種ごとの集合・値) と解決済み記録の `provenance.lj_source`・`lj_resolved` に出る
    (互換性ハッシュには値だけが入り、集合名は入らない)。双極子 (H2O・NH3) は種レベルの値だが、`kinetic` で適用するのは
    `stockmayer` の集合 (`gri30`・`legacy_v1`) で解決した種だけ。`svehla1962` (`lj12-6`) で解決した種には適用しない (起動ログに NOTE; #14-L1b)。
- **化学種の解決済み記録と入力場の照合** (TP `thermalMethod: 2` のみ、2026-09-27、[plan §4.3 #3a](../plans/active/thermophysics-solver-owned-species-db.md))。
  - 起動時に使用した全種の物性 (順序・名前・MW・datum 前の絶対係数と温度区間・LJ・`thermoHrefTemp`・来歴) を run ディレクトリへ
    `resolved_species_<互換ハッシュ16桁>.yaml` として書く (**出力=記録であり入力ではない**; 同名で来歴だけ違えば `_<完全性16桁>` 付きの別名)。
    各 `res_*.h5` (境界出力を含む) のルート属性に `species_hash` (互換性ハッシュ; `source` とパスは含まない)・`species_record_sha256` (記録全文)・
    `species_record_file`・`species_input_unverified` が付く。Python からは `tools/forge_species.py` の `load_record` / `find_record`。
  - `valueFileName` に `species_hash` があり自分と違えば**起動を拒否**し、入力側の記録が見つかれば差のある種・係数を表示する (許可手段なし;
    種を変えるなら `tools/convert_species_field.py`)。**属性が無い場は「照合不能」で拒否**する。確認済みなら**その実行だけ**
    `FORGE_ALLOW_UNVERIFIED_SPECIES=1 forge` で通せる (出力に `species_input_unverified=1` が付き、以後の restart に継承される)。
    **config キーでの恒常的な許可は無い** (生成 config に埋め込まない)。
    **属性なしの場は既定で停止する** (2026-09-27 #3c で過渡期の既定を終了。過渡期用の `FORGE_REQUIRE_VERIFIED_SPECIES` はソルバから撤去し、
    立てても何も変わらない)。停止したら (1) 属性を付ける生成処理で IC を作り直す (IC 生成は `forge --resolve-species` で宛先を解決して
    `species_hash` を付ける; 付いた場から `restart_field.py` / `interp_field.py` で作った場は属性を継承する) か、(2) 同じ種・datum で作った場だと
    確認できているなら**その実行だけ** `FORGE_ALLOW_UNVERIFIED_SPECIES=1 forge` で通す (起動ログの警告は環境変数で許可された旨を示し、
    出力に `species_input_unverified=1` が付く)。種・DB が変わった場は `convert_species_field.py`。
  - **引き継ぎツールも同じ規約** (2026-09-27 #3c 残)。`restart_field.py` / `interp_field.py` / `convert_species_field.py` と設計 runner の
    IC 付与・段間継承 (`runner_axismach` の `_restart_same_mesh`・`runner_sern` の `restart_by_index` / `warm_from_run`) は、
    SRC が未検証 (属性なし / ハッシュが宛先と違う `species_input_unverified=1`; 一致する印付きは下記) で宛先が TP のとき・宛先の `solverConfig.yaml` が無く CPG か TP か判定できないとき・
    宛先を解決できない (`--resolve-species` を持つ forge が無い = 旧バイナリ; `--forge` / `FORGE_BIN` で新しいバイナリを渡す) とき・
    記録が壊れているときに**既定で書き込まずに停止**し、ソルバと同じ 2 通り (IC を属性を付ける処理で作り直す / その実行だけ許可) を案内する。
    許可は**その実行だけ** `FORGE_ALLOW_UNVERIFIED_SPECIES=1 python3 tools/restart_field.py ...` か各ツールの `--force-species`
    (記録の破損・係数不一致を通すのは `--force-species` だけ; 環境変数は未検証・解決不能だけを通す)。**許可して通した宛先には属性を付けない**ので、
    その場から起動するソルバにも同じ許可が要る。検証済みの SRC からの継承は従来どおり (宛先を解決して互換性ハッシュが一致すれば属性と記録を継承)。
    **印付きの場 (`species_input_unverified=1`) はソルバと同じ規則** (2026-10-01 #3d): 場の `species_hash` が宛先のハッシュ (種変換は SRC run の設定の
    ハッシュ) と一致すれば許可なしで通し、宛先に同じハッシュと印を継承する (印は消さない)。一致しなければ停止。
    したがって属性なしの旧場の移行は「その実行だけ `FORGE_ALLOW_UNVERIFIED_SPECIES=1 forge` で 1 回回す (出力が印付きになる) → 以後は許可なし」で済む。
    宛先が CPG なら対象外。過渡期の `FORGE_REQUIRE_VERIFIED_SPECIES` はツールからも撤去した。
  - **凝縮 ON の気液ペア検査** (2026-10-01 #13-3): 外部 DB の気相 H2O が内蔵と違う (例: #9 以前の runner が写した旧 MW 0.0180153 の H2O) と
    起動を拒否し、キーと両値・移行先 (外部 DB から H2O を外す / lump 記法で `speciesDBFile` をやめる) を示す。`FORGE_ALLOW_UNVERIFIED_SPECIES=1`・
    `--force-species` では通らない。その config の旧場からの継続は `convert_species_field.py` (種の物性が変わるため)。
  - `forge --resolve-species`: GPU を使わず `solverConfig.yaml` を解決して記録を書き、互換性ハッシュを**標準出力の最終行**に出して終了
    (CPG は終了コード 2)。記録を既存の場へ貼っても検証済みにはならない。
- **`thermoHrefTemp: 298.15` を必ず指定する** (反応熱は sensible datum の残差項 $\dot Q=-\sum_s h^{abs}_s(T_{ref})\dot\omega_s$ として入る。絶対 datum (0) でも動くが陰解法は不安定)。
- 機構に現れる種は `species` に全て含めること (無ければ起動時エラー)。`species` にだけある種は不活性として扱う。
- 熱力学 DB は `tools/cea_thermo_to_species_db.py thermo.inp --species ...` で CEA から生成する (内蔵 DB は H・O・OH・H2 までは持つが HO2・H2O2 などは無い; 下の「内蔵の化学種」)。
- 出力: `chemQdot` [W/m³], `chemTau` [s] (=1/max|∂ω_s/∂ρY_s|、化学時間の目安。`dt` や `cfl_pseudo` の妥当性判断に使う)。
- 検証: `case/35.uniform_periodic_box/run_0049_node_h2_ignition` (0-D 着火 vs Cantera)。

**TP (`thermalMethod: 2`) × 亜音速 `outlet_statPress` の γ 混用は修正済 (2026-09-07 移植)**: 出口列が指定 Ps を無視して
一様に加圧/逆流する症状 (node では出口壁ノードが常に亜音速なので必ず通る) はこのバグの指紋
([plan §2.12](../plans/active/boundary-node-nozzle-wall-outlet-stability.md))。

## turbulence.sst* — SST 整合オプション (2026-09-08, codex レビュー由来)

```yaml
turbulence: {model: "sst", ..., sstNodeWallKPin: 1, sstOmegaProdFromPk: 1, sstSigmaBlend: 1, sstIsotropicStress: 0, sstEnergyKSource: 0, sstEnergyIncludesK: 0}
```

| キー | 既定 | 意味 |
|---|---|---|
| `sstNodeWallKPin` | 1 | node 低 Re 壁 (`wallTreatmentSST: 0`) で壁ノードの k/roK を 0 にピンし、壁ノードの k/ω 残差を 0 化する。旧 (0) は ghost 反射のみで、node は境界半割面の拡散を skip するため壁 k=0 が実際には効いていなかった。cell には無関係 |
| `sstOmegaProdFromPk` | **1** (2026-09-08 から) | 1: P_ω = α P_k/ν_t (k 側のリミッタ・dilatation・壁関数置換後の P_k と整合。NASA TMR は 2003 論文の αρS² を誤植とし、訂正式をこの形としている)。0: P_ω = α ρ S² (旧既定)。リミッタ・等方項・壁関数置換の非発動域では同一。ν_t は `vis_turb` (μt<1e-6μ では 0 側にフォールバック)。等方項 (dilatation 2) はリミッタ**前**に加算 (2026-09-08 変更, キー非依存) |
| `sstSigmaBlend` | **1** (2026-09-08 から) | 1: 拡散係数 σ_k = F1·0.85 + (1−F1)·1.0, σ_ω = F1·0.5 + (1−F1)·0.856 (正式ブレンド)。F1 は同 step の前処理 `ransBlendF1` が `sstF1` に書く (拡散の前に k/ω 勾配→F1 を評価する順序に変更、ラグ無し)。0: k-ω 側定数 (旧既定)。壁近傍 (F1=1) では同一 |
| `sstIsotropicStress` | 0 | 1: 運動量/エネルギーの応力に等方 Reynolds 応力 −(2/3)ρk δᵢⱼ を加える (内部面のみ; `dilatationCorrection: 2` の k 生産側 −(2/3)ρk∇·u と対をなす)。0: 無し (現行) |
| `sstEnergyKSource` | 0 | 1: エネルギー式に −(P_k − D_k)·V を源として加える (E = e + u²/2 に k を含まない定式化で、k に溜まる分を平均流エネルギーから引く。SU2 型 E に k を含める方式の代替)。0: 無し (現行)。**`sstEnergyIncludesK: 1` のときは無効化される (後継)** |
| `sstEnergyIncludesK` | 0 | 1: **全エネルギーに乱流運動エネルギーを含める** (E_t = ρ(e + u²/2 + k), SU2 形; plan [turbulence-sst-energy-includes-k](../plans/accepted/turbulence-sst-energy-includes-k.md))。実装は「分割保持」: 保存量 `roe` は平均流 E_m のまま (T/P/IC/restart は不変)、エネルギー残差を E_t の流束 (面エンタルピー +(5/3)k、圧力流束 p* = p + (2/3)ρk、k 拡散 (μ+σ_k μt)∇k·S、軸対称 hoop p*) で組み、k 更新後に `roe -= Δ(ρk)` とする。周期箱の一様減衰で ΣV(E_m+ρk) が 1e-7 で保存し、k の散逸がそのまま T に戻る。1 のとき `sstIsotropicStress`/`sstEnergyKSource` は強制 0。SLAU/SLAU2 のみ。後処理の全温は T + u²/2c_p **+ k/c_p** で評価する。0: 従来 (平均流 E_m のみ保存、k の散逸は熱に戻らない) |

検証 (plan [turbulence-sst-consistency-options](../plans/active/turbulence-sst-consistency-options.md)): case/16 2D node SST では全オプションが壁圧比 ≤0.25 % の差、case/26 平板 Cf と case/16 3D の結果は同 plan 参照。
`tools/test_scale_invariance.py` (相似メッシュで無次元残差の一致を見る) は単位付き閾値の検出用。

### 遷移モデル `turbulence.transition` (2026-09-22)

```yaml
turbulence: {model: "sst", ..., wallTreatmentSST: 0, transition: "lm2009"}   # 既定 "none"
```

| 値 | 内容 |
| --- | --- |
| `none` (既定) | 遷移モデルなし。SST は前縁から完全乱流 |
| `lm2009` | Langtry–Menter 2009 の $\gamma$–$\tilde{Re}_{\theta t}$ 2 方程式 (SU2 の `KIND_TRANS_MODEL= LM` と同じソース式・相関)。輸送変数の上下限は SU2 と違う: forge は $\gamma\in[10^{-4},1]$・$\tilde{Re}_{\theta t}\ge20$、SU2 は $\gamma\in[10^{-4},5]$・$\tilde{Re}_{\theta t}\ge10^{-4}$ (T3A の準定常場で下限の作動は 0 %、$\gamma=1$ は自由流の 7 %)。理論は [`methods/turbulence/theory.md`](../methods/turbulence/theory.md) §11、実装は同 `implementation.md` |

- **受け付ける組み合わせは検証したものだけ**: node・`model: sst`・`wallTreatmentSST: 0`・非軸対称・DES なし・`sstEnergyIncludesK: 0`・`scalarDiffusion: 1`・定常陰解法。外れると起動時に止まる。
- 入口は $\gamma=1$、$\tilde{Re}_{\theta t}$ は入口の $k$ と速度から決まる局所 $Tu$ の相関値。**追加の入口キーは無い**。入口の $k$・$\omega$ (乱れの強さと減衰の速さ) が遷移位置を決めるので、
  $\omega$ は「入口粘性比」の感度を必ず見る。
- **`turbulence.transitionRethMin`** (既定 **20** = LM2009/SU2 の推奨値 `Corr_Ret_lim`): 相関値と輸送変数 $\tilde{Re}_{\theta t}$ の下限。上げると遷移が遅れる。
  文献には case ごとに 100〜200 へ上げる調整例がある (Lin ら JGPP 6(3):9-15, 2014) が、**Mark II で 20/130/200 を試したところほぼ空振り**だった
  (下限に触れるのは衝撃下流の $\gamma\to1$ の領域だけで、層流域の輸送値は 130 以上。plan §6.5)。**20 以外にしたら run の README に「実験に合わせた調整」と明記する**。
- 壁は $y_1^+\le1$ が前提 (モデルの定義)。残差は `rms_roGamma` / `rms_roReth` が増える (`check_convergence.py` が読む)。
- 出発場は完全乱流の SST 場でよい (平板 T3A で、$k$/$\omega$ を入口値に戻した場からの解と 4 桁一致)。
- 検証: `case/57.transition_flat_plate` (ERCOFTAC T3A、同一メッシュの SU2 LM と比較)。

## output — 出力する場の量の絞り込みと h0 (2026-09-08)

```yaml
output: {level: 1, extraFields: [ducros, volume]}   # 省略時 level 1
```

| level | 出力 | 目安 (2.08M 節点 SST) |
|---|---|---|
| 0 | リスタート最小: `ro roUx roUy roUz roe roK roOmega roY*` + 凝縮モーメント保存量 | ~190 MB (うちメッシュ CONNE/COORD 98 MB) |
| **1 (既定)** | 0 + 原始量 `P T Ux Uy Uz k omega sonic Y*`、`vis_lam vis_turb wall_dist`、**`h0`** | ~290 MB |
| 2 | 従来どおり全診断量 (勾配 27 本・リミッタ・SST/DES 診断・`cfl dt_local ducros volume` …) | ~680 MB |

- `extraFields` で level 0/1 に個別追加する (名前は `variables.hpp` の `output_cellValNames`)。旧解析スクリプトが読む診断量
  (`volume`: case/09 の保存則解析、`cfl`: `sweep_cfl_implicit.py`、`ducros`: `analyze_ducros.py`) はここで足すか level 2 にする。
- **`h0`** = 全エンタルピー (単位質量) = e + p/ρ + u²/2 で、`sstEnergyIncludesK: 1` のときだけ + k を含む (属性 `h0_includes_k`)。
  **全温・全圧の後処理は `VALUE/h0` から作る** (CPG: T₀ = h₀/c_p, semi-perfect: h(T₀)=h₀ の逆算, P₀ は s°(T₀)−s°(T)=R ln(P₀/P))。
  スクリプトで `T + u²/2c_p` を自前で組まない (k を含めるかどうかの判断がソルバ側に閉じる)。逆算は
  **`solver_density_cuda/tools/total_quantities.py RUN_DIR [--res res_N.h5] [--write] [--Tt ..]`** (CPG 閉形式 / TP は NASA-9 Newton、
  P₀ は s° 差、凝縮は凍結組成の気相のみ・警告) を使う。`--write` で `VALUE/T0`, `VALUE/P0` を res に追記 (属性 includes_k / method)。
  Python からは `total_state(run_dir, res_path)`。
- 原則: **後処理で導出できる量はソルバから出さない**。勾配・リミッタ・診断は必要な run だけ level 2 で取る。

### 本段の監視用の出力と協調停止 (2026-10-15)

長い本段を「全場を出し続けず、境界・選んだ節点・流量で監視し、止める step の全場だけを確実に残す」ための出力。
どれも**出力だけ**で数値 (流束・勾配・リミッタ・境界条件・陰解法・ソース) には触らない。キーを書かなければ従来と同じ出力
(協調停止と随時の書き出しは常に有効だが、要求のファイルが無ければ毎 step `stat` を 2 回するだけ)。設計は
[plans/active/tooling-nozzle-card-local-mesh.md](../plans/active/tooling-nozzle-card-local-mesh.md) §4.4 1′・§5.1 #4a と
[plans/active/tooling-settle-convergence-and-monitor.md](../plans/active/tooling-settle-convergence-and-monitor.md) §4.1 の 6 (随時の書き出し)。

```yaml
time:
  outStepInterval: 2000          # 境界の h5・point probe・節点の出力・流量の行
  fieldOutStepInterval: 20000    # 全場 res_<step>.h5 だけ (省略時 = outStepInterval)
output: {nodeProbeFile: nodes.txt, boundaryMassflux: 1}
```

要求カードの雛形 (`axismach_final_v1`・`axismach_explore_v1`) の NS は、問題 YAML の `evaluate.ns_output`・`evaluate.settle` から
`runner_axismach.prepare_ns` が**本段の設定だけ**にこれらを書く (間隔 2,000 / 全場 20,000、入口・出口・壁・軸の境界の h5、流量、
判定の基準 `settle_criteria.yaml` の節点の表 `settle_nodes.txt`、等温壁なら `interfaceDiag: 1`; 2026-10-16、
plan tooling-settle-convergence-and-monitor §4.5)。前の段 (S1・S2) の設定には入らない。

| キー | 既定 | 効き方 |
| --- | --- | --- |
| `time.fieldOutStepInterval` | `outStepInterval` | 全場 `res_<step>.h5` の間隔。`outStepStart` は共通。**最後の step・協調停止の step・随時の書き出しの step は間隔に関係なく書く**。1 以上。 |
| `output.nodeProbeFile` | (なし) | 節点の番号の表 (1 行 1 つ、`#` 以降と空行は無視)。`outStepInterval` ごと・最後の step・協調停止と随時の書き出しの step に `res_nodes_<step>.h5` を書く。**node の離散化のみ**。 |
| `output.boundaryMassflux` | 0 | 1 で `boundary_massflux.csv` に、節点の出力と同じ step に境界ごとの質量流束の和を 1 行。node・`gpu: 1`・(軸対称なら) `axisymMethod: 0` のみ。 |

**`res_nodes_<step>.h5`**: `MESH/COORD` (その節点の座標、全場と同じ並び x,y,z の 1 次元)、`MESH/NODE_ID` (表の順の節点の番号)、
`VALUE/ro Ux Uy Uz P T sonic` (全場と同じ変数・同じ型の、その step の節点の値)、属性 `step`・`nodeProbeFile`。
**同じ step の全場の同じ節点とビット一致する** (case/71 の Euler で 486 節点 × 10 時点を確認)。表の誤り (数でない行・
1 行に 2 つ・範囲外・重複・節点 0 個・ファイルが無い) は起動時に止まる (終了コード 1、`[output] 起動時の検査で停止: ...`)。
節点の番号は node の CV 番号 = 格子の節点の番号 (要求カードのノズルの構造格子なら `i·nj + j`)。

**`boundary_massflux.csv`**: 先頭の `#` 行に定義、続いて見出し `step,mdot_<名前>_<physID>,raw_<名前>_<physID>,...`
(周期境界は除く)。`raw` はその境界の面 (node: 境界の半割の双対面) の `massflux` = $|S_f|\rho_b(\mathbf U_b\cdot\mathbf n_f)$ の
面の順の double の和で、**ソルバが境界の流束に使った値そのもの** (`convectiveFlux_boundary_d`、$\rho_b,\mathbf U_b$ は境界の状態 bvar)。
その step の最後の流束評価の値 (定常の陰解法では、その step の更新に使った step の頭の状態の流束)。**符号は外向きが正
(流出が正、流入が負)**。軸対称 (`axisymMethod 0`) の面ベクトルは r で重みづけした 1 ラジアンあたりなので `raw` は [kg/s/rad]、
`mdot = 2π·raw` が実寸 [kg/s]。非軸対称は `mdot = raw` (2 次元の格子なら奥行き 1 m あたり)。`step` は全場・境界の h5 と同じ番号
(`residual_history.csv` の `step` 列は step−1)。確認 (case/71 の Euler、466,570 節点): 面ごとの `massflux` の書き出し
(`FORGE_DUMP_MASSFLUX`) を格子ファイルの境界の面の一覧から別に足した値と、入口・出口・壁・軸とも**ビット一致**。
全場の折れ線の積分 (`stop_rule.mass_flows`) との差は相対 $\le 3\times10^{-7}$ (診断値)。軸の境界は 0、slip 壁は ~1e−6 kg/s。

**協調停止** (`STOP_REQUEST` → `STOPPED_AT`):

1. 見張りは run dir (forge の作業ディレクトリ) に `STOP_REQUEST` を置く。中身は停止の理由 (空でもよい、先頭 4 KiB を記録)。
   途中まで書いたファイルを読まれないよう、**一時ファイルに書いてから `mv`** する。
2. forge は毎 step 末尾の出力の前に `stat` で確かめ、見つけた step で全場・境界の h5・節点の出力・流量の行 (ソルバ内 CHT の run は
   固体の状態も、下記) を**間隔に関係なく**書き、log に `[coopStop] step N: STOP_REQUEST を検出 ...` を 1 行出し、その step で
   main ループを抜けて通常の終了処理に進む。
   残差履歴の最後の行は `step` 列 N−1 (全場 `res_N.h5` と同じ step を記録した行; 完走のときと同じ規約)。
3. 終了処理の最後に、残差履歴を flush し、停止の step の書き込みで例外が出ていないこと・残差履歴と流量の CSV の最後の行の
   step・書くべきファイル (全場と `outputHDFflg: 1` の境界の h5/xmf、節点の h5、CSV、CHT の固体の状態) がディスクにあって空でないことを確かめ、
   **すべて通ったときだけ** `STOPPED_AT` (JSON: `status`・`step`・`residual_step`・`request`・`nStepOuter`・`files` [{`path`, `bytes`}])
   を一時ファイル経由で置く。終了コード 0。
4. **確認・通知に失敗したら `STOPPED_AT` を置かず終了コード 3** (理由は `[coopStop] ERROR:` の行)。見張りは `STOPPED_AT` と
   終了コード 0 がそろったときだけ、その step の場で判定・抽出する。

- forge は `STOP_REQUEST` を消さない。**起動時に `STOP_REQUEST`・`STOPPED_AT`・`SNAPSHOT_REQUEST` のどれかが既にあると止まる**
  (終了コード 1、前の run の残りで 1 step 目に止まる・見張りが完了と取り違えるのを防ぐ)。同じ run dir で起動し直すなら消してから。
- 最後の step の出力の後に置かれた `STOP_REQUEST` は見ない (通常の完走、`STOPPED_AT` なし)。
- **ソルバ内 CHT の固体の状態** (`conjugateWall::writeConjugateState`、2026-10-16): 協調停止・随時の書き出しの step でも、通常の出力の
  step と同じもの (壁温 `conjugate_Tw_<physID>.csv`、mode fem2d は固体場 `res_solid_<physID>_<step>.h5/.xmf` と再開用の状態
  `conjugate_state_<physID>.h5` (属性 `step` = 累積 step)) を書き、`STOPPED_AT` / `SNAPSHOT_DONE` の `files` に載せて存在と大きさを確かめる。
  その step の書き込みの失敗は停止なら終了コード 3 (`STOPPED_AT` なし)、随時の書き出しなら `SNAPSHOT_FAILED` (通常の出力の step では
  従来どおり警告だけ)。確認 (case/65 C1 n16 の写し、warmup 100・interval 50): 格子でない step 1267 の停止で壁の h5 の `step`・`step_abs`、
  `conjugate_state_5.h5` の `step` がどれも 1267、`res_solid_5_1267.h5` の `T` = 再開用の状態の `SOLID/T` (ビット同一)、壁の h5 の `Ts` と
  壁温 CSV の差 ≤ 5e−9 K (CSV の 10 桁)。step 724 の随時の書き出しも同じ (書き出しの直後に写した再開用の状態の `step` 724)。
- `tools/check_field_regress.py --boundary` は `res_<名前>_<physID>_<step>.h5` (physID が整数) だけを境界の出力として拾う
  (`res_nodes_<step>.h5`・detectNaN の `res_nan_<step>.h5` は拾わない; 2026-10-16 まではどちらも境界 `nodes`・`nan` として比べていた)。

**随時の書き出し** (`SNAPSHOT_REQUEST` → `SNAPSHOT_DONE` / `SNAPSHOT_FAILED`、2026-10-16): `fieldOutStepInterval` を大きくした本段でも、
見たいときにその step の全場を取る。**計算は止めない**。

1. 見張り (または人) は run dir に `SNAPSHOT_REQUEST` を置く (中身は理由、空でもよい; 一時ファイルに書いてから `mv`)。
2. forge は協調停止と同じ毎 step の `stat` で見つけ、その step の全場・境界の h5 (`outputHDFflg: 1`)・節点の出力・流量の行
   (ソルバ内 CHT の run は固体の状態も) を間隔に関係なく書く。出力の格子の step と重なっても各出力は 1 回だけ (流量の CSV の行も 1 行)。
3. 書いたもの (全場と境界の h5/xmf、節点の h5、流量の CSV の最後の行の step、CHT の固体の状態) をディスク上で確かめ、成功なら `SNAPSHOT_DONE` に
   JSON を 1 行書き足す: `{"step", "pid", "request", "files": [{"path", "bytes"}]}` (log に `[snapshot] step N: ...` を 1 行)。
4. 書き込みに失敗しても計算は止めず、`SNAPSHOT_FAILED` に `{"step", "pid", "request", "errors": [...]}` を 1 行残す。
5. **記録の行を書いた後に `SNAPSHOT_REQUEST` を消す** (失敗でも消す: 同じ失敗を毎 step 繰り返さない)。見張りは「要求が消えた」を
   見たら `SNAPSHOT_DONE` / `SNAPSHOT_FAILED` の最後の行 (自分の `pid`・`step`) を読めばよい。要求を消せなかった run では以後の要求を見ない
   (`SNAPSHOT_FAILED` にその旨)。
6. `SNAPSHOT_DONE` / `SNAPSHOT_FAILED` は書き足しなので、同じ run dir で起動し直すと前の run の行が残る (行の `pid` で区別する)。
7. その step の残差の行 (`step` 列 N−1) を**記録の行の直後に flush する** (2026-10-16; 以前は `monitorInterval` ごとだけ)。
   `SNAPSHOT_DONE` の時点ではまだ書かれていないことがあるので、見張りは残差の行 N−1 が現れるのを待ってから読む
   (`monitorInterval` 1000 の CHT の run で、記録の直後に残差の最後の行が `723,-1,outer_end` (N = 724) になることを確認)。
- 確認 (case/71 の Euler、466,570 節点): 格子でない step 163・401 の要求で全場・境界・節点・流量が書かれ計算は 600 step まで続く
  (終了コード 0、節点の値は全場とビット一致)。格子 (`outStepInterval: 1`) と重なった step でも流量の行は各 step 1 行。
  全場を書けない (同名のディレクトリ) 要求は `SNAPSHOT_FAILED` 1 行・要求は消え、計算は続いて次の要求は成功。要求の後の場と要求なしの run の差は、
  要求なしの run どうしの再実行の差と同じ大きさ (node の残差の atomicAdd で再実行はビット一致しないので、ビットの比較はできない)。

**出力の処理が数値の状態を変えないことの確認** (診断、既定 off; plan tooling-settle-convergence-and-monitor §5.1 #1c・§6 V1 (f′)(g′)):
`FORGE_DIAG_OUTPUT_INVARIANCE=<K>` を付けて起動すると、step K (通常の出力の step でも最後の step でもないこと; 違えば終了コード 2) の
出力の処理の直前の状態 S0 を退避し、S0 に戻しては出力の処理だけを 4 つの腕 (A 要求なし・B `SNAPSHOT_REQUEST` で節点の出力と流量の積算を
止める・C `SNAPSHOT_REQUEST` でキーのまま・A2 要求なし) で実行して、device の全配列 (`var.c_d`・`var.p_d`・境界の `bvar_d`/`bint_d`・
FP64 影アキュムレータ・ライン陰解法の配列) と host の写し (`var.c`・`var.p`・境界の `bvar`/`bint`/`diagVar`・CHT の固体の状態) を S0 と
ビットで比べ、`output_invariance_<K>.txt` (配列の一覧と大きさ・腕ごとの差) を書いて終了する (0 = 腕の間で全配列が同一、4 = 差あり)。
`FORGE_DIAG_OUTPUT_INVARIANCE_CONTINUE=restore` は S0 に戻して計算を続け、`=poison` はさらに「腕の間で違った host の配列」を NaN で
埋めて続ける (数値更新がその写しを読むかの検査; 結果は `check_field_regress.py` で要求なしの反復と比べる)。
- 結果 (2026-10-16、float のビルド): case/71 の Euler (466,570 節点、K 163)・B0 の NS (442,520 節点、ライン陰解法・SST・`interfaceDiag: 1`、
  K 163)・case/65 の CHT fem2d (K 263) のどれも **device は 4 つの腕の間で全配列がビット同一** (307・329・311 配列)。host は全場・境界の
  出力が D2H で写す配列 (43・58・57 個; 写した値は S0 の device と同一) と界面の診断 `diagVar` だけが B・C で変わる (A = A2、B = C)。
  `poison` (それらを NaN で埋めて 137〜1,237 step 続ける) は NaN が出ず、要求なしの反復との差はノイズ床の 2 倍以内
  (Euler は基準 3 本 対 1 本で最大 1.90 倍、NS は 10 本 対 10 本で最大 1.25 倍、CHT は 10 対 10 で丸めでゼロの slip 壁の `roUy` を
  除いて PASS)。
- host の写しを数値の経路が読む箇所 (コードの点検): 遠方境界 `farfield` の流束が自由流を host の `bc.bvar` の先頭から毎回読む
  (`convectiveFlux_d.cu` の `farfieldFlux_d_wrapper`)。device の farfield の `bvar_d` を書くカーネルが無いので出力の D2H は同じ値を
  写すだけだが、その前提が崩れると出力が数値を変える。CHT は使う前に D2H し直す (`updateConjugateWalls`)。

## mesh.wallDistExtraPhysIDs — 壁距離に含める非 wall 境界 (変換時)

```yaml
mesh: {wallDistExtraPhysIDs: [6]}
```

壁距離 (`wall_dist`, SST の F1/F2・ω 壁 BC が読む) は既定で `wall` / `wall_isothermal` の境界だけを壁点集合にする。
出口バッファのように **no-slip 壁の下流を slip 壁で延長**すると、接合 x で wall_dist が 0 → 数 mm に跳び、
BL をそのまま運んでいる slip 壁近傍で SST の ω が爆発する (case/16 3D run_0223: 接合直後 x=100 mm, y=壁, ω 1e21 → NaN)。
このキーで指定した physID (slip 等) を壁点集合に加えると wall_dist が連続になる。対称面 (slip) は含めないこと。既定は空。

## mesh.nodeInletCornerWall — node の入口∩壁コーナー (変換時)

node (median-dual) で入口 (`inlet_*`) と no-slip 壁が角ノードを共有するメッシュでは、角 CV に入口流束が入るのに
壁ピン ($\mathbf u=0$) で流し出せず $P$ が暴走する (1 次で延命・2 次で NaN)。`mesh: {nodeInletCornerWall: 1}` を
**`convertGmshToForge` 実行時の solverConfig に置く**と、入口エッジの壁ノード側半割面を壁 bcond に帰属させて回避する
([methods/discretization.md](../methods/discretization.md) §7.2 (D))。solver 実行時には読むだけで挙動は変わらない (メッシュに焼き込まれる)。
既定 0。入口壁角を slip で逃がす旧回避 (case/47 run_0006–0008) は不要になる。2D/3D とも対応 (3D は 2026-09-07)。

## mesh.axisSegmentRWeight — 軸対称 r 重みを区間ごとに掛ける (2026-10-10)

既定 1: 軸対称 (`axisymMethod: 0`) の面ベクトルを、変換器が書く区間ごとの $W_f=\sum_k r_k S_k$ (`/PLANES/rSurfVect`、**`isAxisymmetric: 1` の solverConfig で 2D の node を変換したときだけ書かれる**) にし、0 で旧来の $\bar r_f S_f$ (A/B・古い run の再現。古い run の再現は保存した旧 HDF5 + 0)。node・平面の 2D・`axisRFloor: 0` 以外とデータセットの無い格子では 1 でも旧来のまま (起動時に `[axisym] axisSegmentRWeight: ON/OFF` を 1 行出す)。plan [axisymmetric-freestream-hoop-gauge](../plans/active/axisymmetric-freestream-hoop-gauge.md) §4.5。

## discretization / bndFirstOrder — 離散化レイアウト (node-centered)

`mesh.discretization` (任意, 既定 `"cell"`)。`"node"` で node-centered (中点双対 median-dual) 化。
詳細は [`methods/discretization/`](../methods/discretization/)。

### `mesh.bndFirstOrder` — **使用禁止 (2026-08-12〜)**

**`bndFirstOrder` は新規・既存を問わず使わないこと。** 設定ファイルに書かない。既存 run の
config に残っている場合も、そこから複製するときに必ず落とす。削除は将来課題
([plans/active/architecture-bndfirstorder-removal.md](../plans/active/architecture-bndfirstorder-removal.md))。

禁止理由 (2026-08-12, case/26 の 2 次精度切り分けで判明):

1. **粘性応力を破壊する**。`zeroBndNodeGradient_d` は再構成勾配として `dUxdx…dUzdz` を 0 にするが、
   同じ配列を `viscousFlux_d.cu:114-118` が Newton 応力の面勾配として読む。「対流再構成のみを
   対象とする」というコメントの意図に反し、**境界隣接 CV の粘性せん断応力が消える**。
2. **適用範囲が意図より遥かに広い**。`bnode_flag` は *いずれかの* bcond の `iCells` を立てるため、
   疑似 2D (1 セル厚) メッシュでは spanwise の `side1`/`side2` が全ノードを覆い
   **flag が全域 1** になる (case/26: 21510/21510)。結果として「境界近傍だけ 1 次化」ではなく
   **全域 1 次化 + 全域の粘性応力喪失**になる。
3. そのため本フラグを付けた A/B は**切り分けとして成立しない** (2 つの副作用が交絡する)。
   過去に「bndFirstOrder で解決/改善しなかった」と記録された判定は再解釈が必要。

過去に本フラグで安定化したと記録されたケース (bump hiM node, 旧 node ノズルレシピ) は、
実際には上記 1.–2. の副作用による安定化であり、正しい対策ではない。node ノズルのレシピからは
すでに撤去済み ([plan §2.7](../plans/active/boundary-node-nozzle-wall-outlet-stability.md))。

```yaml
mesh:
  discretization: "node"   # cell | node
  # bndFirstOrder: 使用禁止 (書かないこと)
```

## node の離散化スキーム (固定・2026-08-16)

`discretization: node` では次が**常時 ON** で、config で切れない (旧キー `nodeAxisDirichlet` /
`nodeMidpointFx` / `nodeValueAtNode` / `nodeReconEdgeMidpoint` / `nodeAxisUrDirichlet` を書くと**起動時にエラー**):

- 値の位置 = ノード座標 (双対 CV 重心は軸対称の回転体積 $V=\bar r A$ にのみ使う)
- 対流 2 次再構成の目標点 = エッジ中点 (SU2 と同じ)
- 勾配 = LSQ (`gradLSQ: 2` 固定。GG は伸縮した壁行で市松を作るため node では禁止)
- 軸対称の軸ノード: 通常 DOF として解き、u_r=0 を状態ピン+残差 0+block-DPLUR 行の三点セットで課す

根拠と検証は [case/43.node_axis_dof](../case/43.node_axis_dof/README.md) と
[plans/active/architecture-node-option-consolidation.md](../plans/active/architecture-node-option-consolidation.md)。

## mesh.scalarGradient — node のスカラー勾配の作用素 (2026-09-27 から node の既定は lsq)

`lsq` (**node の既定、2026-09-27 から**) / `gg`。node で $k,\omega$・化学種 $Y_s$・受動種 $\xi$・凝縮モーメントの勾配を、
NS と同じ事前計算 LSQ 係数 (差分形、境界は内部隣接のみ、並進周期は合併係数を共有) で計算する。`gg` は Green–Gauss で、
2026-09-26 までの node の既定。cell は常に GG (キーは無視)。
起動エコー `'scalarGradient' effective: <値> (default|explicit)`。node で省略すると既定変更の警告が 1 行出る。

- **旧結果の再現**: `mesh.scalarGradient: gg` を明記する。
- **切り替え日をまたぐ run は途中から再開しない** (ユーザ決定 2026-09-27): 旧既定 gg で始めた段階起動 run を省略のまま
  新バイナリで再開すると、`stage_manifest.py` が gg 段と lsq 段を 1 区間につないで収束判定する (起動記録との結び付けは
  plan [`tooling-stage-manifest-launch-binding`](../plans/active/tooling-stage-manifest-launch-binding.md) で後回し)。最初から回し直すか、`gg` を明記して続ける。
- 設計 DB (SERN) は `FLAG_POLICY` 2026-09-27 以降、全作動点の実効値が lsq と確認できた評価だけを学習に使う。
- 検証: plan [`gradient-scalar-lsq-unification.md`](../plans/accepted/gradient-scalar-lsq-unification.md) (S0/S1・S2 [現行 gg 基準]・S3 合格)。

## 文書に無かった設定キー (2026-09-18, plan [config-key-pruning](../plans/accepted/config-key-pruning.md) §5.3 c/g/n)

棚卸しで「使われているのに `procedures/` `methods/` に一度も出てこない」ことが分かったキー。位置づけを 1 行ずつ残す。

| キー | 既定 | 何をするか | 実績 |
| --- | --- | --- | --- |
| `physProp.thermCondMethod` | 0 | 熱伝導率の与え方。0 = `thermCond` の定数、1 = 層流 Prandtl 数 `prandtlLam` から $\lambda=\mu c_p/Pr$ で作る | 797 run が 1 |
| `physProp.prandtlLam` | 0.72 | 上の 1 で使う層流 Prandtl 数。**壁関数の回復係数 $r=Pr^{1/3}$ (断熱壁温) でも読む**ので、`thermCondMethod: 0` でも効く | 797 run が書いているが**全て既定値** |
| `time.deltaT.precondEps` | 0.15 | 低マッハ前処理の停留点フロア $\varepsilon$ ($U_r=\min(c,\max(|u|,\varepsilon c))$)。小さいほど低マッハ振動を減衰するが $\varepsilon\lesssim0.1$ は発散 ($0.05$ で NaN) | 184 run 中 24 run が非既定 |
| `time.deltaT.monitorInterval` | 1 | 残差モニタ行を出す間隔 [step]。長時間 run でログを絞る用途 | 427 run |
| `time.deltaT.axisTimestepBeta` | 0.0 | 軸対称 near-axis 安定化。擬似時間のスペクトル半径に軸項 $\lambda_{axis}=\beta(|u_r|+c)A_{planar}$ を足し、$r\to0$ で $\Delta\tau\propto \mathrm{CFL}\, r/(|u_r|+c)$ を自然に与える | 39 run |
| `time.deltaT.ducrosLimiter` | 0 | KEEP の blend に Ducros センサを掛ける opt-in。**SLAU では効果が場に出ない** (SLAU 自身の散逸が支配) | 56 run 中 1 run が非既定 |
| `turbulence.sstCrossDiffJac` | 0 | SST $\omega$ の交差拡散を point-implicit 対角に入れる。$CD_\omega>0$ のとき対角へ $+CD_\omega/(\rho\omega)$。**収束解は不変**でサブ反復の収束だけ良くなる | 68 run |
| `bodyForceCtrl` / `bodyForceCtrlTarget` / `bodyForceCtrlRelax` | 0 / 0.0 / 1.0 | 質量流量一定制御。体積平均 $\langle\rho u_x\rangle_V$ を `bodyForceCtrlTarget` [kg/m²/s] に保つよう体積力を緩和で追従させる。既定 0 = `bodyForce` の固定値駆動 | 70 run |
