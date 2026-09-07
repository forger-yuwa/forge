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
- `cfg.cfl` は表示用の `cfg.dt` (および `max cfl` 表示) をスケールするだけで、
  定常計算の積分自体には効かない。したがって定常で `cfl` だけを上げても収束は速くならない。

運用ルール:

- **定常 (`unsteady: 0`) の収束を速めたいときは `cfl_pseudo` を上げる** (例: 陰解法 `blockDPLUR: 1`
  なら `cfl_pseudo` を 20〜50 程度まで)。`cfl` を上げても無意味。
- 非定常 (`unsteady: 1`, `dualTime: 0`) の物理時間刻みは `cfl` (または `dt`) で決まる。
- ログの `max cfl` 表示は `cfg.cfl` に追従する値であり、実効積分 CFL ではない点に注意する。

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
- **`lineDtDirectional: 1`**: 方向別 dt — 内部 line 面の λ (音響込み) を擬似 dt の CFL max から除外。
  壁ノードの境界半割面は除外されない (壁 CV 自身の Δτ は境界面律速のまま) 点に注意。
- `lineViscCoupling` / `lineViscousDtRelief`: line 面のスカラー粘性結合と粘性 CFL 割引。
  圧縮性の壁法線 pseudo-dt 律速は音響 (λ_visc/λ_ac=2ν/(Δn·c)≪1) なので通常は効果僅差 — 既定 off で可。

**使い分け (2026-09-03 実測)**: 定常 M6 ノズル (case/45) では cfl 上限を上げない (streamwise
対流律速) — 使わない。**壁解像 DDES/dual-time (case/39) では point の cfl_pseudo 上限 (2 で発散)
を 8 まで引き上げ、`cfl_pseudo 4 + nSub 13-15` で現行比 ~13% 高速・ωバースト低減** — 採用候補
(長時間 run でのバースト余裕検証は未了)。

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
| `2` | S2+S3: species 移流も同一 face 組成 | **使用しない** (cfl≤2 限定の experimental。高 CFL で発散) |

```yaml
time:
  deltaT:
    speciesFaceReconstruction: 1   # 多成分 TP の face thermo 組成を ρ と同次数で整合
```

効果の実測 (case/28 He/空気 coaxial): 圧力振動振幅 cfl2 で −77% / cfl4 で −48%。単成分・CPG では
無効果 (組成を thermo に使わないため)。詳細は [`../methods/convection/theory.md`](../methods/convection/theory.md)
の「多成分 TP の face 組成整合」節。

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

- **`thermoHrefTemp: 298.15` を必ず指定する** (反応熱は sensible datum の残差項 $\dot Q=-\sum_s h^{abs}_s(T_{ref})\dot\omega_s$ として入る。絶対 datum (0) でも動くが陰解法は不安定)。
- 機構に現れる種は `species` に全て含めること (無ければ起動時エラー)。`species` にだけある種は不活性として扱う。
- 熱力学 DB は `tools/cea_thermo_to_species_db.py thermo.inp --species ...` で CEA から生成する (ラジカルは内蔵 DB に無い)。
- 出力: `chemQdot` [W/m³], `chemTau` [s] (=1/max|∂ω_s/∂ρY_s|、化学時間の目安。`dt` や `cfl_pseudo` の妥当性判断に使う)。
- 検証: `case/35.uniform_periodic_box/run_0049_node_h2_ignition` (0-D 着火 vs Cantera)。

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
