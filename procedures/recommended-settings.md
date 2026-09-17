# 推奨解析設定 (現行レシピの正本)

**このファイルが「今どう設定するべきか」の正本** (2026-09-08 起票)。各キーの意味は
[`solver-settings.md`](solver-settings.md)、起動手順は [`divergence-and-startup.md`](divergence-and-startup.md)、
メッシュ/実行手順は [`calculation-workflow.md`](calculation-workflow.md) を参照し、ここには**選ぶべき値と理由・根拠 run**だけを書く。

- 各節の見出しに **現行 (YYYY-MM-DD)** を付ける。設定を変えたら日付を更新し、旧設定は末尾の
  「§9 旧設定 (superseded)」に移す。**日付の無い記述・§9 にある記述を新規 config に使わない**。
- 記憶や過去の run の config をコピーする前に、必ずこのファイルと突き合わせる (古い run の config には
  廃止キー・旧既定が残っている)。
- 新しい config を書くときは `.claude/skills/forge-config` (Claude) の手順どおり、解析種別 → 該当節 → 段階起動 → 検証の順で組む。

## 0. 解析種別の選び方

| 解析対象 | 節 | 主な根拠 case |
|---|---|---|
| 平面/3D ノズル・ダクトの定常 NS (層流・SST) | §1 + §2 (+§3 TP) | case/16, case/26, case/40 |
| 非粘性 (Euler) 設計評価・MOC 比較 | §1 + §4 | case/23, 44, 45, 46 (design chain) |
| 軸対称 | §1 + §5 | case/23, 29, 41, 42 |
| 多成分 semi-perfect (TP)・凝縮・化学 | §3 | case/16 (TP+凝縮), 28, 44, 47 |
| 非定常 (LES/DES, 音響, 周期箱) | §6 | case/09, 18, 38, 39 |
| 変換メッシュ (Fluent/非直交) | §7 | case/24〜27 |

## 1. 共通の基本設定 — 現行 (2026-09-08)

```yaml
mesh: {meshFormat: "hdf5", discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: X.h5, valueFileName: X.h5}
solver: "SLAU"
space: {convMethod: 1, limiter: 2}                      # 本段。起動は convMethod 0 (§1.2)
time:
  unsteady: 0
  last: {control: 0, nStepOuter: N}
  deltaT: {control: 1, dt: 1.0e-5, cfl: C, cfl_pseudo: C, implicitRelax: 0.7, blockDPLUR: 1, dt_min: 1.0e-8, dt_max: 1.0, detectNaN: 1}
  timeIntegration: 11
  nStepInner: 4                                          # 2026-09-12: 5→4 (3D node SST の 12000 step 継続で 3 と 5 の残差経路が一致、余裕で 4)
output: {level: 1}                                      # 保存量 + 原始量 + h0 (2026-09-08〜)
```

- **離散化は node (median-dual) が生産**、cell は回帰対照 ([user-prefers-node-base])。node のメッシュは
  **node 用に変換した h5** (`discretization: node` を書いた config で `convertGmshToForge`) を使い、2D は
  **平面メッシュ** (押し出し 2 ノード spanwise は 2 次 MUSCL の散逸が消えて発散)。
- 対流は SLAU。`convMethod: 1, limiter: 2` が本段の標準、`limiter: 0` は使わない、**`mesh.bndFirstOrder` は禁止**。
- 定常は陰解法 `timeIntegration: 11` + `blockDPLUR: 1`。実効 CFL は `cfl_pseudo` (§solver-settings「CFL の定義」)。
  `cfl` は表示用なので同じ値を入れておく。**`nStepInner` は 4** (node NS の soft/mid 段は 10)。
  根拠 (2026-09-12, 3D node SST TP case/16 run_0410–0412): 本段の `nStepInner: 3` は 5 と残差経路が全列一致
  (到達 step も同じ) で 1 step が 11 % 速く、2 は発散 → 余裕を見て 4 をユーザ決定 (旧レシピ 5 は過剰反復)。
- **TP の温度反転は `physProp.thermoFloat: 1` が既定** (2026-09-12 ユーザ決定): float Newton + double 研磨 (収束まで最大 3 段, 通常 1 段) で
  基準バイナリと同じ (未収束プラトーの) 状態に留まり、12000 step 後の壁圧差 ≤1e-4・温度反転誤差 ≤4e-10·T、1 step −13 %。`thermoHrefTemp: 298.15` が前提で、datum 無し config では
  自動で 0 (double Newton) に落ちて警告が出る (plan performance-3d-node-sst-speedup)。
- **`cfl_pseudo` の目安**: node NS ノズル 4〜6 + `implicitRelax: 0.7` (上限は EOS 圧力床の洗浄で決まり、
  relax のみ有効: cfl 8 + relax 0.7 ≈ 3〜4 倍速 [implicit-cfl-ceiling-eos-floor])、Euler 設計評価 4、
  bump など易しい流れは 50 まで。TP 多成分の陰解法は 0.5〜2 から (§3)。
- `detectNaN: 1` は常時。EOS 床 `pMin/roMin/tMin` は既定 (1 Pa) のままで良いが、**無次元・低圧ケースは下げる**
  (Taylor-Green: 1e-6)。
- 出力は `output: {level: 1}` (既定)。勾配・リミッタ・診断が要る run だけ `level: 2` か `extraFields`。
  全温・全圧は `VALUE/h0` から `tools/total_quantities.py` で作る (AGENTS.md「出力と後処理の原則」)。

### 1.1 境界条件 — 現行 (2026-09-07)

- 入口: `inlet_Pressure` (Pt/Tt/組成/k/ω)。node NS では `mesh.nodeInletCornerWall: 1` を**変換時**に付ける
  (入口∩壁角の質量溜まりを根治)。入口境界層・全温分布・組成 (H2O) 分布は `ints: {inletProfile: 1}` +
  `inlet_profile_<physID>.csv` (`tools/gen_inlet_profile.py` で生成・照合, 手順は
  [`inlet-profile.md`](inlet-profile.md); 超音速入口 `inlet_uniformVelocity` は Tt を ρ/U/Ps に換算)。
- 出口: 亜音速は `outlet_statPress` に **逆流用 Pt/Tt を必ず添える**。超音速出口は `outflow`、または
  `outlet_statPress` の Ps を実出口圧に合わせる (node は壁列が常に亜音速なので Ps ≪ 実圧だと SST が出口列から
  unstart [node-supersonic-exit-outflow])。3D の出口角線で unstart が続く場合は出口バッファ (slip 延長, physID 別)
  + `mesh.wallDistExtraPhysIDs: [そのID]` + Ps 一致 (case/16 run_0226〜0228)。
- 壁: NS は `wall` (断熱 no-slip)、等温は `wall_isothermal` (キーは `Ts`)。Euler は `slip`。
  **SST の壁は変換時に no-slip `wall` でないと wall_dist=0 になり ω が step 0 で発散**。
- 対称面/疑似 2D 側面は `slip`。

### 1.2 段階起動 — 現行 (2026-09-07, case/16 `run_user_profile.py` の手順)

1. **soft**: `convMethod: 0, limiter: 0`, `cfl_pseudo 0.5`, `nStepInner 10` (node NS), 2000 step。
2. **mid**: 同上で `cfl_pseudo 1.0`, 2000 step。
3. **本段**: `convMethod: 1, limiter: 2`, `cfl_pseudo 4〜6` + `implicitRelax 0.7`, 8000〜12000 step。

段間の引き継ぎは**同一メッシュなら index コピー** (`continue_run.py` / `--ic-index-from`)。**3D で `interp_field.py` を
同一メッシュ引き継ぎに使わない** ((x,y) 最近傍で z 列・座標一致ノードを混同 [interp-field-coincident-nodes-trap])。
メッシュを変えたときは `interp_field.py` (cross-mesh、原始変数補間)。一様 IC から超音速/SST を直接始めない。
初期 k/ω は非ゼロ (例 k=1, ω=1000; ω=0 は初手発散)。

### 1.3 検証と報告 (必須)

投入前 `check_mesh_quality.py` (AR ≤ 1000, skew ≤ 0.9)。まとめ前 `check_convergence.py` (全列) と
`check_quasisteady.py` (報告する量) の VERDICT を貼る。NaN チェック、run パスの明示、case README の run 一覧更新。

## 2. SST (RANS) — 現行 (2026-09-08)

```yaml
turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1,
             wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInf: 1.0, omegaInf: 1000.0}
```

- 既定で ON (キー省略時): `sstNodeWallKPin: 1`, `sstOmegaProdFromPk: 1`, `sstSigmaBlend: 1` (2026-09-08〜)。
  旧挙動を再現するときだけ 0 を明記。`sstIsotropicStress` / `sstEnergyKSource` は 0 のまま (分離型、非推奨)。
- **`sstEnergyIncludesK: 0` が既定で確定** (2026-09-08 ユーザ決定: roe に k を含めない E_m 形が安全。
  乱流の跳ねが T に伝播しない)。エネルギー整合の検証で必要なときだけ 1 (opt-in) にし、全温は h0 (k 込み属性) から。
- `dilatationCorrection: 2` はノズル BL を SU2 比 −16 % 薄くする (モデル形式差、ソルバ無罪)。SU2 と揃える A/B は 0。
- `katoLaunder: 1` はノズル喉のよどみ偽生産抑制用。平板など剪断層主体では 0 でも可 (A/B で無影響)。
- 壁処理: y⁺≈1 のメッシュは `wallTreatmentSST: 0` (低 Re)。高 Re で y⁺≈1 と AR ≤ 1000 が両立しないときは
  y⁺ 30〜80 + `wallTreatmentSST: 1`。node 壁関数は Cf −6 % の既知欠損あり ([node-wallfunction-pk-convention-deficit])、
  3D の角線ノードでは代表点なし→u_τ=0 になる (未対応)。断熱壁の壁温出力は `sstThermalWallFunction: 3` (defect-flux)。
- 3D node SST の角部加熱は k/ω 拡散の相対ゼロ割ガード (2026-09-08 修正) で解消済み。**3D SST の壁圧は実験より +9 %
  (側壁合流域の乱流 BL 過厚) が未解決** (case/16 run_0228)。定量比較には 3D 層流 (0.7 %) か 2D SST (+1.5 %) を使う。
- restart で `vis_turb` が再現されない (敏感な擬似衝撃波は位置が動く) [forge-sst-restart-nonfidelity]。

## 3. 多成分 semi-perfect (TP)・凝縮・化学 — 現行 (2026-09-07)

```yaml
physProp: {thermalMethod: 2, species: [MIXDRY, H2O], speciesDBFile: species_db.yaml, thermoHrefTemp: 298.15, ...}
```

- **`thermoHrefTemp: 298.15` は必須** (絶対基準 h では χ_eos が桁違いになり TP×node 軸対称が発散 [isobutane-wt-semiperfect])。
- 種 DB は `forge_design.gas.semiperfect.mixture_pseudo_species_split` (乾き空気を擬似種 MIXDRY にまとめ H2O を残す)。
- TP 陰解法の `cfl_pseudo` は 0.5〜2 から上げる (H2O 生成エンタルピーの増幅で上限が低い)。**`implicitRelax: 0.7` を付ければ 6〜8 まで可**
  (2026-09-16 case/44 va3 M4.19 node Euler 軸対称 TP 2 種 + 非平衡凝縮 `run_0181`–`0189`: cfl 6/8 + relax 0.7 は乾き一様場からの起動でも安定で場・残差床が cfl 2 と同じ;
  relax 無しの cfl 6 は軸列の EOS 床洗浄で発散 [implicit-cfl-ceiling-eos-floor と同じ機構]、relax 無しの cfl 4 は完走するが残差床が 2〜10 倍高い)。
  定常 precond × 多成分は `speciesPrecondDt: 1` (既定)。TP 亜音速 `outlet_statPress` の γ 混用は修正済。
- 凝縮: 平衡凝縮を選ぶなら `condensation: 1, condEquilibrium: 2` (EOS 拘束形、厳密 S=1) を推奨 (設定既定値は 0 = 非平衡)。蒸発は既定 ON。
- 凝縮 (2026-09-15, plan [condensation-source-limiter-steady](../plans/accepted/condensation-source-limiter-steady.md)): 非平衡の θ 律速は
  **`condLimiterMode: 1` (既定) で更新クランプ**になり、ソース残差の明示的な Δτ 依存 (θ×Δτ_loc) を除去した (旧 0 は残差に θ を掛け、大型ノズルで成長を 1/4 に絞っていた; case/44 で cfl 2 の解が旧 cfl 0.5 と一致)。cfl 間の固定点一致の検証状況は plan condensation-source-limiter-steady §9。
  凝縮 ON の定常 run は `output.level 2` の `condLim_<s>` が収束時に全域 ≈1、`condClampCorr_<s>` が 0 であることを確認する。
  RK 陽解法では自動で 0 に降格 (起動ログ `[condensation] condLimiterMode=`)。dual-time + S3 (`speciesFaceReconstruction 2`) では既定 `passiveFct 1` の保存的 FCT 補正が受動種の有界性を担う (`check_passive_budget.py` で収支 PASS を確認)。dual-time は既定 `passiveScalarScheme 1` ならモーメントに BDF 物理時間項が付き
  更新クランプ (mode 1) のまま有効 (2026-09-17, plan species-passive-scalar-unification); `passiveScalarScheme 0` のときだけ 0 に降格。上限は `condDgMaxStep` 5e-3 / `condDTmaxStep` 1 K。
  凝縮 run は h0 保存を確認する (面温度修正済み)。onset は実験より ~5 mm 下流 (case/16 2026-09-08 比較)。
- 受動スカラ (2026-09-17, plan [species-passive-scalar-unification](../plans/accepted/species-passive-scalar-unification.md)): トレーサ・凝縮モーメントは既定で化学種経路
  (`passiveScalarScheme 1`)。2 次面移流 (S3) を使うなら `speciesFaceReconstruction 2` + `speciesImplicitCoupling 1` の組で (coupling 0 + S3 は定常で発散)。
  **`implicitRelax 0.7` は定常 (擬似時間) の安定化として推奨**で、**dual-time の非定常計算には使わない** (§6 の dual-time 節): sub-iter の残差ノルムは下がるのに遅いモードが
  収束せず、同じ物理時刻の解が sub-iter 数に依存する (case/44 `run_0399`–`0404`: nSub 40→80 の ρ 差が緩和なしの 1.3e-5 に対し 6.3e-4; 2026-09-17,
  plan [species-passive-scalar-unification](../plans/accepted/species-passive-scalar-unification.md) §5.1 #26)。dual-time で安定化が要るときは `cfl_pseudo` を下げるか nSub を増やす。
  S3 は凝縮 onset を 0.2 r_t 程度下流に動かす (数値拡散減) ので、実験比較の基準を変えるときは明記する。
- 凝縮 (2026-09-10, plan [condensation-kantrowitz-gamma-twophase-sonic](../plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)):
  Kantrowitz 補正 (`condKantrowitz: 1`) の γ は凝縮種 (蒸気) の γ_v が既定 (`condKantrowitzGammaMode: 0`; 旧=1)。
  凝縮セルの音速は `condSonicModel` (未指定=自動): **TP carrier H2O・`condEquilibrium 0`・境界が
  `inlet_Pressure`/`outflow`/`wall`/`slip`/`periodic` のみ** のとき二相 frozen 音速 (1)、それ以外は旧の全蒸気音速 (0)。
  起動ログ `[condensation] condSonicModel=... (理由)` で解決値を確認する。`outlet_statPress`/`wall_isothermal` を含む構成と
  `condEquilibrium 1/2`・pure N2・CPG は未検証 (明示 1 は警告つきで可)。
- 凝縮 (2026-09-13, branch feature/condensation-air; plan [condensation-kantrowitz-carrier](../plans/accepted/condensation-kantrowitz-carrier.md) /
  [condensation-air](../plans/accepted/condensation-air.md)):
  - `condKantrowitz`: **省略時既定 0** (非等温補正なし; `solverConfig.cpp`) / 1=Kantrowitz (純蒸気形) / **2=Feder carrier 形** (キャリアの冷却込み, $\hat q=b-\tfrac12$) / 3=同+表面項 ($\hat q=b-\tfrac12-\ln S$)。
    Wysłouzil 参照設定 (case/16 run_0335 系) は 1 を明示する。2/3 は θ 167→3–4 で計算上の onset が mode 1 より 8.2 mm 上流に動く (**モデル間差**; 実験との位置差ではない) ため推奨は 1 のまま。診断出力は `condTheta_0`/`condLim_0` (`extraFields`, 種番号付き)。
  - `condFloat` (既定 1, 2026-09-13): 凝縮経路の float 実体 (凝縮 ON の速度; [methods/condensation.md](../methods/condensation.md) 実装 §9)。A/B や旧結果の再現には 0。
  - `condSigmaScale` (既定 1.0): 表面張力の定数倍感度。Wysłouzil mode 3 基準で σ×1.03 は onset +2.5 mm (収束 PASS+STEADY)、σ×0.97 は −2.2 mm だが**未収束の準定常参考値** (残差 plateau, h0 偏差 OSCILLATING; case/16 run_0354)。Tolman 補正の代替ではない。
  - 空気凝縮 (CPG carrier 形): `physProp: {cp: 1008.7, gamma: 1.4}` (二成分空気 R 288.19) + `condensation: {condModel: 0, nCondSpecies: 1, condVaporMassFraction: 0.7671}`。
    受付は SLAU × `thermalMethod 0` × 非平衡 × `condKantrowitz ≤1` × 境界 `inlet_uniformVelocity`/`outlet_statPress`/`slip` のみ (config 検査で拒否)。
  - N2 低温物性: `condN2LatentLowT: 1` / `condN2PsatLowT: 1` (既定, 70 K 未満の潜熱線形外挿と整合する飽和圧), `condN2LiquidCp: 2000` (J/kg/K; 1500/2000/2500 で onset 38.95/39.67/40.25 K = +500 で +0.58 K, −500 で −0.72 K の同符号)。
    0/0 は旧一式 (A/B 用)。ログの `[condensation] WARNING: CPG two-phase temperature inversion failed` が出た run は使わない (保存量を凍結して継続する退避)。
- 化学 (H₂): 定常陰解法は `speciesImplicitCoupling: 2` + `jacobianMode: 2` 必須、Cabra は dual-time 必須
  ([chemistry-finite-rate-direction], branch feature/chemistry-finite-rate)。

## 4. Euler 設計評価 (design chain) — 現行 (2026-09-04)

`design/forge_design/evaluate/runner.py` が生成する config が正本: 全域 2 次 (`convMethod 1, limiter 2`) + 陰解法
`timeIntegration 11, blockDPLUR 1, cfl_pseudo 4`, `nStepInner 5`, `lowMachPrecond 0`, 壁 `slip`。
δ* (排除厚) は CONTUR 初期壁 + 固定 Euler 基準・帯局所抽出 1 pass ([deltastar-production-recipe])、質量流量ゲート必須。

## 5. 軸対称 — 現行 (2026-08)

- `isAxisymmetric: 1`, `axisymMethod: 0` (r 重み方式) が既定。SU2 流 `axisymMethod: 1` は実装済みだが implicit 深収束未達。
- `space.pRef` に動作静圧 (hoop 源とゲージ整合、自由流誤差 1.6e7 倍改善)。`hoopAreaFromClosure` 利用可。
- node の軸ノードは通常 DOF + u_r=0 ピン (固定、キー不要。旧 `nodeAxisDirichlet` 等は起動時エラー)。
- 近軸の block-DPLUR 粘性対角は幾何修正済み (float で解ける)。

## 6. 非定常・LES/DES・周期箱 — 現行 (2026-09-17)

- 過渡 (音響) を含む閉じた系は **`unsteady: 1` 必須** (定常局所 dt は数 step で発散 [steady-localdt-acoustic-transient-instability])。
  陽解法は `timeIntegration: 4` (RK4) か 3 (RK3) を `unsteady: 1` で。dual-time は `timeIntegration 11 + dualTime 1`、
  物理 CFL ≲ 12 ([dualtime-subiter-divergence-fingerprint])。
- **dual-time の内部反復 — 現行 (2026-09-17)**: `cfl_pseudo: 12`〜`20`、`nSubIterDualTime: 10`〜`20`、**`implicitRelax` は書かない** (既定 1.0 = 緩和なし)。
  - 擬似時間の CFL に安定限界が見当たらない (物理時間項が対角に $V/(a\Delta t)$ を足すため)。case/44 で `cfl_pseudo` 2〜**160** を緩和なしで完走、NaN 0。
  - 収束解は `cfl_pseudo` に依存しない (nSub 80 同士で ρ 差 3.1e-4 = 同一 cfl の nSub 40 vs 80 差 2.3e-4 と同水準)。
    一方**必要な nSub は `cfl_pseudo` で決まる**: nSub 80 との ρ 差は `cfl_pseudo` 2 で nSub 10 → 1.2e-1 / 20 → 4.5e-3 / 40 → 2.3e-4 なのに対し、
    **`cfl_pseudo` 12 以上なら nSub 10 で既に 2.2e-4** (float のノイズ床)。壁時間は nSub 10 が 0.045 s/step、40 が 0.14、80 が 0.265 →
    **同じ解が 1.8〜3 倍速い**。根拠は case/44 `run_0459`–`0480` (生産 float, 湿り核生成場, convMethod 1 + S3 + FCT)。
  - **`implicitRelax` (定常の推奨 0.7) を dual-time に持ち込まない**: 安定性は物理時間項が担うので緩和の効果は無く、遅いモードの収束だけ遅らせる。
    しかも残差ノルムには出ない (2.8〜3.5 桁下がって見える) ので、同じ物理時刻の解が nSub 依存になる
    (case/44 `run_0399`–`0404`: nSub 40→80 の ρ 差が緩和なし 1.3e-5 に対し 0.7 で 6.3e-4)。
  - **sub-iteration が足りているかは残差の桁数でなく `nSub` を倍にして解が動かないかで見る** (float では残差床がノイズ床に当たるため)。
    安定化が要るときは緩和ではなく `cfl_pseudo` を下げるか nSub を増やす。
  - 確認は `python3 solver_density_cuda/tools/check_solver_config.py <run_dir>` (投入前) と、上の nSub 倍増比較 (精度が要る run)。
  - 限定: 検証は case/44 の 1 形状・物理 CFL 2・dt 8e-6 の範囲。物理刻みを大きく取ると対角が痩せるので擬似 CFL の余裕も減るはず。
- 解像 LES/DNS は `solver: "KEEP"` + `keepDissType/Coeff` (**トップレベル**キー; `space` 配下は無視される) で
  σ=0.05 (市松抑制) / 0.02 (解像 LES)、`WALE` は off + ES 散逸 ([keep-es-dissipation-status], [wale-inactive-fix])。
  一様流 U∞≠0 は `space.roRef/uRef` (KEEP, CPG) で機械精度保存。fdblend は音響モード成長のため opblend (増分フル c) を使う。
- 低マッハ (backstep 等) の市松は `lowMachPrecond: 2` (1 は不可)。
- 周期境界は node/cell とも実装・検証済 (DOF 同一視)。体積ソースは `volumePartial` (seam 2 倍バグ修正済)。

## 7. 変換メッシュ (Fluent / 非直交) — 現行 (2026-07)

`calculation-workflow.md`「変換メッシュを安定に回す実行レシピ」を厳守: 1 次で投入、初期場を入口流速に合わせる、
陰解法 `cfl_pseudo`≈1 `nStepInner`≈10、非直交は `pRef`、出口は `outlet_statPress` + 逆流 Pt/Tt。

## 8. メッシュ — 現行 (2026-09)

- AR ≤ 1000, skew ≤ 0.9 (`check_mesh_quality.py`)、y⁺≈1 が要るときは第一層と接線長のバランスで AR を守る。**壁法線の構造層に限り AR ≤ 5000 まで緩和可** (2026-09-12; `--ar-max 5000` / 問題 YAML `mesh.ar_max`, 台帳に明記)。冷却壁では y⁺ が ×5〜6 に上がるので `mesh.wall_first_frac_throat` (スロートだけ第一セルを詰める) と併用する。
- node 2D は平面メッシュ、3D は六面体押し出し可。3D の角線ノード (2 壁交線) は内部隣接ゼロだが受動で問題なし
  (2026-09-07 検証)。y₁ 0.6 µm × z₁ 2 µm の極端な異方角セルは旧バイナリで角部加熱を起こした (相対ガードで解消)。
- SST メッシュは壁 bcond を no-slip `wall` で変換 (wall_dist)。slip 延長壁は `wallDistExtraPhysIDs`。

## 9. 旧設定 (superseded) — 新規 config に使わない

| 旧設定 | 状態 | 代替 |
|---|---|---|
| `mesh.bndFirstOrder: 1` | **禁止** (粘性応力破壊・疑似 2D で全域に効く) | 段階起動 (§1.2) |
| `nodeAxisDirichlet` / `nodeMidpointFx` / `nodeValueAtNode` / `nodeReconEdgeMidpoint` / `nodeAxisUrDirichlet` | 廃止 (2026-08-16, 書くと起動エラー) | node は固定スキーム (§1) |
| `turbulence: {LESorRANS: 2, RANSmodel: 1}` 旧キー体系 | 旧 config に残存 | `turbulence: {model: "sst", ...}` |
| `sstOmegaProdFromPk: 0` / `sstSigmaBlend: 0` | 2026-09-08 まで既定 | 既定 1 (旧挙動が要るときだけ 0 明記) |
| `sstEnergyKSource: 1` / `sstIsotropicStress: 1` (分離型) | 非推奨 (境界未完備・離散保存せず) | 必要なら `sstEnergyIncludesK: 1` (opt-in) |
| 定常 + 陽解法 (`timeIntegration: 3`, `unsteady: 0`) | 非推奨 (局所 dt で不安定) | 陰解法 11 + blockDPLUR、または `unsteady: 1` |
| `thermoHrefTemp` 未指定 (絶対 h 基準) の多成分 TP | 発散要因 | `thermoHrefTemp: 298.15` |
| 同一 3D メッシュの段間引き継ぎに `interp_field.py` | 禁止 (z 列混同) | index コピー |
| 超音速 node 出口に `outlet_statPress` (Ps ≪ 実圧) | unstart の原因 | `outflow` か Ps 一致 |
| `output` 未指定 = 全量出力 | 2026-09-08 まで | `output.level` 既定 1 (全量は `level: 2`) |
| k/ω 拡散の絶対ゼロ割ガード 1e-12 [m³] | バグ (2026-09-08 修正) | 相対ガード (コード側、キー無し) |
| `wall_dist` を双対重心から測る変換 | バグ (2026-09-08 修正) | ノード座標 (コード側) |

## 変更ログ

- `2026-09-17` — §6 に dual-time の内部反復レシピを追加 (`cfl_pseudo` 12–20 + `nSubIterDualTime` 10–20 + 緩和なし; 擬似 CFL に安定限界が見つからず、必要な nSub は `cfl_pseudo` で決まる)。定常の `implicitRelax 0.7` は据え置き。投入前チェック `check_solver_config.py` を追加。

- `2026-09-08` — 初稿 (散在していた推奨値を集約。ユーザ要請「既定の解析設定を 1 か所に、最新/旧を明記」)。
