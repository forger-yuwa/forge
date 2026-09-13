# 凝縮 (非平衡・二相) 経路の float 化による凝縮 ON 計算の高速化

## メタ

- **area**: `architecture / performance` (condensation)
- **status**: `in_progress`  <!-- codex plan GO-with-changes 全件採用 → 実装 → codex result 1 回目 NO-GO (M6/m2) 全件採用 (#8–#14) -->
- **related_docs**:
  - [`methods/condensation.md`](../../methods/condensation.md) (実装 §9 混合精度実装)
  - [`methods/architecture/performance.md`](../../methods/architecture/performance.md) (数値精度方針・判定基準)
- **related_plans**:
  - [`performance-3d-node-sst-speedup.md`](../accepted/performance-3d-node-sst-speedup.md) (dry 経路の float 化。判定基準 §4.3・回帰ツール `perf_regress.py` を継承)
  - [`condensation-air.md`](../accepted/condensation-air.md) / [`condensation-kantrowitz-carrier.md`](../accepted/condensation-kantrowitz-carrier.md) (CPG carrier 二相 EOS・Feder carrier 形; 本 plan はモデル式を変えない)
  - [`condensation-evaporation.md`](../accepted/condensation-evaporation.md) / [`condensation-equilibrium-eos.md`](../accepted/condensation-equilibrium-eos.md)
- **created**: `2026-09-13`
- **owner**: `Claude (branch feature/perf-3d-speedup 継続, 計算は AWS A10G)`

## 1. 目的

dry の 3D node SST TP 生産計算は 33.8 ms/step (A10G, run_0400) まで下がったが、**凝縮 ON にすると 89.5 ms/step (2.6 倍)** になる
(§4.1)。凝縮の相変化ソース・二相 EOS・面エンタルピー補正が全て double で書かれており (CC 8.6 の FP64 は FP32 の 1/64)、
しかも凝縮 ON では dry 側で導入したハイブリッド温度反転が全セルで無効化される。凝縮経路を float 化し (モデル式・評価点は不変)、
**凝縮 ON の追加コストを dry +10 ms 以内 (≤ 45 ms/step)** にする。ユーザ指示 (2026-09-13)「凝縮ありの double 計算は無くして速く」。

## 2. スコープ

- **やる**:
  - 物性 ($\ln p_{sat}, L, \sigma, \rho_l, k_{gas}, \mu_{gas}$) の**区分 3 次表** (double で構築、float で評価) と、それを使う `condensation_source_d`
    (核生成 CNT [対数空間・対数上限]・Kantrowitz/Feder θ・成長 HK/Goodheart/Gyarmathy・蒸発・数値ヤコビアン・θ 律速) の float 化、dry セルの早期退出。
  - SLAU 二相面エンタルピー補正の潜熱 $L(T_{cell})$ / CPG 二相 $L(T_f)$ を float 表で評価、$g_L=g_R=0$ 面のスキップ (v1 の `condL` セル配列案は codex M4 で不採用)。
  - `dependentVariables_d`: 凝縮 ON でもハイブリッド温度反転を g≈0 セルで使う。g>0 の一温度二相反転 (TP carrier / pure TP) を float Newton + double 研磨に。
  - `cond_realizability_clamp_d` の float 化、モーメント移流 4 起動の `scalar_advection_multi_d<4>` への融合。
  - `condensation.condFloat` (既定 1)。0 で従来 double 経路 (ビット不変) を残す。
  - 単体試験 `tests/unit/test_cond_float.cpp` (float/double の物性・J・r*・dr/dt・ソースベクトル・温度反転の一致)。
- **やらない**:
  - モデル式・評価点・律速規約の変更 (Kantrowitz/Feder 形、成長則、蒸発規約、θ 律速はそのまま)。
  - 平衡形 (`condEquilibrium` 1/2) の float 化 (double のまま; 風洞ノズルの上限評価用で律速でない)。
  - 二温度 (`condTwoTemp`, 既定 0) の Newton の float 化 (double のまま)。
  - CPG 二相の括弧付き Newton `cond_T_from_e_cpg` (case/34 N2 / 空気; 全セルが二相になり得る 2D 小規模ケースで律速でない) — double のまま。必要なら後続 plan。
  - モーメント保存量の無次元化・格納精度の変更 (現状 float 格納、[`methods/condensation.md`](../../methods/condensation.md) 実装 §8)。

## 3. 関連 docs と前提

- 前 plan の §4.3 判定基準 (絶対基準 A/B/C **または** 基準×基準ノイズの 2 倍 **または** double 高精度参照との距離) と `tools/perf_regress.py` をそのまま使う。
- FP64 律速の同定方法 (nsys/ncu, `cuobjdump -sass`) は [`methods/architecture/performance.md`](../../methods/architecture/performance.md) §2。
- 計算は AWS A10G (`~/forge-perf`, `~/forge-bin/forge_merged` = マージ後 3b77cc4b, sha256 a2b931c4…) で行う (ユーザ指示: ローカル PC はメモリ不足で落ちる)。

## 4. 設計方針

### 4.1 計測 (2026-09-13, A10G, 100 step ×2 交互, `bench_steps.sh`)

- テンプレート `case/16.nozzle_wys/run_0420_perf_cond_baseline` = run_0400 (3D node SST TP 2.37 M 節点, dry) + `condensation: {condensation: 1, nCondSpecies: 1, condModel: 1, condGasSpecies: 1, condKantrowitz: 1, condGrowthModel: 0}` (H2O in MIXDRY, TP carrier, Kantrowitz 原形, Hertz–Knudsen)。IC は dry (モーメント 0 から)。
- **dry 33.83 / 33.81 ms/step, 凝縮 ON 89.49 / 89.46 ms/step** (+55.7 ms)。
- nsys (30 step, kernel 合計 85.6 ms/step) の上位:

| カーネル | ms/step | dry のとき | 内容 |
| --- | --- | --- | --- |
| `condensation_source_d` | **31.7** | — | 全セルで CNT 核生成 (psat/σ/ρ_l/ln/exp)・成長・θ・src_jac 用の摂動評価 2 回 (T+0.1 K, Q1+1 %) を **double** で。dry セル (S≤1, g=0, Q0=0) も同じ経路 |
| `SLAU_d` | **11.5** | 2.6 | TP carrier 二相補正 `h −= g·L(T_cell)` の `cond_latent` (NASA-9 気相/液相の 2 評価, double) を **全面・両側**で評価 (g=0 の面でも) |
| `dependentVariables_d` | **10.2** | 3.2 | `useHybrid = (thermoFloat && condensation == 0)` により凝縮 ON では**全セル**が旧 double 反転。加えて g>0 セルの二相 Newton (double, 潜熱の数値微分込み) |
| block-DPLUR ×5 | 10.1 | 9.3 | 変化なし |
| `cond_realizability_clamp_d` ×2 | 1.6 | — | 蒸発塵判定に psat/ρ_l/cbrt を double で |
| モーメント移流 `scalar_advection_first_order_d` ×4 + `cond_primitive_d` ×8 + RK ×6 | 3.1 | — | 起動数が多い (k/ω は multi 融合済) |

凝縮固有の超過 55.7 ms のうち **~49 ms が上の 5 項目** (FP64 と無駄な評価)。

- **発達した凝縮場での速度と非決定性 (2026-09-13, codex plan m8 採用)**: run_0420 config を dry IC から 3000 step 回した `run_0422_cond3d_spinup`
  (`check_convergence` NOT CONVERGED plateau: rms_roe 5.5e-4 で flat、rms_rog 1.9e-10 flat; `check_quasisteady` pmax/machmax ALL STEADY) の res_3000 を
  IC にした `run_0421_perf_regress_cond3d` (300 step) は **108.3 / 105.8 ms/step** (dry 起動区間の 89.5 より遅い: 湿潤セル 533k / 2.37M = 22.5 %,
  x 0.011–0.115 m)。**同一バイナリ 2 本の差は T 1.8e-2・ρg 0.115・h0 0.14 (正規化最大差)** で、差は x 0.06–0.10 m の湿潤帯 (壁距離 ~1e-3 m,
  T≈212 K, g≈1e-2) に集中する。dry の run_0400 (≤6e-5) と違い、3D 凝縮場は atomicAdd の擾乱が湿潤帯で O(1e-2) まで増幅される
  (プラトーの非定常性)。したがって **3D 凝縮 run は速度と分布統計の確認に使い、場の一致判定は 2D (case/16 node/cell H2O) と case/34 (N2/空気)
  の静かな case で行う** (§6)。

### 4.2 方針 (v2, 2026-09-13 codex plan レビュー M1–M7/m8–m9 を全件採用して改訂)

**codex M1 の実測** (NumPy float32 で現行式を再評価): H2O の潜熱 $L=h_v-h_l$ は液相 NASA-9 (係数 $1.3\times10^9/T^2$ 級) の相殺で
float だと相対 **9.2e-3**、N2 の Jacobsen 飽和圧は **2.9e-4**、係数だけ float に丸めても $L(300\,\mathrm{K})$ が 1.7e-4。式をそのまま float 化する案
(v1) は棄却し、**物性は double で作った区分 3 次テーブルを float で評価する**方式にする。

1. **物性テーブル** (`cuda_forge/condensationProperties_d.cuh`, 新規 `condensationTables_d.{cuh,cu}`)
   - 対象: $\ln p_{sat}(T)$, $L(T)$, $\sigma(T)$, $\rho_l(T)$, $k_{gas}(T)$ (成長則の気相熱伝導率)。それぞれ現行 double 関数 (クランプ・低温外挿込み)
     から**ホストで double で**作る。一様格子 $h$ (N2 0.1 K, H2O 0.25 K)、区分ごとに 4 係数 $(c_0..c_3)$ の 3 次 Hermite (両端の値と**片側**微分;
     接続点 N2 50 K [psat 切替]・70 K [潜熱線形外挿]・45 K [物性床]・$T_c-0.5$、H2O 273.15 K [液相線形外挿] は格子点に置き、区間をまたぐ
     漏れを無くす)。範囲: N2 [20, 125.7] K、H2O [120, 1200] K (case/44 燃焼室 1161 K の乾き高温セルも表で引く; $\ln p_{sat}$ を表にするので
     $S=\exp(\ln p_v-\ln p_{sat})$ は溢れない)。範囲外は端でクランプ (double 関数のクランプと同じ側に倒れる)。
   - 評価: `CondTablesF` (device グローバル配列へのポインタと $T_0, 1/h, n$) を kernel に値渡し。1 回の評価は 4 float 読み + FMA 3 回。
     H2O 表は 5 物性 × 4320 区間 × 4 係数 × 4 B ≈ 345 KB (L2 常駐)。
   - 誤差: Hermite 3 次の打切りは $h^4 f^{(4)}/384$ で 1e-9 以下、float 丸め (6e-8) が支配。**単体試験で 0.01 K 刻みの全域 (接続点・床・臨界直下を含む)
     で double と比較し相対 ≤ 2e-6 を要求** (§4.3)。表が満たせない物性は double に残す (その場合は本表に記録)。
   - $\ln p_{sat}$ の微分 ($T_{sat}$ Newton の $d\ln p/dT=L/RT^2$ の代替) と $dL/dT$ (二相反転の Newton 微分) は表の解析微分 ($c_1..c_3$) から取る。
2. **核生成 (CNT) は float 経路では対数空間で組み、指数化前に上限を掛ける** (`condensationSource_d.cuh`; codex M2)
   - $\ln J=\tfrac12\ln(2\sigma/\pi)-\tfrac32\ln m+2\ln\rho_v-\ln\rho_l-\Delta G^*/(k_BT)+\ln(\mathrm{corr})$; $m^3=2.7\times10^{-77}$ が float 範囲外。
   - **上限は対数で**: $\ln J\leftarrow\min(\ln J,\ln J_{max})$ ($J_{max}=10^{35}$, $\ln J_{max}=80.6<\ln(\mathrm{FLT\_MAX})=88.7$) を**本体・摂動評価の両方**に掛ける
     (codex の実例: N2 60.2 K, S=100 で $\ln J=93.3$ → float で Inf)。下限は $\ln J<-80$ で 0。これは float 経路の**意図的な設計変更**: 上限に張り付いた
     セルでは $\partial S/\partial T$ の摂動側も上限で評価されるため src_jac が 0 になる (double 経路は摂動側が無制限で、上限セルでは
     線形化が不整合だった)。double 経路 (`condFloat: 0`) は現行のまま。
   - 誤差伝播 (「最終 ULP」ではなく連鎖で評価): $\delta\ln J\approx 2(\Delta G^*/k_BT)\,\delta\ln S+\delta\ln K$ で、表の $\delta\ln p_{sat}\le2\times10^{-6}$、
     $\Delta G^*/k_BT\le200$ (それ以上は $J<10^{-80}$) なら $\delta\ln J\le10^{-3}$。**単体試験は連鎖全体 ($J$, $r^*$, $dr/dt$, $S$ ベクトル, src_jac) を
     double と比較**し、$|\Delta\ln J|\le 2(\Delta G^*/k_BT)\cdot2\times10^{-6}+10^{-4}$ を場所ごとに検査する (§4.3)。
   - src_jac (数値ヤコビアン) は float でも有限・非負であることと、double との相対差 ≤ 1e-3 (上限セルを除く; 上限セル数は報告) を単体で見る。
3. **dry セルの早期退出** ($S\le1$, $g=0$, $\rho Q_0=0$): 核生成 ($S\le1$ で $J=0$)・成長 ($Q_0=0$)・蒸発 ($g=0$)・摂動評価 (同条件で 0) が恒等的に 0
   (codex: 非平衡の通常物性域では妥当)。**先頭の sj_*/diag* 初期化は維持し、輸送残差 `res_*` には触れない** (現行と同じ)。$T_{sat}$ 診断は前 step の
   `condTsat` を初期推定にした float Newton (表の $\ln p_{sat}$ とその微分; 収束 |ΔT|<1e-3 K, 最大 25 回)。
4. **SLAU 面潜熱** (codex M4 でセル配列案を棄却): `condL` キャッシュは dependentVariables 後の境界処理 (`nodeWallDirichlet_d` が壁ノードの T を
   書き換える) で評価時点がずれるため**採用しない**。面カーネル内で現行どおり $L(T_{cell})$ を評価するが、float 経路は表 (4 読み + 3 FMA) を使い、
   **$g_L=g_R=0$ の面は評価しない** ($g=0$ で補正は恒等 0)。CPG 二相 `cond_face_h_cpg` の $L(T_f)$ も表で。`condFloat: 0` は現行 double 関数。
5. **温度反転** (`dependentVariables_d.cu`, `condensationEOS_d.cuh`; codex M3)
   - `useHybrid` の `condensation == 0` 条件を `condensation == 0 || condFloat != 0` に変える (組成 $Y$ の float 構築も同じ条件で切り替わる; M5。v2 初稿の `condFloat == 0` は誤記 [codex result m7])。
     $g_{liq}\le10^{-12}$ のセルは dry と同じハイブリッド反転。
   - $g>0$ の一温度二相反転 (`cond_T_from_e_carrier` / `_onetemp`): float Newton (`SpeciesThermoF` の cp/h + 表の $L, dL/dT$、最大 12 回) →
     double 研磨 (現行 double 関数を 1 反復ずつ、最大 3 段) → **成功条件 $|G(T)|=|e_{mix}(T)-e_{in}|\le10^{-9}|e_{in}|+0.05$ J/kg かつ有限** を
     二相エネルギー残差で判定。不成立なら現行の double Newton (30 回) を続行し同じ条件で再判定。**それでも失敗したセルは `roe` を上書きしない**
     (T/P/Ht は反転値で更新、保存量は保護; CPG 経路 `cond_T_from_e_cpg` の `ok` 規約と同じ) — 既存カウンタ `g_condTinvFail` に計上し wrapper が警告。
     現行の TP carrier 反転にも同じ成功判定と保護を入れる (double 経路の挙動変更: 失敗時のみ; 成功時はビット不変)。
   - `thermoFloat: 0` / datum 無し (`thermoHrefTemp` 未設定) では従来どおり全 double (ハイブリッドも二相 float も使わない)。
   - 単体試験に「float 格納 roe の 10 往復ドリフト」(thermo_float と同形式) を二相セル ($g$=1e-3〜$Y_w$) で追加。
6. **`condFloat` の分岐表** (codex M5):

   | 経路 | `condFloat: 1` (既定) | `condFloat: 0` |
   | --- | --- | --- |
   | `condensation_source_d` | float 実体 (表・対数 CNT・早期退出・上限の対数適用) | 現行 double 実体 (コード不変) |
   | `cond_realizability_clamp_d` | float 実体 (表) | 現行 double |
   | SLAU 面潜熱 (TP carrier / CPG 二相) | 表 + g=0 面スキップ | 現行 double 関数 (全面評価) |
   | 温度反転 (g≈0 セル) | ハイブリッド (thermoFloat 経路) | 現行 (凝縮 ON では全 double) |
   | 温度反転 (g>0 セル) | float Newton + double 研磨 + 成功判定 | 現行 double Newton + 成功判定 (新規; 失敗時のみ挙動差) |
   | 平衡形 (`condEquilibrium`≠0) / 二温度のソース kernel と EOS 反転、CPG 二相反転 `cond_T_from_e_cpg`、二相音速 | double (共通; ただし平衡形でも面潜熱・clamp・g≈0 セルの組成/反転は上の float 経路) | double |
   | モーメント移流 4 起動 → `scalar_advection_multi_d<4>` | 共通 (精度でなく起動数の変更) | 共通 |

   **ビット一致の要求は決定的な単体・凍結状態評価にだけ掛ける** (反復 run は atomicAdd で同一バイナリでも `ro` が 2.9 万要素不一致 [codex]):
   `condFloat: 0` の source kernel 出力 (res_rog/Q, src_jac, diag) と clamp 出力は同じ入力状態で変更前バイナリとビット一致、SLAU 残差と
   融合移流の残差は凍結状態 1 step で相対 ≤ 1e-6 (atomic 順序の丸め)。
7. **移流融合の保存性**: 凍結状態で融合前後の 4 モーメント残差の全セル和 (= 境界流束の和) が相対 1e-6 で一致すること、周期境界を含む
   case (case/09 周期 dual-time の凝縮 ON 版 20 step) で全セル和が丸め (1e-6) 以下であることを確認。
8. **FP64 命令の監査** (codex m9): float 実体の kernel について `cuobjdump -sass` で `DADD|DMUL|DFMA|DSETP|DMNMX|F2F.F64|F2F.F32.F64|MUFU.*64` を数え 0
   (呼び出し先を含む; 同名の double 実体は別シンボル)。

### 4.3 判定基準 (v2)

- **短期回帰 (場の一致)**: 前 plan §4.3 (`perf_regress.py cmp --noise … [--truth …]`; 凝縮量 `rog_*/roQ*_*` は B 1e-4; 診断 `condS/condTsat/condDrdt/
  condR30/condTheta/condLim` は除外 — `condTheta`・`condLim` を `SKIP_PREFIX` に追加 [codex m9])。**静かな case だけ**に適用: case/16 2D node H2O
  (`run_0456`, 基準同士 ρ 9.6e-7)、2D cell H2O (新設)、case/34 の cell/node × N2 pure/空気 CPG (`run_0100–0102` + pure N2 の新設 2 本、床は同一バイナリ
  反復 3 本)、case/44 `run_0201` (平衡, double のまま)。3D `run_0421` は基準同士が 1e-2 なので合否には使わず、分布統計 (湿潤帯の T/g の
  パーセンタイル・onset x・湿潤セル数) を並記する。
- **凝縮の律速規約**: `condLim` (θ 律速係数) の「<1 のセル数」と最小値を基準と比較し、核生成域で θ≈1 の割合が変わらないこと (別ゲート, 数を報告)。
- **物理検証 (短期回帰とは別; codex M6)**: 対象 run は基準・新版とも同 IC から**十分長く** (case/34 12000 step, case/16 2D 12000 step) 回し、
  (a) `check_convergence.py` の verdict 区分が同じ、(b) `check_quasisteady.py` (pmax/machmax) が両方 STEADY、(c) **onset**: case/34 は
  `onset_analysis.py --series` が両方 STEADY で onset x の差 ≤ 同一バイナリ反復 3 run の再現幅 (README の onset 表)、(d) **$h_0$ 保存**: `VALUE/h0`
  の軸 (case/16 2D は y=0 列、case/34 は中心線) について $\max_x|h_0(x)-h_{0,in}|$ [J/kg] を新設 `tools/cond_axis_h0.py` で各スナップショットで
  計算し、時系列が頭打ち (末尾 3 枚の変化 ≤ 10 J/kg) で、**新版の値 ≤ 基準の値 + 100 J/kg** ($h_{0,in}\approx3\times10^5$ J/kg の 0.03 %) [非劣化ゲート] かつ
  絶対水準 ≤ 0.1 % × 3e5 J/kg (twophase face-T fix の規約) [絶対ゲート]。**2 つは別判定として残す**: node 場は両方 PASS、cell 場 (case/34 cell 空気 429 J/kg, N2 cell 482 J/kg) は基準バイナリ自体が絶対ゲートを超えており (1.4–1.6e-3), 本 plan では非劣化のみ PASS・絶対ゲートは既存課題として condensation-followups へ転記 (codex result-2 M3)。
- **単体 (host)**: 表 vs double: $\ln p_{sat}$ 絶対 ≤ 2e-6、$L,\sigma,\rho_l,k_{gas}$ 相対 ≤ 2e-6 (0.01 K 刻み全域、接続点・床・臨界直下を含む)。
  連鎖: $|\Delta\ln J|\le 2(\Delta G^*/k_BT)\cdot2\times10^{-6}+10^{-4}$、$r^*$ 相対 ≤ 1e-5、$dr/dt$ 相対 ≤ 1e-4 ($|dr/dt|<10^{-6}\max$ では絶対)、
  ソースベクトル相対 ≤ 1e-3 (連鎖の上限)、蒸発 λ 絶対 ≤ 5e-4 (S=0.999 の (p_v−p_d) 相殺で float は 1.4e-4 まで出る; v2 初稿の 1e-5 を訂正)、src_jac 相対 ≤ 3e-3+tolJ (J 上限セル除く; 蒸発端 0.99<S<1 は |Δsj|·dt ≤ 1e-4)、二相反転 |ΔT| ≤ 1e-8·T と 10 往復ドリフト ≤ 従来 double 反転。
  範囲 (codex M7, 実測を覆う): N2 T 39–126 K (床 45 K 未満を含む; case/34 空気 run の $T_{min}$ 39.2 K)、H2O 150–400 K (+ 表は 1200 K まで)、
  S 1–1000 (run_0456 の $S_{max}$ 171)、g 0–$Y_w$ (case/34 $g_{max}$ 0.060)、$Q_0$ 0–1e19、r 1e-9–1e-6 m、$J_{max}$ 前後、消滅閾値 $r_{30}=2r_{min}$ 前後、
  Kantrowitz 0/1/2/3、成長 0/1、蒸発 on/off、Kelvin on/off。
- **単体 (device)**: `tests/unit/test_cond_float_device.cu` — 合成状態配列 (上の範囲を格子で) に float/double の両 kernel 実体を掛け、出力 (res, src_jac,
  diag) を host 側の許容で比較。`condFloat: 0` 実体は変更前と**ビット一致**。
- **速度**: A10G で (a) 発達場 `run_0421` IC (基準 108/106 ms/step) と (b) dry 起動区間 `run_0420` (89.5) を 100 step ×2 交互。目標 **(a) ≤ 50、(b) ≤ 45 ms/step**、
  dry `run_0400` は 33.8±0.2 で不変。ログ・nsys CSV・バイナリ/入力 sha256 を `case/16.nozzle_wys/_aws_perf_evidence/cond/` に回収し、湿潤セル比率・
  反転失敗数 (`g_condTinvFail`) を併記する。

## 5. 実装ステップ (codex 推奨順: ①物性・CNT → ②反転・互換 → ③device 試験 → ④定常性・保存性 → ⑤発達場の性能)

1. `cuda_forge/condensationTables_d.{cuh,cu}` (新規): 区分 3 次表の構築 (host, double) と device 評価 (`cond_tab_lnpsat_f` 等)、`CondTablesF`。
   `tests/unit/test_cond_float.cpp` の表精度試験を**先に**通す。
2. `cuda_forge/condensationSource_d.cuh`: `template<Real>` の θ/核生成 (float は対数 CNT + 対数上限)/蒸気状態/成長/ソースベクトル/蒸発;
   double 実体は現行式のまま。連鎖の単体試験。
3. `cuda_forge/condensationSource_d.cu`: kernel `template<Real>`、早期退出、$T_{sat}$ warm start、`condFloat` 起動切替。
4. `cuda_forge/condensationEOS_d.cuh` + `dependentVariables_d.cu`: 二相ハイブリッド + 成功判定 + `roe` 保護 + カウンタ; `useHybrid` 条件。
5. `convection/convectiveFlux_slau_d.inc.cuh` / `common_d.cuh` / `convectiveFlux_d.cu`: `CondArgs.tables`/`condFloat`、g=0 面スキップ、CPG 面潜熱の表。
6. `cuda_forge/condensationTransport_d.cu`: clamp の float 実体、`scalar_advection_multi_d<4>`、`cond_primitive_d` 融合。
7. `input/solverConfig.{hpp,cpp}`: `condFloat`; `procedures/solver-settings.md`; `tools/perf_regress.py` の除外リスト; `tools/cond_axis_h0.py`。
8. `tests/unit/test_cond_float_device.cu`; `methods/condensation.md` 実装 §9 (v2 反映); `methods/architecture/performance.md`。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 0 | codex plan レビュー | **済** (2026-09-13, GO-with-changes M7/m2, 全件採用 → §4.2 v2) |
| 1 | ① 物性表 + CNT 対数化 + host 単体試験 (§4.2-1,2; codex M1/M2) | **済** (8b32f3f1): 表は全物性 ≤2e-6 (ln p_sat は float 表現 5 ulp を含めた形), 連鎖 (J, r*, dr/dt, S ベクトル, 蒸発 λ, T_sat) ALL PASS。許容の最終形は methods/condensation.md 実装 §9 |
| 2 | ② 反転の成功判定・`roe` 保護・`condFloat` 分岐表 (§4.2-5,6; codex M3/M5) | **済** (8629f569, b021bca1): 研磨は残差 1e-3 J/kg まで → |ΔT|/T ≤ 7.5e-12、10 往復ドリフト = double 反転と同値、失敗 0/360; 失敗セルは roe 不変 + `g_condTinvFail`。分岐表どおり `condFloat: 0` で source/clamp 診断が変更前とビット一致 (run_0458 凍結 1 step) |
| 3 | ①③ ソース kernel・clamp・SLAU 面潜熱 (表, g=0 スキップ)・移流融合 (§4.2-3,4,7) | **済** (8b32f3f1, 8629f569, 0449725a): float 実体は template でなく別関数 (double はビット不変)。蒸発は δ=λ−1 で組む (float で λ=1 に丸まる問題を device 試験で検出)。移流 4 起動→1、原始量 8→1 |
| 4 | ③ device 単体試験 + `condFloat: 0` のビット一致 + FP64 命令監査 (§4.3; codex M7/m9) | **済**: `test_cond_float_device.cu` (H2O carrier 7426 状態 × Kw1/HK・Kw3/Gyar、N2 CPG 3744 状態 × Kw1/Goodheart・Kw0/Gyar) ALL PASS; `cuobjdump -sass` (8b32f3f1 時点, 退避なし): `condensation_source_f_d` FP64 0/6432 命令, `cond_realizability_clamp_f_d` 0/360 (double 実体 16802/40016, 617/1552)。**最終形 (SHA 8c709692…, 範囲外 double 退避と蒸発 Jacobian の double 摂動を含む)**: float kernel 19943/53600, clamp 931/2496 — いずれも範囲外・蒸発分岐に閉じた命令で、表範囲内の dry/核生成セルは通らない (速度で確認: 起動区間 84→86–88, 発達場 113→110–112 ms/step) |
| 5 | ④ 短期回帰 (静かな case) + 物理検証 (onset / h0 / 定常性; codex M6) | **ローカル予備** (RTX 3060, 2026-09-13): case/16 run_0456 (node 2D H2O) PASS 32/32 (T 9e-7 ≤ ノイズ 2e-6), case/44 run_0201 PASS 30/30, `cond_axis_h0.py` run_0456: 軸 h0 の非保存 271.8 J/kg (9.1e-4 of 3e5) で新旧差 −0.04 J/kg PASS。case/34 (12000 step, 床 = 同一バイナリ反復): cell 空気 PASS 28/28 (T ノイズ ×1.1, ρg ×1.2 = マージ後バイナリと同等), node 空気 PASS 28/28 (T 1.7e-6 ≤ ノイズ 1.8e-6), cell dry PASS 20/20; 反転失敗 0。**pure N2** (新設 `run_0103_condf_regress_n2_cell` [床 run_0045/0046], `run_0104_condf_regress_n2_node` [床 = マージ後バイナリ反復 2 本]): cell PASS 28/28 (T ノイズ ×0.7, ρg ×0.6), node PASS 28/28 (T 1.4e-6 ≤ 床 1.5e-6, ρg 4.7e-5 ≤ 4.2e-5×2)。**onset** (`onset_analysis.py --series`, 全 run STEADY): cell は float 版が 2.112/2.131 in を交互 (1 セル = 0.019 in) するが、基準の反復 run_0045 も同じ 2 値を取る (再現幅内)、T_on 39.52–39.67 K・p/p_dry @3/4/5 in は 4 桁一致; node は 2.102 in / 39.59 K で float・マージ後・基準が同一。 **2D node H2O 12000 step 物理検証** (`case/16.nozzle_wys/run_0460_condf_phys_node2d_cond`, run_0213/res_24000 IC, ラベル merged/merged_r2/condf): `check_convergence` は両方 NOT CONVERGED (plateau; rms_roe 6.0e-4 vs 6.2e-4, rms_rog 2.2e-11 vs 2.3e-11 で同区分), `check_quasisteady` 両方 ALL STEADY, `cond_axis_h0.py`: 軸 h0 の非保存 271.9 J/kg (9.1e-4 of 3e5) で新−基準 +0.07 J/kg (許容 100) PASS・末尾変化 0.1 J/kg STEADY, 12000 step 後の場の差は `cmp --noise merged_r2` で PASS 32/32 (T 5.6e-6 は絶対基準 A、Y1 1.2e-5 はノイズ ×0.7、ρg 1.3e-5 は**絶対基準 B 1e-4 で合格** [ノイズ床 5.7e-6 の 2 倍は超える; codex result m8 で訂正])。**AWS での正式取得 (pure N2 cell/node 新設・12000 step 物理検証・反復床) はインスタンス停止中 (44.211.54.88 応答なし) で未了** |
| 6 | ⑤ 発達場・起動区間の速度 + 証拠回収 (codex m8) | **ローカル RTX 3060 (AWS はキャパ不足で起動不可, 2026-09-13)**: 3D 2.37M `run_0459_perf_cond3d_local` (run_0455 dry テンプレート + 凝縮 H2O carrier, 起動区間 100 step, warm-up 後 3 反復): **凝縮 ON 197.3/197.2/198.1 → 83.9/84.0/84.1 ms/step** (2.35 倍), dry `run_0455` 62.8 → 63.3 (不変)。凝縮固有の超過 135 → 21 ms/step (−84 %)。A10G 換算 (帯域比 ≈0.53): 凝縮 ON ≈ 34 + 11 ≈ 45 ms/step の見込み → AWS 起動後に run_0420/0421/0400 で確定。証拠: `run_0459_perf_cond3d_local_bench/{m,f}_14*_r*_n100/bench_*.log` |
| 7 | codex result レビュー (1 回目) | **済**: NO-GO (M6/m2, 2026-09-13) → 全件採用、#8–#13 に展開。再レビューは #8–#13 完了後 |
| 8 | 表の範囲外は double 関数へ退避 (codex result M1) | **済** (ecd56f6a, 2 回目 M1 の退避漏れを修正): 範囲外セルの dry 判定は旧 double 飽和圧で行い (表の端値だと N2 15 K・p_v=0.9 p_sat(20 K) の過飽和セルを未飽和と誤判定して核生成を消した)、dry 以外は委譲。device 試験にその反例 (15/18 K, 液相なし, 過飽和) を追加し範囲外 28 状態がビット一致。旧 kernel 本体を `condensation_source_cell_d` (`__device__`) にし、float kernel は湿潤セル (g>0 または Q0>0 または S>1) の T が表範囲 [Tmin, TwetMax] (N2 20–125.6 K, H2O 120.15–647 K) の外なら丸ごと委譲 (device 試験: 範囲外 32 状態が両実体でビット一致)。面潜熱 (`cond_latent_tab_or_d`) と clamp も範囲外は旧 double 関数。ln p_sat の dry 診断だけ表範囲 [Tmin, Tmax] (H2O 1200 K) で表 (外は端クランプ, 判定対象外)。双方の単体試験は実温度で比較 (ALL PASS)。double kernel の SASS は変更前バイナリと**同一** (40016 命令, `cuobjdump`) |
| 9 | 二相反転の成功条件を全経路で統一 (codex result M2) | **済** (ecd56f6a): `cond_twophase_polish` (double Newton を tol=1e-9|e|+0.05 J/kg まで最大 10 段) を退避後と `condFloat: 0` 経路の両方に適用、満たさなければ失敗。強制失敗試験: host (e が T_max 外 → ok=0 両経路)、device 凍結 `run_0461_frozen_roe_protect` (湿潤 20 セルの roe を +5e7 J/kg → 20/20 で roe 不変・T=6000 K 張り付き・警告 "inversion failed in 20 cells" [TP でも報告するよう修正]) |
| 10 | `cond_axis_h0.py` の cell 対応と参照側の定常性 (codex result M3) | **済** (4f927e86): cell 場は `CELLS/centCoords`、新旧とも ≥3 枚・有限・末尾定常・共通 step を要求、合否は「参照 STEADY かつ new−ref ≤ 100 J/kg かつ new ≤ max(0.1 % h0scale, ref+10 J/kg)」(参照が既に 0.1 % を超える cell 場では悪化しないこと)。再取得: case/34 cell 空気 428.5 vs 434.1 (−5.6), N2 cell 482.4 vs 482.0 (+0.4), N2 node 128.9 vs 128.9 (−0.02), case/16 2D node 271.9 vs 271.9 (+0.07) J/kg — 非劣化ゲートは全 PASS。**絶対ゲート (≤300 J/kg) は cell 2 件で FAIL** (基準バイナリも同値: 既存の cell 凝縮の h0 非保存; ツールは両判定を別々に出す) → condensation-followups へ転記 (codex result-2 M3) |
| 11 | 単体試験の有限性検査と蒸発端の陰的更新差 (codex result M4) | **済** (ecd56f6a): device 試験は除外判定の前に全出力の有限性を検査 (非有限 0)、蒸発端 0.99<S<1 は |Δsj|·dt を実測 (最大 3.4e-5、許容 1e-4: Δg/ΔT 律速の境界で分岐が ULP で変わる)。蒸発 Jacobian は float kernel でも double で摂動評価 (数値微分の相殺で 1e-4 ずれるため)。蒸発 λ の許容 5e-4 (S=0.999 の (p_v−p_d) 相殺で 1.4e-4 が出る) を §4.3 に明記 |
| 12 | 互換性・保存性・cell H2O 回帰 (codex result M5) | (a) **済**: `condensation_source_d` の SASS が変更前と同一 (40016 命令)、凍結 1 step の凝縮診断 6 量と T/P ビット一致 (`run_0458`); (c) **済**: 2D cell H2O `run_0457_perf_regress_cell2d_cond` (run_0196 IC + 凝縮, 基準 4 本): PASS 30/30 (roUy ノイズ ×1.2–1.6, 基準同士も 1.3e-4); (b) **済**: 周期箱 `case/09.Taylor-Green/run_0053_cond_periodic_conservation` (run_0052 node 周期 dual-time + N2 CPG 凝縮 `condEquilibrium: 2` [Q0–Q2 は移流のみ], 非一様 ρQ の IC): 1 step の 4 モーメント全セル和の変化は融合前 (merged) と融合後 (condf) で ρQ0/Q1/Q2 とも 5.60e-6 で一致、差 2〜5e-10 = 同一バイナリ反復ノイズ (3〜6e-10)、ρg (平衡射影) −1.3e-7 vs ノイズ 2e-8。**保存ゲート (20 step ≤1e-6) は未達**: 1 step でも +5.6e-6、20 step で +1.1e-3 が両バイナリで出る (融合前後で同値 → 融合の差ではないが、node 周期の凝縮モーメント輸送自体が非保存)。原因の切り分け (処理別収支: 移流 / dual-time 部分反復 / clamp / seam) は本 plan のスコープ外として condensation-followups §5.1 に転記 (2026-09-13)。融合の同値性は「同一凍結状態の 1 step 和の一致 (差 ≤5e-10 = 反復ノイズ)」で担保 |
| 13 | 発達場の速度・湿潤統計・A10G (codex result M6) | **ローカル済 / A10G 未達成のまま**: `run_0462_cond3d_local_spinup` (condf, 3000 step, 反転失敗 0) の res_3000 を IC にした `run_0463_perf_cond3d_local_dev` (湿潤セル 542k / 2.37M = 22.9 %) で 100 step ×3 交互 (RTX 3060): **merged 220.6/221.6/220.8 → condf 114.4/112.1/114.6 ms/step** (1.94 倍; 起動区間は 197→84)。dry 63 に対する凝縮固有の超過は発達場 158 → 50 ms/step。反転失敗 0。A10G の ≤50/≤45 ms/step は未計測 (インスタンス起動不可) で、帯域比換算 (発達場 ≈ 34+27 ≈ 60 ms, 起動 ≈ 45 ms) は見積りとしてのみ記載し達成判定には使わない |
| 14 | 文書同期と記録の修正 (codex result m7/m8, 2 回目 m5/m6) | **済** (3 回目 m3 で残っていた λ 5e-4 / 蒸発端 1e-4 / 平衡形の記述 / 変更ログを再同期): performance.md / plans/README / solver-settings.md (`condFloat`) / plan §4.2-5 表記 / 平衡形の扱い (`condEquilibrium 1` の湿潤 TP セルはハイブリッド反転、2 は反転も double) / methods §9 の範囲外退避・蒸発端許容 / §4.3 の λ 5e-4 と h0 二重判定; run_0460 の ρg は絶対基準で PASS (ノイズ床内ではない) と訂正; 各ラベルに VERDICT ファイル保存 (§5.1 #15) |
| 15 | 最終バイナリでの証拠統一 (codex result-2 M4/m6) | **済 (ローカル)**: 最終バイナリ **SHA 8c7096927b03b5d6…** (ecd56f6a + 範囲外 dry 判定修正) で RTX 3060 100 step ×3 交互: dry `run_0455` 62.2–62.7 → 61.8–62.5 (不変), 起動区間 `run_0459` 194.8–196.0 → 85.9–88.1, 発達場 `run_0463` 215.2–220.5 → 110.5–111.9 ms/step (反転失敗 0)。凍結 `run_0458` cf0 ビット一致・`run_0456` / `run_0201` PASS を同 SHA で再取得。VERDICT ファイル (`CONVERGENCE_VERDICT.txt` / `QUASISTEADY_VERDICT.txt` / `CMP_VERDICT.txt` / `H0_VERDICT.txt`) を run_0456/0457/0460・run_0201・case/34 run_0100–0104 の判定ラベルに保存。**A10G は未計測 (§8 に明記)** |

## 6. 検証

- **単体 / ビルド**: `test_cond_float.cpp` (host: 表精度・連鎖・反転ドリフト)、`test_cond_float_device.cu` (device: 両実体の比較・ビット一致)、
  既存 `test_cond_air` / `test_cond_kantrowitz_carrier` / `test_cond_evaporation` / `test_cond_equilibrium_eos` / `test_cond_sonic` ALL PASS (double API 不変)。
  `cuobjdump -sass` の FP64 命令監査 (§4.2-8)。
- **短期回帰 (AWS A10G, `perf_regress.py`, 基準 = 変更前バイナリ, ノイズ床は基準の反復)**:
  - 2D node H2O TP: `case/16.nozzle_wys/run_0456_perf_regress_node2d_cond` (基準 2 本あり); 2D cell H2O TP: `run_0457` (run_0192 config, 新設, 基準 3 本)。
  - case/34: 空気 cell/node・dry cell (`run_0100–0102`, 12000 step, 床 = 反復 3 本) + **pure N2 cell/node** (`run_0103/0104`, run_0029/0032 config, 床 = run_0045/0046 + 新反復)。
  - Euler + 平衡: `case/44.vitiated_air_wt/run_0201` (double 経路のまま; 差はノイズ床内であること)。
  - 3D H2O: `run_0421_perf_regress_cond3d` は分布統計と速度のみ (基準同士 1e-2)。
- **物理検証**: case/34 pure N2 cell / 空気 node の 12000 step を新版で再取得し onset (`onset_analysis.py --series`) と反復幅; case/16 2D node H2O は
  run_0213 IC から 12000 step (基準・新版) で `check_convergence` / `check_quasisteady` / `cond_axis_h0.py` (§4.3)。
- **速度**: `bench_steps.sh run_0421_perf_regress_cond3d 100` (発達場) と `run_0420_perf_cond_baseline 100` (起動区間) を基準/新版で交互 2 回、dry `run_0400` 再測。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result (3 回目) | `2026-09-13` | [`notes/reviews/2026-09-13-condensation-float-speedup-result-3.md`](../../notes/reviews/2026-09-13-condensation-float-speedup-result-3.md) | **NO-GO**, C0/M1/m2 | **全件採用**: M1 (A10G 目標は完了条件のまま → `in_progress` 維持、AWS 復旧後に最終バイナリで dry/起動/発達場を実測: §5.1 #13) / m2 (h0 ツールを `ABSOLUTE_VERDICT` と `NONDEGRADATION_VERDICT` の二重表示に、H0_VERDICT.txt 再生成: cell 2 件は ABSOLUTE FAIL / NONDEGRADATION PASS) / m3 (§4.3 の λ 5e-4、methods §9 の蒸発端 1e-4、solver-settings の平衡形の記述、変更ログの「ノイズ床内」を訂正)。ローカル高速化 (起動 195→87, 発達場 218→111, dry 不変) と範囲外 dry 判定修正は codex 側で再現 |
| result (2 回目) | `2026-09-13` | [`notes/reviews/2026-09-13-condensation-float-speedup-result-2.md`](../../notes/reviews/2026-09-13-condensation-float-speedup-result-2.md) | **NO-GO**, C0/M4/m2 | **全件採用**: M1 (範囲外セルの dry 判定を表の端値でなく旧 double 飽和圧で → 修正, 反例を device 試験に追加: §5.1 #8) / M2 (周期保存ゲートは未達と明記、既存の非保存を followups へ転記、融合の同値性は 1 step 和の一致で: #12) / M3 (h0 は非劣化ゲートと絶対ゲートを別判定に、cell の絶対 FAIL を既存課題として followups へ: #10) / M4 (最終バイナリで dry/起動/発達場を再計測、A10G は未検証のまま: #15) / m5 (methods §9・§4.3 の同期: #14) / m6 (VERDICT ファイル保存、README/変更ログの「ノイズ床内」訂正: #15) |
| result | `2026-09-13` | [`notes/reviews/2026-09-13-condensation-float-speedup-result.md`](../../notes/reviews/2026-09-13-condensation-float-speedup-result.md) | **NO-GO**, C0/M6/m2 | **全件採用** → §5.1 #8 (表範囲外の double 退避) / #9 (反転成功条件の統一) / #10 (h0 ツールの cell 対応・参照定常性) / #11 (単体の有限性・蒸発端) / #12 (SASS 同一・周期保存・cell H2O) / #13 (発達場速度・A10G 未達成の明記) / #14 (文書同期・記録訂正)。短期回帰 PASS と run_0460 の h0 は codex 側で再現 |
| plan | `2026-09-13` | [`notes/reviews/2026-09-13-condensation-float-speedup-plan.md`](../../notes/reviews/2026-09-13-condensation-float-speedup-plan.md) | **GO-with-changes**, C0/M7/m2 | **全件採用**: M1 (物性の float 直接評価は L 9e-3・psat 3e-4 の誤差 → 区分 3 次表方式に変更, §4.2-1) / M2 (CNT の対数上限を本体・摂動の両方に, 連鎖誤差の単体評価, §4.2-2) / M3 (二相ハイブリッド反転に成功判定・double 継続・roe 保護・カウンタ・ドリフト試験, §4.2-5) / M4 (condL キャッシュ不採用: 境界処理で T が変わる; 面内で表評価 + g=0 スキップ, §4.2-4) / M5 (`condFloat` の分岐表とビット一致は凍結状態評価のみ, §4.2-6) / M6 (短期回帰と物理検証の分離: onset・h0 の VERDICT と絶対許容 [J/kg], `cond_axis_h0.py`, §4.3) / M7 (試験範囲を実測 [N2 39 K, g 0.06, S 171] まで拡張, pure N2 cell/node 回帰を新設, device 試験, 移流融合の保存性, §4.3/§6) / m8 (発達場の速度 108 ms/step と非決定性 1e-2 を計測し §4.1 に記録, 証拠回収先を明記) / m9 (`condTheta/condLim` を除外リストへ, FP64 命令監査を DFMA 以外に拡張, §4.2-8)。実装順は codex 推奨 (①物性・CNT → ②反転 → ③device 試験 → ④定常性 → ⑤性能) に揃えた |

## 7. 影響範囲

- `cuda_forge/condensation{Properties,Source,EOS,Transport}_d.cu{,h}`, `dependentVariables_d.cu`, `convection/convectiveFlux_slau_d.inc.cuh`,
  `convection/convectiveFlux_common_d.cuh` (`CondArgs`), `convection/convectiveFlux_d.cu`, `input/solverConfig.*`, `tests/unit/test_cond_float.cpp`。
- 既存ケース: `condensation: 0` はビット不変。凝縮 ON は既定 `condFloat: 1` で float 経路 (`condFloat: 0` で旧値)。
- docs: `methods/condensation.md` (実装 §9 新設), `methods/architecture/performance.md`, `procedures/solver-settings.md` (`condFloat`),
  `plans/README.md`。

## 8. 完了条件

- [ ] `methods/condensation.md` 実装 §9 と `performance.md` を更新済み
- [ ] §5.1 #1–#5 完了、§4.3 の判定 (場・物理量・単体・速度・ビット不変) を満たす — **速度の A10G 目標 (発達場 ≤50 / 起動 ≤45 ms/step) は未計測** (インスタンス起動不可; ローカル RTX 3060 の 1.9–2.3 倍だけ)。周期保存ゲート (≤1e-6) と cell の h0 絶対ゲートは既存挙動で未達 (followups F-cf1/F-cf2)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し Critical / Major の採否を §5.1 に反映
- [ ] `status: done`、§9 変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-13` — 起票。A10G で dry 33.8 / 凝縮 ON 89.5 ms/step、nsys で凝縮ソース 31.7・SLAU 面潜熱 11.5・反転 10.2 ms/step を同定 (§4.1)。
- `2026-09-13` — codex plan レビュー GO-with-changes (M7/m2) を全件採用し §4.2/§4.3/§5/§6 を v2 に改訂 (物性は区分 3 次表、CNT は対数上限、反転に成功判定と保存量保護、condL 不採用、3D は速度/統計のみ)。
- `2026-09-13` — 実装 (8b32f3f1 → b021bca1): 物性表 + 対数 CNT の float ソース kernel、float 面潜熱 (g=0 面スキップ)、float clamp、二相ハイブリッド反転 (残差判定・roe 保護)、
  移流/原始量の融合、`condFloat`、host/device 単体試験 ALL PASS、FP64 命令 0。実装上の差 (template でなく別関数、蒸発 δ 形、研磨 tol 1e-3 J/kg、
  単体許容の最終形) は methods/condensation.md 実装 §9 に記載。ローカル予備回帰 PASS (run_0456 / run_0201 / h0 保存)。AWS (速度・正式回帰・物理検証) はインスタンス停止で未了。
- `2026-09-13` — AWS がキャパ不足で起動できないため速度と物理検証をローカル RTX 3060 で実施: 3D 凝縮 ON 197 → 84 ms/step (dry 63 不変)。pure N2 cell/node 回帰 (`run_0103/0104`) と 2D node H2O 12000 step 物理検証 (`run_0460`) を取得中。
- `2026-09-13` — ローカル検証完了: pure N2 cell/node 回帰 PASS、onset STEADY (1 セル交互は基準反復と同じ)、2D node H2O 12000 step で同 verdict 区分・ALL STEADY・h0 新旧差 +0.07 J/kg・場は 32/32 PASS (ρg は絶対基準 1e-4 で合格; ノイズ床 5.7e-6 の 2 倍は超える — 2026-09-13 訂正)。codex result レビューへ。
- `2026-09-13` — codex result レビュー 1 回目 NO-GO (M6/m2) 全件採用: 表範囲外の double 退避、反転成功条件の統一、h0 ツールの cell 対応、単体の有限性、SASS/周期保存/cell H2O、発達場速度、文書同期 (§5.1 #8–#14)。
- `2026-09-13` — codex result 1 回目の全項目 (#8–#14) を実装・取得: 表範囲外 double 退避 (SASS 同一), 反転成功条件統一 (強制失敗 20/20 で roe 保護), h0 ツール cell 対応 (4 case PASS), 単体の有限性, cell H2O 回帰 PASS, 周期保存 (融合前後 1 step 一致), 発達場 221→113 ms/step (RTX 3060)。A10G は未計測。2 回目の result レビューへ。
- `2026-09-13` — codex result 2 回目 NO-GO (M4/m2) 全件採用: 範囲外セルの dry 誤判定を修正、周期保存と cell h0 絶対ゲートは「未達 (既存挙動)」として followups へ転記、最終バイナリで速度再計測、VERDICT 保存。
- `2026-09-13` — 最終バイナリ SHA 8c7096927b03b5d6 で速度・凍結・回帰を再取得 (起動 195→87, 発達場 218→111, dry 62 不変 ms/step @RTX 3060)、VERDICT を各ラベルに保存。codex result 3 回目へ。
- `2026-09-13` — codex result 3 回目 NO-GO (M1/m2) 全件採用: A10G 実測まで `in_progress` 維持 (AWS 復旧 34.201.22.141 で計測開始)、h0 ツールの二重 VERDICT、閾値の文書同期。
