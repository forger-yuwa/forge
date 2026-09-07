# 有限速度化学 (H₂ 燃焼・ノズル上流の化学非平衡) 導入 — 文献調査と方針

<!-- ファイル名規約: <area>-<short-slug>.md -->

## メタ

- **area**: `thermophysics / chemistry` (新規 area 候補: `methods/chemistry.md`)
- **status**: `draft` (調査・方針。Phase 0 の結果は §9)
- **related_docs**:
  - `methods/thermophysics.md` (多成分 TP gas: NASA-9・種輸送・sensible datum・陰解法結合)
  - `methods/condensation.md` (ソース項 + `src_jac` 配管の既存例)
  - `methods/time_integration/` (block-DPLUR / general-EOS Jacobian)
- **related_plans**:
  - `plans/accepted/thermophysics-multicomponent-tpgas.md` (§10 に「化学反応源項は将来拡張点」と明記)
  - `plans/accepted/convection-multispecies-contact-pressure.md`
- **created**: `2026-09-04`
- **owner**: `CFD Dev`
- **PDF**: [`papers/combustion/`](../../papers/combustion/) (Jachimowski NASA TP-2791, US3D NTRS 2018)

---

## 0. 結論先出し

1. **forge には反応流に必要な土台がほぼ揃っている**。多成分 TP (NASA-9 絶対エンタルピー・エントロピー $a_8$ 込み)、
   化学種輸送 (`speciesTransport_d`)、種ごとのソース項と対角 Jacobian の配管 (`src_jac_Y{s}` → `species_dplur_sweep_d`)、
   種↔エネルギーの陰解法結合 (案C `speciesImplicitCoupling=2`: EOS-JVP $\delta p_Y$)、SST/凝縮で実績のある
   segregated point-implicit ソース項の型。**追加するのは「反応ソース項 $\dot\omega_s$ とその Jacobian」「反応熱のエネルギー式への注入」
   「熱力学 datum の整合」「反応機構/熱力学データの入力配管」の 4 点**で、ソルバ構造の変更は不要。
2. **推奨アーキテクチャ**: 自前 device 実装 (Cantera/Mutation++ 等の CPU ライブラリは GPU セル単位評価に使えない。
   SU2-NEMO は Mutation++ 依存だが host 側評価)。反応機構は Cantera YAML 互換のサブセットを読む自前パーサとし、
   **Cantera (pip, CPU) は参照解 (0-D 着火遅れ・PFR ノズル) と機構ファイルの供給源としてのみ使う**。
3. **剛性処理は point-implicit (種ブロック $n_s\times n_s$ 解析 Jacobian) を主、演算子分離 (Strang + sub-cycling) を陽解法向けの副**とする。
   forge の生産経路は定常陰解法 (block-DPLUR) であり、Candler ら (2013) の *decoupled implicit* (流れ 5 変数 → 化学種を逐次) が
   forge の既存構造 (flow block → species DPLUR) と同型。既存の `src_jac_Y{s}` (対角のみ) を **種ブロック密行列に拡張**するのが最小変更。
4. **熱力学 datum**: forge は陰解法安定化のため sensible-enthalpy datum (`thermoHrefTemp: 298.15`, 係数 $a_7$ に焼き込み) を
   必須にしている。反応ありでは平行移動が相殺しないので、**sensible datum を維持したままエネルギー式に反応熱
   $\dot Q=-\sum_s h^{\circ}_s(T_{\rm ref})\,\dot\omega_s$ を陽に加える** (絶対 datum に戻すと $\chi_{\rm eos}$ 起因の陰解法発散が再発する)。
5. **段階計画**: Phase 0 (CEA による効果スクリーニング + 熱力学 DB 生成 + 機構パーサ) → Phase 1 (反応ソース項 + 0-D 検証) →
   Phase 2 (陰解法結合 + Q1D ノズル再結合検証) → Phase 3 (燃焼器 RANS: laminar finite-rate → PaSR) → Phase 4 (NOx・超音速燃焼)。
   **Phase 0 の CEA スクリーニング (平衡 vs 凍結 の出口 M/T 差) で、対象案件に有限速度化学が本当に必要かを先に確定する**
   ($T_t\approx1160$ K の case/44 では平衡解離がほぼ無く効果は小さい。効くのは $T_t\gtrsim1800$–2000 K)。

---

## 1. 対象アプリケーションと物理

### 1.1 想定する 3 つの用途

| 用途 | 流れ | 化学の役割 | 難所 |
| --- | --- | --- | --- |
| (A) 燃焼加熱器 (vitiated air heater: H₂ + O₂ + air) | 亜音速・低 M・高 $T$ | 着火・燃焼完了・NO 生成 | 低マッハ (`lowMachPrecond=2` 併用)、着火 (IC)、TCI |
| (B) ノズル膨張中の化学非平衡 (再結合・凍結) | スロート付近で凍結、下流は frozen | 試験部組成 (OH/H/O/NO 残留)・$T$・$M$ への影響 | 三体再結合の剛性、凍結点の解像 |
| (C) 超音速燃焼 (将来: 試験部での燃焼器試験) | 超音速・衝撃波誘起着火 | 着火遅れ・保炎 | 剛性最大、TCI |

案件の近さから **(B) を第一目標、(A) を第二目標** とする。(B) だけなら燃焼器は解かず、
**燃焼室出口を CEA 平衡組成の入口 BC** として与え、ノズル内の有限速度再結合だけを解けば済む
(TDK/ODE 型の古典的アプローチ。Bray の sudden-freezing 判定はその簡易版)。

### 1.2 ノズル化学非平衡の物理 (文献の要点)

- 燃焼室 (高 $p$, 高 $T$) では平衡に近い。膨張で $T$ が下がると再結合 (H+OH+M→H₂O+M, H+H+M, O+O+M, H+O₂+M→HO₂+M)
  が発熱するが、三体反応は $\rho^2$ に比例するため低圧化で急減し、**スロート付近〜直後で凍結**する
  (Bray 1959 の sudden freezing; NASA の評価 NTRS 19650044133)。凍結後は組成一定の frozen 膨張。
- 結果として frozen 計算は再結合熱を取り逃がし $T$ を低く見積もり、平衡計算は逆に $T$ を高く見積もる。
  有限速度は両者の間で、**ロケットノズルの $I_{sp}$ で 1–3 % 程度、H₂/O₂ で顕著** (研究例: "Numerical study of a nonequilibrium H₂–O₂ rocket nozzle flow")。
  スクラムジェットノズルでは未燃燃料の後燃え + 再結合で frozen が性能を過小評価する
  (Chemical non-equilibrium flow analysis of H2 fueled scramjet nozzle, 13 種 33 反応)。
- 燃焼加熱風洞では試験部に H₂O・CO₂ に加え **ラジカル (O, H, OH) と NO が残留**し、模型燃焼器の着火・保炎を変える
  (Pellett et al. AIAA 2002-3880 "Review of Air Vitiation Effects on Scramjet Ignition and Flameholding")。
  8-ft HTT は CH₄/air/LOX 燃焼加熱で M7、有限速度で NO を評価 (Langley arc-heated/8-ft HTT の資料)。
- **効果の大きさは $T_t$ と $p_t$ で決まる**: $T_t\approx1160$ K (case/44 va3) では平衡 OH は ppm 未満、解離無視可。
  $T_t\approx1600$ K (case/45) で OH ~$10^{-4}$ 級、$T_t\gtrsim2000$ K (M7 級加熱器) で再結合熱が M・T に 0.5–数 % 効く。
  → **Phase 0 で CEA の `rocket` 問題 (equilibrium と frozen の出口 $M$, $T$, 組成) を回し、案件ごとに効果の上限を先に数値化する**
  (CEA2 は `.venv-cea/nasa_cea/FCEA2` にビルド済)。

---

## 2. 文献調査

### 2.1 H₂ 反応機構 (どれを積むか)

| 機構 | 種 / 反応 | 特徴 | forge での位置付け |
| --- | --- | --- | --- |
| Jachimowski 1988 (NASA TP-2791) | 13 種 / 33 反応 (N 化学込み)、H₂/O₂ 部分は 9 種 / 19 反応 | スクラムジェット CFD の事実上の標準。三体効率あり、**falloff 無し** (単純 Arrhenius+M) | **最初の実装対象** (falloff 不要でカーネルが単純。文献との比較性が高い)。N 化学 (N, NO, NO₂, HNO) は加熱器の NO 評価で使う |
| Evans–Schexnayder 1980 | 7 種 / 8 反応 | 最小構成。HO₂/H₂O₂ 無し | デバッグ・速度検証用の縮約 |
| Ó Conaire et al. 2004 | 10 種 / 21 反応 | 広い $p$ 範囲、Troe falloff あり | 生産候補 |
| Burke et al. 2012 | 9 種 (+bath) / 27 反応 | 高圧 (>10 atm) 検証済、H+O₂(+M) の Troe falloff と bath 依存効率 | **高圧加熱器の生産候補** (Phase 2 以降) |
| Kéromnès 2013 / Konnov 2019 | 〜10 種 / 30 反応前後 | 最新の再校正 | 任意 |

Zhukov 2012 (ISRN "Verification, Validation, and Testing of Kinetic Mechanisms of Hydrogen Combustion in Fluid-Dynamic Computations")
は 1 段/2 段/縮約/詳細 4 種を 1-D 火炎・着火遅れで比較し、**詳細機構でも H₂ は 9–10 種で済み CFD コストの支配要因ではない**、
1 段総括反応は着火遅れを再現できない、と結論している。**総括 1 段反応は採用しない** (ノズル再結合は本質的にラジカル反応)。

実装に必要な反応形式: (i) 修正 Arrhenius $k=AT^b\exp(-E_a/RT)$、(ii) 三体 $M$ + 種別効率、(iii) 逆反応を平衡定数
$K_c=K_p\,(p^\circ/R_uT)^{\Delta\nu}$, $K_p=\exp(-\Delta G^\circ/R_uT)$ から算出 (NASA-9 の $a_7,a_8$ が必要 → forge の DB に既にある)、
(iv) falloff (Lindemann / Troe) は Phase 2 で追加。

### 2.2 剛性 (stiffness) の扱い

化学時間: 加熱器 ($T\sim2500$ K, $p\sim1$ MPa) でラジカル反応 $\tau_c\sim10^{-9}$–$10^{-8}$ s。
forge 定常陰解法の局所 $\Delta t$ は `cfl_pseudo` 2–8 で $\Delta x/(u+c)\sim10^{-7}$–$10^{-6}$ s。**比 $10^1$–$10^3$ で剛性**、
陽解法 RK では化学だけで $\Delta t$ が 2–3 桁縮む。文献上の選択肢:

| 方式 | 内容 | 代表 | forge への適合 |
| --- | --- | --- | --- |
| **Point-implicit (ソース項の陰化)** | セル毎に $(I/\Delta t-\partial\dot\omega/\partial U)\delta U=R$。対流は陽 or 別途陰 | Bussing & Murman 1988 (AIAA J 26:1070), Eberhardt & Imlay 1992 | **既存 `src_jac_Y{s}` 配管の拡張で実装可** (現状は対角のみ → 種ブロック密行列へ) |
| **Fully coupled block implicit** | $(5+n_s)\times(5+n_s)$ ブロック DPLUR | DPLR / LAURA / US3D 既定 | block-DPLUR の 5×5 を可変サイズ化する大改修。コスト $O(n_s^2)$。**採らない** |
| **Decoupled implicit** | 流れ 5 変数を先に解き、化学種 (+内部エネルギー) を後から逐次陰解 | Candler, Subbareddy, Nompelis AIAA J 2013 (doi 10.2514/1.J052070); US3D のオプション、fully coupled 比 ~3 倍速・メモリ小 | **forge の flow block → `speciesImplicitDPLURSolve` の順序そのもの**。案C の EOS-JVP はこの「流れ側が組成変化を見る」補正に相当 |
| **演算子分離 (Strang)** | 対流/拡散 → 化学 ODE (CVODE / backward Euler + sub-cycling) → 対流/拡散 | AMReX 系 (arXiv 2412.00900), XFLUIDS, OpenFOAM 系 | 非定常 (LES/DES, 燃焼振動) 向け。定常擬似時間では 1 次分離誤差が定常解に残る |
| Component-splitting implicit | 種ごとに陰化した半陰解法 | arXiv 2403.03440 | 参考 |

**forge の結論**: 定常は **decoupled point-implicit** (種ブロック解析 Jacobian、流れ側へは案C の EOS-JVP に $\dot Q$ の温度応答を追加)、
非定常は **Strang 分離 + セル内 backward Euler sub-cycling** (同じ Jacobian を再利用) の 2 経路。どちらもカーネルは同じ
$\dot\omega_s$ と $\partial\dot\omega_s/\partial(\rho Y_k), \partial\dot\omega_s/\partial T$ を使う。

### 2.3 化学 Jacobian

- 解析 Jacobian が有限差分より 3–7 倍速く、GPU では疎行列 SIMT 化で更に効く (pyJac, Niemeyer et al. CPC 2017;
  Perini et al. Energy & Fuels 2012 の三体・Troe を含む閉形式)。**H₂ 機構 ($n_s\le13$) は密 $13\times13$ で十分**、
  セル毎に LU (部分ピボット) を device で解く。
- 温度微分: $\partial k_f/\partial T=k_f\,(b/T+E_a/R T^2)$、$\partial K_c/\partial T$ は $\Delta H^\circ/R_uT^2$ (van 't Hoff)。
  $T$ は $(\rho e,\rho Y)$ の従属変数なので $\partial T/\partial(\rho Y_s)=-e_s/(\rho c_v)$ を連鎖 (既存 `species_eos_cross_response` と同式)。
- 変数の選択: 保存量 $\rho Y_s$ 基準の Jacobian にすると流れ密度更新との整合が取りやすい (案C の接空間 $z_s=\rho\,\delta Y_s$ と同じ規約)。

### 2.4 GPU 化学の実装動向

- pyJac / KinetiX / ChemGen (arXiv 2510.10005) / Pyrometheus: 機構 → C/CUDA コード生成。
  forge は機構が小さく固定 (H₂ 系) なので、**コード生成は不要、YAML 読込 → device 定数配列 (反応表) → 汎用ループ**で足りる。
  ただし 3 体効率・falloff パラメータの構造体を最初から用意する。
- preJacGPUFoam (2025) は GPU 1 枚で CPU 1 コア比 75–180 倍。コスト構造は「$\exp/\log$ 評価 × 反応数」が支配。
  forge (consumer GPU, FP64 = FP32/64) では **温度依存項 ($\ln k_f$, $\ln K_c$) を FP32 で評価し、質量作用則の積と Jacobian 解だけ FP64** にする
  (M6 の「桁落ちする量だけ double」方針を踏襲)。$\ln k$ の温度テーブル化は後回し。

### 2.5 乱流‐化学相互作用 (TCI)

- Burrows–Kurkov (1973, NASA TM X-2828, M2.44 vitiated air + H₂ 壁噴射) と DLR スクラムジェットが標準検証。
  NASA GRC の Wind 検証アーカイブに実験データがある (grc.nasa.gov/www/wind/valid/bk/)。
- 高速 H₂ 燃焼の RANS では **laminar finite-rate (No-TCI) が意外に健闘、PaSR が最良、EDC は DLR で保炎失敗**の報告
  (Acta Astronautica 2024 "Mixing time scale analysis of the PaSR model for high-speed turbulent combustion of hydrogen in vitiated air";
  同著者の "Chemical timescale analysis of the PaSR model for a hydrogen-fuelled scramjet")。
- **方針**: Phase 3 の燃焼器は No-TCI から始め、PaSR ($\kappa=\tau_c/(\tau_c+\tau_{mix})$ で $\dot\omega$ をスケール、既存 SST の $k/\omega$ から
  $\tau_{mix}$) を 1 スカラ係数として後付けする。EDC・flamelet は当面採らない。

### 2.6 他ソルバの構造 (参考)

- **SU2-NEMO**: Mutation++ にソース項を委譲、Jacobian は roadmap 段階 (Maier et al. Aerospace 2021)。forge 既存調査
  [`su2-nemo-contact-thermo-investigation.md`](su2-nemo-contact-thermo-investigation.md) 参照。同一メッシュ比較の相手にはなるが、
  化学 Jacobian 未整備で剛性ケースの収束は弱い。
- **Eilmer** (UQ): 有限速度化学を「専用の更新スキームで時間積分に結合」= 演算子分離 + 陰的 ODE。非定常・衝撃波管向け。
- **US3D / DPLR**: 上記 decoupled implicit。定常極超音速の標準。forge の目標像に最も近い。
- **AMReX 系 GPU ソルバ** (arXiv 2412.00900): 演算子分離 + Cantera 生成率。LES 向け。

---

## 3. forge の現状棚卸し (コードパス)

| 機能 | 所在 | 反応流での再利用 |
| --- | --- | --- |
| NASA-9 熱力学 ($c_p,h,s$, Newton 温度反転, 混合則) | `cuda_forge/thermo_d.{cuh,cu}` (`SpeciesThermo`, `THERMO_HD` host/device 共通) | $K_c$ 計算に $s^\circ$ ($a_8$) をそのまま使える。radical 種 (H, O, OH, HO₂, H₂O₂, N, NO, ...) は内蔵 DB に無く **`speciesDBFile` で供給** |
| sensible datum (`thermoHrefTemp`) | `thermo_d.cu` L199– ($a_7$ に $\Delta a_7$ 焼き込み) | 反応ありでは $h^\circ_s(T_{\rm ref})$ を別配列に保持し反応熱を陽に加える (§4.1) |
| 化学種輸送・拡散・エネルギー結合 | `cuda_forge/speciesTransport_d.cu` (`species_diffusion_d`, `species_energy_correction_kernel`) | そのまま。ラジカルの LJ パラメータを DB に追加 |
| 種ソース項 + 対角 Jacobian 配管 | `res_roY{s}`, `transport_diag_Y{s}`, `src_jac_Y{s}` → `species_dplur_sweep_d` (L270) / `speciesTimeIntegration_d_wrapper` | `src_jac` を対角 → 種ブロックへ拡張 (§4.3) |
| 種↔エネルギー陰解法結合 (案C) | `species_eos_coupling_d.cuh`, `speciesEOSCrossPredictInject_d_wrapper` (main.cpp L1128–) | 予測 $\delta(\rho Y)^*$ に反応寄与を含めれば流れ側が反応熱を同一 block 解で見る |
| ソース項 point-implicit の既存例 | `ransSource_d` (SST), `condensationSource_d` (`src_jac_g`: 潜熱自己抑制) | 実装パターンの雛型 |
| 実現可能性 | `species_renormalize_d` ($\sum Y=1$, クリップ) | そのまま |
| EOS 床 (`pMin`/`tMin`) | `dependentVariables_d` | 反応で $T$ が上限 (6000 K) を超えないよう $T$ 上限クランプも検討 |
| CEA 熱力学データ | `.venv-cea/nasa_cea/thermo.inp` (全種 NASA-9) | `thermo.inp → species_db.yaml` 生成ツール (未作成。`design/forge_design/gas/semiperfect.py` は手書き係数) |
| CEA 実行 | `.venv-cea/nasa_cea/FCEA2` | Phase 0 スクリーニング・平衡入口組成・平衡/凍結の参照解 |

**既知の制約**: 多成分 TP + 陰解法の安定 `cfl_pseudo` は 1–2 と低い (H₂O の生成エンタルピー結合、EOS 床律速 —
[`plans/accepted/thermophysics-multicomponent-tpgas.md`](../../plans/accepted/thermophysics-multicomponent-tpgas.md) §10、memory
`wys-tp-divergence-is-cold-not-multispecies`, `implicit-cfl-ceiling-eos-floor`)。反応流ではこれに化学剛性が上乗せされるため、
**種ブロック陰化を対角 Jacobian で済ませると CFL 上限が更に落ちる**。最初から $n_s\times n_s$ ブロックで組む理由。

---

## 4. 設計方針

### 4.1 熱力学 datum と反応熱

保存エネルギー $\rho e$ は sensible datum ($h_s(T_{\rm ref})=0$) のまま。種 $s$ の絶対エンタルピー基準との差 $c_s=h^{\rm abs}_s(T_{\rm ref})$
(生成エンタルピー相当) を `thermo_init_db` 時に別配列 `h_form_ref[s]` に保持し、エネルギー式に

$$
\dot Q = -\sum_s c_s\,\dot\omega_s \qquad [\mathrm{W/m^3}]
$$

を加える ($\sum_s\dot\omega_s=0$ なので datum 共通定数は落ちる。$c_s$ は $T_{\rm ref}$ での絶対 $h$ なので生成エンタルピー + 顕熱 $T_{\rm ref}$ 分)。
`thermoHrefTemp: 0` (絶対 datum) なら $c_s=0$ で $\dot Q=0$ となり、同じコードで両基準が整合する。
IC 生成器・BC の $h$ 評価は既に datum を共有しているので変更不要。

### 4.2 反応ソース項カーネル `chemistrySource_d` (新規 `cuda_forge/chemistrySource_d.{cuh,cu}`)

- 入力: セルの $(\rho, T, Y_s)$、反応表 (device 定数: $A,b,E_a$、反応物/生成物の化学量論、三体効率、falloff パラメータ、可逆フラグ)。
- 出力: `src_roY{s}` ($\dot\omega_s$ [kg/m³/s]) を `res_roY{s}` に加算、`src_roe` に $\dot Q$、
  Jacobian `src_jac_Y{s,k}` ($n_s\times n_s$, セル毎に連続配置) と $\partial\dot\omega_s/\partial T$。
- 数値: 濃度 $[X_k]=\rho Y_k/W_k$ は FP64、$\ln k_f$, $\ln K_c$ は FP32 で評価して $\exp$ を 1 回。
  $K_c$ は NASA-9 の $\Delta G^\circ$ から (絶対 datum の $a_7$ で評価 — sensible 焼き込み後の $a_7$ では $\Delta G$ が狂うため、
  **反転前の絶対係数を別に保持する**か $c_s$ で補正する)。
- ホスト単体テスト (`tools/test_chemistry.cpp`): `THERMO_HD` 共有コードで 0-D 定積/定圧反応器を積分し Cantera と比較 (§5)。

### 4.3 陰解法結合 (定常主経路)

1. **種ブロック point-implicit**: `species_dplur_sweep_d` の対角 $\texttt{transport\_diag}+\texttt{src\_jac}$ を
   $\big(\tfrac{V}{\Delta t}+D_{\rm conv}\big)I-\tfrac{\partial\dot\omega}{\partial(\rho Y)}$ の $n_s\times n_s$ 密ブロックに置き換え、
   セル毎 LU で $\delta(\rho Y)$ を得る。非対角 (隣接) は既存の scalar-DPLUR sweep のまま (Jacobi 反復)。
   Patankar 型の負対角 (消費項) が正値性を担保、生成項は陽側に置く。
2. **エネルギー結合**: 案C の予測 $\delta(\rho Y)^*$ に反応寄与が入るので、EOS-JVP $\delta p_Y$ に反応熱 $\dot Q\,\Delta t$ の分が自動で載る。
   更に $\partial\dot Q/\partial(\rho e)$ (温度応答) を flow block の対角 (5,5) に加える (発熱の自己増幅を陰化。既存 `src_jac_g` の潜熱項と同型)。
3. `speciesImplicitCoupling=0/1` (segregated / 対角) でも動くようにするが、生産推奨は 2 (案C) + 種ブロック。

### 4.4 陽解法・非定常 (副経路)

Strang 分離: RK ステージの前後で半ステップずつセル内 ODE を backward Euler (同じ Jacobian・$n_s\times n_s$ LU) で sub-cycle。
sub-cycle 数は $\Delta t/\tau_c$ から適応 ($\tau_c=1/\max|\lambda(J)|$ の Gershgorin 近似)。dual-time 非定常でも内側反復は 4.3 を使う。

### 4.5 入力配管

- `solverConfig.yaml`: `physProp.chemistry: {mechanismFile: h2_jachimowski.yaml, enabled: 1, tci: 0|1(PaSR), tMaxReaction: 6000, freezeBelowT: 0}`。
- 機構ファイル: **Cantera YAML のサブセット** (`reactions:` の `equation`, `rate-constant {A,b,Ea}`, `efficiencies`, `type: three-body|falloff`, `Troe`)
  をそのまま読む。単位は Cantera 既定 (cm³/mol/s, cal/mol) を内部 SI へ変換。これで Cantera との機構同一性が保証される。
- 熱力学 DB: `tools/cea_thermo_to_species_db.py thermo.inp --species H2 O2 H O OH HO2 H2O2 H2O N2 ... > species_db.yaml`
  (LJ パラメータは Cantera `transport:` 節または内蔵表から)。
- BC: `inlet_Pressure` の指定組成 `Yb` 配線が未 (§10 既知課題) → 燃焼室出口組成を与えるために **先に解消**する。

### 4.6 実現可能性・安全策

- $T$ 上限クランプ (反応熱の暴走), $Y_s$ 下限 0 → 再正規化 (既存), 反応の凍結温度 `freezeBelowT` (低温外挿域でのゴミ反応抑止),
  NASA-9 範囲外 (6000 K 超) は既存の線形外挿。
- 初期場: 燃焼器は「高温平衡組成の IC」から始める (着火は解かない)。ノズルは既存 warm-start (`interp_field.py`) + 平衡組成。

---

## 5. 段階計画と検証

| Phase | 内容 | 検証 (合否) | 成果物 |
| --- | --- | --- | --- |
| **0** 前提整備 (1 週) | CEA スクリーニング (case/44, 45 + 想定 $T_t$ 2000 K の平衡/凍結出口差)、`thermo.inp→species_db.yaml` ツール、機構 YAML パーサ、`inlet_Pressure` 組成 BC | 平衡/凍結の $\Delta M$, $\Delta T$ 表で「要否」を確定; radical 種の $c_p,h,s$ が CEA と一致 | `notes/`, `tools/` |
| **1** ソース項 (2 週) | `chemistrySource_d` (Arrhenius・三体・逆反応 $K_c$)、$\dot Q$、解析 Jacobian、ホスト 0-D テスト | 定積着火遅れ・定圧平衡到達が Cantera (同一 YAML) と $\pm2$ % / $\pm$ 数 K; Jacobian の有限差分照合 $<10^{-6}$ | `methods/chemistry.md`, plan |
| **2** 陰解法結合 (2–3 週) | 種ブロック point-implicit、案C 結合、Strang 副経路 | Q1D ノズル (H₂/O₂ 平衡入口 → 膨張) を forge 2D/軸対称 vs Cantera PFR (面積則) で $T$, $Y_{OH}$, 凍結位置; frozen / equilibrium / finite-rate の三者図; `cfl_pseudo` が非反応と同等 | case 新設 (`case/46.nozzle_h2o2_kinetics`) |
| **3** 燃焼器 RANS (3 週) | 低マッハ + SST + No-TCI → PaSR; NO を含む 13 種 | Burrows–Kurkov (壁噴射 H₂, M2.44 vitiated) の出口 $Y_{H_2O}$/$T$ 分布 (NASA GRC データ) と文献 CFD 帯域内 | case 新設 |
| **4** 応用 | 加熱器→ノズル→試験部 end-to-end、falloff (Burke 2012)、TCI 感度 | 試験部の残留 OH/NO と $M$ 偏差の定量 | 設計チェーン (`methods/design/`) へ反映 |

検証ツール: Cantera を `.venv-chem` に pip 導入 (CPU 参照解専用、ソルバ依存にはしない)。CEA は既存。

---

## 6. リスク・未決事項

1. **剛性 × 既存 CFL 上限**: 多成分 TP 陰解法の `cfl_pseudo` 上限 1–2 に化学剛性が重なると実用速度が出ない可能性。
   対策は種ブロック陰化 (4.3) と $\partial\dot Q/\partial(\rho e)$ の陰化。効果は Phase 2 で `cfl_pseudo` スイープで測る。
2. **FP32 桁落ち**: ラジカル $Y\sim10^{-6}$–$10^{-4}$ を `flow_float` (FP32) で輸送すると相対精度 $10^{-7}$ で足りるが、
   $\rho Y$ の増分 $\delta(\rho Y)$ が主値の $10^{-3}$ 以下なら丸めが効く。濃度・Jacobian は FP64 に留め、輸送量は FP32 のまま様子を見る。
3. **$K_c$ と datum**: sensible 焼き込み後の $a_7$ で $\Delta G$ を組むと $K_c$ が狂う。絶対係数の保持を Phase 1 の設計に入れる (4.2)。
4. **凍結温度域の熱力学**: 試験部 $T\sim250$ K は NASA-9 下限近傍。ラジカルは実質 0 になるので問題ないが `freezeBelowT` で明示的に止める。
5. **TCI の妥当性**: 加熱器は乱流火炎。No-TCI では火炎長が短くなる傾向 (文献)。Phase 3 で PaSR まで進めて感度を見る。
6. **GPU コスト**: 13 種 33 反応で $\exp$ 33 回 + 13×13 LU/セル/step。粘性流束と同程度と見込むが、FP64 律速なら FP32 化 (2.4)。
7. **case/44 (va3) では効果が小さい見込み**: 投資判断は Phase 0 の CEA 表で行う。高 $T_t$ 案件 (M6–7 加熱器) が無ければ Phase 2 で止める選択もある。

---

## 7. 次アクション

1. Phase 0 を開始: CEA で case/44 (Tt 1161 K) / case/45 (Tt 1600 K) / 仮想 Tt 2000・2500 K の `rocket equilibrium` vs `frozen` 出口 $M$, $T$, $Y_{OH}$ を表にする。
2. `plans/active/chemistry-finite-rate-h2.md` を `_template.md` から起草し、`methods/chemistry.md` の骨子 (4.1–4.4) を先に書く (AGENTS.md 開発フロー)。
3. 機構ファイルは Cantera 同梱 `h2o2.yaml` と Jachimowski 9 種 19 反応の YAML 化から着手。

---

## 9. Phase 0 結果 (2026-09-04)

方針 §5 の Phase 0 を実施した。数値は [`plans/active/chemistry-finite-rate-h2.md`](../../plans/active/chemistry-finite-rate-h2.md) §9 に転記。

- **CEA スクリーニング** ([`scripts/cea_kinetics_screen/`](scripts/cea_kinetics_screen/), 表 `screen_table.md`): 燃焼室 $T_c$ 1161 / 1600 / 2000 / 2400 K
  (va3 組成) と H₂ 加熱器組成 (H₂O 9–26 %) について平衡 vs 凍結 (燃焼室凍結 / スロート凍結) の出口 $M$, $T$ 差。
  **$T_c\le1600$ K で差 0.1 % 未満、2000 K で 0.3–0.4 %、2400 K で 1.1 %** (M4.2 出口; M6 出口は 1.5 倍)。
  → case/44 (1161 K)・case/45 (1600 K) は frozen 設計で十分。本機能の投資対象は $T_c\gtrsim2000$ K の加熱器案件。
- **熱力学 DB 生成ツール** `solver_density_cuda/tools/cea_thermo_to_species_db.py` (CEA `thermo.inp` → `species_db.yaml`)。内蔵 N2 と相対差 0。
- **機構ファイル** `solver_density_cuda/tools/mechanisms/` (Jachimowski 1988 13 種 33 反応 / 9 種 20 反応、CEA NASA-9 熱力学)。
  Cantera で読めることと着火遅れの妥当性 ([`scripts/chem_ignition_check.py`](scripts/chem_ignition_check.py)) を確認。
- **参照解環境**: `.venv-chem` (Cantera 3.2, pypdf)。ソルバ依存にはしない。
- 罠: 種名 `NO`/`N` は YAML 1.1 で真偽値になる (DB 生成は引用符付きに変更)。

方針文書 (`methods/chemistry.md`) と plan (`plans/active/chemistry-finite-rate-h2.md`) を起草済。次は Phase 0b (`inlet_Pressure` 組成 BC) と Phase 1 (ソース項)。

---

## 8. 参考文献・URL

- Jachimowski, C. J., "An Analytical Study of the Hydrogen-Air Reaction Mechanism with Application to Scramjet Combustion," NASA TP-2791 (1988). [NTRS](https://ntrs.nasa.gov/api/citations/19880006464/downloads/19880006464.pdf) (PDF 保存済)
- Zhukov, V. P., "Verification, Validation, and Testing of Kinetic Mechanisms of Hydrogen Combustion in Fluid-Dynamic Computations," ISRN Mech. Eng. (2012). [Wiley](https://onlinelibrary.wiley.com/doi/10.5402/2012/475607)
- Burke, Chaos, Ju, Dryer, Klippenstein, "Comprehensive H2/O2 kinetic model for high-pressure combustion," Int. J. Chem. Kinet. 44 (2012).
- Ó Conaire et al., "A comprehensive modeling study of hydrogen oxidation," Int. J. Chem. Kinet. 36 (2004).
- Bussing & Murman, "Finite-volume method for the calculation of compressible chemically reacting flows," AIAA J. 26(4) 1070 (1988). [ADS](https://ui.adsabs.harvard.edu/abs/1988AIAAJ..26.1070B/abstract)
- Eberhardt & Imlay, "Diagonal implicit scheme for computing flows with finite rate chemistry," J. Thermophys. Heat Transfer 6, 208 (1992).
- Candler, Subbareddy, Nompelis, "Decoupled Implicit Method for Aerothermodynamics and Reacting Flows," AIAA J. (2013). [AIAA](https://arc.aiaa.org/doi/10.2514/1.J052070)
- Candler et al., "Development of the US3D Code for Advanced Compressible and Reacting Flow Simulations," NTRS 20180002133. [NTRS](https://ntrs.nasa.gov/api/citations/20180002133/downloads/20180002133.pdf) (PDF 保存済)
- Niemeyer, Curtis, Sung, "pyJac: Analytical Jacobian generator for chemical kinetics," Comput. Phys. Commun. 215 (2017). [arXiv](https://arxiv.org/pdf/1605.03262)
- Perini, Galligani, Reitz, "An Analytical Jacobian Approach to Sparse Reaction Kinetics...," Energy & Fuels (2012). [ACS](https://pubs.acs.org/doi/10.1021/ef300747n)
- Curtis, Niemeyer, Sung, "An investigation of GPU-based stiff chemical kinetics integration methods," Combust. Flame (2017). [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0010218017300354)
- "A GPU-accelerated stiff ODE solver for combustion CFD with fully analytical chemical kinetics..." (preJacGPUFoam), Fuel (2025). [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0016236125033538)
- ChemGen: Code Generation for Multispecies Chemical Kinetics. [arXiv 2510.10005](https://arxiv.org/html/2510.10005v1)
- "An AMReX-based Compressible Reacting Flow Solver for High-speed Reacting Flows relevant to Hypersonic Propulsion." [arXiv 2412.00900](https://arxiv.org/html/2412.00900v2)
- "A component-splitting implicit time integration for multicomponent reacting flows simulations." [arXiv 2403.03440](https://arxiv.org/pdf/2403.03440)
- Maier et al., "SU2-NEMO: An Open-Source Framework for High-Mach Nonequilibrium Multi-Species Flows," Aerospace 8(7) 193 (2021). [MDPI](https://www.mdpi.com/2226-4310/8/7/193)
- Gibbons, Damm, Gollan, Jacobs, "Eilmer: An open-source multi-physics hypersonic flow solver," CPC (2023). [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0010465522002703)
- NASA NTRS 19650044133, "Evaluation of the Bray sudden-freezing criterion for predicting nonequilibrium performance in multireaction rocket nozzle expansions." [NTRS](https://ntrs.nasa.gov/citations/19650044133)
- "Chemical non-equilibrium flow analysis of H2 fueled scramjet nozzle," Case Studies in Thermal Eng. (2015). [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2214157X15000052)
- "Numerical study of a nonequilibrium H2–O2 rocket nozzle flow." [ResearchGate](https://www.researchgate.net/publication/331408649_Numerical_study_of_a_nonequilibrium_H2_-_O2_rocket_nozzle_flow)
- Pellett et al., "Review of Air Vitiation Effects on Scramjet Ignition and Flameholding Combustion Processes," AIAA 2002-3880. [NTRS](https://ntrs.nasa.gov/api/citations/20030010281/downloads/20030010281.pdf)
- "Test Capabilities and Recent Experiences in the NASA Langley 8-Foot High Temperature Tunnel." [NTRS](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20000058170.pdf)
- Burrows & Kurkov, NASA TM X-2828 (1973); NASA GRC Wind validation archive. [GRC](https://www.grc.nasa.gov/www/wind/valid/bk/bk.html)
- "Mixing time scale analysis of the Partially Stirred Reactor model for high-speed turbulent combustion of hydrogen in vitiated air," Acta Astronautica 218 (2024). [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0094576524000675)
- "Turbulence–chemistry interaction models with finite-rate chemistry and compressibility correction for simulation of supersonic turbulent combustion," Eng. Appl. CFD (2020). [T&F](https://www.tandfonline.com/doi/full/10.1080/19942060.2020.1842248)
- Cantera reactor / PFR examples (参照解用). [Cantera](https://cantera.org/3.1/examples/python/reactors/pfr.html)
