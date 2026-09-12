# 乱流‐化学相互作用 (TCI): 1st-order CMC の導入

## メタ

- **area**: `thermophysics / chemistry / turbulence`
- **status**: `in_progress` (2026-09-05 起票・設計確定、Phase A 着手)
- **related_docs**:
  - `methods/chemistry_cmc.md` (現在仕様: 理論・実装設計)
  - `methods/chemistry.md` (有限速度化学の基盤、§5 TCI)
  - `notes/investigations/cabra-liftoff-model-fidelity-survey.md` (方式選定の文献根拠)
- **related_plans**:
  - `plans/active/chemistry-finite-rate-h2.md` (親。§5.1 P0-1/P0-4 の決定で本 plan を起票)
- **created**: `2026-09-05`
- **owner**: `CFD Dev`

## 1. 目的

laminar chemistry (セル平均で Arrhenius を評価) では自着火安定化火炎の着火遅れを表現できず、Cabra (case/48) で
リップ付着・BK (case/47) で早着火になる (機構交換・PaSR では直らない — 親 plan §9 2026-09-05)。
混合分率に条件付けた組成分布で反応を評価する 1st-order CMC を導入し、**Cabra の浮き上がり高さ H/d≈10 とその
コフロー温度応答、BK の着火位置**を再現できる状態にする。加熱器設計 (火炎位置・壁熱負荷まで) の前提。

## 2. スコープ

- **やる**: Bilger 混合分率の診断、分散輸送 + $\tilde\chi$ (SST)、$\eta$ 格子上の条件付きスカラー ($n_s$ 種 + $h$) の輸送
  (物理空間 1 次風上 + 乱流拡散、$\eta$ 拡散陰的 (AMC)、化学は各 $\eta$ 点で点陰解)、β-PDF 積分によるソース結合
  (`tci: 2`)、Cabra $H(T_c)$ 応答曲線と BK 着火位置での検証。定常 (擬似時間) と dual-time の両方で動く。
- **やらない**: 2nd-order CMC (条件付き分散)、密度重み付き CMC 形、多流 (3 流以上) 混合分率、LES-CMC、
  transported PDF (粒子/stochastic fields)、条件付き速度の高次モデル、$\tilde Y$ を PDF 積分で上書きする full coupling。

## 3. 関連 docs と前提

- 理論・設計: `methods/chemistry_cmc.md` (本 plan と同時に起草)。
- 既存資産: `scalarTransport_d` (1 変数の移流+拡散残差と時間積分)、`chemistry_d` (`chem_source`, 種ブロック LU,
  `jacobianInterval`)、`ransTransport_d` (SST の $k,\omega$ → $\varepsilon/k=\beta^*\omega$)、`thermo_d` (NASA-9, $h\leftrightarrow T$)。
- 前提: 全化学種の拡散係数を同一 ($Sc_t$ 共通) とする現行仕様のもとで元素質量分率は保存スカラー → $\tilde\xi$ は診断量。
- 制約: 条件付きスカラーは node 当たり $(n_s+1)N_\eta$ 個 (9 種 × 41 = 410 double)。Cabra 60k ノードで 25M 値 (200 MB) は可。
  化学評価は node 数 × $N_\eta$ 倍になる ($\approx40\times$) ので `jacobianInterval` と「$\tilde P(\eta_k)$ が閾値未満の $\eta$ 点は
  化学をスキップ」で抑える。

## 4. 設計方針 (差分のみ。本文は `methods/chemistry_cmc.md`)

1. **source coupling**: 流れ側は $\tilde Y,\tilde h$ を従来どおり輸送し、`chemistry_source_d` の $\dot\omega(\tilde Y,\tilde T)$ を
   $\sum_k w_k\tilde P(\eta_k)\dot\omega(Q(\eta_k))$ に置換する (`tci: 2`)。既存の陰解法配管を変えない。
2. **$\tilde\xi$ は診断** (Bilger、元素組成は機構 YAML の `composition`)、**分散だけ輸送** (`scalarTransport_d` 汎用拡散、
   生成 $2\mu_t/Sc_t|\nabla\tilde\xi|^2$、散逸 $\bar\rho\tilde\chi$ は `src_jac` で点陰解)。
3. **$\eta$ 格子**: $N_\eta$=41 (config)、$\xi_{MR}$ と $\xi_{st}$ 付近を tanh で密に。$\eta$ 端は Dirichlet (燃料/酸化剤流の状態)。
4. **$\eta$ 拡散 + 化学の陰解**: node 毎に「三重対角 ($\chi$ 項) + 各 $\eta$ 点の種ブロック LU」を分離反復 (先に $\eta$ 拡散を陰的に、
   次に化学を点陰解)。分離誤差は擬似時間収束で消える。
5. **PDF**: β 分布を $\eta$ 格子上で解析正規化して離散化 (台形則)。$\widetilde{\xi''^2}<\epsilon$ はデルタ扱い。
6. **AMC** で $\langle\chi|\eta\rangle$。$C_\chi=2$ (config)。
7. **精度**: $Q$・化学は double、輸送残差は flow_float。
8. **出力**: `xi`, `xiVar`, `chi`, `cmc_dY`。$Q_T(\eta)$ は指定 node 群のみ (診断)。

## 5. 実装ステップ

| Phase | 内容 | 主要ファイル | 状態 |
| --- | --- | --- | --- |
| **A** 混合分率インフラ | 元素組成読込 (`composition`)、Bilger $\tilde\xi$ 診断カーネル、分散 `roXiVar` の登録・輸送・ソース、$\tilde\chi$、出力。Cabra 混合場で検証 | `chemistry_mech_io.hpp`, `chemistry_d.cuh`, `cuda_forge/cmc_d.{cuh,cu}` (新規), `variables.cpp`, `main.cpp`, `input/solverConfig.*` | **実装済 (2026-09-05)**, 検証中 (§9) |
| **B** 条件付きスカラー | `Q` ストレージ・$\eta$ 格子・初期化 (混合線)・物理空間輸送 (scalarTransport ループ)・$\eta$ 拡散陰解 (AMC)・各 $\eta$ 点の化学点陰解・境界条件 | `cmc_d.{cuh,cu}`, `scalarTransport_d`, `chemistry_d` | **実装済 (2026-09-05)**, 凍結検証 PASS (§9 (3)) |
| **C** ソース結合 | β-PDF 積分で $\bar{\dot\omega},\bar{\dot Q}$, ブロック Jacobian を置換 (`cmc.couple: 1`)、`cmc_dY` 診断、PDF 閾値スキップ | `chemistry_d.cu` (`cmcOmega/cmcQdot/cmcJac` 引数) | **実装済 (2026-09-05)**, 反応 ON 検証中 (`run_0082`) |
| **D** 検証 | 凍結混合整合 → Cabra $H(T_c)$ 応答曲線 (5 点) → BK 着火位置 → 回帰 (`tci 0/1` ビット不変) | `case/48`, `case/47`, `tests/` | todo |

### 5.1 残作業表 (優先順)

| 優先 | 項目 | 状態・次アクション |
| --- | --- | --- |
| 1 | Phase A: 元素組成読込 + Bilger $\tilde\xi$ 診断 | **done 2026-09-05** (カーネルは numpy 再計算と 3e-8 一致; Y_H₂ 由来 ξ との差 ≤0.068 は差分拡散の物理) |
| 2 | Phase A: 分散輸送 + $\tilde\chi$ + 出力、Cabra 混合場で分散の妥当性 | 実装済・300 step smoke OK (実現可能性 1e-12, χ 10²–10³ 1/s)。発達場 `run_0080` (5000 step) で確認中 |
| 3 | Phase B: $Q$ ストレージ・$\eta$ 格子・混合線初期化・凍結整合 | **done 2026-09-05**: `run_0081` で条件付き T が全 node 1045.00 K (混合線保存)、`cmc_dY` ≤0.017 (差分拡散分のみ) |
| 4 | Phase B: $\eta$ 拡散陰解 (AMC) + 化学点陰解、境界条件 | 実装済 (凍結で η 拡散の線形不変を確認; 化学は `run_0082` で検証中) |
| 5 | Phase C: 結合 | couple 1 (PDF 平均ソース) は平均場が燃えず、couple 2 (緩和ソース) は差分拡散差を反応熱に換算して NaN、couple 3 全置換は陰解法と非整合で NaN → couple 3 (Y,h ブレンド) はリップの熱伝導/Le≠1 で非物理 → **couple 4 (Y ブレンド + ブレンドの反応熱) を検証中 (`run_0089_cmc_c4`)**。性能: 化学 ON 652 ms/step (単独) → **2026-09-05 (12) で CMC 部分 ~175 ms/step に短縮** (fp32 化学ほか; 残りは輸送残差の融合) |
| 6 | Phase D: Cabra $T_c$ 1015/1030/1045/1060/1075 K の $H$ 応答曲線 vs 実験・文献、中心軸/半径 T ±30 K | **1045 K は達成** (run_0107: 半径 T 差 +20〜+42 K, 基部 x/d 8.4〜9.7 vs 実験 ≈10)。run_0110 (1030 K) は ξ シンク (§9 (17)) で無効。**修正後スイープを AWS で実行中: `run_0112`–`run_0116` (1045/1030/1060/1015/1075 K, 各 3000 step, 直列)**。旧メモ: 残り 1060/1015/1075 K は run_0110 完了後に 3000 step で投入 (`cd case/48.cabra_h2n2 && DT=1e-5 nohup bash run_tc_sweep.sh ../run_0107_dt_c7_from0077/res_2000.h5 ../run_0107_dt_c7_from0077/cmc_Q_2000.bin 111 3000 "1060 1015 1075" > tc_sweep_2026-09-06b.log 2>&1 &`; GPU を 2 本以上共用すると 0.2 step/s まで落ちるので直列で)。基準 1045 K は `run_0108` (40 ms 継続, 基部 x/d 8.1〜8.4 で定常) |
| 7 | Phase D: BK 着火位置、`tci 0/1` 回帰 | todo |
| 7b | 性能: 物理空間輸送残差 (410 スライス × 2 カーネル, 41 ms/step) を面幾何共有の融合カーネルに。PDF 閾値スキップは効果小 (Ω<1e-4 は活性対の 17 %) | todo (低優先) |
| 8 | ドキュメント: `procedures/solver-settings.md` に `cmc` キー、`methods/index.md` | todo |

## 6. 検証

- **単体 / ビルド**: Bilger $\xi$ が純燃料/純酸化剤で 1/0、混合線上で線形 (ホストテスト)。β-PDF の正規化・モーメント
  ($\int P=1$, $\int\eta P=\tilde\xi$, $\int(\eta-\tilde\xi)^2P=\widetilde{\xi''^2}$) を $\eta$ 格子上で 1e-3 以内。
- **検証ケース**:
  1. 凍結混合 (case/48 化学 OFF): $Q_\alpha(\eta)$ = 混合線、`cmc_dY` < 1 %。
  2. Cabra H₂/N₂ (case/48): $H(T_c)$ 応答曲線 (5 点)。合格: 1045 K で H/d が 10 ± 3、傾きが実験 (10 K で約 2 倍) と同符号・同桁、
     $T_c\pm30$ K のどこかで H/d=10 を通る。中心軸・半径 T の平均差 ±30 K。判定は `check_quasisteady --quantity ignx` STEADY の場で。
  3. Burrows–Kurkov (case/47): 着火位置 18–25 cm (現状 4–6 cm)。
  4. 回帰: `tci: 0` / `tci: 1` の既存 run が ビット不変。
- **判定基準**: 上記。加えて質量・元素・全エンタルピー収支 (`check_quasisteady` の exit 量) が反応 ON でも閉じること。

## 7. 影響範囲

- 新規: `cuda_forge/cmc_d.{cuh,cu}`, `tools/test_cmc.cpp` (ホストテスト), `methods/chemistry_cmc.md`。
- 変更: `chemistry_mech_io.hpp` (composition), `chemistry_d.{cuh,cu}` (`tci: 2` 分岐), `variables.cpp` (`xi`, `roXiVar`, `Q`),
  `main.cpp` (呼び出し), `input/solverConfig.*` (`cmc` キー), `procedures/solver-settings.md`。
- 既定 (`tci: 0/1`) では全経路ビット不変。

## 8. 完了条件

- [x] `methods/chemistry_cmc.md` 起草
- [ ] Phase A–C 実装
- [ ] Phase D 検証 (§6) — Cabra 応答曲線と BK 着火位置
- [ ] `status: done` + §9 変更ログ、`plans/active/` → `plans/accepted/`、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-05` — 初稿。ユーザ決定「分布で考える TCI を入れて検証する」(親 plan §5.1 P0-1/P0-4) を受け、方式を 1st-order CMC
  (source coupling, Bilger 診断 + 分散輸送, AMC, β-PDF) に確定。根拠: 文献で Cabra $H(T_c)$ を再現した実績が最も強いのが
  transported PDF / CMC で、決定論的な CMC は定常 (擬似時間) ソルバと GPU の既存スカラー輸送・種ブロック陰解に乗る
  (stochastic fields は擬似時間定常と相性が悪い)。Phase A に着手。
- `2026-09-05 (2)` — **Phase A 実装** (commit 本): `ReactionTable` に元素組成 (`nElem/elemW/atoms/elemH,O,C`)、`chemistry_mech_io` が `species[].composition` を読む
  (`chemistry_init` は mixfrac のみでも機構を読み、反応ソースは `chemEnabled` のみ)、config `physProp.chemistry.mixfrac {enabled, cChi, fuelX, oxidizerX}`、
  `variables::registerMixfrac` (`xi, chi, roXiVar/xiVar/N/M/res/res_m/src_jac_xiVar/transport_diag_xiVar, dXi*, dXiVar*`; リスタート復元)、
  新規 `cuda_forge/cmc_d.{cuh,cu}` (Bilger カーネル、1 スカラ GG 勾配、分散ソース (生成/散逸=src_jac)、原始量+実現可能性クランプ、境界、N/M)、
  `main.cpp` に凝縮モーメントと同じ 8 箇所の呼び出し (init/restart/assembleResidual/implicit inner/explicit RK/dual-time)。
  **落とし穴**: (i) `cuda_forge/CMakeLists.txt` は明示リストなので新規 .cu の追加が要る、(ii) 変数名は `condMomentCellVarNames` の
  「先頭 2 文字を落とす」規約で `XiVar` になるため明示登録した (不一致は `c_d[]` が空エントリを作り GPU 不正アクセス→ `.at()` に統一)。
  **smoke `case/48 run_0079_mixfrac_smoke`** (run_0067 混合場, 化学 OFF, 300 step): NaN なし、`xi` は numpy 再計算と 3e-8 一致、
  `|xi − Y_H2/Y_H2,fuel|` 最大 0.068 はノズル出口リップの層流層 (差分拡散、`|ξ_H2−ξ_N2|` 0.148) で物理、分散の実現可能性違反 9e-13、
  χ̃ はせん断層で 10²–10³ 1/s、sqrt(var)/ξ ≈ 0.35–0.41 (文献 RMS/mean と同桁)。発達場は `run_0080` (5000 step)。
- `2026-09-05 (3)` — **Phase B/C 実装**: `cmc_d.cu` に Q(η) の packed ストレージ ([slice][node], flow_float ×6 = 100 MB×6 @ Cabra 60k)、
  混合線初期化 (流れの顕エンタルピーは device thermo で評価)、各スライスを `scalarTransport_d` に渡す物理空間輸送、
  node 毎カーネル `cmc_step_d` (β-PDF 重み Ω_k: 端ビン解析処理 + 数値正規化、デルタは隣接 2 点配分; AMC の erfinv は host Newton で表化;
  η 拡散は Thomas 陰解; 各 η 点で `chem_source` (jacMode 2) + 9×9 Gauss 消去の点陰解; PDF 平均 ω̄/Q̇̄/J̄ と Ỹ_pdf, T_Q max)、
  `chemistry_source_d` に `cmcOmega/cmcQdot/cmcJac` を渡して平均場の Arrhenius を置換 (`cmc.couple: 1`)。config `physProp.chemistry.cmc`
  (`nEta, etaPow, pdfFloor, couple, chem, fuelT, oxidizerT, dtScale, interval`; timeIntegration 11 のみ)。
  **凍結検証 `run_0081_cmc_frozen`** (化学 OFF, couple 0): 初版は条件付き T が 998–1087 K と混合線上限 1045 K を超えた →
  原因は保存形輸送のドリフト (未収束の平均流では Σṁ_f≠0 で一様 Q でも −Q Σṁ_f が残り、流れ側の Δρ と点陰解の Δρ が一致しないため
  Q が (ρ+Δρ_s)/(ρ+Δρ_f) 倍ずつ動く)。`res_Q -= Q·res_ro` の非保存補正で **全 node 1045.00 K (厳密)** に。ρ̄Q の再同期も追加。
  `cmc_dY` ≤0.017 は出口リップの差分拡散分 (Phase A と同じ)。反応 ON・結合 ON は `run_0082` (run_0080 混合場から 6000 step)。
- `2026-09-05 (4)` — **反応 ON 初回 `run_0082` (couple 1 = source coupling) は結合方式の欠陥を露呈**: 条件付き空間は着火 (T_Q max 1546→1615 K,
  前線 x/d 18.8→6.4 と上流へ) するが、平均場は T_max 1045→1123 K しか上がらず (軸 z/d 26: 840 K, 実験 1410)、`cmc_dY` が 0.07 で定常化。
  条件付き状態が燃え切ると ω(Q)→0 で、未燃の平均場に渡る ω̄ が消える (両者の時間履歴が違うので積分値が一致しない)。
  → **couple 2 (PDF 積分状態への緩和) を実装**: ω̄_s = ρ̄(Ỹ_pdf,s − Ỹ_s)/τ, Q̇̄ = ρ̄(h̃_pdf − h̃)/τ, τ = relax·Δτ_local, J = −I/τ (点陰解)。
  Σω̄=0 で質量保存、定常で Ỹ → Ỹ_pdf (文献 RANS-CMC と同じ)。診断 `cmc_TQst` (η≈ξ_st の条件付き T) を追加。検証 `run_0083`。
- `2026-09-05 (5)` — **結合方式の試行 (run_0083 系)**: couple 2 (緩和ソース ρ̄(Ỹ_pdf−Ỹ)/τ) は (i) 初回残差時に Ỹ_pdf が未計算 (=0) で
  全種剥離 → 初期化時に PDF 積分を実行、(ii) 独立なエンタルピー緩和 ρ̄(h̃_pdf−h̃)/τ が 案C 予測子の ∂Q̇/∂ρY と非整合 → 反応熱を標準規約
  −Σc_sω̄_s に、(iii) それでも差分拡散由来の組成差 (リップ 0.017) を τ=Δτ で反応熱に換算し Q̇ ±1e11 W/m³ で step 2 NaN → **couple 2 は不採用**。
  couple 3 (平均 Y,T を PDF 積分値で毎ステップ全置換) は密度ベース陰解法と整合せず step 14 で NaN (dt 1.5e-9)。
  **couple 3 + α ブレンド (Y ← Y + α(Y_pdf−Y), h も同様, α=0.05/step, 平均方程式の化学ソースは 0)** で NaN なく進行 (`run_0083_cmc_react_couple3`)。
  併せて **リスタートで `roXiVar` が復元されていなかった** (host に読むだけで H2D 転送なし + 凝縮分岐の内側に置いていた) を修正 —
  run_0081/0082 は分散 0 (δ-PDF) から始まっていたので結果は再解釈が必要。別セッションが `run_0079_react_li_dualtime_cont` を作成しており
  `run_0079_mixfrac_smoke` と番号衝突 (README に注記)。GPU 共有で 852 ms/step。
- `2026-09-05 (6)` — **couple 3 の往復バグと差分拡散の非整合**: (i) 上書きカーネルが流れ更新後の保存量 ρ と更新前の原始量 T, u を混ぜて
  roe を組み直していた → α=0 の往復だけで発散 (`run_0084_cmc_ab_a0`)。保存量のみから T (thermo_T_from_e)・KE を再構成して受動 run と一致
  (残差 3 桁一致)。(ii) それでも α=0.05 で rms_roe が 25 step で 300 倍成長: 診断 `run_0085` でノズルリップ直後 (x 0.2 mm, r 2.4–2.6 mm) の
  平均 T が 884→510 K と引き下げられ P 124 kPa まで励起。**差分拡散 (`speciesDiffusionMethod: 1`) の平均場が Bilger ξ の混合線から
  130 K 以上ずれ**、ブレンドと平均場の差分拡散が綱引きになる。→ CMC 計算は `speciesDiffusionMethod: 0` を前提にする (methods §4)。
  混合場を定数 Sc で作り直し (`run_0086_mix_sc0`)、そこから CMC couple 3 (`run_0087_cmc_sc0`)。GPU は他 2 セッションと共有中。
- `2026-09-05 (7)` — **等拡散でも couple 3 は step 322 で NaN** (`run_0087_cmc_sc0`): 未着火のまま、リップ (x 0.4 mm, r 2.7 mm, ξ 0.31) の
  平均 T が 1034→400 K に引き下げられ P 95 MPa。ここは**管壁の熱伝導と Le≠1 で温度が混合線 T(ξ) から切り離される**場所で、
  エンタルピーのブレンドが原理的に非物理。→ **couple 4: 組成だけを PDF 積分状態へ α ブレンドし、その組成変化の反応熱 −Σ c_s Δ(ρY_s) を
  エネルギーに加える** (平均エネルギー式が伝導・境界を担い、CMC は反応進行だけを渡す)。400 step 診断 (`run_0088`) で NaN なし、
  P 100–102 kPa、rms_roe 横ばい、条件付き空間で着火開始 (T_Q 1356 K)。本番 `run_0089_cmc_c4` (6000 step)、受動対照 `run_0090`。
- `2026-09-05 (8)` — **重大バグ 2 件を Q(η) プローブ (`CMC_PROBE_NODE`) で特定**: (i) `cmc_step_d` の `doChem` 引数を初期化呼び出しと通常更新で
  **逆に**渡していた (初期化時に未設定 dt_local で η 拡散が走り内点 Q が 0/NaN に、通常更新では η 拡散も化学も動かず) — run_0083〜0091 の
  「反応 ON」結果はすべて無効 (run_0088 の T_Q 1356 K 等)。fmax が NaN を無視するため `cmc_TQmax` が端点値 1045 を返して偽合格になっていた
  (NaN を伝播させる判定に変更)。(ii) couple 4 の反応熱が混合線に沿った組成変化 (coflow の H2O) まで生成エンタルピー換算していた → 混合線からの
  離脱分だけに限定。修正後 `run_0091_cmc_c4_diag2` (200 step): ξ_pdf=ξ (差 0.003)、P 101–102 kPa、NaN なし、リップの受動対照との差 −89 K
  (初期組成補正の一過性、Qdot 1e9→3e6 と収束)。本番 `run_0093_cmc_c4`。
- `2026-09-05 (9)` — **couple 4 の熱源を確定**: 「混合線離脱分の反応熱」でも `run_0093` は平均 T 3477 K (T_Q 1392 K) で NaN —
  平均場の数値的な混合線からのずれ (種ごとのリミッタ差) を毎ステップ補正した分まで反応熱に計上していた。「D_pdf の時間差分」は PDF 重みの
  移動で ±2e9 W/m³ の偽吸発熱。**最終形: 条件付き点陰解で実際に進んだ組成変化の反応熱 q_rel = Σ_k Ω_k (−Σ_s c_s δY_s,k) [J/kg] を
  PDF 平均して ρ̄ q_rel を平均エネルギーに加える** (組成のブレンドは熱を伴わない)。600 step 診断で熱源が非負、P 101→110 kPa の緩やかな
  上昇、条件付き着火 x/d 39–41 (希薄側) → 平均場 OH>2e-4 が x/d 38 に出現。本番 `run_0094_cmc_c4` (8000 step)。
- `2026-09-05 (10)` — `run_0094` (最終形の熱源) は step 1000 で平均 T_max 1506 K・平均 OH>2e-4 が x/d 15.9 (条件付き着火は x/d 24.8) と
  燃え始めたが、**P が 96–186 kPa に振れて step 1291 で NaN**: 擬似時間では条件付き化学が数百 step で燃え切り、流れが音響的に追従できない速さで
  熱が入る。→ **1 step の発熱を ΔT ≤ `cmc.dTmax` (既定 10 K) に制限し残りを持ち越す** (総発熱は保存)。`run_0095` は dTmax 10 + dtScale 0.5。
- `2026-09-05 (11)` — **結合方式の決着: couple 5 採用**。A/B: couple 4 (発熱を roe に直接加算) は上限 2 K/step・dtScale 0.25 でも step 4300 で
  NaN (`run_0096`)。**couple 5 (同じ発熱を陰的化学ソース経路 `res_roe` に Q̇ として渡す, 上限 10 K/step, dtScale 0.5) は 6000 step 完走・
  圧力 100–102 kPa 安定 (`run_0097_cmc_c5`)**。平均場の火炎基部 (Y_OH>2e-4) は x/d 22.7 (2k) → 18.4 → 16.6 → 14.7 → 13.8 (6k) と減速しつつ
  上流へ移動、条件付き着火前線 (T_Q>1300 K) も 23.6 → 14.1 と追従。laminar 化学の付着 (0) とは異なり、実験 H/d≈10 に向かって収束しつつある形。
  平均 T_max は 1063 K と低く (軸 z/d 26: 883 K) 平均場の加熱は遅い (上限持ち越し + 陰的経路の減衰) → 継続計算で確認。
  Q(η) の restart 永続化 (`cmc_Q.bin`, `cmc.restartQ`) を実装。
- `2026-09-05 (12)` — **性能: CMC 部分 613→133 (+41) ms/step** (methods/chemistry_cmc.md 実装 §4 に表)。`CMC_TIMING` 区間計測と ncu で
  「node 毎 1 カーネルの化学が FP64 パイプ 87 % で律速 (RTX 3060 は fp64 = fp32 の 1/64)」と特定。対策の順: (node,η) 分割 → 活性対の圧縮リスト
  (`thrust::copy_if`) → T 反転のウォームスタート → **fp32 化学 `chemistry_f32_d.cuh` (`cmc.fp32: 1` 既定; 点陰解は double)** → η 拡散 float +
  時間積分の一括化。各段階で 40 step A/B (run_tmp_prof_* , 削除済) が double 参照と相対 1e-6〜1e-5 で一致 (OH の 4e-4 は Y~1e-13 のノイズ)。
  couple 1/2 の ω̄/J̄ 集計は double atomicAdd (順序非決定) だが生産の couple 5 は使わない。長期 run_0098 は旧バイナリで続行中。
- `2026-09-06 (13)` — **長期 run と結合方式の再判断**: `run_0098` (couple 5, double, 20000 step) と `run_0099` (同 step 4000 から fp32 で継続,
  Q restart 初使用) は同一 step で全量 ~1 K 一致 → **fp32 化学と `cmc.restartQ` を長期で検証**。結果: 火炎基部 (Y_OH>2e-4) は x/d 7.8〜14.5 を
  coflow 柱モード (出口 ṁ ±97 %) と共に往復 (`check_quasisteady` DRIFTING)、平均 T_max ~1150 K、軸 z/d 26 は 1075 K (実験 1410)。
  新診断 `cmc_pdf_mean.py` で **PDF 積分の条件付き T は実験と整合** (z/d 9/11/14 の半径分布はほぼ重なる、軸 z/d 26 1328 K) と判明し、
  平均場の不足は couple 5 の構造欠陥 (熱は燃えた node の PDF 重みでしか渡らず、上流で燃えて運ばれた Q の熱が届かない) と特定。
  → **couple 6** (Y ブレンド + 燃焼領域ゲート付き h̃→h̃_pdf 緩和を陰的経路で; `cmc.dTgate`) を実装、`run_0101` (run_0099 の 16000 step 状態から
  8000 step) で couple 5 継続 `run_0100` と A/B 中。setup の `--cmcq` 文字列の波括弧バグ (yaml "illegal flow end") で両 run が一度起動失敗
  → 修正済 (dry run は `yaml.safe_load` で検証すること)。
- `2026-09-06 (14)` — **couple 6 成功 + 定常反復の柱モードは無害でないと判明**: `run_0101` (couple 6, 8000 step) で平均 T が PDF 診断 T に数 K で一致し
  (z/d 9/11/14 の半径分布は実験とほぼ重なる、火炎基部 x/d 9.4〜10.9 vs 実験 ≈10)、**couple 6 を採用**。一方、軸 ξ が z/d 20→30 で 0.51→0.02 に
  崩落し z/d ≥ 26 で T が実験より −220 K。原因調査 (`xi_flux_check.py`, 全 H/O 元素の Bilger ξ の断面フラックス ∫ρuξdA): 定常反復の run
  (`run_0067` 混合, `run_0099`, `run_0101`) では **ξ フラックスが z/d 5 で +60 %・z/d 30 で −70 %、全質量流量も +30〜45 %** と保存されず、
  dual-time の `run_0077/0079` では ±8 % で保存され軸 ξ も実験に近い (z/d 20: 0.57 vs 0.62, z/d 30: 0.32 vs 0.41)。すなわち coflow 柱の
  擬似時間リミットサイクル (親 plan の「定常 vs dual-time は近傍場で ≤20 K」) は**下流の混合・巻き込みを壊す**ので、CMC 検証と T_c スイープは
  dual-time で行う。ξ_pdf と ξ の偏り (ブレンドによる ξ の漏れ) は <0.003 で無罪。setup に `--dualtime <dt>` を追加、`run_0102` (couple 6,
  dt 1e-5 s, 4000 step) を run_0101 の状態から投入。`run_0100` (couple 5 継続) は step 10000 で中断 (不要)。
- `2026-09-06 (15)` — **dual-time で couple 6 は不可 → couple 7 を実装**: `run_0102` (couple 6, dual-time, run_0101 の状態から) は 10 ms で
  条件付き場が消炎・平均 ξ が z/d 5 以降で 0 に崩落・ξ=0 の coflow が 770 K に過冷却。切り分け: 純混合 dual-time (`run_0103/0104`, sdm 0/1 はビット同一)
  と CMC 受動 (`run_0105`, couple 0) は ξ フラックス保存 (1.51〜1.60 g/s) で無罪 → 原因は結合。`run_0106` (couple 6, run_0077 の燃焼場から) で
  質量流量が 172〜233 g/s (真値 37)・出口逆流 −63 g/s・P 97〜110 kPa: **擬似反復ごとの α ブレンド (10 K/反復 × 5 反復/物理 step) は
  物理時間で 5e6 K/s の加熱/冷却になり音響過渡で流れを壊す** (定常反復は音響を解かないので無害だった)。couple 6 の熱持ち越し qdebt も
  目標到達後の過冷却を生むので廃止。→ **couple 7**: 物理緩和時間 τ_c (`cmc.tauRelax` 1e-4 s) の通常ソース ω_s=ρ(Ỹ_pdf−Y)/τ_c,
  Q̇=ρ g (h_pdf−h)/τ_c (Jacobian 対角 −1/τ_c, ∂Q̇/∂ρe=−g c_p/(c_v τ_c)) を `chemistry_source_d` の通常経路に追加 (ゲート g は cmc 側で配列化)。
  `run_0107` (run_0077 から 2000 step = 20 ms) で ξ フラックス保存と着火/火炎基部を確認中。
- `2026-09-06 (16)` — **couple 7 + dual-time で Cabra が合った** (`run_0107`, run_0077 の燃焼場から 20 ms): 圧力・ξ フラックスとも健全、条件付き
  火炎基部 (T_Q>1300 K) が x/d 22→8.4 に落ち着き (平均 T>1200 K は x/d 9.7; 実験 H/d≈10)、半径 T の実験差は z/d 9/11/14/26 で
  +27/+32/+20/+42 K (読取誤差 ±30 K 級)、軸 T は z/d 26 で 1466 vs 1410 K。PDF 診断 T と平均 T は数 K で一致 (結合が整合)。
  残る差は z/d 26 の半径方向の肩 (r 10〜20 mm で +100 K; RANS 丸噴流の拡がり)。欠陥 1 件: `chemistry_source_d` の Tfreeze 早期 return で
  冷たい噴流コアに緩和ソースが入らず初期場の OH が残る → couple 7 は Tfreeze を通すよう修正。**couple 7 + dual-time を生産構成に確定**
  (定常反復は柱モードで不可、couple 5/6 は dual-time で音響過渡)。`run_0108` (続き 40 ms) と T_c スイープ `run_0110–0113` を投入。
- `2026-09-12 (17)` — **ξ 崩落の真因 = 離散 β-PDF の 1 次モーメント偏り × 緩和ソース**: T_c 1030 K の `run_0110` は 15 ms で吹き飛んだ後、
  ξ が z/d 10 以降で 0.008 に崩落 (couple 6 の `run_0102` と同じ症状)。`cmc_xiOm−xi` を ξ ビン別に見ると希薄側で −0.5〜−2.3 %、濃厚側で
  +0.1 % の系統偏り (端ビンの質量を η=0/1 節点に載せる + 台形則)。couple 7 の ω_s=ρ(Ỹ_pdf−Y)/τ_c はこの偏りを 1/τ_c=1e4/s で刷り込むので
  希薄側の燃料元素シンク (dξ/dt ≈ −50ξ /s) になり、定常反復の α/Δτ (~300/s) ではさらに速い → `run_0099/0101` の下流 ξ 崩落もこれ
  (柱モードは別件で実在)。**修正 71daa48e**: 目標 Ỹ_pdf, h̃_pdf を混合線に沿って (ξ̃−ξ_Ω) だけ平行移動し Bilger ξ(Ỹ_pdf)=ξ̃ を厳密化
  (反応離脱は不変)。CUDA 13 (CCCL 3) で `thrust::identity` が消えたので自前述語に (da44f74c)。**AWS g5 (44.222.71.7, `~/forge-chem`,
  native build `.build-native/relwithdebinfo`, `build` はシンボリックリンク; `/home/sano/work/...` の絶対パスもリンクで解決)** で
  `run_0111_dt_c7_xifix` (1045 K, 1500 step) を検証中。OK なら T_c スイープを AWS で回す。
- `2026-09-12 (18)` — **スイープ script の couple 5 残り**: `run_tc_sweep.sh` は couple 5 時代に書いたまま `--couple 5` を渡していた (run_0110 と AWS の旧 run_0112/0113 は couple 5 で無効; run_0112 は平均 T が PDF 診断 T より −200 K で気付いた)。`--couple 7` に修正 (commit 本)、AWS の polite ランチャー (30 分毎に GPU 空きを確認して 1 本ずつ) で 1045/1030/1060/1015/1075 K を再投入 (run_0112–0116)。教訓: 投入直後に `cmcQInit_d: ... couple=N` のログ行を確認する。

