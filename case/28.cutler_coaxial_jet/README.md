# 28.cutler_coaxial_jet — Cutler 超音速同軸 He/空気ジェット (多成分 TP 検証)

NASA Langley の超音速同軸ジェット (Cutler 他, AIAA-99-3588, NTRS 20040086859) を、多成分
thermally-perfect gas + 組成依存入口 BC + 乱流化学種拡散の検証ケースとして構築する。

## 計算 run 一覧

過去の run (run_0001〜0032) の詳細は `## 進捗` を参照。以下は DPLUR×(CPG/TP)×SST 検証 (2026-06-20)。

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0033_he_o2_n2_cpg_dplur` | CPG 層流 block-DPLUR 2次, cfl=0.5 (cfl=5 は初手発散) | 残差プラトー (限界周期, NOT CONVERGED)。`res_20000.h5` を SST のリスタート元に使用 | ref |
| `run_0034_he_o2_n2_tp_dplur`  | TP 層流 block-DPLUR 2次, cfl=0.5 (cfl=5 は step67 発散) | 残差プラトー (NOT CONVERGED)。SST リスタート元 | ref |
| `run_0035_he_o2_n2_cpg_dplur_cfl2` | CPG 層流 cfl=2 継続 (0033 から) | cfl=2 で安定 (NaN なし)。SST 要請で中断 | 削除済み(2026-08-31) |
| `run_0036_he_o2_n2_tp_dplur_cfl2`  | TP 層流 cfl=2 継続 (0034 から) | **cfl=2 で発散** (TP は cfl≤0.5 必須)。SST 要請で中断 | 削除済み(2026-08-31) |
| `run_0037_he_o2_n2_cpg_sst` | **CPG + RANS-SST**, 2次, block-DPLUR, cfl=2 (0033 から k/ω seed) | 安定・残差 2.5 桁低下 (`residual_history.png`)。但し **He コア T≈1007K (非物理)**。`cutler_profiles.png` 核長 x/D≈17 | active |
| `run_0038_he_o2_n2_tp_sst`  | **TP + RANS-SST**, 2次, block-DPLUR, cfl=0.5 (cfl=2 は step1 発散) | 安定・残差 2.6〜3.3 桁低下 (まだ降下中)。**He コア T≈300K (物理的)**。核長 x/D≈21 | active |

| `run_0039_tp_sst_cfl1` | TP-SST 継続 (0038 から) cfl=1.0 | **安定** (8000步 NaN なし)。残差は 0038 のフロア(rms_ro~3e-7)に張り付き横ばい | active |
| `run_0040..0044_tp_sst_cfl{2,4,8,1p5,1p2}` | TP-SST 継続 CFL 上限探索 (0038 から) | **全て発散** (cfl2=step50, 4=step1, 8=step9, 1.5=step116, 1.2=step247)。res 出力なし | 削除済み(2026-08-31) |
| `run_0049..0057_rycl_{A,B,C}_{cfl1,2,4}` | **rho-Y 共通リミタ診断** (`multispeciesRhoYCommonLimiter`)。A=S2 / B=S3 / C=S3+共通 min を同一 restart で比較。解析 `analyze_limitcycle.py` | cfl2 settled: A 最良(A_P=382), C は B 改善(2473→1585)も S2 に届かず, B 最悪。cfl4: B step433/C step550 発散, A 安定。**S2 維持・S3 棚上げ** ([plan §14](../../plans/active/convection-multispecies-contact-pressure.md)) | 削除済み(2026-08-31) |
| `run_0058..0063_thermoY_{S0,S2c,S2}_{cfl2,cfl4}` | **「界面Yは中心補間で十分か」検証**。S0(cell R_mix)/S2c(中心 face Y=`FORGE_FACE_THERMOY`)/S2(MUSCL) 比較 | cfl2: S2c が S0→S2 改善の 73-76% を回収=「界面Yで評価」が主因(高次でない)。但し cfl4 で **S2c 発散(step638)**・S2 安定 → ρ-Y 同一リミタが安定の鍵。**S2 が最小の正しい形** ([plan §15](../../plans/active/convection-multispecies-contact-pressure.md)) | 削除済み(2026-08-31) |

### node (median-dual) 版: S3 発散の原因確認 (plan species-passive-scalar-unification §4.2 ステップ 2, 2026-09-16)

同じ `mesh/coaxial.msh` を `discretization: node` で変換 (76800 ノード, wall_dist≡0 = cell と同じ)、IC は `run_0052/restart.h5` (cell 発達場) を `interp_field.py` で最近傍移植。
species/DB/BC は `run_0052` と同一 (旧 `LESorRANS` 系キーは `turbulence.model: sst` に置換)。TP 絶対基準 (thermoHrefTemp なし, cell と同じ)。全 run: block-DPLUR, `implicitRelax 0.7`, 2 次+Venkat, `detectNaN 1`。
「onset」= `rms_ro` がそれまでの最小値の 3 倍を初めて超えた step、「先行列」= step 0 の 3 倍を最初に超えた残差列。
**注意: ベースライン設定自体 (SFR 0 / coupling 1) が cfl 2 でも restart 後 ~2900 step で top 境界 (y≈0.06, outlet_statPress) または軸 (x≈0.08–0.09) から発散する遅い不安定を持つ** (run_0076) ので、
2000 step 以降の NaN は化学種設定と無関係 (ベースライン帯)。S3+coupling 0 だけがそれより 1 桁早く、He/air または coflow/ambient の**組成せん断層**で、**化学種残差列が先行**して発散する。

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0064_node_base_cfl05` | node ベースライン第 1 段: SFR 0 / coupling 1 / cfl 0.5, 2000 step (cross-mesh IC) | NaN 0, 5 ms/step, NOT CONVERGED (plateau; rms_ro 7.8e-5→2.3e-6) | ref |
| `run_0065_node_base_cfl2` | 第 2 段: 同設定 cfl 2, 4000 step (0064 res_2000 から) | NaN 0, NOT CONVERGED (plateau, rms_ro 1.7e-6, rms_roe 0.41 限界周期)。**マトリクスの共通 IC = `res_4000.h5`** | ref |
| `run_0076_node_s0_c1_cfl2_cont` | 0065 の同設定継続 (ベースライン自体の長期安定性対照) | **DIVERGED step 2947** (onset 1834, 先行列 rms_roe@1755), NaN 58 ノード top 境界 x 0.05–0.06 / y 0.056–0.060 | 診断 (ref) |
| `run_0072_node_s0_c1_cfl4` | SFR 0 / c1 / cfl 4 (cfl 4 対照) | DIVERGED step 2208 (onset 1445), top 境界 x 0.014–0.023 | 診断 (ref) |
| `run_0066_node_s2_c1_cfl4` | **S2 / coupling 1 / cfl 4** | DIVERGED step 2221 (onset 1453), top 境界 x 0.025–0.034 (ベースライン帯) | 診断 (ref) |
| `run_0067_node_s3_c1_cfl4` | **S3 / coupling 1 / cfl 4** (cell run_0052 の node 版) | DIVERGED step 2399 (onset 1478, 先行列 rms_roe@1422; **rms_roY は減少 ×0.7**), top 境界 x 0.000–0.009 (ベースライン帯 = **S3 固有の発散なし**) | 診断 (ref) |
| `run_0068_node_s2_c0_cfl4` | **S2 / coupling 0 / cfl 4** | DIVERGED step 2077 (onset 1471, 先行列 rms_roY0@1301 [緩やか ×2.3@1000]), 軸 x 0.080–0.089 (ベースライン帯) | 診断 (ref) |
| `run_0069_node_s3_c0_cfl4` | **S3 / coupling 0 / cfl 4** | **DIVERGED step 402** (onset 101, **先行列 rms_roY0/roY1@78**; s100: roY1 ×20.6, ro ×4.4; s200: roY1 ×716), NaN 106 ノード x 0.099–0.109 / y 0.021–0.024 = coflow/ambient 組成せん断層 | 診断 (ref) |
| `run_0070_node_s3_c0_cfl2` | S3 / coupling 0 / cfl 2 | DIVERGED step 3140 (ベースライン帯だが onset 334, 先行列 rms_roY1@219; s400: roY1 ×80, ro ×7), top 境界 x 0.063–0.072 | 診断 (ref) |
| `run_0071_node_s3_c0_cfl6` | S3 / coupling 0 / cfl 6 | **DIVERGED step 165** (onset 101, 先行列 rms_roY0@44; s100: roY1 ×198), NaN 61 ノード x 0.040–0.049 / y 0.0084–0.0105 = He/air せん断層 | 診断 (ref) |
| `run_0073_node_s0_c0_cfl2` | SFR 0 / coupling 0 / cfl 2 (c0 対照) | DIVERGED step 2548 (onset 1757), 軸 x 0.085–0.094 (ベースライン帯) | 診断 (ref) |
| `run_0074_node_s3_c1_cfl2` | S3 / coupling 1 / cfl 2 | DIVERGED step 2813 (onset 1681), top 境界 x 0.002–0.011 (ベースライン帯) | 診断 (ref) |
| `run_0075_node_s2_c0_cfl2` | S2 / coupling 0 / cfl 2 | DIVERGED step 2460 (onset 1730), 軸 x 0.084–0.093 (ベースライン帯) | 診断 (ref) |
| `run_0077_node_s3_c1_cfl6` | S3 / coupling 1 / cfl 6 | DIVERGED step 1712 (先行列 rms_roe@1154; rms_roY ×0.7), top 境界 x 0.019–0.028 (ベースライン帯; S3+c0 cfl 6 の 165 より 10 倍遅い) | 診断 (ref) |
| `run_0078_node_s2_c0_cfl6` | S2 / coupling 0 / cfl 6 | DIVERGED step 1571 (先行列 rms_roY0@1013 [×3@1000]), top 境界 x −0.002–0.007 (ベースライン帯) | 診断 (ref) |

結論 (2026-09-16): node では **S3 の発散は coupling 0 (segregated point-implicit; LHS が 1 次風上の `transport_diag`, 緩和なし) に固有**で、cfl 4 で step 402 / cfl 6 で step 165 (化学種残差が真っ先に成長し組成せん断層で NaN)。
**S3 + coupling 1 (scalar-DPLUR, relax 0.7) は cfl 4/6 でも S3 固有の発散を示さず** (化学種残差は減少)、cell `run_0052` (S3+c1, step 433) とは逆。化学種の擬似 CFL だけを絞る config キーは存在しない (未試験)。
ベースラインの遅い不安定 (top `outlet_statPress` / 軸, ~1500 step で rms_roe 先行) は本 case の node TP-SST 設定の別問題 (thermoHrefTemp 未設定・top BC 等が候補、未切り分け)。

**TP-SST 継続の CFL 上限 (run_0039〜0044, 発達場 0038 から)**: **cfl=1.0 が安定上限** (1.2 以上は発散)。ただし場は cfl=0.5/20000步 で既に残差フロア(限界周期)に達しており、cfl を上げても収束は速くならない (律速はステップ数でなく**安定性**)。冷流 TP ジェットは剛性が高く、層流時の cfl≤0.5 に対し発達 SST 場でも継続上限は ~1.0。

**結論 (CPG vs TP, DPLUR, SST)**:
- **TP は剛性が高く cfl=2 で発散** (層流・SST とも), cfl≤0.5 が必要。**CPG は cfl=2 で安定**。
- **SST 投入で収束が劇的に改善** (層流はプラトー/限界周期 → SST で残差 2.5〜3.3 桁低下)。せん断層の非定常を渦粘性が抑える。
- **CPG は He ジェットを物理的に表現できない**: 単一 γ/cp (空気値) では He の高 R を表せず、He コア温度が ~1007K に誤熱化する (本来 ~300K の冷流)。**TP は He コア ~300K を正しく再現**。多成分混合では TP が必須であることの明確な実例。

## 実験条件 (AIAA-99-3588)

| | 中心ジェット | coflow |
| --- | --- | --- |
| ガス | He 95% + O2 5% (体積) | 空気 |
| Mach | 1.8 (nominal) | 1.8 (nominal) |
| 出口圧 | 1 atm (両者圧力整合) | 1 atm |
| 形状 | 軸対称、中心ジェット出口径 ~10 mm (r≈5.12 mm)、center body lip r≈1.74 mm | |

- 対流マッハ数 $M_c=|u_j-u_\infty|/(c_j+c_\infty)=0.7$。中心ジェットは音速が大きく速度は coflow の 2 倍超。
- 全温 ~300 K (冷流)。検証データ: He モル分率分布、pitot 圧、全温の probe survey (RELIEF/GSAS)。
- **非反応** (He は燃料 simulant)。化学はスコープ外なので好都合。

## 進捗

- **メッシュ** (`mesh/coaxial.geo`): 軸対称 2D 同軸 (中心ジェット/lip/coflow/ambient の 4 バンド構造格子)。
  境界 physID = axis(1)/jet_in(2)/lip(3)/coflow_in(4)/amb_in(5)/top(6)/outflow(7)。38240 セル。
  → **`convertGmshToForge` で一発変換成功、forge も NaN なしで稼働** (`run_0001_meshcheck`, CPG 単成分で
  流れ場が立つことを確認。残差 3.4e-5→1.8e-5、P/T/Ux 有限)。
- **組成依存超音速入口 BC** (plan §10 / M5, 実装・検証完了): `inlet_uniformVelocity` に入口組成 `Y_s`
  を追加 (config `floats: Y0,Y1,...`)。ghost は $\rho=P/(R_{mix}(Y)T)$, $roe=\rho(e_{mix}(Y,T)+e_k)$,
  $\rho Y_s=\rho Y_s^{in}$ (Dirichlet)。あわせて SLAU 対流流束の被移流エンタルピーをセル組成で再構成
  (He/N2 で cp 質量基準が ~5×差のため `sp[0]` 固定だと混合 contact で圧力発散していた)。
  - `run_0003_he_n2_inviscid/` — binary [He,N2], 1次風上・非粘性, 2000 步。NaN/発散なし、$\sum Y=1$、
    中心 He コア (Y_He=1, Ux~1280 m/s = Mach1.8 He)、coflow N2。
  - `run_0008_he_n2_2nd/` — 上を 2次風上でリスタート (計 4000 步)。He ポテンシャルコアが ~0.07m へ伸長。
    成果物: `field_YHe_Mach.png`, `residual_history.png`。
  - `run_0002_cpg_regr/` — CPG 単成分回帰 (修正が単成分/CPG を壊さないことの確認)。
- **startup 注意**: 静止 IC (u=0) + `outlet_statPress` は出口の backflow 分岐で `Pt`/`Tt` 未指定だと NaN。
  IC に下流速度 (u=100) を与え、出口 bcond に `Pt`/`Tt` を指定すること (`gen_binary_ic.py` 参照)。
- **次 (後続)**: 定量照合 (He コア長・濃度減衰 vs RELIEF survey) は粘性 + RANS-SST + 乱流化学種拡散が必要。
  その後 `[He,O2,N2]` 3 成分へ拡張。

## ガス物性 (TP 用)

- 中心ジェット He-O2 (95/5 vol): MW≈5.40 g/mol, γ≈1.65 (He 単原子 + O2 二原子の混合)。
- coflow 空気: N2/O2。簡易には binary He/N2 から始め、後で He-O2/air へ拡張。
