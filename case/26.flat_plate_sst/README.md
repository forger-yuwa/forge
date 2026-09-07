# 26. 乱流平板 (SST 壁法則検証)

Menter SST k-ω モデルが乱流平板の**壁法則 (law of the wall)** を再現できるかを確認するケース。
陰解法 (block DPLUR, `timeIntegration: 11`) + SLAU + MUSCL で定常計算する。

## 流れ条件

- 入口: 全圧 `Pt=100000 Pa`, 全温 `Tt=288.15 K` (`inlet_Pressure`)
- 出口: 背圧 `Ps=97250 Pa` (`outlet_statPress`)、`M≈0.2` (U∞≈67.8 m/s)
- 空気 `mu=1.8e-5 Pa·s`, `Re/m≈4.46e6` → 平板長 1 m で `Re_L≈4.5e6` (乱流)
- 流入乱流: `k=0.3, omega=300` (freestream `mu_t/mu≈60`, naca_ml 相当)

## 領域・メッシュ (`mesh/flat_plate.geo`)

- `x∈[-0.1,1.0]`, `y∈[0,0.2]`, z 厚さ 0.01 (1 層押し出し)
- 底面: `x∈[-0.1,0]` slip (対称), `x∈[0,1.0]` no-slip wall (平板)
- 上面 slip, 入口/出口は上記 BC、z 両面 slip
- 壁垂直: 90 点, 第一セル高さ **4.14e-6 m**, 等比 1.10 → **y⁺₁≈0.35〜0.4 (<1)**
- streamwise: 前縁 (x=0) へ上流側・平板側の両方を寄せ、前縁でセルサイズ連続 (≈0.4 mm)

## 実行手順 (段階法)

**重要**: 定常計算の実効 CFL は `cfl` ではなく `cfl_pseudo`
(`procedures/solver-settings.md` の「CFL の定義」を参照)。
本ケースは壁セルが極めて細く陰解法でも CFL 制限が低いため、段階的に進める。

1. **stage A (`run_0005_slau_spinup`)**: 一様初期値からの cold start。
   `convMethod: 0` (1次風上) + `cfl_pseudo: 2.0` で 25000 step 回し、平滑な発達場を作る。
2. **stage B (`run_0006_slau_muscl`)**: stage A の `res_25000.h5` から restart。
   `convMethod: 2` (MUSCL) + `cfl_pseudo: 3.0` で 40000 step 回し、2次精度で収束させる。

> cold start で `cfl_pseudo≳5` や、いきなり MUSCL + 高 `cfl_pseudo` にすると発散する。
> 発達場からの restart でも本ケースの安定上限は `cfl_pseudo≈3` 程度。

メッシュ生成・変換・実行は `procedures/calculation-workflow.md` の標準手順に従う
(Docker, `convertGmshToForge`, `forge`)。

## ポスト処理

```bash
python3 tools/postprocess_wall_law.py run_0006_slau_muscl 0.3 0.6 0.89
```

`wall_law.png` (u⁺-y⁺ と Cf-Re_x) を出力する。

## 収束について

陰解法 `blockDPLUR` は減衰付き point-Jacobi の defect correction で、近似ヤコビアン+高アスペクト比
壁セル+陽的 SST 生産項のため安定 `cfl_pseudo` は MUSCL で ≈5(cold start は 2〜3)が上限。
そのため収束が遅く、`run_0006`(40000 step)時点では外層がまだ発達途中(Cf が 10000 step あたり
3〜4% ドリフト)だった。

- **stage C (`run_0007_slau_muscl_innersweep`)**: `run_0006/res_40000.h5` から restart し、
  `cfl_pseudo: 5.0`, `nStepInner: 20`, MUSCL で **120000 step** 継続。残差 3.8e-7、
  最後の 20000 step での Cf ドリフト 0.1% 以下まで収束。**最終結果はこの run を使う。**

> 安定 CFL が低い構造的原因と改良案(SST 生産項の陰的化, 真の粘性ヤコビアン)は
> `.github/plans/implicit-acceleration-session-prompt.md` に整理(別セッションで実装予定)。

## 結果 (run_0007, 収束済み MUSCL)

| station | Re_x | y⁺₁ | u_τ | Cf | Cf/Schlichting |
|---|---|---|---|---|---|
| x=0.30 | 1.34×10⁶ | 0.36 | 2.707 | 3.12×10⁻³ | 0.885 |
| x=0.60 | 2.66×10⁶ | 0.34 | 2.580 | 2.82×10⁻³ | 0.918 |
| x=0.89 | 3.99×10⁶ | 0.33 | 2.511 | 2.67×10⁻³ | 0.944 |

- u⁺ vs y⁺: 3 ステーションが普遍プロファイルに重なる。粘性低層 `u⁺=y⁺` (y⁺<5)、
  対数則 `u⁺=(1/0.41)ln y⁺+5.0` (y⁺=30〜300)、外層 wake までよく再現。
- Cf vs Re_x: 発達領域で Schlichting 1/5 乗則と 1/7 乗則の間に乗る。1/5 乗則より数% 低いのは、
  前縁からの層流→遷移区間ぶんの virtual-origin オフセット(Re_x を前縁起点で測るため)による系統差。
- y⁺₁≈0.33〜0.36 (<1) を全域で満たす。

## 計算 run 一覧

| run_* | 目的・設定 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_rans_implicit` | 標準: SLAU + SST + block-DPLUR 陰解 (2D) | rms_ro 2.08e-9 収束、壁法則再現 | active |
| `run_regr_cf` | 回帰: 閉形式 FVS 既定化 (`implicitSolvePrecision=0`) の確認。run_0001 と同条件 | 残差収束が legacy と一致 (rms_ro 2.0e-9)、NaN なし → **閉形式は 2D 陰解に無害**。plan [precision-mixed-axisym.md](../../plans/archived/precision-mixed-axisym.md) | active |
| `run_0009_ewt_fine_mode1` | EWT 回帰: 細メッシュ (y⁺₁≈0.35) で `wallTreatmentSST:1`。run_0007 収束場から restart 20000 step | mode1 ≈ mode0(run_0007): Cf/Schl 0.93/0.96/0.99、残差 floor 全列一致 → 壁解像に縮退 | active |
| `run_0010_ewt_yp30_mode0` | EWT 検証: 粗メッシュ y⁺₁≈30, mode0(low-Re)。cold start→MUSCL | **Cf 崩壊** Cf_model/Schl 0.13/0.13/0.14 (low-Re BC が粗メッシュで破綻) | active |
| `run_0011_ewt_yp30_mode1` | EWT 検証: 粗メッシュ y⁺₁≈30, mode1(enhanced) | **Cf 回復** Cf_model/Schl 0.89/0.91/0.93、u_τ≈細メッシュ一致 → y⁺ 非依存 | active |
| `run_0012_ewt_yp10_mode0` | EWT 検証: バッファ帯 y⁺₁≈10, mode0(low-Re) | Cf_model/Schl 0.51/0.56/0.60 (過小) | active |
| `run_0013_ewt_yp10_mode1` | EWT 検証: バッファ帯 y⁺₁≈10, mode1(ω-blend+τ_w, k=0 Dir) | Cf_model/Schl 1.05/1.10/1.14 (過大) → full WF で改善 | ref |
| `run_0016_ewt_fine_wf`  | **full WF** (ω pin + k zero-grad + P_k log則化): 細 y⁺₁≈0.35 回帰 | Cf/Schl 0.89/0.92/0.95 (=壁解像), Cf_molec≈Cf_model → 壁解像極限を保持 | active |
| `run_0017_ewt_yp30_wf`  | **full WF**: 対数帯 y⁺₁≈30 | Cf/Schl 0.92/0.94/0.97, u_τ≈2.5–2.6 | active |
| `run_0018_ewt_yp10_wf`  | **full WF**: バッファ帯 y⁺₁≈10 | Cf/Schl 0.95/0.98/1.01 (ω-blend単独の過大を解消) | active |
| `run_0019_des_regr_mode0` | **SST-DES T1-A 回帰**: 発達場(run_0007/res_120000)から `DESmode:0` 1000 step。新旧バイナリ比較 | DESmode:0 は baseline と atomicAdd 床内で一致 (1 step base–base ≡ base–des: ro 1.19e-7/roe 1.56e-2/roK 3.73e-9) → **ビット不変**。plan [turbulence-iddes-sst.md](../../plans/active/turbulence-iddes-sst.md) | active |
| `run_0020_des_mode1` | **SST-DES T1-A シールド (wall-resolved y+1)**: 同場から `DESmode:1` 1000 step | mode1 vs mode0 `roUx` relL2 **5.7e-7**≪Cf 0.1%, `roK` 1.8e-5。付着 BL シールド (DES 発火は外縁 3%)、NaN なし | active |
| `run_0021_des_ewt_yp30_mode0` | SST-DES T1-A: y+30 wall-modeled 場 (run_0011) から `DESmode:0` 1000 step (mode1 比較基準) | NaN なし | ref |
| `run_0022_des_ewt_yp30_mode1` | **SST-DES T1-A シールド (wall-modeled y+30)**: 同場 `DESmode:1` 1000 step | mode1 vs mode0 `roUx` relL2 4.7e-5, `roK` 9.3e-4、付着 BL シールド (frac f_d<0.05=0.92)、NaN なし | active |

### node-centered (median-dual) 検証 run

| run_* | 目的・設定 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_node_lam_cont` | node モード層流平板 (visc=4e-4, Re_L≈1e5, M=0.2) | Blasius BL 再現 (δ99=0.99×, 形状誤差 3.3%)、図 `blasius_validation.png` | active |
| `run_node_sst_fine` | node モード SST (cell run_0001 から cross-mesh restart, 20000 step) | omega/k 6 桁収束だが u_τ 過小 (node 1.24 vs cell 1.97)。残差プラトー (rms_roUx≈0.23)。図 `uplus_yplus_node_vs_cell.png` | active |
| `run_node_lam_5e_long` | **node 弱形式境界 (Phase 2: 5a+5b+5e)** 層流 (40000 step, run_node_lam_cont 崩壊場から restart) | **出口 BL 崩壊解消** (δ99(x) 単調・Blasius 一致, x=1.0 δ99=0.0114 vs 旧 1e-5)。**残差プラトー打破** rms_roUx 0.214→3.08e-5 (3.8桁), rms_ro 1.04e-7。NaN なし | active |
| `run_node_sst_5e_long` | **node 弱形式境界 (5a+5b+5e)** SST (40000 step) — 過剰乱流の「before」 | プラトー打破 (rms_roUx 0.597→1.3e-3) だが **過剰乱流**: peak mut/μ=375, Cf≈3×Schl。真因=壁 ω 非ピン + wall_dist バグ (下行で解消) | ref(before) |
| `run_node_sst_final` | **node SST 完成**: massflux 書込 + 壁 ω Dirichlet ピン (omega/roOmega ピン+ω decouple+res_roOmega 壁ゼロ) + **wall_dist 修正** (node で壁点に壁ノード座標)。過剰乱流場から 60000 step。**convMethod 0 (1次風上) のまま** | u⁺-y⁺ collapse は cell 一致だが **Cf/δ99 は cell(MUSCL)と乖離**: Cf −8%, δ99 +13〜20% (1次の数値拡散で BL 過厚)。かつ res_60000 で Cf 未収束 (step30k→60k で +20% ドリフト)。massflux+ω pin+wall_dist 修正で過剰乱流は解消済 | ref(1次) |
| `run_node_sst_muscl_cont` | **node SST を MUSCL(convMethod 2)に切替**え run_node_sst_final/res_60000 から +90000 step。cell 基準 (run_0009, MUSCL) とスキームを揃えた公平比較 | **node≈cell に一致**: Cf 差 **−3.4〜−3.9%**, δ99 差 **−1.0〜−1.5%**, u_τ −2.0〜−2.3% (3 station)。Cf@0.6 ドリフト +0.20%→+0.06%→0.00% で定常化、NaN なし。残差は block-DPLUR 構造プラトー (cell run_0007 と同様)。図 `cf_bl_cell_node_muscl.png` | active |
| `run_node_skip_verify` | **node 境界半割面拡散 skip 検証 (新)**: run_node_sst_muscl_cont/res_90000 から +5000 step、k/ω 境界半割面拡散を ∇·S 弱形式→skip 化したバイナリ。plan [diffusion-node-boundary-real-distance.md](../../plans/accepted/diffusion-node-boundary-real-distance.md) §3(c) | **Cf/u_τ/δ99 は ∇·S と完全一致** (3 station 差 0.00%, ref90k と 0.01% 以内)。k 場差 (relL2 6.5%) は **全て前縁上流 x<0 のスリップ域**に局在 (∇·S が Neumann 漏れで k を ~10 に過剰生成→skip は ∂k/∂n=0 で k≈3=内部値、skip がより正)。平板 BL は無影響。NaN なし | active |
| `run_node_gradS_verify` | 上記の比較基準: 同じ restart から **旧 ∇·S 弱形式バイナリ**で +5000 step | skip と Cf 完全一致を確認するための apples-to-apples ペア (k のみ前縁スリップ域で skip と差) | ref |
| `run_0004_keep_es_rans_m02` / `run_0005_keep_es_rans_m02_cont` | **単一スキーム KEEP+ES の RANS 検証 (第1試行)**: run_0003 (回帰テンプレート) の場から KEEP+matrix ES (σ0.05 jump2 precond) + implicit cfl20 で 10k+40k | KEEP+implicit 初の組で**発散なし**。ただし run_0003 の土俵は Cf 未検証 (理論+30%・x 無依存) で A/B 判定に不適と判明 → run_0012 へ移行。10k 時点は Ux ドリフト 1.1%/5k で未収束だった (準定常確認の教訓) | 記録 (土俵不適) |
| `run_0023_isoT_fine_ref` | **等温壁 (Tw=320K) 壁解像リファレンス**: run_0009 (y+0.35 EWT 細) 収束場から `wall_isothermal` 化、通算 70k step (q_w(x) ドリフト 0.2%/10k まで静定) | 発達域 q_w ≈ 4.0-4.7 kW/m²。**熱壁関数ギャップ実測の基準** | ref (等温 ref) |
| `run_0024_isoT_yp30_wf` | **等温壁 y+30 + `wallTreatmentSST:1`** (熱壁関数なしの現行 EWT): run_0021 場から同 Tw、20k step (ドリフト 0.6% 静定) | **q_w を +39〜+49% (平均 +45%) 過大予測** (x=0.3-0.9, vs run_0023)。原因: 第 1 セル (y+30) の大きな μt を線形 (Tw−T1)/y1 勾配に掛けるため sublayer 熱抵抗を無視。**等温壁×粗 y+ EWT は熱壁関数 (Kader 相当) が必須** — 実装は [wmles plan](../../plans/active/turbulence-wmles-wall-stress.md) の Kader 部品流用で可能 (未着手) | ref (熱壁関数ギャップ実測) |
| `run_0025_sstopt_{kpin0,base,all}` | **SST 整合オプション A/B** (run_node_sst_muscl_cont/res_90000 から 5000 step, config は新キー体系に更新: `model: "sst"`, dilatation 0): kpin0=`sstNodeWallKPin:0` (旧挙動) / base=既定 / all=`sstOmegaProdFromPk,sstSigmaBlend,sstIsotropicStress,sstEnergyKSource` 全 ON。Cf は `tools/cf_node.py` (node 用: 第一内点勾配) | **Cf/Schlichting @x=0.3/0.6/0.9: kpin0 0.8899/0.9292/0.9553 = base (完全一致) / all 0.8897/0.9289/0.9552 (差 ≤0.03 %)**。平板 (F1=1, リミッタ非発動, ρk≪P) では全オプション無影響。plan [turbulence-sst-consistency-options](../../plans/active/turbulence-sst-consistency-options.md) | active (A/B 記録) |
| `run_0026_sstdef_{ref,keys0,new}` | **SST 整合オプション既定 ON の回帰** (run_0025_sstopt_base/res_5000 から +5000, index コピー): ref=旧バイナリ b12d74a4 / keys0=新バイナリで `sstOmegaProdFromPk:0, sstSigmaBlend:0` 明記 / new=新既定 (`sstOmegaProdFromPk:1, sstSigmaBlend:1`, F1 は同 step 前処理) | keys0−ref ≤6e-5 (ビット同等)。**Cf/Schlichting @0.3/0.6/0.9: 0.8920/0.9308/0.9561 (ref) → 0.8920/0.9307/0.9562 (new), 差 ≤0.01 %**; k 2 %・μt 7 % (下流端)。`check_convergence` rms_roK 5e-5 falling / roOmega plateau (参照どおり), NaN 0。plan [turbulence-sst-consistency-options](../../plans/active/turbulence-sst-consistency-options.md) §3.6 | active (A/B 記録) |
| `run_0026_walldist_node` | node 壁距離のノード座標基準化 (2026-09-08) の回帰: res_90000 から 5000 step | wall_dist 不変 (均一メッシュで重心=節点)、Cf/Schl 0.8899/0.9292/0.9553 = run_0025_sstopt_base と同一 | active (回帰記録) |
| `run_0025_isoT_yp30_qwwf` | **`sstEnergyWallFunction:1` 初版** (Kader q_w 流束置換, [plan](../../plans/active/turbulence-sst-thermal-flux-model.md)) を y+30 で。warm from run_0024 | q_w +42% 過大 → **−12.9%** (vs run_0023)。残った系統過小の切り分けで **Kader T⁺ の式混用バグ発見** (対数部が Pr_t(u⁺+β_Kader) — Jayatilleke 形と Kader β の混用で T⁺ +30% 過大) | ref (混用式の記録) |
| `run_0026_isoT_yp30_qwwf_kader` | 同上 + **Kader 原式修正** (`wallLaw_kader_tplus` = 2.12·ln(1+y⁺)+β 形へ)。warm from run_0025 | **q_w = 真値 +0.1%** (range −0.7〜+1.4%, x=0.3-0.9)。qwall drift 0.17%/10k で静定 | active (**等温壁 y+30 正**) |
| `run_0027_isoT_fine_qwwf` | Kader 原式を細メッシュ y+0.35 で (low-Re 極限回帰)。warm from run_0023 | **q_w = 真値 −0.1%** — 置換が解像伝導へ正しく縮退 | active (low-Re 極限) |
| `run_0028_isoT_yp10_qwwf` | Kader 原式をバッファ帯 y+10 で。warm from run_0018 (断熱 wf 場) | **q_w = 真値 +6.2%** (range +5.7〜+7.3%) — ブレンド帯の残差で許容 | active (バッファ帯) |
| `run_0029_node_outletic_subsonic` | **亜音速出口の背圧アンカー検証** (node outlet の p_tilde/勾配閉包を内部値参照に改訂した最終形, [plan §2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)): run_node_sst_final 収束場から 12000 step | **出口 P 平均 97250.6 Pa (規定 97250, ±3 Pa) を保持** = p_tilde=P[ic] でも特性構成 bvar の質量流束が背圧を係留。残差悪化なし・NaN なし | active (亜音速アンカー検証) |
| `run_0030_node_yp30_wf_omgoff` | **node × y+30 壁関数の平板** (ω irep ピン A/B 用に新規 node 変換 `flat_plate_yp30_node.h5`)。段階起動: 1次 cfl2 (清浄収束) → MUSCL cfl3 | 1 次段は健全 (3.9 桁)。**MUSCL 段は残差 2-4 桁上昇のプラトーで Cf が 10%/10k step ドリフト = 準定常に達せず**。convMethod 1 は清浄場からも発散。**node × 亜音速平板 y+30 の数値レシピは未確立** (別課題)。ω ピン判定は case/40 の y+1 真値比較で代替済 ([plan §10c](../../plans/accepted/turbulence-sst-su2-taw-coupling.md)) | 削除済み(2026-08-31) (レシピ未確立の記録) |
| `run_0031_node_yp30_stageA` | **node × y+30 壁関数の段階起動 stage A** (1次 `convMethod:0` + cfl_pseudo 2 + nStepInner 20, cell `run_0011` からcross-mesh interp) | 健全 (roK/roOmega 3.6-3.9 桁低下、falling)。**node 亜音速平板の唯一の安定な起点** | active (node 平板 seed) |
| `run_0032_node_yp30_2nd_base` | 2次 (MUSCL, limiter venkata) を stage A 場から 3000 step | **残差上昇 (発散性)**。原因は**前縁 x≈19mm の局所暴走** (Ux 199 m/s = 自由流 67.8 の 3 倍・P 90.5 kPa) | 削除済み(2026-08-31) (2次失敗の記録) |
| `run_0033_node_yp30_2nd_lowmach` | 同 + `lowMachPrecond: 2` (低マッハ市松対策の転用) | **流れ場が完全凍結** (`ro`/`roe` がビット不変、`roUx` の変化 4.7e-34 = 非正規数)。乱流のみ進行 → 本構成では使用不可 | 削除済み(2026-08-31) (要調査の記録) |
| `run_0034_node_yp30_2nd_lsq` | 同 + `gradLSQ: 2` (LSQ 勾配の検証も兼ねる) | **NaN 発散** (GG より悪化)。LSQ は本不安定の対策にならない — [LSQ plan](../../plans/active/discretization-lsq-gradient.md) 検証4 | 削除済み(2026-08-31) (LSQ A/B 記録) |
| `run_0035_node_yp30_2nd_barth` | 同 + `limiter: 1` (barth) | **NaN 発散** | 削除済み(2026-08-31) |
| `run_0036_node_yp30_2nd_cfl05` | 同 + `cfl_pseudo: 0.5` | 残差上昇 (cfl 2 と同率) → **時間積分でなく空間離散の問題**と確定 | 削除済み(2026-08-31) |
| `run_0037_node_yp30_2nd_bfo` | 同 + `bndFirstOrder: 1` | 残差上昇 (改善せず) | 削除済み(2026-08-31) |
| `run_0038_node_yp30_1st_long` | **node y+30 の生産検証** (1次, stage A から通算 50000 step) | 残差床到達 (rms_ro 1.1e-9)。**Cf は完全定常** (20k/30k/40k で 0.00242 不変)。**Cf/真値 (run_0007 壁解像) = 0.835-0.840 で x 依存なく一様**。1次風上の散逸が支配的とみられ、node 壁関数自体の評価には 2 次が必要 | active (**node 平板の現状ベスト**) |
| `run_0039_node_yp30_planar_2nd` | **真 2D (平面) メッシュへの切替 + 2 次精度 (3000 step)**。メッシュ = 新規 `mesh/flat_plate_yp30_planar.geo` (押し出しなし, 10755 ノード, 品質 PASS: AR max 58.8 / skew 0.000)。IC = run_0038 res_40000 から cross-mesh interp | **2 次が初めて安定**。全残差 falling・NaN なし・前縁暴走なし。旧 2 次 (run_0032-0037) の発散が**メッシュ側の欠陥**だったことの直接証拠 | active (2 次成立の初回証拠) |
| `run_0040_node_yp30_planar_2nd_long` | 同構成の本収束 (run_0039 から 40000 step) | NaN なし・Ux max 71.3 m/s (run_0032 の 199-249 の暴走は消滅)・Uz ~1e-14。残差は block-DPLUR 構造プラトー (本 case 既知, run_0007 も同様) だが **quasisteady `ALL STEADY`**、**Cf は 20k/30k/40k で 4 桁一致**。**Cf/真値 = 0.848/0.847/0.848** | active (**node 平板の生産基準**) |
| `run_0042_node_yp30_planar_2nd` | **壁関数 y+ 非依存性テスト**: 第一層 1.768e-4 (`mesh/flat_plate_ny52_planar.geo`, node 第一 DOF y+27.1 = cell 第一 DOF と同高さ)、run_0040 から cross-mesh interp、20000 step | $u_\tau$=2.343 (ny=45 と同値)。**cell (同高さ 1.700e-4) の 2.586 に対し 9.4% 低いまま** = 第一層高さ仮説の棄却 | active (y+ 掃引) |
| `run_0043_node_yp112_planar_2nd` | 同掃引の粗い側: 第一層 6.588e-4 (`mesh/flat_plate_ny38_planar.geo`, y+101.5) | $u_\tau$=2.350。**y+ 27→102 (3.7 倍) で $u_\tau$ の変化 0.3%** = node 壁関数は y+ 非依存 | active (y+ 掃引) |
| `run_0044_cell_yp30_wf_regr` | cell y+30 壁関数を**現行バイナリで再走** (run_0017 は旧バイナリで `wf_pk`/`Pk_diag` 未出力のため対照にならない)。run_0017 res_10000 から restart、5000 step | $u_\tau$=2.586 (run_0017 と一致)。生産収支比較の cell 基準 | active (cell 対照) |
| `run_0045_node_yp30_kwfdir0` | **`nodeKwfDirichlet: 0` A/B** (既定 1 = 第一内層 k を Dirichlet 固定)。run_0042 と同一 IC/設定 | 実効生産が 3.1 倍になるが $u_\tau$ は +0.55% のみ → **「生産規約が真因」を反証** | active (反証の証拠) |
| `run_0046_node_yp30_omgwfdir1` | **`nodeOmegaWfDirichlet: 1` A/B** (第一内層の $\omega$ も固定)。run_0042 と同一 IC/設定 | $\omega$@第一内層 2.645e5→1.213e5。ただしこのスイッチは $\omega$ ピンと limiter 迂回を**束ねている**ので単独因子の証拠にならない (§3 の 2×2 で分離済: 平均主効果 pin +0.108 / bypass +0.016 / 相互作用 +0.021) | active (2×2 の $Y_{11}$ 相当) |
| `run_0050_node_lowre_fine_regr` | **node の壁解像 (壁関数なし) 検証**: planar 細メッシュ (`flat_plate_planar.h5`, y+0.7) + `wallTreatmentSST: 0` + MUSCL、現行 binary、`run_node_sst_muscl_cont` res_90000 から restart 20000 step | **$C_f/C_{f,\mathrm{KS}}$ = 0.943** (cell 壁解像 0.954 と 1.1% 差)、収支 1.036、quasisteady **ALL STEADY** → **node 離散化そのものは健全**と確定 | active (**node 無罪の証拠**) |
| `run_0047_su2_sst_ny52` | SU2 v8.5.0 low-Re (壁関数なし)、`run_0042` と同一メッシュ | $C_f/C_{f,\mathrm{KS}}$=0.365、収支 0.181、$\mu_t/\mu$ 340–452 で**破綻** → low-Re を y+27 に当てる誤りは SU2 でも同じ | 削除済み(2026-08-31) (設定誤りの記録) |
| `run_0048_su2_sst_wf_ny52` | **SU2 + `STANDARD_WALL_FUNCTION`**、`run_0042` と**完全同一メッシュ・BC・物性・流入乱流** | 壁 $C_f$/KS=0.793–0.796 だが **運動量収支 0.790 = 21% 破綻**、残差も rms −4.4 止まり (壁解像版は −12.4) → **解として棄却**。「SU2 も 0.75」という当初結論は撤回 | 削除済み(2026-08-31) (未収束・収支破綻の記録) |
| `run_0049_su2_sst_lowre_fine` | **SU2 壁解像** (planar 細メッシュ, 第一 DOF y+0.70, low-Re) — node 壁解像の独立検証 | **rms −12.4 まで収束**。壁 $C_f$/KS=**0.938–0.944**、運動量収支 **0.999/0.995 (ほぼ完璧)** → forge node 壁解像 0.943 と **0.5% 以内で一致** | active (**node 離散化検証の正**) |
| `run_0051_node_wf_widiag` | **W–I 実力診断の初回** (`wi_ftan`/`wi_fnrm`/`wi_ftan_res` 追加後、`run_0042` 場から 500 step) | $\Sigma$`wi_ftan`/$\int\rho u_\tau^2dx$=1.016–1.043、再スケール係数 0.640 | active (§4.2 初回) |
| `run_0052_regr_bitinv` | `wf_irep_flag` 追加後の既定経路回帰 (`run_0051` と同一入力) | 差は rel 1e-6 = node の run-to-run 床 → 既定経路不変 | 削除済み(2026-08-31) (回帰記録) |
| `run_0053_2x2_Y00` | 2×2 factorial $Y_{00}$: $\omega$ pin=0 / limiter bypass=0 (= 既定、`run_0042` 相当を新 binary で) | (w) W–I 実力 $C_f$/KS = **0.7615**、(c) 積分 0.784 | active (factorial 基準) |
| `run_0054_2x2_Y10` | 同 $Y_{10}$: pin=1 / bypass=0 (`FORGE_WF_LIMITER_BYPASS=0`) | (w) **0.8594**、(c) 0.840 → 単純効果 pin +0.0979 | active (factorial) |
| `run_0055_2x2_Y01` | 同 $Y_{01}$: pin=0 / bypass=1 (`FORGE_WF_LIMITER_BYPASS=1`) | (w) **0.7670**、(c) 0.794 → 単純効果 bypass +0.0055 | active (factorial) |
| `run_0056_2x2_Y11` | 同 $Y_{11}$: pin=1 / bypass=1 (= 従来の `nodeOmegaWfDirichlet:1`) | (w) **0.8864**、(c) 0.846。平均主効果 pin +0.1087 / bypass +0.0163 / 相互作用 +0.0215 | active (factorial) |
| `run_0057_widiag2` | 法線力を**絶対値**で再測定 (`wi_fnrm_abs` 追加後) | $\Sigma\lvert F_n\rvert/\Sigma\lvert F_t\rvert$=0.0027%、局所最大 6.1e-5、再スケール係数 min 0.6400/median 0.6405/max 0.6410 | active (法線力無害の正) |
| `run_0053_2x2_Y00_widiag` | $Y_{00}$ に `FORGE_WI_FORCE_DIAG=1` を付けた 500 step 再走 (W–I 実力取得用) | (w) 列の出所 | active (診断) |
| `run_0054_2x2_Y10_widiag` | 同 $Y_{10}$ | (w) 列の出所 | active (診断) |
| `run_0055_2x2_Y01_widiag` | 同 $Y_{01}$ | (w) 列の出所 | active (診断) |
| `run_0056_2x2_Y11_widiag` | 同 $Y_{11}$ | (w) 列の出所 | active (診断) |
| `run_0063_B_nasa_wf` | 候補 B: **TMR 推奨範囲内で選んだ条件** ($k$=8.4215e-5, $\omega$=616.2 → $\mu_t/\mu$=0.009) × y+27 壁関数 | $C_f$/KS=**0.764** (現行 65.9 と同値) → **今回の準定常評価量では検出可能な影響なし**。全残差は `NOT CONVERGED`、$\theta$/$C_f(Re_\theta)$ は `STEADY` | active (候補 B 棄却の正) |
| `run_0064_B_nasa_ref` | 同 TMR 範囲内条件 × y+0.7 壁解像 | $C_f$/KS=0.943 (現行 0.944)、実測/Reichardt=**0.9494** (現行 0.9492)。全残差 `NOT CONVERGED` / 判定量 `STEADY` | active (同上) |
| `run_0065_B_mid_wf` | 中間条件 ($\mu_t/\mu$=1.0) × y+27 壁関数 (単調性確認) | $C_f$/KS=0.764。全残差 `NOT CONVERGED` / 判定量 `STEADY` | active (単調性確認) |
| `run_0066_B_mid_ref` | 同中間条件 × y+0.7 壁解像 | $C_f$/KS=0.943、実測/Reichardt=0.9494。全残差 `NOT CONVERGED` / 判定量 `STEADY` | active (単調性確認) |
| `run_0073_diag_reichardt` | §5.2 ② 診断実験の対照 (env OFF)。`run_0063` 収束場から 20 step | `utau` 基準値。変更前バイナリとの差は run-to-run ノイズと区別不能 | active (② の対照) |
| `run_0073b_baseline` | 同上を**変更前バイナリ**で実行 (env OFF 不変性の検証用) | 変更前 vs 変更後の差 ≤1.5e-6 = ノイズ同オーダー | active (検証) |
| `run_0073c_repeat` | `run_0073_diag_reichardt` の**同一バイナリ再実行** (実行間ノイズ測定) | node は run-to-run 非決定的 (~1e-6) と確定 | active (検証) |
| `run_0074_diag_spalding` | §5.2 ② **診断実験**: `FORGE_WF_INV_LAW=spalding` (逆解き則のみ Spalding、`wf_pk` の $g$ は Reichardt 固定、E3 OFF) | $u_\tau$ 比 **1.0597** / $\tau_w$ 比 **1.1230** (凍結場予測 1.0599/1.1234 と 0.02–0.04%)。デバイス実装 = CPU 後処理 **1.0000** | active (② の正) |
| `run_0075_coupled_reichardt` | §5.2 連成 A/B の**対照**: `run_0063` 場から 20000 step、env OFF | (c) 運動量積分/KS=0.754、(d)=0.764、$\tau_w$=6.80 Pa。全残差 `NOT CONVERGED` / $\theta$,$C_f$ `STEADY` | active (A/B 対照) |
| `run_0076_coupled_spalding` | §5.2 連成 A/B の**本命**: 同 restart + `FORGE_WF_INV_LAW=spalding` | (c)=**0.811**、(d)[整合後処理]=**0.822**、$\tau_w$=**7.24 Pa** (+6.5%)。$\theta$ drift 2.0% で `STEADY` の縁 | active (§5.2 の正) |
| `run_0077_ext_reichardt` | 延長対照 (`run_0075` から 60000 step、env OFF) | (c)=0.7557、(d)=0.7643、$\tau_w$=6.80 Pa。**3 判定量 `STEADY`** (全残差は `NOT CONVERGED`)。**同一 run_dir へ二重投入し履歴末尾に 1 行誤追記** → 隔離・復元し `_audit/` に記録 | superseded (→ `run_0079`) |
| `run_0079_ext_reichardt_clean` | 上の**清潔な再実行** (監査用) | (c)=0.7555、(d)=0.7643、履歴 180000 行・巻き戻り 0。**3 判定量 `STEADY`** | active (**対照の正**) |
| `run_0078_ext_spalding` | 延長本命 (`run_0076` から 60000 step、`FORGE_WF_INV_LAW=spalding`) | (c)=**0.8168** (+8.1%)、(d)=0.8199、$\tau_w$=7.23 Pa。**3 判定量 `STEADY`** (全残差は `NOT CONVERGED`) | active (確定値の正) |
| `run_0080_closure_spalding` | §5.2 ③ **整合実験**: `FORGE_WF_INV_LAW=spalding` + **`FORGE_WF_CLOSURE_LAW=spalding`** ($g$→`wf_pk`,`roK_wf`。$\omega$ 生産・E3・熱側は不変)。`run_0078` 場から 60000 step | (c)=**0.8282**、(d)=0.8367、$\tau_w$=7.35 Pa。3 判定量 `STEADY`。②/③ 比: `wf_pk` −7.8% / `roK_wf` +15.0% | active (§5.2 ③ の正) |
| `run_0062_repdiag` | **代表点幾何診断** (`FORGE_WF_REP_DIAG=1`、`run_0053` 場から 100 step) | `best_cos`=1.0000・接線オフセット/y=9.3e-5 → **平板では代表点が数値許容差内で法線上**。§3.3 の粗場側入力 | active (§3 の正) |
| `run_0060_E3` | **E3: $\omega$ 生産を壁法則整合 $S_{\rm wf}^2$ へ** (`FORGE_WF_OMEGA_SOURCE=1`、`run_0053` 場から 20000 step) | (w) $C_f$/KS = **0.8445** (E3 前 0.7615、pin 0.8594)、$\omega$ 2.62e5→1.35e5、quasisteady `ALL STEADY` | active (**E3 の正**) |
| `run_0060_E3_widiag` | 同 E3 の W–I 実力 + $\omega$ 収支取得 (500 step 再走) | $D/P$=0.895、$\lvert T\rvert/P$=0.105 (E3 前 0.078) | active (診断) |
| `run_0059_omgbudget` | **$\omega$ 項別収支診断** (`FORGE_OMEGA_BUDGET=1`、`run_0053` 場から 500 step) | 第一内層で $D/P$=0.922・$\lvert T\rvert/P$=0.0783・交差拡散 0.0%。全残差 `NOT CONVERGED`/snapshot 2 個だが着目 $\omega$ の変化 0.0048% で瞬時サンプルとして妥当 | active (§4.1 の正) |
| `run_0058_diagoff_regr` | 診断 OFF 経路の回帰 (`FORGE_WI_FORCE_DIAG` 未設定、500 step) | `wi_*` が確保・出力とも除外されることを確認 (通常経路のメモリ/D2H/HDF5 不変) | 削除済み(2026-08-31) (回帰記録) |
| `run_0041_node_yp30_z4_2nd` | **機構対照**: 押し出しを 1 層 → **4 層** (z ノード 2 枚 → 5 枚) にした同一 2 次設定 (`mesh/flat_plate_yp30_z4.geo`, 53775 ノード, 3000 step) | スパン方向モードの成長が **350 分の 1** (step3000 で spread 0.598 m/s vs run_0032 の 211.9)。= 発散は「node の slip BC 欠陥」でも「前縁」でもなく **spanwise が 2 ノードしかないこと**が原因と確定 | 削除済み(2026-08-31) (機構確定の記録) |
| `run_0012_keep_es_ewt_fine` / `run_0013_keep_es_ewt_cont` | **単一スキーム KEEP+ES の RANS 検証 (正)**: 確立済み Cf 基準 run_0009 (EWT 細メッシュ y+0.35) と同一設定で flux のみ KEEP+ES 化、40k+20k ([iddes plan §4.8 設計更新](../../plans/active/turbulence-iddes-sst.md)) | **合格**: implicit cfl20 安定・quasisteady `ALL STEADY`・**Cf ドリフト +0.006〜0.041%/20k で完全収束**。Cf/Schl = **0.88/0.91/0.94** (SLAU 0.91/0.95/0.97 比 **−3.6〜−3.9%** = 本 case の確立済みスキーム間ばらつき node vs cell MUSCL −3.4〜−3.9% と同帯)。**「KEEP+ES で SST-RANS」は成立** | active (KEEP+ES RANS 検証) |

> **node 弱形式境界 (Phase 2)**: node モードで inlet/outlet/slip も ghostless 弱形式化。(5a) 主対流ループを内部+periodic
> のみに、(5b) 全境界を `convectiveFlux_boundary_d` の bvar 弱形式に、(5e) **block-DPLUR Jacobian の境界半割面で
> 退化粘性対角 `viscous_diag` をスキップ** (境界ノードが境界上で `dcc·S≈0`→delta 爆発→dq≈0 で凍結し出口 BL 崩壊・
> 残差プラトーを生んでいた真因)。**出口 BL 崩壊と残差プラトーの双方を解消**。plan
> [discretization-node-boundary-ghostless.md](../../plans/active/discretization-node-boundary-ghostless.md)。コーナー BC
> (壁∩出口) は `ow=ib` マルチマーカ emit + `wall_flag` Dirichlet で別途解決済み (commit 5ce92dc)。
> **過剰乱流を解消** (commit 2b19d1d/0f6d53d): 真因は node で**壁 ω が Dirichlet ピンされていない**こと (ghost BC が
> 壁ノード CV 中心・dcc≈0 退化で omega[ic] を固定できず ω 過小→k 消散不足→過生成)。修正: (i) 境界 massflux 書込
> (出口で k が移流流出できず蓄積するバグ)、(ii) `wall_y_eff` MEAN→MIN (正しい ω_w)、(iii) 壁ノードで omega/roOmega を
> 直接ピン (保存変数 Dirichlet)、(iv) point-implicit で ω 行 decouple (dω=0)、(v) res_roOmega 壁ゼロ化。
> （commit 2b19d1d/0f6d53d で peak mut/μ 367 まで改善するも x=0.89 が ~1.5× 残った）
>
> **さらに wall_dist バグを修正 (commit 8b60f2c) → cell 基準と完全一致**: `calcWallDistance_kdtree` が node で壁点に
> 半割面重心 (壁ノードから x に ~dx/8 ずれ) を使い、近壁 wall_dist が法線距離 y でなく x ずれ (≈1e-4·x, 下流で増大)
> になっていた → ω_w 過小 → 過剰乱流 (下流ほど)。node では壁点に**壁ノード座標**を使うよう修正。結果
> [run_node_sst_final](run_node_sst_final/): peak mut/μ **217** (cell 199), k_wall=0, rms_roUx 0.597→1.5e-4 (3.6桁)。
> **massflux + omega pin + wall_dist の 3 点で過剰乱流 (旧 peak mut/μ 375, Cf≈3×) を解消**し乱流レベルは cell 並みに。
>
> **ただし Cf・境界層厚さの cell 一致には convection スキームを揃える必要がある (追検証)**: `run_node_sst_final` は
> `convMethod 0` (1次風上) のままで、同一測定 (分子勾配 τ_w, δ99=0.99U_e) で cell(MUSCL, run_0009) と比べると
> **Cf −8% / δ99 +13〜20%** と乖離していた (1次の数値拡散で BL が厚くなる + res_60000 で Cf がまだ +20%/30k step
> ドリフト中=未収束)。node を **MUSCL に切替えて Cf 収束まで回した [run_node_sst_muscl_cont](run_node_sst_muscl_cont/)**
> では **Cf 差 −3.4〜−3.9% / δ99 差 −1.0〜−1.5% / u_τ −2.0〜−2.3%** と cell に一致する。u⁺-y⁺ collapse (各 run 自身の
> u_τ で正規化) は 1次でも揃って見えるため、**絶対量 (Cf, δ99) は必ず同一スキームで比較する**こと。比較図:
> `cf_bl_cell_node_muscl.png` (再生成 `tools/compare_cf_bl_cell_node.py`)。残る ~3.5% の Cf 差は median-dual 壁勾配と
> cell の離散差由来で、全 station で一定符号・小さい。

> **SST-DES (DDES) T1-A** の合否: `DESmode:0` ビット不変 + `DESmode:1` で付着 BL が RANS から
> 不変 (Cf 差≪0.1%)。**発達乱流場 (nut/nu≫1) から restart すること** (`run_regr_cf` 等 nut/nu≤1 の
> 未発達場では νt≈0 で rd 小→シールドが効かない)。詳細 plan [turbulence-iddes-sst.md](../../plans/active/turbulence-iddes-sst.md)。

> EWT (enhanced wall treatment) 検証の詳細は plan [turbulence-enhanced-wall-treatment.md](../../plans/accepted/turbulence-enhanced-wall-treatment.md)。
> まとめ図は `ewt_comparison.png` (壁処理進化の Cf 比較) と `ewt_uplus_yplus.png` (u⁺-y⁺ collapse)。
> 再生成: `python3 tools/plot_ewt_summary.py`。

## EWT 検証結果 (`wallTreatmentSST`, 3 帯の y⁺ 非依存性)

`wallTreatmentSST:1` を第一セル高さの異なる 3 メッシュ (y⁺₁≈0.35/10/30) で検証。壁関数 run の Cf は
分子勾配ではなく壁関数整合の **modeled τ_w=ρu_τ²** で測る (`tools/wall_law_modeled.py`、x=0.6)。
実装は 3 段階で進化した:

1. **mode0 (low-Re, 既定 0)**: 粗メッシュで Cf 崩壊 (y⁺30 で相関の 13%)。
2. **ω-blend のみ (k=0 Dirichlet)**: ω は Menter ブレンド+modeled τ_w。改善するが y⁺10 で 1.10 と過大、
   y⁺30 で 0.91。k=0 が近壁 k を潰し非整合 ([architecture-rans-sst] / theory §6.5(b'))。
3. **full WF (本番)**: ω を wall-adjacent セルにピン留め + k zero-gradient + 近壁 P_k を wall-function 値
   `ρu_τ⁴/ν·g(1-g)` に置換。3 点セットで k が平衡値 u_τ²/√β* に自律収束し全帯で相関に乗る。

| 帯 | y⁺₁ | mode0 | ω-blend(k=0) | **full WF** | full WF u_τ |
|---|---|---|---|---|---|
| 細 (壁解像) | ≈0.35 | 0.92 | — | **0.89–0.95** | 2.50–2.69 |
| バッファ | ≈10 | 0.51–0.60 | 1.05–1.14 | **0.95–1.01** | 2.57–2.79 |
| 対数 | ≈30 | 0.13–0.14 | 0.89–0.93 | **0.92–0.97** | 2.51–2.73 |

(数値は x=0.30/0.60/0.89 の Cf/Schlichting。full WF は全帯で相関 ±8% 以内、u_τ≈2.5–2.8 に揃い
**第一セルの y⁺ 位置に依存しない**。細メッシュでは Cf_molec≈Cf_model で wall-resolved 極限に縮退し、
既定 `wallTreatmentSST:0` は従来と完全不変。)

### 収束判定 (VERDICT 引用)

本ケースは block-DPLUR の近似ヤコビアン由来で残差が**構造的にプラトー**する (受理済み本番
`run_0007` も同様)。`tools/check_convergence.py` の VERDICT は full WF 3 run とも:

```
run_0016_ewt_fine_wf / run_0017_ewt_yp30_wf / run_0018_ewt_yp10_wf
  -> NOT CONVERGED (stalled/plateau)   # 各 run の CONVERGENCE_VERDICT.txt
```

で、参照の `run_0007` (120000 step) も `NOT CONVERGED (plateau)`。**運用判定 A** に従い、
「全残差列が低レベルで横ばい (NaN/上昇なし) **かつ** 着目積分量 Cf のドリフト <0.3% (res_15000→20000)」
を満たすため **実用上 OK (発達・定常)** と判断する。リミットサイクルでないことは Cf 定常で担保。
(`rms_roOmega` は ω ピン留めセルの残差を 0 化したため 1e-3 台で正常。)

## node × 2 次精度の発散 — 真因は「spanwise が 2 ノードしかないメッシュ」 (2026-08-12 解決)

`run_0032`–`run_0037` で node × y+30 × 2 次精度がどの設定でも発散した件の決着。
**前縁 (slip→wall 遷移) は無関係**で、真因は**メッシュの作り方**だった。

### 症状の実体 (run_0032 の時系列を再解析)

前セッションが記録した「前縁 x≈19mm の局所暴走 (Ux 199 m/s)」は**最終スナップショットだけ**を
見た記述で、発生源ではない。時系列を追うと:

| step | 上部 slip 行 | 中間場 | 近壁 |
| --- | --- | --- | --- |
| 0 | 0.0023 | 0.0026 | 0.0045 |
| 500 | **0.86** | 0.05 | 0.11 |
| 1000 | **119.9** | 5.9 | 2.1 |
| 2000 | 73.1 | 114.1 | **203.2** |

(数値 = スパン方向非対称 $|U_x(z{=}0)-U_x(z{=}0.01)|$ [m/s])

このメッシュは z 方向 1 層押し出しの**疑似 2D スラブ**なので厳密解は z 非依存でなければならない。
実際には**上部 slip 境界の 1 行**で先に z 対称性が破れ (1 つ下の行は 17 倍小さい)、そこから
中間場 → 近壁へ流れ込む。前縁の overspeed は最後に現れる症状。

### 機構: 2 ノード方向では 2 次 MUSCL の上流差分散逸が厳密に消える

node-centered (median-dual) では 1 層押し出しメッシュの spanwise に**ノードが 2 枚**できる。
z=0 のノード $i$ と z=dz のノード $j$、その間の双対面 (z=dz/2) を考える。GG 勾配は
境界半割面 (slip, owner-state) と内部 z 面の寄与で

$$\left(\frac{\partial\phi}{\partial z}\right)_i=\left(\frac{\partial\phi}{\partial z}\right)_j=\frac{\phi_j-\phi_i}{dz}$$

となり、面への MUSCL 再構成は

$$\phi_L=\phi_i+\frac{\partial\phi}{\partial z}\Big|_i\frac{dz}{2}=\frac{\phi_i+\phi_j}{2},\qquad
\phi_R=\phi_j-\frac{\partial\phi}{\partial z}\Big|_j\frac{dz}{2}=\frac{\phi_i+\phi_j}{2}$$

で **$\phi_L=\phi_R$ が任意の場について厳密に成立**する。つまり z 面のジャンプが恒等的に 0 になり、
**SLAU/上流差分の散逸がこの方向で完全に消える**。2 ノードしかない方向で表現できる唯一のモードは
spanwise 市松 ($\phi_i\neq\phi_j$) なので、**そのモードだけが無減衰**になり丸め誤差から指数成長する。

これが過去の全観測を一度に説明する:

- **1 次では安定** — $\phi_L=\phi_i,\ \phi_R=\phi_j$ でジャンプ $\neq 0$、強く減衰する (`run_0038`)。
- **リミッタで直らない** (`run_0035` barth / `run_0032` venkata) — 再構成値が隣接 2 点の**算術平均**で
  常に min/max 内に入るため、どのリミッタも作動しない (limiter = 1)。
- **CFL で直らない** (`run_0036` cfl_pseudo 0.5) — 散逸の欠落は空間離散の性質で、時間刻みと無関係。
- **LSQ で悪化** (`run_0034`) — 2 点ステンシルでは LSQ も同じ相殺を起こし、GG の fx 平均による
  わずかな平滑化すら失う。
- **cell モードでは起きない** — primal は z 方向に**セル 1 個**しかなく、この方向のモードが存在しない。
- **case/40 ノズルの node 2 次は成功している** — あちらは**真 2D メッシュ (z レベル 1)**。

### 対照実験 (`run_0041`)

押し出しを 1 層 → 4 層 (z ノード 2 枚 → 5 枚) にすると、5 点ステンシルでは市松の
$\phi_L\neq\phi_R$ が回復して散逸が戻る。実測でスパン方向モードの成長は **350 分の 1**
(step 3000 で spread 0.598 vs 211.9 m/s)。→ 「node の slip BC が壊れている」仮説は棄却され、
**spanwise ノード数**が支配因子と確定。

### 対策 (恒久)

**node モードで 2D 問題を解くときは、押し出しメッシュを使わず必ず平面 2D メッシュ
(`Physical Curve` タグ付け, z レベル 1) を作ること。** 本 case では
`mesh/flat_plate_yp30_planar.geo` を新設した (`flat_plate_planar.geo` の y+30 版)。
3D で node を使う場合も、均質方向は**最低 3 ノード (2 層以上)** 確保する。

### 副産物: 2 次化しても Cf 欠損は残る (前セッションの仮説を反証)

| run | 離散化 | 次数 | Cf/真値 @ x=0.30 / 0.60 / 0.89 |
| --- | --- | --- | --- |
| `run_0007` | cell (壁解像 y+0.35) | 2 次 | 1.000 / 1.000 / 1.000 (**真値の定義**) |
| `run_0017_ewt_yp30_wf` | cell (y+30 壁関数) | 2 次 | **1.034 / 1.027 / 1.024** |
| `run_0038` | node (y+30 壁関数) | 1 次 | 0.840 / 0.844 / 0.848 |
| `run_0040` | node (y+30 壁関数) | **2 次** | **0.848 / 0.847 / 0.848** |

前セッションは「Cf の 16% 欠損はおそらく 1 次風上の散逸」と推定していたが、**2 次にしても
Cf は 1% 未満しか動かない**ので**この帰属は誤り**。同じ y+30 メッシュ・同じ壁関数で
**cell は真値の +2〜3%** に入るのに **node は −15%** であり、欠損は
**node 側の壁関数/離散に固有**。$u_\tau$ で見ると node 2.34 / cell 2.586 / 真値 2.580 m/s
(x=0.6) で、node の $u_\tau$ が 9.3% 低い。**次の調査対象はここ** (壁関数の代表点速度の取り方が
有力: node では壁ノードが `nodeWallDirichlet=1` で u=0 に固定されるため、代表点の選択と
壁距離の組が cell と一致しているかの監査が要る)。

### 「真値」の定義と Cf 欠損の出所 (図: `uplus_truth_compare.png`)

再生成: `python3 tools/plot_uplus_truth_compare.py 0.6`

**本 case に実験値はない。**「真値」と呼んでいるものは 2 種類あるので混同しないこと。

1. **内部基準** = forge 自身の cell 壁解像 run (`run_0007`, y+₁≈0.34) の**分子勾配**
   $\tau_w=\mu U_1/y_1$。CFD の自己参照であって外部真値ではない。
2. **経験式** = Schlichting $C_f=0.0592\,Re_x^{-0.2}$。内部基準自身がこれの 0.92 倍
   (前縁からの遷移区間ぶんの virtual-origin オフセット)。

壁関数 run の $C_f$ は分子勾配ではなく**ソルバが課している modeled $\tau_w=\rho u_\tau^2$**
で測る (Reichardt 則の逆解き。`tools/wall_law_modeled.py` と同じ)。

| run | y+₁ | $u_\tau$ | $C_f$ | /内部基準 | /Schlichting |
| --- | --- | --- | --- | --- | --- |
| `run_0007` cell y+0.35 壁解像 | 0.34 | 2.580 | 2.821e-3 | 1.000 | 0.920 |
| `run_0017` cell y+30 壁関数 | 28.80 | 2.586 | 2.896e-3 | 1.027 | 0.944 |
| `run_0038` node y+30 壁関数 1 次 | 52.09 | 2.338 | 2.368e-3 | 0.839 | 0.772 |
| `run_0040` node y+30 壁関数 2 次 | 52.20 | 2.343 | 2.378e-3 | 0.843 | 0.775 |

**(a) 同じ primal メッシュでも node と cell では y+₁ が 2 倍違う**。node の第一 DOF は
**節点そのもの** (= 第一セル高さ $h=3.40\times10^{-4}$ m) にあるが、cell の第一 DOF は
**セル重心** ($h/2$) にある。$y^+$ 比 = $2\times(2.338/2.586)=1.81$ で実測 28.80→52.09 と一致。
median-dual の構造上の性質であり欠陥ではないが、**「同じメッシュ」は「同じ y+」を意味しない**。

**(b) 欠損は壁近傍 2 点に局在する** (図 C, 内部基準 $u_\tau$ で正規化した対数則からのずれ):

| $y^+$ | node $u^+$ | 壁解像 $u^+$ | 差 |
| --- | --- | --- | --- |
| 52 (node の第一 DOF) | 13.91 | 14.79 | **−0.88 (−2.26 m/s)** |
| 80 | 15.15 | 16.07 | **−0.92 (−2.37 m/s)** |
| 130 | 17.41 | 17.46 | −0.05 |
| 200 | 18.90 | 18.67 | +0.23 |
| 400 | 21.11 | 20.70 | +0.41 |

node の第一 DOF では速度が**対数則を 0.92 $u^+$ 下回る** (cell の第一 DOF は逆に +0.34 上回る)。
壁関数はまさにこの点の速度を Reichardt 則に食わせて $u_\tau$ を逆算するので、
**低い速度 → 低い $u_\tau$ (−9%) → 低い $C_f$ (−15%)** が直接の因果。$y^+\gtrsim130$ では
プロファイルは壁解像と一致し (むしろ僅かに fuller = 壁抵抗が小さいことの帰結)、
**外層は健全**。図 A (各 run 自身の $u_\tau$ で正規化) では collapse して見えるので、
**A だけを見て「node は壁法則に乗っている」と判断してはならない**。

**(c) 次に当たる候補 (未検証)**: `viscousFlux_d.cu:174-191` の AddTauWall 再スケールは、
W–I 双対面の**接線** traction を modeled $\tau_w$ に合わせる意図だが、実装は
`tau_x *= scale` と **traction ベクトル全体**を掛けており、法線成分 $T_n$ も同じ倍率
(本 case では $C_{f,\rm model}/C_{f,\rm molec}\approx2.1$) で増幅される。コメントの意図と
コードが食い違っている。壁法線方向に余分な運動量流束を注入するので、第一内点の速度欠損の
候補として最初に潰す価値がある。正しい形は接線成分だけを再スケールして
$\boldsymbol\tau = \hat{\mathbf t}\,\tau_w A + T_n\hat{\mathbf n}$ とすること。

### 追試: 第一層高さは原因ではない — node の乱流量が低いことが真因 (2026-08-12)

上の図の読み方を**訂正**する。**パネル A は collapse しない** (node が上へずれる)。**パネル B が一致する**。
A のずれ量は $u^+_\infty=U_\infty/u_\tau$ が 26.1→28.8 になったもので、**$u_\tau$ 欠損 9% の可視化そのもの**。
したがって「A で揃うから壁法則に乗っている」ではなく、**「B で揃う = 速度場自体は cell と近い、
A でずれる = 報告される $u_\tau$ が低い」**が正しい。

「第一層高さの扱い (node は節点 = $h$、cell は重心 = $h/2$) の違いが分布差を生んでいるのでは」という
仮説を、**node の第一層高さを振って直接検証**した (planar メッシュ 3 種、`mesh/flat_plate_ny{38,45,52}_planar.geo`)。

| run | $y_1$ [m] | y+₁ | $u_\tau$ | /真値 | $C_f$ | /真値 |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0017` cell y+30 壁関数 | 1.700e-4 | 28.8 | 2.586 | 1.002 | 2.896e-3 | 1.027 |
| `run_0042` node ny=52 | **1.768e-4** | 27.1 | 2.343 | 0.908 | 2.377e-3 | 0.843 |
| `run_0040` node ny=45 | 3.400e-4 | 52.2 | 2.343 | 0.908 | 2.378e-3 | 0.843 |
| `run_0043` node ny=38 | 6.588e-4 | 101.5 | 2.350 | 0.911 | 2.393e-3 | 0.848 |

**結論 1: node の壁関数は y+ 非依存である** — y+ を 27→102 (3.7 倍) 振っても $u_\tau$ は 0.3% しか動かない。
壁関数としては正しく機能している。

**結論 2: 第一層高さは原因ではない** — cell (1.700e-4) と node ny=52 (1.768e-4) は**ほぼ同じ物理高さ**なのに
$u_\tau$ が 2.586 vs 2.343 で **9.4% の差がそのまま残る**。仮説は棄却。

**結論 3: 各 run は運動量積分と整合している** (von Kármán $C_f=2\,d\theta/dx$):

| run | $\theta$(0.6) | $C_f$ (運動量積分) | $C_f$ (壁で課した値) | 比 |
| --- | --- | --- | --- | --- |
| cell y+30 壁関数 | 1.033e-3 | 2.813e-3 | 2.896e-3 | 1.030 |
| node ny=52 | 8.745e-4 | 2.403e-3 | 2.377e-3 | **0.989** |
| node ny=45 | 8.967e-4 | 2.410e-3 | 2.378e-3 | **0.987** |
| cell y+0.35 壁解像 | 9.978e-4 | 2.757e-3 | 2.821e-3 | 1.023 |

→ node 側に**偽の運動量ソース/シンクは無い**。node は「$\theta$ が 15% 薄い BL」という
**自己整合な別解**に収束している。壁関数はその状態を正しく報告しているだけ。

**結論 4: 真因は乱流量が低いこと** (x=0.6, BL 内):

| run | peak $\mu_t/\mu$ | $\mu_t/\mu$ @y+100 | peak $k$ |
| --- | --- | --- | --- |
| cell y+0.35 壁解像 | 125.7 | 33.4 | 20.7 |
| cell y+30 壁関数 | 132.7 | 37.9 | 27.9 |
| node ny=52 | **103.5** | **26.7** | **18.3** |
| node ny=45 | 104.8 | 21.4 | 19.1 |

因果連鎖は **$\mu_t$ −20% → 乱流運動量輸送が小さい → BL が薄い ($\theta$ −15%) →
$\tau_w$ が小さい ($u_\tau$ −9% / $C_f$ −15%)**。壁関数は犯人ではなく忠実な報告者である。

**次の調査対象**: node の SST 生産項、特に壁関数の $P_k$ 置換 $\rho u_\tau^4/\nu\cdot g(1-g)$ が
**どの CV 体積に対して積分されているか**。node の第一内層 CV は cell の第一セルと体積が違い、
壁ノード側は半 CV なので、体積重み付けが cell と整合していないと $k$ が系統的に低く出る。
cell の壁関数は自身の壁解像解に対し $k$ を +35% 過大に出す一方、node は −12% 過小である点も
生産項側を示唆する。関連: [`turbulence-node-wall-function-coverage.md`](../../plans/accepted/turbulence-node-wall-function-coverage.md)。

### 【撤回済み】node の壁関数 $P_k$ は構造的に cell の 78% しか注入していない (2026-08-12 起票 → 同日撤回)

> **この節の「真因」「78%」「閉じた因果」は撤回する** (Codex レビュー指摘 + 実験で反証、下の
> 「反証と真の主因」節を参照)。**測った量が間違っていた**: `nodeKwfDirichlet` の既定は `1` で
> ([`solverConfig.cpp:460`](../../solver_density_cuda/input/solverConfig.cpp))、第一内層ノードは
> `roK_wf` に Dirichlet 固定され `res_roK=0` にされる。実測で `roK == roK_wf` が厳密一致し、
> **$\Sigma P_kV$ の 67.7% は k 収支から消される量**だった。以下は記録として残す。

### 【撤回された記述】node の壁関数 $P_k$ の規約差 (2026-08-12)

**先の「B が一致 = 速度場はほぼ同じ」という記述を撤回する。** $\theta$ が 15% 違う以上、速度場は
違う。B が揃って見えたのは、共通 $u_\tau$ で正規化すると $u^+_\infty=U_\infty/u_\tau$ が
**$U_\infty$ 共通ゆえ必ず一致する**うえ、対数軸と 0–32 のレンジで差が潰れているだけである。
実差はパネル C (= B から対数則を引いた拡大図) と $\theta$ に出る。
**パネル A も B も「一致」ではない**というのが正しい。

`wf_pk` と `volume` を出力から直接積分して cell と比較した (x=0.6, spanwise 単位厚, `run_0044` は
現行バイナリで cell y+30 を再走したもの — `run_0017` は旧バイナリで `wf_pk` を出力していない)。

| | CV 数 | $P_k$ 評価距離 $y_1$ | `wf_pk` | 被覆体積 | $\Sigma P_k V$ |
| --- | --- | --- | --- | --- | --- |
| cell y+30 (`run_0044`) | 1 | 1.700e-4 (セル重心) | 3.913e5 | 4.144e-6 | **1.622** |
| node ny=52 (`run_0042`) | 2 | 1.768e-4 (壁ノード→第一内層) | 2.857e5 | 3.370e-6 | **0.963** |

**規約 factor = (壁関数が覆う層厚)/($P_k$ 評価距離 $y_1$)** で規格化すると構造が見える:

| run | 覆う層厚 | $y_1$ | factor |
| --- | --- | --- | --- |
| cell y+30 | 3.400e-4 | 1.700e-4 | **2.000** |
| node ny=52 | 2.765e-4 | 1.768e-4 | **1.564** |
| node ny=45 | 5.317e-4 | 3.400e-4 | **1.564** |
| node ny=38 | 1.030e-3 | 6.589e-4 | **1.564** |

- **cell**: $P_k$ を**セル重心** $y_1=h/2$ で評価し、**第一セル全体** (厚さ $h$) に適用 → factor $h/(h/2)=2.000$。
- **node**: $P_k$ を**壁ノード→第一内層ノード距離** $y_1=h$ で評価し、
  **壁ノードの半 CV** ($h/2$) と**第一内層ノードの CV** ($(h+rh)/2$) の**両方**に適用
  ([`ransWallFunction_d.cu:209,213`](../../solver_density_cuda/cuda_forge/ransWallFunction_d.cu) の
  `wf_pk[ic]` と `wf_pk[irep]`) → 層厚 $(2+r)h/2$ で factor $(2+r)/2\approx1.55$。

**この factor は 3 つの node メッシュすべてで厳密に 1.564 = メッシュ非依存**であり、
$r\to1$ でも 1.5 にしかならない。つまり **node の壁関数生産は構造的に cell の
$1.564/2.000=0.782$ しかない**。

因果の分解 (実測と 2% で一致):

$$\underbrace{0.782}_{\text{規約の差 (}u_\tau\text{ を揃えた場合)}}\times\underbrace{0.744}_{(u_{\tau,\rm node}/u_{\tau,\rm cell})^3\ \text{= 自己増幅}}=0.582\quad\text{vs}\quad\text{実測}\ \Sigma P_kV\ \text{比}=0.594$$

**閉じた因果連鎖**: 生産規約が 22% 小さい → $k$・$\mu_t$ が低い (−20%) → 乱流輸送が小さい →
BL が薄い ($\theta$ −15%) → $\tau_w$ が小さい ($u_\tau$ −9%) → $P_k\propto u_\tau^3$ が更に下がる
→ 合計 41% 減で釣り合う。

**修正方針の候補** (未実装・要判断):

1. **CV 平均生産に統一** — 各 CV の実 $y$ 区間 $[y_a,y_b]$ に対し
   $\bar P_k=\frac{1}{y_b-y_a}\int_{y_a}^{y_b}\rho u_\tau^3/(\kappa y)\,dy$ を使う。離散化非依存に
   なるが cell の既定挙動も変わるため cell の再検証が要る。
2. **node の factor を cell の 2.000 に合わせる** — 壁半 CV と第一内層 CV をそれぞれの重心距離で
   評価する、または `wf_pk` を $2.000/1.564=1.279$ 倍する。cell はビット不変で済む。
3. まず **(2) の 1.279 倍を env ゲートで入れて $C_f$ が cell に乗るかを見る**のが、上記の
   因果連鎖全体を一発で検証する最小実験。

なお cell の factor 2.000 が物理的に厳密という保証はない (Launder–Spalding 系の慣用)。ただし
cell は壁解像解に対し $C_f$ を +2.7% で再現しており**較正済みの規約**なので、当面は
これを基準に node を合わせるのが妥当。

### 反証と、第一内層ノード $\omega$ への局在化 (2026-08-12) — **主因は未分離 (仮説段階)**

> **本節の「真の主因」表現は 2 度目のレビューで再度差し戻した**。$\omega$ 周辺への局在化は妥当だが、
> ①機構の説明が現行実装と違い、②A/B が limiter bypass と交絡している。§2'/§3' に訂正を置く。

前節の「生産規約 78% が真因」を**実験で反証**し、主因を $\omega$ に特定した。

#### 1. 生産規約は主因ではない (反証)

`nodeKwfDirichlet: 0` (第一内層 k のピン解除) の A/B で、**k 収支に実際に入る生産量**を測り直した:

| run | $\Sigma P_kV$ 全体 | うち k がピンされた分 | **実効生産** | 対 cell |
| --- | --- | --- | --- | --- |
| cell y+30 (`run_0044`) | 1.622 | 0 | **1.622** | 1.000 |
| node kwfDir=1 (`run_0042`, 既定) | 0.963 | 0.652 (67.7%) | **0.311** | **0.192** |
| node kwfDir=0 (`run_0045`) | 0.977 | 0 | **0.977** | **0.602** |

**実効生産が 3.1 倍 (0.311→0.977) になっても $u_\tau$ は 2.343→2.356 (+0.55%) しか動かない。**
$C_f$/真値 0.843→0.852。→ **壁関数の $P_k$ 積分量は $\tau_w$ をほとんど支配していない**。
前節の「規約 factor 0.782 × $(u_\tau$比$)^3$ 0.744 = 0.582 が実測 0.594 と 2% 一致」は、
$u_\tau$ 比を説明対象の解から取っている**循環論**であり、因果の証明ではなかった。
`wf_pk` × 1.279 の実験は**やらない** (物理補正でなく cell への経験的チューニングであり、
かつ第一内層が k ピンされている既定では仮説の検証にすらならない)。

#### 2. 主因は第一内層ノードの $\omega$ スパイク

内部スケーリング ($\omega_{eq}=u_\tau^2/\nu/(\sqrt{\beta^*}\kappa y^+)$) で見ると node だけ近壁で $\omega$ 過大:

| $y^+$ | cell $\omega/\omega_{eq}$ | node $\omega/\omega_{eq}$ |
| --- | --- | --- |
| 30 | 1.09 | **2.55** |
| 60 | 1.49 | 1.89 |
| 100 | 1.28 | 1.55 |
| 400 | 1.25 | 1.38 |

壁法線プロファイルを見ると異常が直接見える (`wall_dist` は両者とも座標と厳密一致でバグではない):

| | 壁ノード $y=0$ | 第一内層 $y=1.768$e-4 | 第二 $3.707$e-4 |
| --- | --- | --- | --- |
| node $\omega$ | 1.146e5 (ピン=Menter blend) | **2.645e5** | 9.635e4 |
| node $k$ | 20.72 | **10.46** | 18.27 |

**$\omega$ が壁の値より高く跳ねている** (単調減少すべきところが逆転)。その値 2.645e5 は
**第一層距離の半分**での Menter blend (2.663e5) にほぼ一致する。高い $\omega$ は
$\mu_t=\rho k/\omega$ を下げ、同時に散逸 $\beta^*\rho k\omega$ を上げて $k$ も下げる = 両症状を同時に説明する。

**構造的な非対称が原因と考えられる (仮説)**:

- **cell**: `wf_pk` と $\omega$ ピンが**同じ DOF** (第一セル) にある → $\omega$ は固定なので暴走しない。
- **node**: $\omega$ ピンは**壁ノード**だが `wf_pk` は**壁ノードと第一内層の両方**に入る
  ([`ransWallFunction_d.cu:209,213`](../../solver_density_cuda/cuda_forge/ransWallFunction_d.cu))。
  第一内層は $\omega$ を**解きながら**壁関数生産を受け、そこでは $\nu_t$ が小さい ($\mu_t/\mu\approx2$)
  ため $\omega$ 生産 $\propto P_k/\nu_t$ が過大になる。

#### 3. $\omega$ が主因であることの A/B (`nodeOmegaWfDirichlet`)

| run | $u_\tau$ | /真値 | $C_f$/真値 | $\omega$@第一内層 | $k$@第一内層 | $\mu_t/\mu$@y+100 |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0042` node 既定 | 2.343 | 0.908 | 0.843 | 2.645e5 | 10.46 | 26.7 |
| `run_0045` kwfDir=0 | 2.356 | 0.913 | 0.852 | 2.642e5 | 11.74 | 26.8 |
| `run_0046` **omgWfDir=1** | **2.497** | **0.968** | **0.957** | 1.213e5 | 12.46 | 32.8 |
| `run_0044` cell y+30 (対照) | 2.586 | 1.002 | 1.027 | 1.307e5 | 27.93 | 37.9 |

$\omega$ を固定すると $\omega$ が cell 並み (1.213e5 vs 1.307e5) に戻り、**$C_f$ 欠損の 73% が回復する**
(0.843→0.957, 目標 1.027)。→ **主因は $\omega$**。

ただし **`nodeOmegaWfDirichlet: 1` は解ではない**: case/40 ベルノズルでは同じスイッチが
$\tau_w$ を y+1 真値の **1.237 倍に過大化**して正式化が見送られている
([plan §10c](../../plans/accepted/turbulence-sst-su2-taw-coupling.md))。本ケースでは 0.957 と不足、
ノズルでは 1.237 と過大 = Dirichlet で殴るのは対症療法。**正しい修正は第一内層の $\omega$ ソース項を
cell と整合させること** (壁関数生産を受ける DOF と $\omega$ を解く DOF の関係を cell に合わせる) で、
これは未着手。

#### 4. 準定常確認

残差は全 run で block-DPLUR の構造プラトー (`NOT CONVERGED`) だが、**報告している派生量 $u_\tau$ は
定常**: 最終 2 スナップショット間のドリフトは run_0040 0.001% / run_0042 0.075% / run_0043 0.025% /
run_0044 0.001% / run_0045 0.013% / run_0046 0.070%。

#### 2'. 機構の訂正: $P_\omega$ は `wf_pk` と結合していない

「第一内層が `wf_pk` を受け、$\nu_t$ が小さいので $P_\omega\propto P_k/\nu_t$ が過大」という説明は
**現行実装として誤り**。[`ransSource_d.cu:178`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) は

$$P_\omega=\alpha\rho S^2$$

であり `wf_pk` と独立。`wf_pk` を変えても $P_\omega$ は動かない (それが §1 の反証とも整合する)。

**実測による正しい説明**: スパイクは **$\alpha\rho S^2$ と $\beta\rho\omega^2$ の局所源項平衡そのもの**。
`run_0042` 第一内層 (x=0.6) を再構成すると

| 量 | 値 |
| --- | --- |
| 解像ひずみ $S$ | 1.0750e5 s⁻¹ |
| $F_1$ | 0.621 → $\alpha$=0.5118, $\beta$=0.07796 |
| 局所平衡 $\sqrt{\alpha/\beta}\,S$ | 2.754e5 |
| **実測 $\omega$** | **2.645e5** (比 0.960) |

で、ほぼ完全に源項平衡で決まっている。「第一層距離の半分での Menter blend と一致」は結果的な近さで、
こちらが直接の説明。**cell では同じ DOF の $\omega$ がピンされていて $S$ に応答しない**のが差の出所。

さらに、壁法則整合のひずみ $S_{\rm wf}=u_\tau^2 g/\nu$ と比べると

| $y^+$ | 解像 $S$ | $S_{\rm wf}$ | $S/S_{\rm wf}$ | $(S/S_{\rm wf})^2$ |
| --- | --- | --- | --- | --- |
| 27.1 | 1.075e5 | 5.146e4 | 2.09 | **4.36** |
| 56.9 | 3.192e4 | 1.623e4 | 1.97 | 3.87 |
| 89.6 | 1.531e4 | 9.611e3 | 1.59 | 2.54 |

→ y+30 メッシュの**未解像な離散勾配**が壁法則の 2 倍あり、$S^2$ で 4.4 倍の $P_\omega$ を生んでいる。

#### 3'. A/B の交絡: `nodeOmegaWfDirichlet` は $\omega$ だけを変えていない

[`turbulent_viscosity_d.cu:222`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu):

```cpp
const bool wfPin = (roOmega_wf != nullptr && roOmega_wf[ic] >= 0);
vis_turb[ic] = wfPin ? rho * k_c / w_c
                     : rho * a1 * k_c / max(a1 * w_c, S_mag * F2);
```

このスイッチは (a) 第一内層 $\omega$ のピン、(b) **SST shear limiter $\max(a_1\omega, SF_2)$ の迂回**、
(c) $\mu_t=\rho k/\omega$ への切替を**同時に**行う。したがって $C_f$ 0.843→0.957 は
「$\omega$ 低下」と「limiter bypass」の**合成効果**であり、$\omega$ 単独の効果ではない。
case/40 で $\tau_w$ が 1.237 倍に過大化したのも、既存調査では主に limiter bypass による $\mu_t$ 過大とされている
([plan §10c](../../plans/accepted/turbulence-sst-su2-taw-coupling.md))。

#### 4'. 現時点の正確な結論

**第一内層ノードの「k 半拘束ピン + $\omega$ free + 解像 $S$ による $P_\omega$」の周辺に問題があり、
`nodeOmegaWfDirichlet` 経路を通すと大幅改善する。ただし効いているのが $\omega$ 状態なのか
strain limiter なのか両方なのかは未分離。** 「真因確定」とは言わない。

分離実験と恒久修正候補は [plan: turbulence-node-wf-omega-source.md](../../plans/accepted/turbulence-node-wf-omega-source.md) に起票した。

## 検証基準の改訂: 外部相関 $C_f(Re_\theta)$ (Kármán–Schoenherr) を主基準に (2026-08-12)

**結論から**: 外部相関に対し **cell は 4–5% 低い / node は 22–23% 低い**。
これまで「cell = 真値」として node を 0.843 と評価していたが、**cell 自身が外部相関を外している**ので、
その枠組みでは node の欠損量を正しく測れていなかった。

$$C_{f,\mathrm{KS}}=\left[17.08(\log_{10}Re_\theta)^2+25.11\log_{10}Re_\theta+6.012\right]^{-1},\quad
Re_\theta=\frac{\rho_\infty U_\infty\theta}{\mu_\infty}$$

(NASA TMR flat-plate validation と同じ。<https://tmbwg.github.io/turbmodels/flatplate_val.html>)

**呼称**: cell 壁解像解は今後「**内部基準**」と呼び「真値」と呼ばない。K–S も厳密解ではないので
「**外部相関**」と呼ぶ。

### 前提ゲート (これを満たす station だけで判定する)

| ゲート | 結果 | 判定 |
| --- | --- | --- |
| **ZPG** | 加速パラメータ $K=(\nu/U_e^2)dU_e/dx=2.1\times10^{-9}$ ($\ll3\times10^{-6}$)。$\Delta U_e/U_e=+0.6\%$, $\Delta p/q=-1.4\%$ (x=0.2→0.9) | **合格 (実質 ZPG)** |
| **発達乱流** | $\theta(x)$ 単調増加は $x\gtrsim0.009$ (cell 壁解像) / $0.024$ (cell wf) / $0.052$ (node)。ただし **BL 内 peak $\mu_t/\mu$ が自由流値 65.9 を超えるのは $x\approx0.25$–$0.37$** | **$x\ge0.3$ を条件とする** |
| **K–S 有効域** ($4000<Re_\theta<13000$) | $Re_\theta$ は $x=0.90$ でも **5571 (node) – 6427 (cell)**。**上限 13000 に全く届かない**。$Re_\theta>4000$ は $x\gtrsim0.60$ (cell) / $x\gtrsim0.75$ (node) | **下端のみカバー** |
| **収支の信頼域** | $C_f/(2d\theta/dx)$ は $x=0.45$–$0.75$ で 0.99–1.05、$x\ge0.9$ で 0.49–1.10 に崩れる (局所 fit の station 不足+出口影響) | **$x\le0.8$ に限る** |

→ **判定帯は $x\in[0.60,0.80]$** ($Re_\theta\approx3900$–$5500$)。**K–S を唯一の合否基準にはできない**。

### $C_f$ の 4 定義 (混ぜないこと)

1. **wall-resolved**: $\tau_w=\mu\,\partial u/\partial y$ (第一 DOF の分子勾配)
2. **wall-function**: Reichardt 逆解きの $\tau_w=\rho u_\tau^2$ (ソルバが運動量式に課す値と整合)
3. **`twall` 出力**: AddTauWall で 2. の目標値へ**再スケール済み**なので**独立検証にならない** (使わない)
4. **運動量積分**: $C_f=2\,d\theta/dx$

本表の $C_f$ は **2.** で統一 (壁解像 run では 1. と 2. がほぼ一致する: x=0.6 で 2.9149e-3 vs 2.8411e-3)。
**node の W–I 双対面に実際に加わった接線力からの $C_f$ は保存出力から取得できず「未確認」**
(取得には plan §4.2 の診断追加が必要)。

### 主表 (判定帯を含む 5 station、$C_f$ は wall-function 定義)

$d\theta/dx$ は $x\pm0.08$ m の局所 2 次 fit (使用 station 7–26 個)。$\theta$ は全域積分。前縁 $x<0.10$・出口 $x>0.95$ 除外。

| run | $x$ | $Re_\theta$ | K–S域 | $C_f$ | $C_f/C_{f,\mathrm{KS}}$ | $C_f/(2d\theta/dx)$ |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0007` cell 壁解像 (内部基準) | 0.596 | 4614 | ○ | 2.915e-3 | **0.954** | 1.102 |
| | 0.745 | 5499 | ○ | 2.820e-3 | **0.956** | 1.069 |
| `run_0044` cell y+30 wf | 0.596 | 4599 | ○ | 2.933e-3 | **0.960** | 1.026 |
| | 0.745 | 5541 | ○ | 2.834e-3 | **0.962** | 1.035 |
| `run_0042` node h=1.77e-4 (y+27) | 0.602 | 3905 | △ | 2.406e-3 | **0.762** | 1.023 |
| | 0.753 | 4696 | ○ | 2.332e-3 | **0.766** | 0.993 |
| `run_0040` node h=3.40e-4 (y+52) | 0.602 | 4009 | ○ | 2.408e-3 | **0.766** | 1.043 |
| `run_0043` node h=6.59e-4 (y+102) | 0.602 | 4184 | ○ | 2.426e-3 | **0.779** | 1.049 |
| `run_0045` node kwfDir=0 | 0.753 | 4686 | ○ | 2.361e-3 | **0.775** | 0.999 |
| `run_0046` node omgWfDir=1 | 0.602 | 4464 | ○ | 2.734e-3 | **0.889** | 1.044 |

- **cell 壁解像と cell y+30 wf は外部相関に対し 0.954–0.962** で互いに 0.6% 以内。
- **node は 0.762–0.779** で、第一 DOF の y+ を 27→102 と 3.7 倍振っても **±0.011 に収まる**
  (壁関数としての y+ 非依存は成立している)。
- `nodeOmegaWfDirichlet:1` で 0.889 まで戻るが、**これは $\omega$ ピンと shear limiter 迂回の合成効果**
  (plan §1.3) なので単独因子の証拠にならない。

### $\theta$ 積分の規約と感度

$\delta_{99}$ は BL 厚さの定義であって $\theta$ の打切り位置ではないので、**積分核が十分ゼロになる
外部流まで (領域上端 $y=0.2$) 積分**するのを主値とする。node は壁点 ($y=0,u=0$) をそのまま含め、
cell は $(y,u)=(0,0)$ を補う。上端感度 (x=0.6):

| 上端 | cell 壁解像 | cell y+30 | node h=1.77e-4 | node omgWfDir=1 |
| --- | --- | --- | --- | --- |
| $y=0.2$ (全域) | 1.0373e-3 | 1.0341e-3 | 8.7775e-4 | 1.0036e-3 |
| $y=0.1$ | 1.0147e-3 | 1.0183e-3 | 8.7507e-4 | 9.9387e-4 |
| $2\delta_{99}$ | 9.9807e-4 | 1.0171e-3 | 8.7426e-4 | 9.8944e-4 |
| **最大偏差** | **3.79%** | 1.64% | 0.40% | 1.41% |

**「上端を下げても値が変わらない」は厳密には成立しない** (自由流の微小非一様性を長い距離で拾うため。
壁解像 run で最大 3.8%)。ただし $C_{f,\mathrm{KS}}$ の $Re_\theta$ 感度は
$\partial\ln C_f/\partial\ln Re_\theta\approx-0.2$ なので、**比への影響は $\pm0.8\%$** に留まる。

### Schlichting 係数 0.0592 と 0.0576 の整理

リポジトリ内で混在していた 2 つは**同じ式ではない**。1/7 乗則から自己整合に導出すると:

$u/U=(y/\delta)^{1/7}$ ($\theta/\delta=7/72$) + Blasius $\tau_w=0.0225\rho U^2(\nu/U\delta)^{1/4}$
+ 運動量積分 → $\delta/x=0.3707\,Re_x^{-1/5}$、$C_f=\mathbf{0.0577}\,Re_x^{-1/5}$

- **0.0576** = $\delta/x=0.37Re_x^{-1/5}$ と**対になる 1/7 乗則の自己整合値** (`tools/compare_cf_bl_cell_node.py`)
- **0.0592** = 同じ 1/7 乗則族の**別の経験フィット** (Schlichting, $5\times10^5<Re_x<10^7$)
  (`tools/postprocess_wall_law.py`, `wall_law_modeled.py`, `plot_uplus_truth_compare.py`)
- 差は **2.8%** で、数 % を論じる本 case では**混用不可**。
- **提案 (未実施)**: 補助基準は 0.0592 に統一し、0.0576 を使う箇所はラベルを
  「1/7 乗則自己整合値」に変える。ツール修正は別途。

### 収束・準定常 (使用 snapshot と VERDICT)

`check_convergence.py` は **8 run すべて `NOT CONVERGED (stalled/plateau)`** — これは本 case の
block-DPLUR 近似ヤコビアン由来の構造プラトーで既知 (受理済み `run_0007` も同じ)。したがって
**派生量の準定常で判定する**。`check_quasisteady.py` に `theta` / `cf_retheta` を追加した
(2026-08-12。これが唯一のツール変更)。

```bash
python3 solver_density_cuda/tools/check_quasisteady.py <run_dir> \
    --quantity theta,cf_retheta --cf-x 0.6 --tail 0.5 --drift 0.02 --osc 0.02
```

| run | 使用 snapshot | VERDICT | $\theta$ drift/tail | $C_f/C_{f,\mathrm{KS}}$ drift/tail |
| --- | --- | --- | --- | --- |
| `run_0007` | `res_120000.h5` | TRANSIENT-UNSETTLED | 1.6% | 0.1% |
| `run_0017` | `res_10000.h5` | STEADY | 0.1% | 0.0% |
| `run_0044` | `res_5000.h5` | STEADY | 0.0% | 0.0% |
| `run_0040` | `res_40000.h5` | STEADY | 0.0% | 0.0% |
| `run_0042` | `res_20000.h5` | TRANSIENT-UNSETTLED | 1.4% | 0.1% |
| `run_0043` | `res_20000.h5` | STEADY | 1.3% | 0.5% |
| `run_0045` | `res_20000.h5` | STEADY | 0.4% | 0.0% |
| `run_0046` | `res_20000.h5` | STEADY | 1.0% | 0.8% |

**報告する比 $C_f/C_{f,\mathrm{KS}}$ は全 run で drift ≤0.8%/tail** (窓 = 末尾 50%)。
`run_0007`/`run_0042` の UNSETTLED は **$\theta$ 側**の残ドリフト (1.4–1.6%) が
系列の極値を末尾に持つため。$\theta$ の 1.6% は比へ 0.3% しか効かないので、
**比の判定には影響しない**が、$\theta$ 単独を報告するときは注記が要る。

### NASA TMR ケースとの条件差 (直接比較は不可)

| 項目 | case/26 | NASA TMR flat plate |
| --- | --- | --- |
| $M$ | 0.2 | 0.2 |
| $Re$ | $Re/m=4.464\times10^6$, $Re_{x,\max}=4.46\times10^6$ | $Re_L=5\times10^6$ |
| SST variant | forge 実装 ($P_\omega=\alpha\rho S^2$ = **SST-2003 の誤記形**) | SST-Vm 等 |
| **freestream $\mu_t/\mu$** | **65.9** (`k=0.3, omega=300`) | **≈0.009** |
| transition model | 無し (モデル自身の activation 位置あり) | 無し |
| 前縁 | $x<0$ に slip 助走区間 0.1 m | 同様 (symmetry) |
| BC | 入口 全圧全温 / 出口 静圧 / 上端 slip | riemann / 上端 farfield 等 |
| $Re_\theta$ 範囲 | ~2200–6400 | K–S 主比較域を広くカバー |

**freestream $\mu_t/\mu$ が 4 桁違う**ため、**NASA の SST 数値解との直接比較は apples-to-apples ではない**。
K–S 相関との物理比較 (本節) と、NASA SST 数値解による実装 verification は分けて扱う。
後者には流入乱流量を合わせた別ケースが要る (未実施)。

## 【重要】壁処理 × 離散化 × ソルバの分解 — node 離散化は無罪、SU2 も同じ 0.75 (2026-08-12)

「node 壁関数に欠陥がある」という作業仮説を、**壁関数を外した node** と **同一メッシュの SU2**
で切り分けた。結論は 2 つとも仮説に不利である。

### 決定表 (外部相関比 $C_f/C_{f,\mathrm{KS}}$、有効帯)

| 構成 | 壁処理 | 離散化 | $C_f/C_{f,\mathrm{KS}}$ (x=0.60 / 0.75) | 収支 $C_f/(2d\theta/dx)$ |
| --- | --- | --- | --- | --- |
| `run_0007` forge cell y+0.34 | **壁解像** | cell | 0.954 / 0.956 | 1.069 |
| **`run_0050` forge node y+0.7 (現行 binary)** | **壁解像** | **node** | **0.943 / 0.944** | **1.036** |
| `run_node_sst_muscl_cont` forge node y+0.7 (旧 binary) | 壁解像 | node | 0.943 / 0.947 | 1.033 |
| `run_0044` forge cell y+29 | 壁関数 | cell | 0.960 / 0.962 | 1.035 |
| `run_0042` forge node y+27 | 壁関数 | node | **0.762 / 0.766** | 0.993 |
| **`run_0048` SU2 y+27 (forge と同一メッシュ)** | **壁関数** | **node** | **0.754 / 0.755** | 0.775 |
| `run_0047` SU2 y+27 | low-Re (不適) | node | 0.365 / 0.365 | 0.181 |

### 結論 1: node 離散化そのものは健全

**壁関数を外して壁解像で解くと node は 0.943–0.944** で、cell 壁解像 (0.954–0.956) と **1.1% しか違わない**。
しかも**運動量収支は node の方が良い** (1.036 vs 1.069)。旧 binary の既存 run とも 0.943 で一致するので、
本セッションの一連の修正とも独立している。
→ **median-dual 離散化が $C_f$ を落としているのではない。** 欠損は**壁関数経路に限定**される。

### 【撤回】結論 2: SU2 も同じ 0.75 を出す — forge node 壁関数は「異常」ではない

> **この結論は撤回する (2026-08-12、同日の追検証)。** SU2 壁関数版 (`run_0048`) は
> **運動量積分を満たしていない** (壁 $C_f$ が BL 発達が要求する値の 0.79 倍) ため、
> forge node の値を裏づける証拠として使えない。残差も rms −4.4 止まり (壁解像版は −12.4)。
> 詳細は下記「運動量収支による再評価」節。以下は撤回された記述として残す。

### 【撤回された記述】SU2 も同じ 0.75 を出す

SU2 v8.5.0 は **vertex-centered median-dual = forge node と同じ DOF 配置**である。
`STANDARD_WALL_FUNCTION` を有効にして **forge `run_0042` と完全同一のメッシュ・同一 BC・
同一物性・同一流入乱流**で回すと **0.754–0.755** で、**forge node の 0.762–0.766 と 1% 以内で一致**する
(SU2 自身の `Skin_Friction_Coefficient` 出力で測ると 0.791–0.799)。

→ **「forge の node 壁関数に固有の欠陥がある」という仮説は成立しにくい。**
独立した検証済みソルバが同じ離散化・同じ条件で同じ答えを出す。
むしろ**外れているのは forge の cell 壁関数 (0.960)** の方で、これが壁解像値に近いのは
偶然の可能性がある (cell の第一 DOF はセル重心なので、vertex-centered とは壁関数の当て方が構造的に違う)。

**未解明として残るもの** (どちらも 0.75 を「正しい」とする根拠にはならない):

- 本ケースの **freestream $\mu_t/\mu=65.9$** は NASA 標準ケース (~0.009) より 4 桁大きい。
  両ソルバで共通の条件なので相対比較は成立するが、**壁関数がこの高い自由流乱流下で
  どう振る舞うべきか**の外部基準は無い。
- SU2 壁関数版の**運動量収支が閉じない** (0.775)。forge node は 0.993 で閉じている。
  同じ $C_f$ に至りながら収支の質が違うので、両者が同じ理由で 0.75 なのかは未確認。

### 結論 3: SU2 low-Re も y+27 で崩壊する (forge 固有ではない)

`run_0047` (壁関数なし、y+27) は $C_f/C_{f,\mathrm{KS}}=0.365$、運動量収支 **0.181**、
$\mu_t/\mu$ ピーク 340–452、$\theta$ が forge の 2.6 倍。**解として使えない**。
forge cell の y+30 low-Re が崩壊した記録 (`run_0010`, Cf/Schl 0.13) と同種で、
**「low-Re 壁処理を y+≫1 メッシュに当てると壊れる」は forge 固有ではない**ことの確認になる。

### SU2 クロスチェックの設定 (再現用)

`run_0047`–`run_0049`。手順は [`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md)。
メッシュは **forge と同一の `.geo`** から `gmsh -2 ... -format su2` で生成 (交絡を避ける)。

| forge 設定 | SU2 設定 |
| --- | --- |
| `Pt=100000, Tt=288.15` | `MARKER_INLET= ( inlet, 288.15, 100000.0, 1,0,0 )` + `INLET_TYPE= TOTAL_CONDITIONS` |
| 背圧 `Ps=97250` | `MARKER_OUTLET= ( outlet, 97249.67 )` |
| `M=0.2` | `MACH_NUMBER= 0.2`, `FREESTREAM_PRESSURE= 97249.67`, `FREESTREAM_TEMPERATURE= 285.8631` |
| `visc=1.8e-5` (定数) | `VISCOSITY_MODEL= CONSTANT_VISCOSITY`, `MU_CONSTANT= 1.8E-5` |
| `thermCond=0.0251` (定数) | `CONDUCTIVITY_MODEL= CONSTANT_CONDUCTIVITY`, `THERMAL_CONDUCTIVITY_CONSTANT= 0.0251` |
| 流入 `k=0.3, omega=300` | `FREESTREAM_TURBULENCEINTENSITY= 0.006598`, `FREESTREAM_TURB2LAMVISCRATIO= 65.853` (等価換算) |
| `wall` (断熱 no-slip) | `MARKER_HEATFLUX= ( wall, 0.0 )` |
| `sym`/`top` (slip) | `MARKER_SYM= ( sym, top )` |
| `wallTreatmentSST: 1` | `MARKER_WALL_FUNCTIONS= ( wall, STANDARD_WALL_FUNCTION )` |

評価は forge と**同一の後処理**(同じ $\theta$ 積分規約・同じ Reichardt 逆解き $C_f$・同じ K–S 比)で行う。

### 速度 (参考、公平比較ではない)

| | 反復 / 時間 | 1 反復 |
| --- | --- | --- |
| forge `run_0042` (GPU RTX 3060) | 20000 step / **38.8 s** | **1.94 ms** (内部反復 20 回込み) |
| SU2 `run_0047` (CPU 8 スレッド) | 40000 iter / 314 s | 7.86 ms |

同一メッシュ (12428 ノード) で forge が約 4 倍速いが、**GPU vs CPU なので効率の比較にはならない**。

## 運動量収支による再評価 — SU2 壁関数版は棄却、node 壁解像は SU2 と一致 (2026-08-12)

### 収支の取り方 (簡易形は不十分だった)

当初は ZPG 極限 $C_f/2=d\theta/dx$ で評価していたが、本ケースは弱い加速があるので**圧力勾配項を
落とせない**。von Kármán 運動量積分の完全形を使う:

$$\frac{C_f}{2}=\frac{d\theta}{dx}+\frac{\theta}{U_e}\left(H+2-M_e^2\right)\frac{dU_e}{dx},
\qquad H=\frac{\delta^*}{\theta},\quad
\delta^*=\int_0^{y_{top}}\left(1-\frac{\rho u}{\rho_e U_e}\right)dy$$

圧力勾配項は $d\theta/dx$ の **2.0–3.1%** あり、収支比を系統的に ~2.5% 押し下げる
(例: forge cell 壁解像 1.102→1.076)。$d\theta/dx$ と $dU_e/dx$ はいずれも $x\pm0.08$ m の局所 2 次 fit。

さらに、**壁 $C_f$ と「BL 発達が要求する $C_f$」($=2(d\theta/dx+\text{PG})$) の両方を外部相関で
正規化**すると、どちらが外れているかが分離できる。

### 決定表

| 構成 | 残差 | 壁 $C_f$/KS | 積分 $C_f$/KS | 壁/積分 |
| --- | --- | --- | --- | --- |
| `run_0007` forge cell 壁解像 | plateau | 0.954 / 0.956 | 0.886 / 0.919 | 1.076 / 1.039 |
| **`run_0050` forge node 壁解像** | plateau | **0.943 / 0.944** | 0.936 / 0.937 | **1.008 / 1.008** |
| **`run_0049` SU2 壁解像 (y+0.7)** | **rms −12.4** | **0.938 / 0.944** | 0.939 / 0.950 | **0.999 / 0.995** |
| `run_0044` forge cell 壁関数 | plateau | 0.960 / 0.962 | 0.958 / 0.956 | 1.001 / 1.006 |
| `run_0042` forge node 壁関数 | plateau | 0.762 / 0.766 | 0.759 / 0.790 | 1.003 / 0.970 |
| **`run_0048` SU2 壁関数** | **rms −4.4** | 0.793 / 0.796 | **1.003 / 1.005** | **0.790 / 0.792** |

### 結論 A: SU2 壁関数版は棄却する

**壁 $C_f$ が BL 発達の要求値の 0.790 倍** = 定常 2D 解なら恒等的に閉じるはずの運動量積分を
**21% 破っている**。$\theta$ の積分上端感度は 0.67%、自由流 $u/u_e=0.99998$、$U_e(x)$ は forge と
ほぼ同一なので、**後処理側の問題ではない**。残差も rms −4.4 止まり (壁解像版は −12.4) で、
**`STANDARD_WALL_FUNCTION` 経路が収束しきっていない可能性が高い**。
なお `CD` は iter 3333 以降 0.002611 で 6 桁不変なので「積分量が動いている」という意味では収束済み。

→ **「SU2 も 0.75 を出すから forge node は異常でない」という結論は撤回する。**

### 結論 B: node 壁解像は SU2 とほぼ完全に一致 (node 離散化の検証)

**`run_0050` forge node 壁解像 0.943/0.944 と `run_0049` SU2 壁解像 0.938/0.944 が 0.5% 以内で一致**し、
**両者とも運動量収支がほぼ完璧** (1.008 / 0.999)。SU2 側は rms −12.4 まで落ちた素性の良い解である。

→ **median-dual (vertex-centered) 離散化は検証済みソルバと一致する。** node 離散化は無罪、という
結論は**むしろ強化された** (今度は良く収束した SU2 解が根拠)。

### 結論 C: 本ケースの「正しい」壁解像値は $C_f/C_{f,\mathrm{KS}}\approx0.94$

forge node 0.943 / SU2 0.938–0.944 / forge cell 0.954。**3 者が 0.94–0.95 に収まる**。
K–S より 5–6% 低いのは、本ケースの **freestream $\mu_t/\mu=65.9$** (NASA 標準の 4 桁大) が
効いている可能性が高い (未検証)。
**したがって node 壁関数の 0.762 は、この 0.94 に対して −19% の欠損**である。

### 結論 D: forge cell 壁解像の収支が 3 者中で最も悪い

`run_0007` の壁/積分は 1.076 / 1.039 で、node (1.008) や SU2 (0.999) より劣る。
壁 $C_f$ が 0.954 と高めなのは**収支が閉じていないため**の可能性がある
(積分 $C_f$ で見ると 0.886–0.919 で node/SU2 の 0.936–0.950 より低い)。
**内部基準として使ってきた `run_0007` 自身に疑いがある** — 別途要確認。

## 追検証: SU2 の壁 $C_f$ 出力は何か / どの定義を信じるか (2026-08-12)

「SU2 の壁 $C_f$ は壁関数を反映していないのでは。積分 $C_f$ が正しいなら流れ場は正しいのでは」
という指摘を検証した。**指摘は概ね当たっている**が、単純に「出力だけの問題」でもない。

### 4 定義の突き合わせ (x=0.6、すべて $C_f/C_{f,\mathrm{KS}}$)

| 構成 | (a) 解像勾配 | (b) ソルバ出力 | (c) 運動量積分 | (d) Reichardt 逆解き | 最大差 |
| --- | --- | --- | --- | --- | --- |
| **`run_0049` SU2 壁解像 (rms −12.4)** | **0.945** | **0.938** | **0.939** | **0.947** | **1.0%** |
| `run_0050` forge node 壁解像 | — | — | 0.936 | 0.943 | 0.7% |
| `run_0007` forge cell 壁解像 | — | — | 0.886 | 0.954 | **7.4%** |
| `run_0048` SU2 壁関数 | 0.373 | **0.793** | **1.003** | 0.754 | **86%** |
| `run_0042` forge node 壁関数 | 0.369 | — | **0.759** | **0.762** | 62% |
| `run_0044` forge cell 壁関数 | — | — | 0.958 | 0.960 | 0.2% |

### 検証 1: 後処理パイプラインは正しい

**良く収束した壁解像解 (`run_0049`, rms −12.4) では 4 定義が 1.0% 以内で一致**する。
$\theta$ 積分・$d\theta/dx$ の局所 fit・圧力勾配項・Reichardt 逆解きのどれにも系統誤差が無い実証。

### 検証 2: SU2 の壁 $C_f$ 出力は解像勾配ではない (指摘は正しい)

SU2 壁関数版で (a) 解像勾配 0.373 に対し (b) 出力 0.793 = **2.1 倍**。
→ **SU2 の `Skin_Friction_Coefficient` は壁モデルを反映している**。
「(b) が低いから SU2 の壁応力が低い」と読んだのは誤りだった。

### 【撤回】検証 3: SU2 の外側の場は正しい

> **撤回する (2026-08-12, レビュー指摘)。** `run_0048` は**全残差で見ると強い未収束**である:
> `rms[Rho]` −4.43 だけを見ていたが、**`rms[RhoE]` = +1.03 (残差 10)**、`rms[w]` −0.29、
> `rms[RhoV]` −2.13 ([history.csv](run_0048_su2_sst_wf_ny52/history.csv) 最終行)。
> 対して壁解像 `run_0049` は全列 −5.0〜−12.6。
> **AGENTS.md の「`rms_ro` だけで判断しない」を SU2 に対して破っていた。**
> 未収束場では定常運動量積分式に残差項が欠けるので、(c)=1.003 が K–S に近いことは
> 「場が正しい」根拠にならない。正確な表現は
> **「現スナップショットの境界層発達率はたまたま K–S に近いが、壁応力との収支が 21% 閉じず、
> 定常解として使用不能」**。`CD` の停止だけでは $k,\omega,\rho E$ の未収束を覆せない。
>
> **また SU2 に Reichardt 逆解き (d) を当てた比較も成立しない。** SU2 の
> `STANDARD_WALL_FUNCTION` は **Nichols–Nelson (2004) の圧縮性 Spalding/White–Christoph 型**
> (`CNSSolver.cpp` の Newton 反復、Crocco–Busemann 込み) で Reichardt ではない。
> 「SU2 の第一ノード速度と壁応力が不整合」という記述は**別の壁法則を比べていたための誤り**。
>
> 以下は撤回された記述として残す。

### 【撤回された記述】SU2 の外側の場は正しい / forge node は場そのものが薄い

**SU2 の (c) 運動量積分 = 1.003** (x=0.3/0.6/0.9 で 0.961/1.003/1.034) →
**SU2 の境界層は外部相関どおりに育っている**。
一方 **forge node 壁関数の (c) は 0.759** → **BL 発達そのものが 24% 遅い**。

→ **この 2 つは「同じ 0.75」ではない。** SU2 は場が正しく $C_f$ の定義がずれているだけ、
forge node は**場そのものが薄い**。前節の「SU2 も 0.75 だから forge node は異常でない」は
**二重に誤り**だった (SU2 側の未収束と、$C_f$ 定義の取り違え)。

### ただし SU2 壁関数版も内部整合していない

(a) 0.373 / (b) 0.793 / (c) 1.003 / (d) 0.754 と **4 定義が 86% 開いている**
(壁解像解では 1.0%)。特に **(d) 第一 DOF 速度から逆解きした値が (c) と整合しない** =
SU2 の第一ノード速度は、(c) が示す壁応力を対数則に当てはめた値より低い。
残差 rms −4.4 止まりと合わせ、**「場は正しいが壁近傍が未整合」**。
定量比較の相手にするには収束させ直すか壁応力の実装確認が要る。

### 【訂正】(c) は「定義非依存の壁摩擦」ではない

$2(d\theta/dx+\mathrm{PG})$ は壁出力に依存しない有用な推定だが、**定常性・境界層近似・
$U_e,\rho_e$ の取り方・積分上端・微分 fit に依存する**。主判定は次の 3 点セットにする:

1. **実際に運動量残差へ入った W–I 接線力** ← **取得済 (2026-08-12, §4.2)。`tools/cf_retheta_analysis.py --wi-force` で再生成可**
2. それと運動量積分値の収支
3. 両者と外部相関 K–S との比較

### いま信じるべき数字 (収支診断として)

| | (c) $C_f/C_{f,\mathrm{KS}}$ |
| --- | --- |
| SU2 壁解像 (rms −12.4) | **0.939** |
| forge node 壁解像 | **0.936** |
| forge cell 壁関数 | 0.958 |
| forge cell 壁解像 | 0.886 ← 3 者中の外れ値。要確認 |
| **forge node 壁関数** | **0.759** ← 明確な外れ値 |
| SU2 壁関数 | 1.003 (ただし内部不整合) |

**壁解像基準は 0.94 前後 (SU2 0.939 / forge node 0.936)。forge node 壁関数の 0.759 は −19%、
しかも $\theta$ 発達そのものが遅い**ので、壁応力の出力だけの問題ではない。

## 現状認識の確定版 (2026-08-12 レビュー反映)

> **case/26 では forge node 壁解像が SU2 壁解像と 0.5% 程度で一致するため、node 離散化一般は
> 主因ではない。forge node 壁関数では、課した低い壁応力と境界層の遅い発達が互いに整合しており、
> 欠損は壁関数経路に局在する。ただし、狙った壁応力が W–I 残差へ正しく作用しているか、また
> $\omega$ 状態・limiter・源項のどれが低摩擦解を作るかは未分離である。SU2 壁関数 run は
> 全残差未収束のため比較証拠に使用しない。**

### 確定してよいこと

| # | 内容 | 根拠 |
| --- | --- | --- |
| 1 | forge node 壁解像 (`run_0050`) と SU2 壁解像 (`run_0049`) が $C_f/C_{f,\mathrm{KS}}\approx0.94$ でほぼ一致 | (a)(b)(c)(d) が 1% 以内、SU2 は全残差 −5〜−12.6 |
| 2 | したがって −19% 欠損は **median-dual/node 離散化一般ではなく forge の node 壁関数経路に局在** | 1 の対偶 |
| 3 | SU2 low-Re を y+27 に適用した `run_0047` は比較対象にならない | 収支 0.18、$\mu_t/\mu$ 340–452 |
| 4 | SU2 の `Skin_Friction_Coefficient` は**壁関数応力を反映している** | `CFVMFlowSolverBase.inl` で `scale = Tau_Wall / WallShearStress` により接線応力を再スケールしてから出力 |

### 使ってはいけないこと

- **`run_0048` (SU2 壁関数) を定量比較の相手にしない** — 全残差未収束 (`rms[RhoE]`=+1.03)。
- **SU2 に Reichardt 逆解きを当てない** — SU2 は Nichols–Nelson。整合確認には SU2 自身の
  `Y_Plus` (壁ノードに出力。x=0.6 で **27.21**)・$U_\tau$・Spalding 則を使う。
- **(c) 運動量積分を「定義非依存の $C_f$」と呼ばない** — 収支診断である。

### 次にやること (この順)

1. **実 W–I 接線力の診断** — 「狙った摩擦力が W–I 間で本当に作用しているか」。plan §4.2 の未完タスク。
   **E1/E2/E3 より先**。
2. 2×2 factorial ($\omega$ pin × limiter bypass) の分離
3. E3 ($P_\omega$ の壁法則整合形)

### 数値の再生成

README の $C_f(Re_\theta)$ 系の数値はすべて [`tools/cf_retheta_analysis.py`](tools/cf_retheta_analysis.py) で再生成できる:

```bash
python3 tools/cf_retheta_analysis.py run_0042_node_yp30_planar_2nd
python3 tools/cf_retheta_analysis.py run_0049_su2_sst_lowre_fine --su2 --wall-law spalding
python3 tools/cf_retheta_analysis.py <run> --ytop 0.05 --fit-window 0.12 --fit-order 3 --edge global
```

fit 窓幅/次数・積分上端・`local`/`global` の edge 基準・ソルバ固有の壁法則をすべてパラメータ化してある。

## §4.2 実 W–I 接線力の診断 — 狙った壁摩擦は**正しく作用している** (2026-08-12)

「狙った摩擦力が W–I 間で本当に作用しているか」(当初からの本題) に答えるため、
**残差へ加わる直前の traction** を壁ノードに集計する診断を実装した
(`wi_ftan` / `wi_fnrm` / `wi_ftan_res`、[`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) の
AddTauWall 再スケール直後で `atomicAdd`、毎ステップ 0 クリア)。
`twall` 出力はモデル目標値へ再スケール済みで独立検証にならないため、これが必要だった。

`run_0051_node_wf_widiag` (= `run_0042` 収束場から 500 step、設定同一):

| 帯 | 実際に残差へ入った $\Sigma$`wi_ftan` [N] | 狙い値 $\int\rho u_\tau^2 dx$ [N] | **比** | 法線力/接線力 |
| --- | --- | --- | --- | --- |
| x=[0.30, 0.90] | 3.8796e+00 | 3.7980e+00 | **1.0215** | 0.003% |
| x=[0.20, 0.95] | 4.9573e+00 | 4.8780e+00 | **1.0163** | 0.003% |
| x=[0.45, 0.75] | 1.9399e+00 | 1.8595e+00 | **1.0432** | 0.002% |

### 結論 1: 配線は健全 (仮説を 1 つ消せた)

**【表現の訂正】** 比 1.016–1.043 を「伝達精度 2–4%」と書いたのは不正確だった。
[`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) の再スケールは
**局所的に $|F_t|=\tau_w A_{\rm WI}$ になるようコード自身が強制している**ので、
比の 2–4% は主に **W–I 双対面積による和と壁面側 $\int\tau_w dx$ の幾何・積分差、
および station の切り方/端点処理**である。正確な結論は:

> **AddTauWall 後の接線力は W–I 双対面上で意図した $\tau_w A_{\rm WI}$ になり、
> その同じ力が内部ノード残差へ加算される。したがって本ケースの低摩擦は
> 「壁応力が内部残差へ届かない」ことでは説明できない。**

壁ノード側の運動量残差は [`nodeWallDirichlet_d.cu`](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu)
でゼロ化されるので、**内部側に反力が残る構造**もコード上確認できる。

### 結論 2: 余計な法線力は実害なし (絶対値・局所最大で再確認)

初回は**符号付き和** `wi_fnrm` で見ており相殺の可能性があったため、
**絶対値和 `wi_fnrm_abs`** を追加して再測定した (帯 x=[0.3,0.9]):

| 量 | 値 |
| --- | --- |
| $\Sigma\lvert F_t\rvert$ | 3.87959e+00 N |
| $\Sigma F_n$ (符号付き) | 1.06627e-04 N (0.0027%) |
| **$\Sigma\lvert F_n\rvert$ (絶対値)** | **1.06627e-04 N (0.0027%)** |
| 壁ノードごとの $\lvert F_n\rvert/\lvert F_t\rvert$ | max 6.13e-5 / median 2.90e-5 / min 8.62e-6 |

符号付き和と絶対値和が一致 = **相殺は起きていない**。

**ただし結論は発達域に限定する**:

| 範囲 | $\Sigma\lvert F_n\rvert/\Sigma\lvert F_t\rvert$ | 局所最大 |
| --- | --- | --- |
| **発達域 x=[0.3,0.9]** | **2.75e-5** | 6.13e-5 @x=0.31 |
| 壁全体 (x≥0) | 6.87e-4 | **0.343 @x=0** (壁∩sym コーナー) |
| 前縁近傍 x=[0,0.05] | 7.54e-3 | 0.343 |

→ **発達域では法線力は無視できる。前縁コーナーには大きい局所比があるが、
壁全体の積分比は 0.069% で判定帯への影響は小さい。**
「本ケース全体で無害」とは一般化しない。コード/コメントの不一致は残る。

### 結論 3: 再スケールの向きは「増幅」ではなく「低減」だった (前提の訂正)

$\Sigma$`wi_ftan_res` (再スケール前) = 6.0581 N → 再スケール後 3.8796 N で **集約比 0.6404**。
**壁ノードごとの分布も min 0.6400 / median 0.6405 / max 0.6410 でほぼ一様**
(帯域集約だけで語らないための確認)。
$C_{f,\rm model}/C_{f,\rm molec}\approx2.1$ から「~2.1 倍に増幅される」と推定していたのは**誤り**。
W–I 双対面の解像 traction は**壁の分子勾配ではなく $\mu_{\rm total}$ (乱流粘性込み) と双対面距離**で
評価されるので、モデル $\tau_w$ より**大きい**。AddTauWall はそれを**下げて**モデル値に合わせている。

### 残る問い

伝達も法線力も無罪なので、**なぜ解が低摩擦で釣り合うのか**が残る。
次は §3 の 2×2 factorial ($\omega$ pin × strain limiter bypass) で
$\omega$ 状態 / limiter / 源項のどれが効いているかを分離する。

## §3 の 2×2 factorial — $\omega$ pin が最大因子 (2026-08-12、判定は 3 点セット)

`nodeOmegaWfDirichlet` が $\omega$ ピンと SST shear limiter 迂回を**束ねていた**問題を 2 ビットへ分解。

### 実装

- **`wf_irep_flag`**: node 壁関数の**第一内層ノード (irep) だけ**を 1 でマーク
  (`wf_pk>=0` は壁ノードにも立つのでマスクに使えない)。**cell では常に 0** = cell ビット不変が構造的に保証。
  **node SST 壁関数が有効なときだけ**カーネルへ渡す (未初期化値を読ませない)。
- **`FORGE_WF_LIMITER_BYPASS`** (env, 既定 −1=従来 / 0=OFF / 1=ON@irep)。壁関数非有効時は常に −1 に固定。
- **`FORGE_WI_FORCE_DIAG`** (env): W–I 実力診断の 0 クリア/集計を**診断時のみ**有効化 (通常経路の性能を変えない)。

### 判定 3 点セット (帯 x=[0.60, 0.80]、全 run 同一 IC・20000 step・quasisteady `STEADY`)

**帯の決め方**: 全 4 run で共通に $Re_\theta\ge4000$ となる下限を取ると
$x_{\min}$ = 0.640 (Y00/Y01 が律速、Y10/Y11 は 0.55 前後)。上限は収支の信頼域から 0.80。
→ **判定帯 x=[0.640, 0.784]** (端点は双対幅が片側になる壁ノードを除外)。

| run | (1) **W–I 実力** $C_f$/KS | (2) 運動量積分 $C_f$/KS | $\Sigma\lvert F_n\rvert/\Sigma\lvert F_t\rvert$ |
| --- | --- | --- | --- |
| $Y_{00}$ pin0 byp0 (`run_0053`) | **0.7615** | 0.784 | 1.97e-5 |
| $Y_{10}$ pin1 byp0 (`run_0054`) | **0.8594** | 0.840 | 1.96e-5 |
| $Y_{01}$ pin0 byp1 (`run_0055`) | **0.7670** | 0.794 | 1.92e-5 |
| $Y_{11}$ pin1 byp1 (`run_0056`) | **0.8864** | 0.846 | 1.87e-5 |

**再生成** (正式ツール):

```bash
python3 tools/cf_retheta_analysis.py run_0053_2x2_Y00 \
    --wi-force run_0053_2x2_Y00_widiag --band 0.64 0.80
```

W–I 実力の規約はツールが出力する: `wi_ftan` を壁ノードの双対幅で割って**局所 $\tau_{w,\rm applied}$** に直し、
$C_{f,\rm KS}\cdot q$ を同じ壁ノード列へ内挿して**同一の台形則**で両辺を積分する
(点和 vs 台形積分の端点ずれを避けるため。これを怠ると 10% 級の偏りが出る)。

### factorial (応答量 = (1) W–I 実力 $C_f$/KS)

| | 値 |
| --- | --- |
| 単純効果 pin (byp=0) $Y_{10}-Y_{00}$ | **+0.0979** |
| 単純効果 pin (byp=1) $Y_{11}-Y_{01}$ | **+0.1194** |
| 単純効果 bypass (pin=0) $Y_{01}-Y_{00}$ | +0.0055 |
| 単純効果 bypass (pin=1) $Y_{11}-Y_{10}$ | **+0.0270** |
| **平均主効果 pin** | **+0.1087** |
| **平均主効果 bypass** | **+0.0163** |
| **相互作用** $\Delta_{\rm int}$ | **+0.0215** |

**用語の訂正**: 以前「主効果」と書いていた $Y_{10}-Y_{00}$ / $Y_{01}-Y_{00}$ は
**相手因子が OFF のときの単純効果**であり、factorial の平均主効果ではない。

### 結論

1. **$\omega$ pin が最大因子** (平均主効果 +0.1075 vs bypass +0.0160)。
   **「全効果の 79% が $\omega$ 状態」という寄与率は撤回** — 相互作用の配分方法に依存するので一意に決まらない。
2. **limiter 迂回は「ほぼ無効」ではない**。pin=0 では +0.0054 だが **pin=1 では +0.0267** で
   事前に定めた有意差 $\pm0.02$ を超える。正しくは
   **「単独効果は小さいが $\omega$ pin 時に正の相互作用を持つ。それでも pin 効果より小さい」**。
3. **回復率**: $Y_{11}$=0.8864 は目標帯 0.94±0.02 に対し
   **中心 0.94 まで 70% / 下限 0.92 まで 79%** の回復 (どちらを基準にしたか必ず明記する)。
4. **「残り ~29% は $\omega$ 状態でも limiter でもない別要因」は撤回**。
   $Y_{11}$ が目標に届かないことから言えるのは
   **「今回試した特定の $\omega$ Dirichlet 値と bypass の組合せでは目標を再現しない」**までである。
   残差分は ①pin 値そのものが不適切 ②pin する DOF/範囲が不適切
   ③$P_\omega$・拡散・交差拡散など**別の $\omega$ 機構** ④$k$–$\omega$ 結合、のいずれもありうる。
   E3 はまさに $\omega$ 方程式への介入なので、「$\omega$ 状態でも limiter でもない」と断定すると自己矛盾する。
5. 収束は**全 run `NOT CONVERGED`**、派生量が `STEADY`。「収束した」とは言わない。

### 位置づけと次

`nodeOmegaWfDirichlet: 1` は依然**採用しない** (case/40 で $\tau_w$ 過大化)。本実験は因子分離が目的。
**次は §4.1 の $\omega$ 方程式の項別収支** ($P_\omega$/$D_\omega$/交差拡散/輸送) を出し、
**それを根拠に E3 を設計**する。E3 へ直行しない。

## §4.1 $\omega$ 方程式の項別収支 — 第一内層は**局所 P–D 平衡が 96%** (2026-08-13)

E3 を設計する前に、$\omega$ の平衡がどの項で決まっているかを直接測った。
診断 `omg_prod` / `omg_dest` / `omg_cross` / `omg_trans` を
[`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) に追加
(`omg_trans` はソース加算**前**の `res_roOmega` = 対流+拡散の寄与)。
`FORGE_OMEGA_BUDGET=1` のときだけ確保・出力される。

`run_0059_omgbudget` (= `run_0053` 場から 500 step)、x=0.694:

| $y$ | $y^+$ | irep | $P_\omega$ | $-D_\omega$ | $CD_\omega$ | 輸送 | 和 | $\omega$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 (壁) | 0 | 0 | 2.451e4 | −1.427e3 | 0 | 4.472e2 | 2.353e4 | 1.137e5 |
| 1.768e-4 | 27.2 | **1** | **1.793e4** | **−1.652e4** | 2.36 | −1.403e3 | **3.08e-1** | 2.621e5 |
| 3.707e-4 | 57.1 | 0 | 1.880e3 | −2.315e3 | −2.20 | 4.377e2 | 2.15e-2 | 9.539e4 |
| 5.833e-4 | 89.9 | 0 | 4.729e2 | −6.996e2 | 1.16e-4 | 2.267e2 | 3.32e-2 | 5.007e4 |

(単位は `res_roOmega` と同じ = 体積込み。壁ノードの「和」が非零なのは $\omega$ が
Dirichlet ピンで残差を後段でゼロ化するため = 設計どおり)

### 第一内層ノードの内訳 (絶対値比)

| 項 | 値 | 比 |
| --- | --- | --- |
| 生産 $P_\omega$ | +1.793e4 | **50.0%** |
| 消滅 $-D_\omega$ | −1.652e4 | **46.1%** |
| 交差拡散 $CD_\omega$ | +2.36 | 0.0% |
| 輸送 (対流+拡散) | −1.403e3 | 3.9% |
| **残差和** | **3.08e-1** | **0.00%** |

**【表現の訂正】** 「$P$–$D$ だけで 96% 平衡」は言い過ぎ。正確には:

$$\frac{D_\omega}{P_\omega}=0.922,\qquad \frac{|T_\omega|}{P_\omega}=0.0783$$

→ **現在の解では局所源項が支配的で、消しきれない生産の 7.8% を輸送が釣り合わせている**。
交差拡散は無視できる。

**収支が 0.00% で閉じるのは代数的にほぼ保証されている** — `omg_trans` は
「ソース加算**前**の残差」なので、$P-D+CD+$`omg_trans` は定義上ゼロに近い。
**符号と帳尻の確認にはなるが、独立した物理検証ではない**。

### E3 の事前予測 (実装前に立てておく)

局所 $P$–$D$ 平衡なら $\omega\propto\sqrt{P_\omega}$、$P_\omega=\alpha\rho S^2$ なので
入力ひずみを解像 $S$ → 壁法則 $S_{\rm wf}=u_\tau^2g/\nu$ に替えると $\omega$ は $S_{\rm wf}/S$ 倍になる:

| $x$ | $\omega$ 実測 | $S$ | $S_{\rm wf}$ | $S_{\rm wf}/S$ | **$\omega$ 予測** | (参考) Menter ピン値 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.653 | 2.631e5 | 1.068e5 | 5.143e4 | 0.482 | **1.267e5** | 1.145e5 |
| 0.694 | 2.621e5 | 1.063e5 | 5.143e4 | 0.484 | **1.269e5** | 1.145e5 |
| 0.753 | 2.608e5 | 1.056e5 | 5.143e4 | 0.487 | **1.270e5** | 1.145e5 |

**これは「凍結場での初期見積り」であって最終解の予測ではない。**
E3 では生産が $(0.484)^2\simeq0.234$ 倍まで落ちるので、**現在の輸送量がそのままなら
新しい生産に対して輸送は約 33% を占める**ことになる (現在は 7.8%)。実際にはさらに
輸送項・$F_1$ と $\alpha/\beta$・$u_\tau$ と $S_{\rm wf}$・$k,\mu_t$・速度場がすべて再平衡する。

したがって **「E3 は $Y_{10}$=0.8594 に近い / 目標に届かない」は事前仮説として持つに留め、
確定的には書かない**。同様に **「残り 0.86→0.94 は $\omega$ 以外」もまだ言えない** —
E3 後も $\omega$ の輸送・ブレンド・適用範囲・$k$ との結合が候補に残る。

**E3 の位置づけ**: 目標達成を狙う実験ではなく、
**「free な $\omega$ が源項整合によってどこへ再平衡するか」を検証する**もの。

### この収支サンプルの限界

`run_0059_omgbudget` は**全残差 `NOT CONVERGED`**、snapshot が 2 個しかないので
`check_quasisteady` は `TRANSIENT-UNSETTLED`。ただし着目点の $\omega$ は 500 step で
**0.0048% しか変化していない**ので、**瞬時収支のサンプルとしては十分**と判断した。
定常収支の主張には使わない。

### case/40 (軸対称) で使うときの注意

`omg_trans` は**軸対称 $1/r$ ソースが加わる前**に取得している
([`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) の `axisymMethod==1` 分岐は
`omg_trans` 記録より後)。**case/40 で収支を見るなら `omg_axisym` を別途出力する必要がある。**

## §5 E3 実施 — $\omega$ 生産を壁法則整合に (2026-08-13)

$k$ 側だけ壁関数値 (`wf_pk`) で $\omega$ 側は解像ひずみ、という**非対称の解消**。
第一内層ノード**のみ**、$P_\omega$ の入力ひずみを

$$S_{\rm wf}=\frac{u_\tau^2}{\nu}g,\qquad P_{\omega,\rm wf}=\alpha\rho S_{\rm wf}^2$$

に差し替える。$\omega$ の**状態は free のまま**、strain limiter も**通常式のまま**
(= Dirichlet で殴らない)。

**実装**: `ransWallFunction_d.cu` で $S_{\rm wf}^2$ を計算し新配列 **`wf_sprod`** へ格納
($\alpha P_{k,\rm wf}/\nu_{t,\rm wf}$ の形は低 y+ で 0/0 なので $S_{\rm wf}^2$ を直接使う)。
`ransSource_d.cu` は `wf_sprod[ic]>=0` のときだけ $\omega$ 生産の入力に使う
($k$ 生産は既に `wf_pk` で置換済みなので触らない)。マスクは第一内層のみ
(cell は `wf_sprod`≡−1 でビット不変)。env `FORGE_WF_OMEGA_SOURCE=1` ゲート。
あわせて軸対称 $1/y$ ソースの診断 **`omg_axisym`** を追加 (`omg_trans` はその加算**前**の残差なので、
軸対称ケースの収支にはこれが要る)。

### case/26 の結果 (`run_0060_E3`, quasisteady `ALL STEADY`)

| 量 | E3 前 ($Y_{00}$) | **E3 後** | (参考) $Y_{10}$ pin=1 |
| --- | --- | --- | --- |
| $C_f/C_{f,\mathrm{KS}}$ (W–I 実力, 帯 x=[0.64,0.80]) | 0.7615 | **0.8445** | 0.8594 |
| 第一内層 $\omega$ (x=0.694) | 2.621e5 | **1.345e5** | (ピン値 1.145e5) |
| 第一内層 $k$ | 10.193 | 11.483 | — |
| 第一内層 $\mu_t/\mu$ | 1.98 | 2.17 | — |

**$\omega$ 項別収支の変化**:

| | E3 前 | E3 後 |
| --- | --- | --- |
| $D_\omega/P_\omega$ | 0.9218 | 0.8952 |
| $\lvert T_\omega\rvert/P_\omega$ | 0.0783 | **0.1048** |
| 生産 $P_\omega$ | 1.793e4 | 4.683e3 (0.261 倍) |
| 輸送 $T_\omega$ | −1.403e3 | −4.906e2 |

**事前仮説の検証**: 凍結場の見積り「$\omega$ は $S_{\rm wf}/S$=0.484 倍」に対し**実測 0.513 倍**
(6% 以内)。ただし**輸送も再平衡した** (1403→491) ため、「輸送比が 33% まで上がる」という
素朴な外挿は起きず 7.8%→10.5% に留まった。**場が再平衡するという指摘のとおり**だった。
$C_f$ も $Y_{10}$=0.8594 の近く (0.8445) に着地し、**目標帯 0.94±0.02 には届かない**。

### case/40 ベルノズルの結果 (`run_0071_E3`, NaN なし)

**定常性**: 当初「`ALL STEADY`」と書いたが、それは既定量 (shock/machmax/pmax) の判定であって
**$\tau_w$ を評価したものではなかった**。`check_quasisteady.py` に **`wall_tau` を正式対応**させ
(未知 quantity は**エラー**にするよう修正 — 従来は黙って無視して `ALL STEADY` を返していた)、
改めて判定: `run_0066` / `run_0071` とも **`wall_tau: STEADY`** (tail mean 1224 / 1374、drift 0.0%)。

**採否の分かれ目**は「case/26 の改善が case/40 の過大化を伴わないこと」:

**統一評価** (`tools/wall_tau_eval.py`、`twall` 接線成分・弧長×周長重み・範囲 x=[0.0102,0.0685]、
y+1 内部基準 `run_0052` = **1453.98 Pa**):

| | bell $\tau_w$ 平均 | 壁温平均 | $\tau_w$/生産比 | y+1 内部基準比 |
| --- | --- | --- | --- | --- |
| 生産 mode3 (`run_0066`) | **1388.35 Pa** | 1414.0 K | 1.000 | **0.9549** |
| **E3 (`run_0071`)** | **1659.59 Pa** | 1421.1 K | **1.1954** | **1.1414** |
| (参考) $\omega$ ピン (`run_0068`) | **1795.89 Pa** | — | 1.2936 | **1.2352** (棄却済) |

→ **E3 も過大化する (1.1414)。$\omega$ ピンの 1.2352 より軽いが 1.0 を超える。**
内部基準比は**統一評価で直接測定した値**である (旧記載の「$\tau_w$ 比からの推定」は誤り)。

> **旧評価系による履歴 (再利用しないこと)**: 別の評価関数・別範囲では
> 生産 1224.2 Pa → E3 1374.3 Pa = **応答比 1.1226**、y+1 真値比の推定「≈1.061」だった。
> **これらは上表の統一評価値と等号で結べない**。数値を引用するときは必ず上表を使う。

### 判定

**E3 はこのままでは採用できない。** ただし性質は良い:

1. **予測どおり機能した** — $\omega$ を源項の整合だけで狙った水準へ下げ、Dirichlet を使わずに
   $Y_{10}$ 相当の効果を得た。$k$/$\omega$ の非対称は解消される。
2. **応答の向きは両ケースとも同じ (増加)** だが**大きさは同程度ではない**。
   **各ケース内で同一定義の before/after** を取ると
   case/26 は $C_f(Re_\theta)$/KS 系で $0.8445/0.7615-1$=**+10.9%**、
   case/40 は面積重み $\tau_w$ 系で $1659.59/1388.35-1$=**+19.5%**。
   **両者は同一の汎関数ではない**ので (case/26 は $C_f$/KS、case/40 は壁面範囲の $\tau_w$)、
   大小は参考比較にとどめる。
   **旧記載「両ケースに約 11–12% の同方向応答」は評価系を混ぜたもので撤回**。
   基準に対する絶対誤差は**統一評価**で **case/26 −23.9%→−15.6%** (不足のまま) /
   **case/40 −4.51%→+14.14%** (不足→過剰)。
   (旧記載「case/40 −5.5%→+6.1% でほぼ対称」は旧評価系の値で**撤回**。
   統一評価では対称ではなく case/40 の過大化の方が大きい。)
   case/40 は**事前ゲート (過大化しない) に照らして不合格**である。
3. → **今回試した 3 つの $\omega$ 介入 (Dirichlet ピン / limiter 迂回 / E3) では
   両ケースのゲートを同時に満たせなかった**。case/26 が不足 (−23.9%→−15.6%) で
   case/40 が不足→過剰 (統一評価 −4.51%→+14.14%) という**符号の違い**が次の手がかりであり、
   **独立性の高い代表点経路を次に診断する**。
   **「$\omega$ 側を出し切った」とは言わない** — $\omega$ の輸送・ブレンド・適用範囲・
   $k$ との結合は未検証。また代表点アルゴリズムという単一の構造変更が
   平板と曲面で異なる方向に働く可能性も残る。

**「残り 0.86→0.94 は $\omega$ 以外」とはまだ言わない** — E3 後も $\omega$ の輸送・ブレンド・
適用範囲・$k$ との結合が候補に残る。**「単一ノブでは閉じない」も断定しない** —
代表点アルゴリズムのような単一の構造変更が両ケースで逆向きに働く可能性がある。

## §3 代表点 (Normal_Neighbor) の 3 分離 — 完了 (2026-08-13)

診断 `rep_id`/`rep_y`/`rep_dist`/`rep_cos`/`rep_toff`/`rep_wdratio`/`rep_nx,ny,nz` を
`ransWallFunction_d.cu` に追加 (`FORGE_WF_REP_DIAG=1` ゲート、OFF なら確保も出力もしない)。
**アルゴリズムは一切変更していない純粋な診断**。

### §3.1 選択の幾何 — 平板は完璧、曲面はずれる

| 量 | case/26 平板 (x=0.10–0.95) | case/40 ベル (x=0.010–0.070) |
| --- | --- | --- |
| `best_cos` | **1.0000** (min=med=max) | 0.941 / **0.956** / 0.990 |
| 接線オフセット / $y$ | 0 / 9.3e-5 / 2.2e-4 | 0.143 / **0.306** / 0.359 |
| `wall_dist[irep]`/$y$ | **1.0000** | 1.010 / **1.046** / 1.063 |

→ **平板では代表点が数値許容差内で法線上にある** (接線オフセット/$y$ が 1e-4 級) ので、
**case/26 の欠損は選択の幾何では説明できない**。
**曲面ベルでは代表点が法線から 30% ずれ**、`wall_dist[irep]` は法線射影 $y$ より 4.6% 大きい
(`wall_dist` 自身も離散量なので「真の壁距離」とは断定しない)。

### §3.2 速度場と壁法則の感度

| | case/26 (x=0.694) | case/40 (x=0.030) |
| --- | --- | --- |
| $U_t$ / $u_\tau$ / $y^+$ | 30.61 / 2.3213 / 26.9 | 966.9 / 57.08 / 101.0 |
| 壁法則残差 $\lvert U_t/u_\tau-u^+(y^+)\rvert$ | **0.000** | **0.000** |
| 感度 $du_\tau/dU_t$ | 0.0585 ($U_t$ +1% → $u_\tau$ +0.77%) | 0.0517 (+0.88%) |

→ Newton は厳密に収束している (実装バグではない)。$u_\tau$ の −9% を説明するには
$U_t$ が約 −12% 必要。

### §3.3 法線補間による 3 分離

**case/26 平板** (幾何が完璧なので (1)=(2)):

| $x$ | $U_t$ 粗場 | $U_t$ 壁解像(法線補間) | 比 | $u_\tau$ 粗場 | $u_\tau$ 基準速度から | $u_\tau$ 壁解像真値 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.602 | 30.977 | 33.100 | 0.936 | 2.3429 | 2.4667 | 2.5716 |
| 0.694 | 30.606 | 32.675 | 0.937 | 2.3213 | 2.4420 | 2.5455 |
| 0.800 | 30.246 | 32.287 | 0.937 | 2.3003 | 2.4194 | 2.5223 |

**$u_\tau$ 欠損 (2.3429→2.5716 = −8.9%) の端点比分解**:

| 段 | 比 |
| --- | --- |
| **代表点選択** | **1.000** (数値許容差内で法線上) |
| **代表点速度** ($U_t$ が 6.4% 低い) | **1.0528** |
| **壁法則段** (基準速度を入れても真値に届かない) | **1.0425** |

**★ 解釈の注意**: 積 $1.0528\times1.0425=1.0976=2.5716/2.3429$ が一致するのは
**中間状態を定義した比の積だから代数的に自明**であり、**連成解に対する因果寄与の
一意な分解ではない**。「端点比を 3 段に分けて眺めるための分解」として読むこと。

**物性の規約** (これを固定しないと「速度だけの寄与」と言えない):
`u_tau(粗場)` と `u_tau(基準速度から)` は**どちらも粗場側の $\nu$ と $y_{\rm rep}$** を使い、
**差が速度だけ**になるようにしている。`u_tau(壁解像真値)` は基準 run 自身の分子勾配
$\tau_w=\mu U_1/y_1$ から求めるので、**壁法則段には物性・$y$ の定義差と基準解自身の誤差も混ざる**。

**壁法則段の中身**: 壁解像プロファイルを $y^+\approx29$–30 で Reichardt 相関と比べると

| $x$ | 真値 $u_\tau$ | $y^+$ | 実測 $U_t$ | Reichardt 予測 | 実測/則 |
| --- | --- | --- | --- | --- | --- |
| 0.451 | 2.6222 | 30.24 | 33.967 | 35.739 | **0.9504** |
| 0.602 | 2.5716 | 29.64 | 33.100 | 34.861 | **0.9495** |
| 0.694 | 2.5455 | 29.37 | 32.675 | 34.422 | **0.9492** |
| 0.800 | 2.5223 | 29.08 | 32.287 | 34.016 | **0.9492** |

→ **壁解像基準のプロファイルが、採用した Reichardt 相関を $y^+\approx30$ で約 5.1% 下回る**
($x$ によらずほぼ一定)。

**★ これは「壁法則自体が間違っている」ことの証明ではない**。この差には
**壁解像解自身の誤差・有限 Re・入口条件・自由流乱流・物性定義**が全部含まれる。
自由流 $\mu_t/\mu$=65.9 (NASA 標準の 4 桁大) は**一候補**にすぎない。
「普遍則」「壁関数の前提そのものがずれた」という言い方はしない。

**case/40 ベル** ((1)/(2) が**純粋なサンプリング誤差** — 同じ場・2 つの位置):

| $x$ | toff/$y$ | (1) 粗場 irep | (2) 粗場法線 | (3) 基準法線 | **(1)/(2)** | (2)/(3) |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0151 | 0.354 | 939.9 | 739.4 | 1046.4 | **1.271** | 0.707 |
| 0.0249 | 0.336 | 960.0 | 768.6 | 1088.1 | **1.249** | 0.706 |
| 0.0349 | 0.308 | 972.0 | 789.2 | 1104.0 | **1.232** | 0.715 |
| 0.0498 | 0.251 | 979.4 | 811.0 | 1118.8 | **1.208** | 0.725 |

→ **曲面では代表点が読む速度が、同じ $y$ の法線上より 21–27% 高い**。
case/26 は 0% なので、**サンプリング誤差は case/40 固有**。

### 両ケースの符号が逆だった理由 (仮説)

case/40 では**選択誤差が速度を 21–27% 過大に読む**方向に働き、
粗場そのものの不足 ((2)/(3)=0.71) を**部分的に相殺**している (正味 (1)/(3)≈0.90)。
case/26 は幾何が完璧なので相殺がなく、欠損がそのまま出る。
→ **同じ $\omega$ 介入で case/26 は不足のまま・case/40 は過剰になった**のはこれで説明がつく。

**ただし仮説**である: (2)/(3) は**別々の収束解を同じ $y$ で比べている**ので、
サンプリングと解の違いが混ざる (親 plan §2 の注意)。純粋に確定しているのは
**(1)/(2) = 同一場・2 位置の比較 = 選択誤差**の部分だけ。

### 【撤回】凍結場カウンターファクト: 代表点を法線化すると case/40 は悪化する

> **IDW 版の数値 (0.679) は撤回する。** 近傍数 $k$ に強く依存し、**k=1/4/8/16 で
> 1.000 / 0.679 / 0.735 / 0.825** と動いた = 物理量ではなく補間の産物だった。
> 原因は 4 つ: ①k-NN IDW は凸包外でも値を返すので「外挿なし」は誤り ②壁ノードの
> ゼロ速度を近傍平均に含めて法線上速度を人為的に下げていた ③`hypot(ux,uy)` を補間しており
> 本体の $U_t=|\mathbf U-(\mathbf U\cdot\mathbf n)\mathbf n|$ と定義が違う
> ④「面平均」が壁ノードの算術平均で面積重みが無い。
> 要素ベース版は下記。以下は撤回された記述として残す。

### 【撤回された記述】IDW 版カウンターファクト

「曲面の代表点を法線上サンプリングへ直す」(候補 A) を**実計算の前に**凍結場で見積もった。
既存場のまま、法線補間速度を壁法則へ入れて $u_\tau,\tau_w$ を再計算する
(`tools/rep_point_analysis.py cf`)。

case/40 `run_0072_repdiag` (x=0.010–0.070):

| 量 | min | median | max |
| --- | --- | --- | --- |
| $U$(法線)/$U$(irep) | 0.7746 | **0.8128** | 0.8532 |
| $u_\tau$(法線)/$u_\tau$(現行) | 0.7984 | **0.8346** | 0.8718 |
| $\tau_w$(法線)/$\tau_w$(現行) | 0.6374 | **0.6965** | 0.7600 |

面平均 $\tau_w$: 1719.94 → 1168.53 Pa (**比 0.679**)

→ **法線化すると $\tau_w$ が約 32% 下がる**。case/40 は現状 y+1 基準比 0.945 なので、
**法線化は 0.6 台まで悪化させる**見込み。→ **候補 A は単独では逆効果**であり、
「選択誤差が粗場の不足を相殺している」という仮説と整合する。

**ただし凍結場の即時効果**であり、実際に選択を変えれば場が再平衡する
(E3 では凍結場見積り 0.484 に対し実測 0.513 と 6% 以内だった前例がある)。

**再生成**:

```bash
python3 tools/rep_point_analysis.py geom  <diag_run> --xmin .. --xmax ..
python3 tools/rep_point_analysis.py split <diag_run> --ref <壁解像 run> --stations ..
python3 tools/rep_point_analysis.py cf    <diag_run> --xmin .. --xmax ..
```

補間法 (逆距離重み `--k`、外挿なし)・使用した $\rho,\nu,y$ の規約はツールが出力に明記する。

### 凍結場カウンターファクト (要素ベース・検証付き) — 2026-08-13

IDW 版の $k$ 依存を解消するため [`tools/rep_normal_interp.py`](tools/rep_normal_interp.py) を新設:

1. **メッシュ接続から包含要素を特定** (VIZMESH の primal quad → 2 三角形、重心座標)。
   **包含要素が無ければ除外** = 真に外挿しない。
2. **速度ベクトルを補間**してから本体と同じ $U_t=|\mathbf U-(\mathbf U\cdot\mathbf n)\mathbf n|$ を取る。
3. **壁面積重み付き平均** (軸対称なら周長 $2\pi r$ を掛ける)。
4. **先に本体の `utau` を再構成できるか検証**してから反実値を出す。

case/40 `run_0072_repdiag` (x=0.010–0.070, axisym):

| | 値 |
| --- | --- |
| **★ 再構成検証** (本体 `utau` との相対差) | **median 1.46e-7 / max 1.96e-4** |
| 包含要素なしで除外した壁ノード | 1 個 / 94 |
| $u_\tau$(法線)/$u_\tau$(現行) | 0.7791 / **0.7810** / 0.7833 (min/med/max) |
| $\tau_w$ 比 | 0.6070 / **0.6099** / 0.6136 |
| **壁面積重み付き平均 $\tau_w$** | 1378.43 → **840.91 Pa (比 0.6100)** |

**$k$ 依存が消え** (IDW は 1.000–0.825 で振れた)、**再構成検証も通っている**ので、
**凍結場での即時効果としては 0.610 が妥当**と考える。

### まだ言えないこと

- **「case/40 が 0.6 台まで悪化する」とは言えない**。既往の **0.945** (y+1 基準比) は
  **範囲・重み・$\tau_w$ 定義が本ツールと同一であることを未確認**なので、
  $0.945\times0.610$ のような単純な掛け算はできない。
  → **基準側を同じ評価関数に通す**のが先 (代表点 plan の完了条件に追加済み)。
- **「単独では逆効果」という判定も未確定**。上記の同一定義比較が要る。
- **「再平衡しても方向は覆らない」も撤回**。E3 で凍結場見積りが 6% 以内だったのは
  **別の介入**であり、代表点変更の再平衡量を保証しない。

**確定しているのは**: case/40 の代表点に大きな接線オフセット (0.25–0.36) があること、
および**同じ粗格子場の中で**斜め隣接点と法線上で速度が 21–27% 違うこと。

### 単一評価関数による決着 (2026-08-13)

これまで「y+1 基準比 0.945」と「カウンターファクト 0.610」は**別の評価関数**で出しており
掛け算できなかった。[`tools/wall_tau_eval.py`](tools/wall_tau_eval.py) を新設して**同一定義**に揃えた:

- **定義**: `res_wall_*.h5` の `twall_x/y/z` (壁半割面 traction) の**接線成分**
  $\tau_w=\lvert\mathbf t-(\mathbf t\cdot\mathbf n)\mathbf n\rvert$。法線は壁面メッシュの隣接節点から作る
  (診断出力に依存しないので任意の run に使える)。
  **壁解像 run では解像値 (=「壁解像内部基準が与える数値的な壁応力」。基準解自身も
  離散誤差を持つので「真の壁応力」とは呼ばない)、壁関数 run ではモデル値へ再スケール済み** =
  「そのソルバ構成が壁に与えている応力」なので両者の比較に使える (`utau` は壁解像で 0 なので不可)。
- **重み**: 壁面弧長 $ds$ × 周長 $2\pi r$ (軸対称)。壁ノードの算術平均は使わない。
- **範囲**: x=[0.010, 0.070]、全 run で壁ノード 93 個・x=[0.0102, 0.0700] と一致。

| run | $\tau_w$ [Pa] | **/y+1 内部基準** |
| --- | --- | --- |
| `run_0052_node_yp1_dofonly_regr` (y+1 内部基準) | 1453.98 | **1.0000** |
| `run_0066` 生産 mode3 | 1388.35 | **0.9549** |
| `run_0068` $\omega$ pin (棄却済) | 1795.89 | **1.2352** |
| `run_0071` **E3** | 1659.59 | **1.1414** |
| `run_0072_repdiag` (= run_0066 相当) | 1388.33 | 0.9549 |
| **反実 (代表点法線化, 凍結場)** | 846.88 | **0.5825** |

**範囲は x=[0.010, 0.069]** — 上端を 0.070 から下げたのは、**出口端ノード (x=0.07000) の
法線点が領域外**になり要素内補間が定義できないため。これで**両ツールが同じ 91 点・
同じ x=[0.0102, 0.0685]** を使い、現行値も **1388.41 vs 1388.33 Pa = 差 0.006%** で一致する
(以前は 93 点 vs 92 点で 0.72% ずれていた)。

**★ E3 の過大化を訂正**: 従来「≈1.061」と書いたのは $\tau_w$ 比を 0.945 に掛けた**推定値**だった。
同一定義で評価し直すと **1.1414** で、**推定より大きく外れている**。
生産構成の 0.9549 は従来の 0.945 とよく一致するので、**ずれていたのは E3 の推定側**。
**$\omega$ pin も同一評価で 1.2352** (従来記録の 1.237 とほぼ一致)。

**代表点法線化のカウンターファクトを接続** (同一 91 点):

| | 値 |
| --- | --- |
| `run_0072` の $\tau_w$ (`wall_tau_eval`) | 1388.33 Pa |
| `rep_normal_interp` の「現行」平均 | 1388.41 Pa (**差 0.006%**) |
| 反実 $\tau_w$ = 1388.33 × 0.6100 | 846.88 Pa |
| **→ y+1 内部基準比** | **0.5825** (現行 0.9549 から) |

→ **代表点を法線上へ直すと、case/40 は凍結場で 0.955 → 0.58 に低下する。**

**残る留保**: これは**凍結場**の即時効果であり、実際に選択を変えれば場が再平衡する。
再平衡量は未知 (E3 の 6% 一致は別介入なので保証にならない)。したがって
**「候補 A は単独では採らない」までは言えるが、「0.58 になる」と断定はしない**。

## 候補 B: 自由流乱流条件 — **仮説は棄却** (2026-08-13)

> **収束の位置づけ**: 新 4 run は**全残差では `NOT CONVERGED`** (例 `run_0064` は `rms_roOmega`=2.77 で停滞)。
> 一方**判定対象の $\theta$ と $C_f(Re_\theta)$ は準定常 `STEADY`** (drift $\theta$ 0.0–0.3% /
> $C_f/C_{f,\mathrm{KS}}$ 0.1%、fluctuation 最大 0.4%)。結果は使えるが「全 run STEADY」では不足なので
> **この 2 文で表現する**。

「壁解像プロファイルが Reichardt 相関を 5.1% 下回るのは自由流 $\mu_t/\mu$=65.9 が原因では」
という仮説を検証した。

### 境界条件の決め方 (NASA TMR 一次記述)

$\mu_{t,\infty}/\mu_\infty$ だけでは BC は決まらない ($k_\infty/\omega_\infty$ の比しか決まらない) ので、
[NASA TMR (SST)](https://tmbwg.github.io/turbmodels/sst.html) の推奨レンジから決めた:

$$\frac{U_\infty}{L}<\omega_{\rm ff}<\frac{10U_\infty}{L},\qquad
10^{-5}\frac{U_\infty^2}{Re_L}<k_{\rm ff}<0.1\frac{U_\infty^2}{Re_L},\qquad
10^{-5}<\frac{\mu_t}{\mu}<10^{-2}$$

本ケース ($L$=1.1 m, $Re_L$=4.91e6) では $\omega_{\rm ff}\in[61.6,\,616.2]$、
$k_{\rm ff}\in[9.36\times10^{-9},\,9.36\times10^{-5}]$。

| 条件 | $k_\infty$ | $\omega_\infty$ | 目標 $\mu_t/\mu$ | NASA レンジ |
| --- | --- | --- | --- | --- |
| **現行** | 0.3 | 300 | 65.9 | $\omega$ は内、**$k$ は上限の 3206 倍** |
| 中間 (単調性確認用) | 9.3572e-3 | 616.2 | 1.0 | $k$ はレンジ外 (中間プローブ) |
| **TMR 範囲内条件** | **8.4215e-5** | **616.2** | **0.009** | **$k,\omega$ とも内** |

### 結果 — 4 桁振っても準定常量に検出可能な変化がない

y+27 壁関数と y+0.7 壁解像の**同一条件ペア**を 3 条件で回した (全 run quasisteady `STEADY`)。
$x$ でなく $Re_\theta$ も併記する:

| run | $\mu_t/\mu$ (入口指定) | **実測 主流 $\mu_t/\mu$** (x=0.05/0.3/0.7) | $C_f/C_{f,\mathrm{KS}}$ | $Re_\theta$ | 壁解像 実測/Reichardt |
| --- | --- | --- | --- | --- | --- |
| `run_0042` 現行 wf | 65.9 | 65.58 / 65.09 / 64.39 | **0.764** | 4390 | — |
| `run_0065` 中間 wf | 1.0 | 0.9913 / 0.9777 / 0.9598 | **0.764** | 4413 | — |
| `run_0063` TMR範囲内 wf | 0.009 | 0.008921 / 0.008799 / 0.008638 | **0.764** | 4417 | — |
| `run_0050` 現行 壁解像 | 65.9 | 65.58 / 65.09 / 64.39 | **0.944** | 5132 | **0.9492** |
| `run_0066` 中間 壁解像 | 1.0 | — | **0.943** | 5103 | **0.9494** |
| `run_0064` TMR範囲内 壁解像 | 0.009 | 0.008922 / 0.008799 / 0.008638 | **0.943** | 5110 | **0.9494** |

**BC 変更は確実に場へ効いている** (実測主流 $\mu_t/\mu$ が 4 桁にわたって指定どおり、
流下減衰も 2% 程度)。にもかかわらず:

- 壁関数の $C_f/C_{f,\mathrm{KS}}$ は **0.764 のまま**
- 壁解像は **0.944→0.943**
- Reichardt との差は **0.9492→0.9494**

**「一切動かない」とは言わない** — 0.944→0.943 は小さいがゼロではなく、全残差も未収束である。
正確には **「今回の準定常量では検出可能な変化がない」**。

### 結論

$\mu_t/\mu$ を **65.9 → 0.009 と 4 桁**下げても、壁解像プロファイルの Reichardt 差も
node 壁関数の欠損も**準定常量として検出可能な変化を示さない**。
→ **5.1% 差の主因候補から自由流乱流を外す**。

**副産物**: 本ケースの現行 $k_\infty$ が TMR 推奨上限の **3206 倍**だと判明した。
今回の量では影響が検出されなかったが、**TMR ケースとの比較を「同一条件」と称する場合は
$k_\infty,\omega_\infty$ を合わせる必要がある**ことは変わらない。

### 数値の再生成

```bash
# 収束 (全残差)
python3 ../../solver_density_cuda/tools/check_convergence.py run_0063_B_nasa_wf
# 準定常 (判定対象量)
python3 ../../solver_density_cuda/tools/check_quasisteady.py run_0063_B_nasa_wf \
    --quantity theta,cf_retheta --cf-x 0.7 --tail 0.5 --drift 0.02 --osc 0.02
# Cf(Re_theta)
python3 tools/cf_retheta_analysis.py run_0063_B_nasa_wf --stations 0.70
```

**主流 $\mu_t/\mu$ の採取**: `res_*.h5` の `vis_turb/vis_lam` を
$|x-x_0|<0.01$ かつ $|y-0.15|<0.01$ (BL 外) の節点で平均。$x_0$=0.05/0.3/0.7。

**Reichardt 比の補間**: 壁法線カラム (構造格子なので 1 次元) 上で $y_{\rm rep}$=1.7681e-4 へ
線形補間。case/26 では IDW と要素補間の差は 0.115% なので 5.1% 差の主因ではないが、
**正式値は要素補間 (`tools/rep_normal_interp.py` 系) に統一する**。

### SU2 でも同じ 5% — forge 固有の壁解像誤差はほぼ除外 (2026-08-13)

5.1% 差が forge 側の壁解像誤差かを切り分けるため、**全残差で収束済み**の SU2 壁解像
(`run_0049_su2_sst_lowre_fine`, rms −5.0〜−12.6) について**要素補間**で同じ量を出した:

| $x$ | 真値 $u_\tau$ | $y^+$ | 実測 $U_t$ | Reichardt 予測 | **実測/則** |
| --- | --- | --- | --- | --- | --- |
| 0.451 | 2.6280 | 30.39 | 34.035 | 35.865 | **0.9490** |
| 0.602 | 2.5761 | 29.79 | 33.161 | 34.969 | **0.9483** |
| 0.694 | 2.5516 | 29.50 | 32.748 | 34.545 | **0.9480** |
| 0.800 | 2.5278 | 29.23 | 32.352 | 34.136 | **0.9477** |

**forge 壁解像の 0.9492–0.9504 と 0.1–0.2% で一致**。

→ **forge 固有の壁解像誤差・物性定義・壁応力定義はほぼ除外できる。**
独立した 2 コード (一方は全残差収束) が、この境界層は $y^+\approx30$ で
Reichardt 相関を約 5% 下回ると一致して示している。

**残る候補**: 有限 $Re$ (この $Re_\theta$ 帯)、入口/前縁の扱い、
Reichardt 相関そのものの適合限界。→ 別 plan へ切り出す。

#### 壁法則スイープ — 差は Reichardt 固有で全法則共通ではない

同一 SU2 解・同一 $y_{\rm rep}$ で**壁法則だけ**差し替えた (`tools/su2_wall_law_ratio.py --law ...`):

| 壁法則 | 定数 | 実測/則 ($x$=0.451) | ($x$=0.800) |
| --- | --- | --- | --- |
| Reichardt | $\kappa$=0.41, 7.8 項 | 0.9490 | 0.9477 |
| **Spalding** | **$\kappa$=0.41, $B$=5.0** | **1.0218** | **1.0198** |
| 純対数則 | $\kappa$=0.41, $B$=5.0 | 0.9718 | 0.9672 |

**全法則共通ではない** — Spalding では逆に 2% 超過する。
Spalding は定数流儀の違いがあるので**定数を必ず併記する**こと。

**言えること**: 固定した SST プロファイルに対し、**壁法則選択によるモデル形式差が約 7%**
(0.948–1.022) ある = **壁法則の選択は観測された 5% 差と同程度の不確かさを持つ**。

**言えないこと**: 「差の相当部分は相関側」「Spalding が正しい」「原因が Reichardt に確定した」は
いずれも**未確定**。forge と SU2 は独立コードだが**どちらも SST** であり、
**共通の SST モデル誤差は除外できていない**。
言えるのは「**採用した Reichardt 式・定数に対する −5% の符号付き差は Spalding では再現しない**」まで。

**次の A/B には交絡がある**: 壁法則は $u_\tau$ 逆解きだけでなく
$g=du^+/dy^+$ 経由で `wf_pk`・E3 の `wf_sprod`、$u_\tau$ 経由で `roK_wf`・**Kader 壁熱流束**にも
波及する (`ransWallFunction_d.cu` L226/L237/L251/L259–294)。
一括置換の A/B は結論に使えないので、**逆解き則のみ差し替える診断実験**と
**全閉包を同一法則で揃える整合実験**に分けること。詳細は
[plan §5](../../plans/active/turbulence-reichardt-gap-residual.md)。

なお 4 station の $Re_\theta$ は **3632 / 4592 / 5167 / 5817**。

**位置づけ (重要)**: この差は **forge 固有の壁関数実装バグではない**が、
**壁関数と無関係でもない**。壁解像プロファイル自身が Reichardt と 5% ずれている以上、
その速度を Reichardt 壁関数へ入力すれば $u_\tau$ は低く推定される
(端点比分解の**壁法則段 1.0425**)。→ **低摩擦バイアスに寄与する要因**である。

再生成:

```bash
python3 tools/su2_wall_law_ratio.py run_0049_su2_sst_lowre_fine/flow.vtu \
    --yrep 1.7681e-4 --stations 0.451,0.602,0.694,0.800 [--law spalding|log]
```

### 凍結場の壁法則逆解き比較 (plan §5.2 ①、2026-08-13)

連成 A/B の**前**にモデル形式感度だけを取り出す。ツール: `tools/wall_law_inverse.py`
(定数: Reichardt $\kappa$=0.41+7.8 項 / Spalding・log $\kappa$=0.41, $B$=5.0。
$y^+$ 帯 25–35 で絞り、カラム中央値)。

**A. 壁解像場に各法則を当てる** (`run_0064_B_nasa_ref`、内部基準 $u_\tau$ = 壁第一節点の分子勾配):

| 壁法則 | $u_\tau$/内部基準 | $\tau_w$/内部基準 |
| --- | --- | --- |
| Reichardt | **0.9616** | 0.9246 |
| **Spalding** | **1.0201** | 1.0405 |
| 純対数則 | 0.9805 | 0.9614 |

→ **Reichardt は内部基準より $u_\tau$ を中央値 3.8% 低く出す**、
**採用した Spalding 式・定数 ($\kappa$=0.41, $B$=5.0) は中央値 2.0% 高く出す**。
**どちらも内部基準を厳密には再現しない**。
なお**内部基準は実験的真値ではない** (壁解像 SST 解自身) ので、
「Spalding が物理的に不正解」とまでは言えない。

**B. 壁関数場の凍結逆解き** (`run_0063_B_nasa_wf`、現行 Reichardt 比。
本体 `utau` の再構成検証 ($y^+$ フィルタ後) median **4.2e-8** / max **1.9e-7**):

| 壁法則 | $u_\tau$ 比 | $\tau_w$ 比 |
| --- | --- | --- |
| Reichardt (現行) | 1.0000 | 1.0000 |
| **Spalding** | **1.0599** | **1.1234** |
| 純対数則 | 1.0163 | 1.0328 |

**A と B の法則間比が一致する**: Spalding/Reichardt は
壁解像場で **1.0608**、壁関数場で **1.0599** (**0.08% 差**)。

ただし**両者は同じ平板・ほぼ同じ $y^+$ 帯**であり、
**壁法則間比は主に $y^+$ で決まる**ので近い値になること自体は自然である。
確定範囲は「**case/26 の $y^+$=25–35 では、異なる 2 つの場に対して
凍結場感度が 0.1% 以内で再現した**」まで。
**case/40 や圧縮性流れでも頑健とは言えない**。

#### 何が言えるか

現行ギャップは $C_f$/KS 0.764 (壁関数) vs 0.943 (壁解像) = $\tau_w$ で **1.234 倍**。
凍結場で Spalding にすると $\tau_w$ ×1.1234 なので

$$\frac{1.1234-1}{1.234-1} = 0.527$$

= **必要な応力差の約 53%** を埋める計算になる (0.764 → 約 0.86 相当)。
**残りの必要応答には粗場速度などが含まれる** (どれがどれだけ寄与するかは未確定)。

**ただし連成すれば場が再平衡するので 0.86 になるとは断定しない。**
また A が示すとおり **Spalding は内部基準を 2% 超過する**ので、
壁関数側の他の不足を Spalding の過大側で相殺する形になりうる。
**この相殺は case/26 と case/40 で符号が変わりうる** (E3 の前例) ため、
採否は §5.2 ③ の整合実験 + 両ケース連成 A/B まで保留する。

#### 位置づけ

**「Spalding が有望な候補」ではなく、
「逆解き則のモデル形式だけで摩擦が約 12% 動くことが確認できた」**が最も堅い。
「Reichardt だけが全原因」でも「Spalding に替えれば解決」でもない。

再生成:

```bash
python3 tools/wall_law_inverse.py run_0064_B_nasa_ref --mode resolved \
    --yrep 1.7681e-4 --ypmin 25 --ypmax 35
python3 tools/wall_law_inverse.py run_0063_B_nasa_wf --mode wf --ypmin 25 --ypmax 35
```


### §5.2 ② 診断実験の実装・検証 (2026-08-13)

`FORGE_WF_INV_LAW=spalding` で**逆解き則だけ**を Spalding へ替える診断経路を実装した
(既定 `reichardt`)。**`wf_pk` の $g$ は Reichardt 固定・E3 は同時指定でエラー終了**という
交絡分離の仕様は plan §6 の表に記載。

| 合格条件 | 結果 |
| --- | --- |
| デバイス Spalding == CPU 後処理 (`wall_law_inverse.py`) | **1.0000** (min/max とも) |
| $u_\tau^{Sp}/u_\tau^{R}$ (凍結場予測 1.0599) | **1.0597** |
| $\tau_w^{Sp}/\tau_w^{R}$ (凍結場予測 1.1234) | **1.1230** |
| E3 同時指定 | `[FATAL]` exit=1 |
| env OFF の不変性 | ビット同一は不可 (node は run-to-run 非決定的)。**変更の寄与は実行間ノイズと区別不能** |

**★ 実装上の落とし穴**: 既存の `kNewtonIt=5` では Spalding 直接 Newton が
層流初期値から収束せず、高 $u_\tau$ 側でデバイス値が CPU と **0.9%** ずれた。
**float32 は無関係**で 10 回で機械精度。`kNewtonItSpalding=30` で解消した。

再生成:

```bash
FORGE_WF_INV_LAW=spalding ../../solver_density_cuda/tools/run_case.sh <run_dir>
python3 tools/wall_law_inverse.py <run_dir> --mode wf --xmin 0.05 --ypmin 25 --ypmax 35
```


### §5.2 連成 A/B — 両ケースが同方向に改善 (2026-08-13)

診断構成 (逆解き則のみ Spalding、$g$ は Reichardt 固定、E3 OFF) を**準定常まで**回した。
対照は**同じ restart から**走らせた Reichardt 版。

| 定義 | Reichardt (`run_0075`) | Spalding (`run_0076`) | 変化 |
| --- | --- | --- | --- |
| (c) 運動量積分/KS [壁出力非依存] | 0.754 | 0.811 **(未確定)** | +7.6% |
| (d) 壁法則逆解き/KS [整合後処理] | 0.764 | **0.822** | +7.6% |
| 生 $\tau_w$ (x=[0.05,1.0]) | 6.80 Pa | **7.24 Pa** | +6.5% |

**★ 落とし穴**: `cf_retheta_analysis.py` の (d) は既定が **Reichardt 後処理**。
Spalding run をそのまま通すと 0.764→0.735 と**悪化して見える**が、これは後処理定義の
不整合であって物理ではない。**`--wall-law spalding` を明示すること**。

**【確定 (2026-08-13)】** 延長 run (各 60000 step) で**3 判定量が `STEADY`** になった (全残差は `NOT CONVERGED` のまま):

| 量 | Reichardt (`run_0077`) | Spalding (`run_0078`) | 変化 |
| --- | --- | --- | --- |
| **(c) 運動量積分/KS** | **0.7557** | **0.8168** | **+8.1%** |
| (d) 壁法則逆解き/KS [整合後処理] | 0.7643 | 0.8199 | +7.3% |
| 生 $\tau_w$ | 6.80 Pa | 7.23 Pa | +6.3% |
| $\theta$ | 9.919e-4 | 1.049e-3 | +5.8% |

**20000 step の (c)=0.811 は未収束値**だった。**目標 0.94 までまだ +15.1% 必要**。

`cf_momentum` は正式後処理と実装を一本化済み
(`solver_density_cuda/tools/flatplate_bl.py`)。以前は fit 窓が実質 ±0.12 m (正式 ±0.08 m) で
別の値 (0.7511/0.8136) を報告していた。

(以下は 20000 step 時点の記録) **(c) はまだ動いていた**: `check_quasisteady.py --quantity cf_momentum` で
**drift 6.6%/tail `DRIFTING`** (step 10000/15000/20000 = 0.761/0.798/0.811 と増加中)。
一方 Spalding 整合 (d) は 0.824/0.823/0.822、生 $\tau_w$ は 7.274/7.251/7.240 Pa で落ち着いている。
→ **改善方向は信用できるが (c)=0.811 を最終値としない**。
延長 run `run_0077_ext_reichardt` / `run_0078_ext_spalding` (各 60000 step) で確定させる。

case/40 側は $\tau_w$ が y+1 内部基準比 **0.9548 → 1.0063** (統一評価の従来と同じ引数
x=[0.010,0.069]、91 点。壁温 +0.9 K、断熱壁なので $q_w$=0)。

**これまでの介入で初めて両ケースが同方向に動き、かつ両方が改善した**
($\omega$ pin は case/40 を 1.2352、E3 は 1.1414 に過大化させたが Spalding は 1.0063)。

**採用はまだ決めない**: (1) これは**診断構成で正式仕様ではない** ($g$ が Reichardt のままで
物理的に不整合)、(2) case/26 は 0.822 で**目標 0.94±0.02 に未達**、
(3) 凍結場見積り +12.3% に対し連成は +6.5% = **再平衡が約半分を打ち消した**。
詳細と次段は [plan §7](../../plans/active/turbulence-reichardt-gap-residual.md)。


### §5.2 ③ 整合実験 — ② の改善は偶然の相殺ではなかった (2026-08-13)

$u_\tau$ の逆解きに加え **$k$ 側の壁閉包** ($g$ → `wf_pk`, `roK_wf`) も Spalding へ揃えた。
**$\omega$ 生産 (解像ひずみ)・E3 (OFF)・熱側 (Kader)・圧縮性補正 (未導入) は変えていない**ので、
「壁法則を全部揃えた」わけではない。

| 構成 | case/26 (c) | case/40 $\tau_w$/y+1 基準 |
| --- | --- | --- |
| ① 現行 Reichardt | 0.7557 | 0.9548 |
| ② 逆解きのみ (混成) | 0.8168 | 1.0063 |
| **③ 逆解き+$k$側閉包** | **0.8282** | **1.0027** |
| 目標 | 0.94±0.02 | ≈1.0 |

**両ケースともさらに改善**。→ ② の好結果は混成による相殺ではなかった。

**★ $g$ 変更の符号は $y^+$ で反転する** ($g_{\rm Sp}/g_{\rm R}=1$ の交差点は $y^+\approx41.8$)。
本 case は $y^+$~28、case/40 は $y^+$~92 で**交差点の反対側**にあり、

| ケース | `wf_pk` (③/②) | `roK_wf` (③/②) |
| --- | --- | --- |
| case/26 ($y^+$~28) | **−7.8%** (凍結場予測 −12.3%) | **+15.0%** (+19.1%) |
| case/40 ($y^+$~92) | **+7.3%** (+10.8%) | **−7.8%** (−10.2%) |

と**逆向きに動いたのに、どちらも改善した**。$k$ 側閉包の符号ではなく
**$u_\tau$ と $k$ 側を整合させること自体**が効いていることを示唆する。

**採用はまだ決めない**: (1) case/40 の**等温壁で $q_w$** を確認していない (断熱壁は $q_w$=0)、
(2) **圧縮性補正**未検証、(3) case/26 は目標まで残り +13.5%。
