# 49. 超音速平板の環状深キャビティ (等温壁) — 内部温度場

計画: [`plans/active/case-plate-annular-cavity-m5.md`](../../plans/active/case-plate-annular-cavity-m5.md)

超音速流 (既定 M5) が流れる**断熱平板**に開いた円柱キャビティ (既定 φ50 × 深さ 50 mm) の中に、
**φ45 × 深さ 50 mm の円柱が上面面一で埋まっている**。流体側は**幅 2.5 mm・深さ 50 mm の
底閉じ環状スリット**で、キャビティ壁は等温 (既定 500 K)。**その内部温度場**を求めるのが目的。

対称面 $y=0$ の**半割** ($y\ge0$) で解く。形状・条件は下記のパラメータで振れる
(ユーザ要件 2026-09-19: **内部円柱の下流偏心・すきま量・マッハ数が変わりうる**)。

![形状と計算領域](geom_layout.png)

## パラメータ (単一ソース)

| ファイル | 役割 |
| --- | --- |
| **[`case.json`](case.json)** | **最上位入力** (条件 + 形状 + メッシュ + 評価)。マッハ数・すきま・偏心を変えるのはここだけ |
| [`setup.py`](setup.py) `--resolve` | `case.json` → **`manifest.json`** (自由流・BC 値・グループ別面積・VL 層数・physID・CV 定義まで解決) |
| `manifest.json` | **CAD・メッシュ・IC・BC・評価が読む唯一の解決済み設定** (転記事故を防ぐ) |
| [`tools/geom_common.py`](tools/geom_common.py) | 偏心対応の CV マスク・すきま中央線・周方向平均の重み |

```bash
python3 setup.py --resolve                          # case.json -> manifest.json
cd cad && QT_QPA_PLATFORM=offscreen /home/sano/opt/squashfs-root/usr/bin/freecadcmd build_geom.py
python3 tools/plot_geom.py                          # geom_layout.png を再生成
python3 tools/geom_common.py                        # CV マスクの自己検査 (偏心時の取りこぼし)
```

### 既定の作動条件 (`setup.py` 出力, ISA 20 km · M5)

| 量 | 値 | | 量 | 値 |
| --- | --- | --- | --- | --- |
| $T_\infty$ | 216.65 K | | $Re/m$ | 9.1357e6 /m |
| $P_\infty$ | 5474.72 Pa | | $T_t$ (CPG / TP) | 1299.90 / **1222.73** K |
| $\rho_\infty$ | 0.0880483 kg/m³ | | $T_{aw}$ (CPG / TP) | 1187.55 / **1126.33** K |
| $U_\infty$ | 1475.21 m/s | | $T_w/T_{aw}$ (TP) | 0.444 (強冷却) |
| $\mu_\infty$ | 1.42178e-5 Pa·s | | $k_\infty$ / $\omega_\infty$ | 81.6 m²/s² / 5.05e4 1/s |

> **CPG は $T_{aw}$ を +61 K (+5.4 %) 過大に出す** (高温で $c_p$ が +17 % 上がるため)。成果物が温度そのものなので
> **生産 EOS は semi-perfect (TP)** とする (計画 §4.10)。CPG はパイプライン疎通と段階起動レシピ確立にのみ使う。

## 形状と境界 (既定値)

すきま = $R_o-R_i$ = 2.50 mm (偏心 `x_off` を入れると周方向に $2.5\mp|x_{off}|$ mm)。
開口幅/深さ = 0.05 の**極めて深いキャビティ**。

| physID | グループ | 位置 | BC |
| --- | --- | --- | --- |
| 1 | `inlet` | $x=-80$ mm | `inlet_uniformVelocity` + `inletProfile` (2D 前駆の乱流 BL) |
| 2 | `outlet` | $x=+120$ mm | `outlet_statPress` ($P_s=P_\infty$ + 逆流 Pt/Tt) |
| 3 / 4 / 5 | `top` / `side` / `sym` | $z=60$ / $y=70$ / $y=0$ | `slip` |
| 6 / 7 | `plate` / `plate_in` | $z=0$, $r>55$ / $25<r<55$ | `wall` (断熱 no-slip)。`plate_in` は開口周りの細分用 |
| 8 / 9 / 10 | `cav_outer` / `cav_floor` / `cyl_side` | キャビティ 3 壁 | `wall_isothermal` Ts=500 |
| 11 | `cyl_top` | 内円柱上面 | `wall_isothermal` Ts=500 (断熱版も比較) |

**キャビティの検査体積は `cav_outer` + `cav_floor` + `cyl_side` + 開口面**。`cyl_top` は外部流に
面していて CV を囲まないので、熱収支には入れず**別枠で報告**する (計画 §4.8)。

## ワークフロー

```
case.json --setup.py --resolve--> manifest.json
cad/build_geom.py (FreeCAD)  --> cavity_fluid_half.step + geom_used.json
precursor/ (2D 平板 M5 断熱) --> inlet_profile_1.csv (delta≈5mm の z 分布)
cad/mesh_salome.py (SALOME NETGEN + ViscousLayers) --> .med
cad/med_to_msh41.py --> forge.msh --convertGmshToForge--> forge.h5
   -> check_mesh_quality.py + tools/check_mesh_extra.py (体積/sliver/dual)
gen_runs.py --> S0..S6 段階起動 (AWS A10G) --> res_*.h5
tools/cavity_eval.py --> 温度/熱流束/侵入深さ/開口流束 + 時系列 CSV
   -> check_convergence.py (全列) + check_quasisteady.py --series-csv
```

- **メッシュ生成はローカル** (Salome は AWS 未導入)、**計算は AWS g5.xlarge (A10G 23 GB)**。
- 変換は**最終形の壁タグ**で行う (SST の `wall_dist` のため)。段階起動 S0 の全面 slip は
  実行時に `bcondConfig.yaml` を差し替えて実現する。
- **3D の cross-mesh restart は使わない** — `interp_field.py` は 2 次元最近傍で、深さ 50 mm の
  キャビティでは z が無視される (計画 §3.2)。Stage A / Stage B はそれぞれ S0 から立てる。

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `precursor/run_0001_precursor_m5` | 2D 前駆 第 1 版 (H=0.05 m)。M5 断熱平板 | 前縁圧縮波が上面 slip で反射し x≳220 mm の外縁判定が破綻 (δ\* が乱れる) | 破棄 (H 不足) |
| `precursor/run_0002_precursor_m5_H100` | **2D 前駆 正本** (H=0.10 m, y₁=6 µm, 段階起動 + cfl 2 × 24000) | x=313.45 mm で **δ₉₉ 5.056 / δ\* 2.899 / θ 0.2243 mm, Re_θ 2085, μt/μ 39 (乱流), T_w 1166.5 K (CPG T_aw 1187.5 の −1.8 %)**。ここから入口 CSV を作成 | active (**入口分布の供給元**) |
| `run_0001_stageA_cpg` | Stage A (69k 節点) の疎通・段階起動レシピ確立。CPG, 旧バイナリ | 全段 NaN 0 完走 (6.24 ms/step)。`qwall` は全点 0 (診断未実装) | ref (レシピ確立) |
| `run_0002_stageA_qwall` | 同上を **qwall 診断入りバイナリ**で再実行 | **壁熱流束が取得可能に**。断熱 `plate` は厳密 0、等温壁 q'' 平均 8.3–8.5 kW/m²、**壁温 (断熱平板) 最大 1185.4 K = CPG T_aw の −0.18 %**。y⁺ 0.66–0.73 (mean) | active (**Stage A 基準**) |
| `run_0003_stageB_cpg` | **Stage B** (330k 節点, y₁=10 µm, VL 11 層) | 全段 NaN 0 (21 ms/step)。**キャビティ 3 壁 総入熱 74.8 W (全周)**、q'' 平均 4.89 kW/m²・最大 242 kW/m²、**h_ref 93 W/m²K** (基準温度基準)。開口ガス 696 K / 中央 515 K / 底 535 K。侵入 (25 K) 49.5 mm。**報告量は全て STEADY** (drift ≤0.6 %)。y⁺ 0.12–0.13。図 `cavity_fields.png` / `cavity_profile.png` | active (**現行の主結果**) |
| `run_0004_stageC_cpg` | **Stage C 格子収束 3 点目** (589k 節点, y₁=15 µm, VL 7 層, 内面 0.6 mm) | 計算中 | active |
| `run_0005_stageD_cpg` | **Stage D** (741k 節点, 内面 0.45 mm, 継ぎ目比 0.239, growth 0.15) — ユーザ要望 (内面細分・継ぎ目の段差解消・全体細分) を反映 | 待機中 | active |

> メッシュ本体 (`cad/*.med`, `cad/*.msh`, `mesh/*.h5`) と run 成果物は git に入れない。
> run を作成・破棄したらこの表を必ず同期する (命名 `run_NNNN_<slug>`)。

## メッシュ一覧

| メッシュ | 節点 / 要素 | 品質 | 備考 |
| --- | --- | --- | --- |
| `mesh/stageA.h5` | 69,005 / tet 105,981 + prism 91,399 | AR 184.7 / skew 0.821 **PASS** | VL 7 層・第一層 40 µm。疎通と URANS 用 |
| `mesh/stageB.h5` | 329,868 / tet 367,755 + prism 501,017 | AR 382.0 / skew 0.881 **PASS** | VL 11 層・第一層 10 µm・リップ 0.45 mm。dual 体積 relErr 5.4e-11, 閉性 1.2e-6 |
| `mesh/stageC.h5` | 589,437 / tet 862,563 + prism 814,282 | AR 206.8 / skew 0.874 **PASS** | VL 7 層・第一層 15 µm・内面 0.6 mm |
| `mesh/stageD.h5` | 740,763 / tet 1,514,642 + prism 878,730 | AR **127.6** / skew 0.849 **PASS** (平均 AR 6.1) | VL 6 層・第一層 20 µm・**内面 0.45 mm**・継ぎ目比 0.239・growth 0.15 |

**メッシュ生成の要点 (実測 2026-09-19, 4 例)**: **凸角 (開口リップ) まわりの接線セルサイズが VL 総厚を
下回ると、層が自己交差して NETGEN の tet 充填が無言で失敗する** (`Compute: True` なのに tet 0 で
体積の 2.2 % しか埋まらない)。エッジだけに下限を課しても不十分で (リップ 0.9 mm × 面 0.5 mm × VL 0.9 mm で失敗)、
**壁面サイズ全部に下限 = VL 総厚**が要る。`setup.py` がこれを課し、引き上げた面を manifest の
`size_bumped_by_vl` に残す。**帰結**: 内面を細かくするほど VL も薄くせざるを得ず、
継ぎ目比 (最終層厚/壁面サイズ) は (s−1)/s で頭打ちになる。第一層 y₁ は独立なので y⁺ は保てる。

**熱伝達率の基準温度**: `h_ref = q''/(T0_ref − T_w)` の **T0_ref = すきま中央面の最近傍点の総温**
(中央面の底は「すきま幅の半分」で切り取り) を主とする (ユーザ指定 2026-09-19)。深部では回復温度基準の
`h_aw` が見かけ上 0 に落ちるため。総温は `VALUE/h0` から `total_quantities.py` で逆算。

## EOS (CPG / TP)

`case.json` の `conditions.gas` で切替 (`gen_runs.py run --gas TP` でも上書き可)。

- **CPG** — パイプライン疎通と段階起動レシピ確立用。速い。$T_{aw}$ を +61 K 過大に出す。
- **TP (semi-perfect, 生産)** — 乾燥空気を **1 擬似種 MIXDRY** にまとめた NASA-9
  (`species_db.yaml`, `forge_design.gas.semiperfect.mixture_pseudo_species` で生成, R=287.048)。
  `thermalMethod: 2` + `thermoHrefTemp: 298.15`。多成分 TP の implicit 不安定を避けるため 1 種に畳む。
  **IC の内部エネルギーは TP の datum で組み直す** (CPG の $p/(\gamma-1)$ は datum が違う)。
  自由流密度も TP の R で引き直す。

## 偏心スイープ

```bash
python3 sweep_offset.py --offsets 0 0.5 1.0 1.5 --stage D
```

CAD → メッシュ → 変換 → 段階起動 → 評価を各偏心量で通す。VL 層数は最小すきまから自動で引き直される
(例: すきま 2.5/2.0/1.5 mm → 6/6/5 層)。`--dry` で設定だけ確認できる。

## ソルバ側の変更 (本ケースで入れたもの)

- **低 Re 壁の `qwall` / `utau` 診断出力** (`solver_density_cuda/cuda_forge/viscousFlux_d.cu`):
  従来 `wallTreatment==0` では `qwall_b`/`utau_b` に誰も書かず、壁面ダンプの `qwall`/`utau`/`ypls` が
  **全点 0** で熱流束・熱伝達率・$C_f$ が取れなかった。残差に入れたのと同じ解像熱流束を同じ符号規約
  (壁→流体が正) で格納するようにした。モデル経路 (`wallTreatment` 1/2, `sstEnergyWallFunction`) では
  `qwall_b` は入力側なので触らない → **解はビット不変**。検算: 断熱壁で厳密 0、等温壁で負。
