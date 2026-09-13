# case/34 Arthur (1952) 2D source-flow 極超音速ノズル — dry 超音速膨張

P.D. Arthur の博士論文 (Caltech/GALCIT, 1952) で用いられた **2 次元 source-flow
ノズル**を再現し、**乾いた (dry) 非粘性の超音速膨張**を forge で計算するケース。
窒素凝縮を扱った論文
[`papers/on nitrogen condensation in hypersonic nozzle flows_summary.md`](../../papers/on%20nitrogen%20condensation%20in%20hypersonic%20nozzle%20flows_summary.md)
の検証ノズル (Arthur 1952) に対応する。

> **注意**: forge には核生成・液滴成長など**相変化 (凝縮) モデルは無い**。本ケースが
> 再現するのは「凝縮なし」の等エントロピー的膨張 (Arthur Fig.4/11 の dry 1D 理論線) まで。
> 凝縮による静圧上昇そのものは対象外。

## ノズル形状 (出典: Arthur 1952, II.A Apparatus)

- 2 次元 source-flow ノズル: スロート全高 **0.010 in = 0.254 mm**、出口全高
  **~1 in = 25.4 mm**、幅 1 in 一定、**全開き角 11° (half-angle 5.5°)**、面積比 **A/A*≈100**。
- 壁は仮想ソース点から出る直線 (source flow)。**スロート曲率は原典に記載が無い**
  (製作詳細は Ref.5 = Nagamatsu & Willmarth, GALCIT Memo No.6 / JAP 23(10) 1089, 1952 に委ねられ、
  入手できず)。source-flow は本来直線壁なので曲率は未定義。
- 本ケースは **発散部のみ・上半分**をモデル化 (中心線 y=0 を対称 slip)。スロート直下から
  始め、**入口に超音速インレット (M=1.05) を置いて choking を回避**。
- 壁形状はメッシュ上 **双曲線 `y(x)=√(yt²+(x·tanα)²)`** で表現:スロートで水平接線
  (dy/dx=0) となり鋭い凸コーナーの膨張特異点を回避、下流は 5.5° 直線 (Arthur の source-flow)
  に漸近。実効スロート曲率 R_t = yt/tan²α ≈ 13.7 mm (滑らかな代表値)。

メッシュ生成: [`mesh/arthur_nozzle.geo`](mesh/arthur_nozzle.geo) (gmsh, 400×60 の 1 層押し出し,
24,000 hex)。境界 physID: 1=inlet, 2=outlet, 3=wall(上壁+中心線+z両面, 全て slip)。

## 計算条件

- 貯気: P0 = 844,037 Pa (8.33 atm), T0 = 290 K (Arthur Run 34-1 / Fig.4,11 相当)。窒素 γ=1.4。
- 非粘性 Euler (viscMethod 0)、SLAU、1 次精度、陰解法 block-DPLUR。
- 初期場 `arthur_n2` ([`solver_density_cuda/input/setInitial.hpp`](../../solver_density_cuda/input/setInitial.hpp)):
  スロート超音速状態 (M=1.05) の一様場。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_dry` | 直線壁(スロート鋭頂点) + cfl_pseudo=5 | step 467 で発散 (スロート上端の凸コーナー膨張特異点; Ux~1.5e8)。`res_nan_467.h5` | 破棄 (コーナー特異点の記録) |
| `run_0002_slau_dry_cfl1` | 双曲線スロート(滑らか) + cfl_pseudo=1 + **1次精度** (convMethod 0) | **収束** (rms 6.7–6.9 桁低下)。出口 **M=6.75** (等エントロピー 6.93)、A/A*=99.6、P/P0=2.79e-4 (理論 2.58e-4)。`postproc_arthur.png`, `residual_history.png` | active (1次基準) |
| `run_0003_slau_2nd` | run_0002 の収束場を引き継ぎ **2次精度** (convMethod 1, limiter 2 Venkatakrishnan) | 出口 **M=6.87** (等エントロピー 6.93、誤差 2.5%→**0.8%**)、P/P0=2.71e-4。Mach が A/A*=100 まで等エントロピー線にほぼ完全一致。残差は **~3e-4 でプラトー** (リミタ起因のリミットサイクル; 場は準定常)。`postproc_arthur.png`, `residual_history.png` | active (高精度) |
| `run_0004_cond_skeleton_dry` | 非平衡凝縮 Phase 1 (`condensation:1, nCondSpecies:1`)。run_0003 と同条件で **4 モーメント (ρg,ρQ2,ρQ1,ρQ0) を受動スカラー (ソース=0) として輸送**する骨格の回帰確認 ([plan](../../plans/active/condensation-nonequilibrium.md)) | **dry 回帰一致**: 出口 **M=6.874**, P/P0=2.71e-4 で run_0003 と同一。液相モーメントは全セル厳密 0、NaN/Inf なし。cond ON vs OFF の場差 (ro maxrel 9.5e-4) は cond OFF を 2 回回した run-to-run ノイズ (1.0e-3, リミットサイクル+GPU atomic 非決定) **以下**=凝縮スカラーは気相に不干渉。`postproc_arthur.png`, `residual_history.png` | active (Phase 1 骨格) |
| `run_0005_cond_off_ref` | 上の回帰用 cond OFF 参照 (現行バイナリ, `condensation:0`)。run_0004 とのビット差が plateau ノイズ由来であることの対照 | run_0004 との場差は run-to-run ノイズ内 (上記) | ref (対照, 成果物は破棄可) |
| `run_0006_cond_n2` | 非平衡凝縮 Phase 2 (`condensation:1`, CPG 二相 EOS + 核生成/成長 + **二相エネルギー流束 (全エンタルピー保存)**)。run_0003 収束 dry 場から restart、`restart_dry.h5` | **論文 Fig.2 検証成功**: 凝縮で壁面静圧が dry より上振れし、cond/dry 比が論文と 1–2% 以内一致 (3in 1.21x/4in 1.33x/5in 1.43x vs 論文 1.20/1.31/1.45x)、onset バンプも再現。全エンタルピー保存 (dev −0.0%)、g max 0.08・mean 1.0% (論文 ~0.75–2.3% 範囲)、NaN なし。`compare_fig2_htfixed.png`, `residual_history.png` | **active (N2 検証済)** |
| `run_0007_cond_develop` | 同上を**一様初期から膨張+凝縮を同時 develop** (1次精度, restart でなく `arthur_nozzle.h5`)。over-condensation が restart アーティファクトか確認 | g max 0.21・mean 0.023 で run_0006 とほぼ同じ → **over-condensation は restart 由来でなく定常解そのものの過大予測** (定常で S~100 のまま=潜熱昇温が膨張冷却に追いつかず、レートが過大)。定量一致 (paper の壁圧トレース照合・レート較正) は後続課題 | ref (診断, 成果物は破棄可) |

> `run_0002` の VERDICT は形式上 NOT CONVERGED と出るが、これは**面外 z 方向の `rms_roUz`
> (≈1e-14 のノイズ; 単層押し出しで物理が無い) のみ**による偽陽性。物理 4 残差
> (`rms_ro`,`rms_roUx`,`rms_roUy`,`rms_roe`) は 6.7–6.9 桁低下し平坦=実質完全収束。

| **空気凝縮 (CPG carrier 形) と N2 低温物性の整合 (2026-09-12/13, branch feature/condensation-air, plan [condensation-air](../../plans/active/condensation-air.md))** — run_0006 のプロトコル (dry 収束場 restart, Euler, `convMethod 1`, cfl_pseudo 1) を **12000 step・1000 step 毎保存**にし、cell (押し出し 1 層) と node (平面メッシュ `mesh_node/`, 24000 双対 CV, QC PASS AR 21.9 / skew 0.07) の両方で。判定は `onset_analysis.py --series` (壁セル列の p/p_dry @3/4/5 in, 中心線 onset (Δp/p_dry>1 %), 出口 u_n/c>1)。全 run NaN 0・series STEADY・出口超音速。残差は 2 次のリミットサイクル plateau (既知) で `check_convergence` は NOT CONVERGED → 準定常比較。図 `compare_air_n2_wall.png` | | | | |
| `run_0008_ref_n2` | 旧バイナリ (0512823d) で run_0006 再現 (8000 step, 4000 毎保存) | N2 旧物性 | onset 2.37 in / 678 Pa / 37.85 K, Ṗ 1.7e4 /s, cond/dry @3/4/5 in 1.240/1.357/1.451, g_exit 0.081 | ref (旧バイナリ参照) |
| `run_0013_air_dry` | 空気 dry (**入口速度を N2 のまま = M_in 1.068 で誤り**, codex 指摘) | dry | 破棄予定 (run_0023 に置換) | 破棄予定 |
| `run_0014_n2_ref_cell` | 新バイナリ + 旧物性 (`condN2LatentLowT 0, condN2PsatLowT 0`) — SLAU 面状態一貫化・slip 状態保持の影響 | N2 旧物性 | onset 2.367 in / 678 Pa / 37.85 K (run_0008 と同一), cond/dry 1.241/1.357/1.451, 場差 ≤1e-3 (case/34 の反復ノイズ水準) | active (R0 cell) |
| `run_0015_dry_slip_cell` / `run_0015o_dry_slip_cell_oldbin` | 凝縮 OFF: slip ghost の状態保持 (新) vs 旧バイナリ | dry | 場差 ρ 3e-4, U_y 2e-3 = 既知の run-to-run ノイズ (~1e-3) 水準 → slip 変更は dry に影響なし | 破棄予定 (回帰記録) |
| `run_0016_n2_latent_only` | 潜熱だけ新 (飽和圧は旧 C–C; `condN2PsatLowT 0`) 診断 | N2 | onset 不変 (2.367 in, 37.85 K; S は飽和圧で決まる), **g_exit 0.081→0.059** (潜熱 +25 % で放出熱増), cond/dry **1.158**/1.345/1.446 (3 in が実験 1.133 に近づく) | active (R2) |
| **`run_0017_n2_new`** | 潜熱+飽和圧 新 (既定, c_l 2000) | N2 新物性 | **onset 2.112 in / 800 Pa / 39.67 K** (飽和圧 0.52 倍→S 増で 0.25 in 上流), Ṗ 1.9e4, cond/dry 1.242/1.371/1.476, g_exit 0.065。理論 onset 線 (Ṗ=20000, 800 Pa: 37.5 K) に対し **+2.1 K** (旧 +1.0 K) | **active (R3, N2 新既定)** |
| `run_0018_n2_cl15` / `run_0019_n2_cl25` | c_l 1500 / 2500 J/kg/K 感度 | N2 | onset 38.95 / 40.25 K (2000: 39.67) → c_l ±500 で ∓0.7 K; cond/dry @3 in 1.2375 / 1.246 | active (R4) |
| `run_0023_air_dry` | 空気 dry 基準 (二成分 0.79/0.21: R 288.19, cp 1008.7; IC = N2 dry を ρ×R_N2/R_air, u×√(R_air/R_N2); 入口 ρ 6.161, U 324.48) | dry (空気) | 12 スナップショット STEADY (壁圧比 1.000 ±0.1 %) | active (E2) |
| **`run_0024_air_cpgcarrier`** | **空気 CPG carrier 形** (`condVaporMassFraction 0.7671`: N2 選択凝縮 + O2 キャリア, 新物性) | 空気 | **onset 2.207 in / 750 Pa / 38.94 K**, Ṗ 1.8e4 → 理論線 (37.3 K) +1.7 K。N2 新物性 (39.67 K) より 0.7 K 低温 (S=y_N2 p/p_sat が小さい分; Daum & Gyarmathy の「空気 ≈ N2」と整合)。cond/dry 1.214/1.351/1.455, g_exit 0.059 (≤Y_w), T_min 39.2 K | **active (E1, 空気凝縮の代表)** |
| `run_0021_dry_slip_node` / `run_0021b_dry_slip_node` | node 平面メッシュ dry N2 (IC = cell dry を interp_field) 旧バイナリ / 新バイナリ | dry | 旧: PASS (出口 M 6.93 = 等エントロピー)。新 vs 旧の場差 ≤8.6e-6 (node ノイズ水準) → slip 変更は dry に影響なし | active (node dry 参照 / 回帰) |
| `run_0020_n2_ref_node` / **`run_0022_n2_new_node`** | node: 旧物性 / 新物性 (IC = run_0021 res_12000) | N2 | onset 37.89 K@680 Pa / **39.59 K@793 Pa** (cell 37.85 / 39.67 と 0.1 K 以内), cond/dry 1.245/1.364/1.456 / 1.248/1.379/1.481 (cell と ≤0.6 %) | active (node R0/R3) |
| `run_0025_air_dry_node` / **`run_0026_air_cpgcarrier_node`** | node: 空気 dry / 空気 CPG carrier | 空気 | onset **38.88 K@744 Pa** (cell 38.94@750), cond/dry 1.218/1.357/1.459 (cell と ≤0.4 %) → node/cell 一致 | active (node E2/E1) |

**Arthur 実験との比較 (壁 cond/dry @3/4/5 in; 実験 1.133/1.250/1.500, Lin 2014 計算曲線 1.200/1.308/1.447)**: 旧物性 +9.5/+8.6/−3.2 %、新物性 +9.6/+9.7/−1.6 %、潜熱のみ新 +2.2/+7.6/−3.6 %、空気 +7.1/+8.1/−3.0 % (いずれも実験比)。
**従来「1–2 % 一致」と書いていたのは Lin の計算曲線との差** (旧 +3.4/+3.8/+0.3 %) で、Arthur 実験記号とは 3–4 in で ~9 % 過大 (codex 指摘 2026-09-12)。
潜熱・飽和圧の整合修正は onset を 0.25 in 上流 (+1.8 K) に動かし、実験との差は 3 in で不変・4 in で +1 pt・5 in で −1.6 pt 改善。

## 結果まとめ

- 出力先: `case/34.arthur_n2_nozzle/run_0002_slau_dry_cfl1/`
- 主要図:
  - [`run_0002_slau_dry_cfl1/postproc_arthur.png`](run_0002_slau_dry_cfl1/) — 左: Mach vs A/A* (Arthur Fig.4 相当),
    右: P/P0 vs x (Arthur Fig.11 相当)。forge の中心線・壁が等エントロピー線に一致し、
    かつ中心線と壁が互いに一致 (Arthur Fig.11 の知見を再現)。
  - [`run_0002_slau_dry_cfl1/residual_history.png`](run_0002_slau_dry_cfl1/)
- 後処理: [`postproc.py`](postproc.py) (`python3 postproc.py run_0002_slau_dry_cfl1`)。
- forge 出口 Mach の等エントロピー (6.93) との差は **1次 6.75 (2.5%) → 2次 6.87 (0.8%)** と縮小。
  1次の不足は数値拡散によるもので、2次精度化 (`run_0003`) でほぼ理論線に一致する。
  → A/A*≈100 の Arthur 形状の dry 出口 Mach ~6.9 は妥当。
- **マッハ数の解釈** (重要): A/A*≈100 → 等エントロピー dry で M≈6.93 は正しい。これは「凝縮しなければ」
  の値で、実機・論文では**凝縮の潜熱で M はこれより低下**する (論文 parametric study で凝縮あり
  Ma 5.40–6.67)。論文の「設計 Ma=6.0」ノズルは A/A*≈53 の**別形状**で、本ケース (Arthur 検証ノズル,
  A/A*=100) とは異なる。
- 2次精度の残差プラトー (~3e-4) は Venkatakrishnan リミタのリミットサイクルで、forge で散見される
  挙動。出口諸量は準定常で安定しており Mach 比較には使える。

## 実行手順

```bash
cd case/34.arthur_n2_nozzle/mesh
gmsh -3 arthur_nozzle.geo -o arthur_nozzle.msh -format msh4
cd ../run_0002_slau_dry_cfl1
cp ../mesh/arthur_nozzle.msh .
/home/sano/work/forge/solver_density_cuda/build-native/convertGmshToForge arthur_nozzle.msh arthur_nozzle.h5
/home/sano/work/forge/solver_density_cuda/build-native/forge
python3 ../postproc.py .
```
