# case/48.cabra_h2n2 — Cabra H₂/N₂ 浮き上がり火炎 (vitiated coflow, 燃焼加熱器の検証ケース)

Cabra et al. (UC Berkeley; NASA/CR-2004-212887 Table 6.1, [`papers/combustion/`](../../papers/combustion/)): 中心ジェット d=4.57 mm
(管 OD 6.35 mm) から H₂ 25 %/N₂ 75 % (体積)・305 K・バルク 107 m/s (Re 23,600) を、希薄 H₂/空気予混合火炎の燃焼ガス coflow
(D=210 mm, 1045 K, 3.5 m/s, X O₂ 0.15 / H₂O 0.099 / N₂ 0.751) へ噴射。自着火で浮き上がり高さ H/d≈10、火炎長 30 d。
公開データ (TNF): 中心軸 (z/d=1–34) と半径分布 (z/d=1,8,9,10,11,14,26) の Favre 平均 T・主要種・OH・NO。ここでは NASA 報告書の図を
目視で読み取った `exp_centerline.csv` / `exp_radial_T.csv` を参照にしている (±3 %)。加熱器の「高温 vitiated 流への H₂ 噴射→自着火」と同じ物理で、
**低マッハ前処理 + 反応** の経路と着火位置予測を検証する。計画: [`plans/active/chemistry-finite-rate-h2.md`](../../plans/active/chemistry-finite-rate-h2.md)。

## セットアップ

- `mesh/make_cabra_mesh.py`: 軸対称平面構造 quad (node)。ジェット管を上流 20 mm 延長 (管内 no-slip、入口は 1/7 乗則プロファイル)、
  リップ 0.89 mm、外周 R=100 mm slip、x=0–250 mm。60k ノード、品質 PASS (AR max 70, skew 0)。
- `setup_cabra_case.py run_dir [--chem 0|1] [--ji 5] [--tci 0|1] [--cfl] [--conv] [--relax] [--nstep] [--restart] [--mesh mesh/cabra_ext.msh]`: Cantera で入口状態、
  BC (`inlet_uniformVelocity`×2 + profile, `outlet_statPress`, `wall`, `slip`, `axis`)、config (node 軸対称, `nodeInletCornerWall: 1`,
  `lowMachPrecond: 2`, implicit, SST, TP 9 種, chemistry `jacobianInterval 5`)、IC (coflow 一様 + ジェット柱)。
  切り分け用オプション: `--single 1` (N₂ 単一種), `--jetcof 1` (ジェット組成=coflow), `--jetn2 1` (純 N₂ ジェット), `--iccol 0` (ジェット柱 IC なし=ピストン起動),
  `--planar 1` (平面・軸→slip), `--precond 0|1|2|3`, `--coupling 0|1|2`, `--far slip|outlet_statPress`, `--sfr`, `--nojet`。
- `compare_cabra.py run_dir step out.png`: 中心軸 T/Y と半径 T を実験と重ねる。

## 計算 run 一覧

このケースで **forge の低マッハ node 経路のバグを 3 件特定・修正** した (2026-09-05, 詳細は
[plan lowmach 変更ログ](../../plans/accepted/time_integration-lowmach-preconditioning.md) と
[plan outlet §2.12](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)):
(A) node × `lowMachPrecond>=2` の前処理 block カーネルが境界ノード (軸・壁・出口・遠方) を凍結、
(B) TP 亜音速 `outlet_statPress` の γ 混用で低マッハ出口が一様逆流、
(C) 多成分 × 定常擬似時間 × precond で組成前線と密度前線の擬似速度が不整合 → `speciesPrecondDt` (既定 1)。
`run_0001`–`run_0035` は (A)(B)(C) すべてを含む旧バイナリ、`run_0036`–`run_0046` は (A) 修正後、`run_0047`–`run_0061` は
(A)(B) 修正後、`run_0062` 以降が (A)(B)(C) 修正後。全て化学 OFF (混合のみ)、3000 step 以下の切り分け run は破棄予定。

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_mix_start` | 段階起動 1 次 cfl_pseudo 0.5 (旧バイナリ) | 完走するが後述の凍結境界・偽逆流を含む | 破棄予定 |
| `run_0002`–`run_0013` | 2 次化 / cfl 2 / precond 0 / coflow のみ / 単一種 / ジェット柱なし / coupling 1,0 / 平面 の切り分け (旧バイナリ) | 多成分は全て P 70–230 kPa・T>coflow の非物理場 (`run_0003` outlet∩farfield で NaN)。coflow のみ (`run_0008`) と単一種 N₂ (`run_0009`, 柱 IC) は健全 → 「多成分固有」と誤認した (実は境界凍結 (A) が単一種では見えなかった) | 破棄予定 |
| `run_0014`–`run_0029` | 400–800 step 短期診断: 有界性 / 1 次 / sfr / 純 N₂ ジェット / precond 0,1 / relax 1 / RK3 (局所 dt・時間精度) / cfl 0.1 / 一様組成 | 化学種は有界 (Y_H2 ≤ 入口値, ΣY=1)。時間精度 RK3 (`run_0024`/`run_0025`, 6000 step 120 µs) は健全 = 空間スキーム無罪。implicit は初期 400 step の圧縮過渡は物理 (ピストン) | 破棄予定 |
| `run_0030`–`run_0035` | 単一種ピストン起動で軸ノードの追従を確認 (precond 2 / axisymMethod 1 / centroidShift 0 / precond 0) | **precond 2 では軸ノードが IC (1045 K, 4 m/s) のまま隣接 (305 K, 80 m/s) に追従せず、precond 0 (`run_0035`) と explicit (`run_0026`) は追従** → 真因 (A) を特定 | 破棄予定 (エビデンス) |
| `run_0036`–`run_0046` | (A) 修正後: 単一種/多成分, 柱 IC あり/なし, precond 3, cfl 0.1, farfield=outlet, 出口列の早期追跡 | 軸は追従するが**出口列が全体で一様に逆流** (−50 m/s @100 step, P 単調上昇, 平面・precond 0 でも同一 `run_0045`/`run_0046`) → 真因 (B) を特定 (`run_0044` 出口列時系列) | 破棄予定 (エビデンス) |
| `run_0047_fix2_outlet_early` | (B) 修正後の出口列 200 step | 出口列 3.5 m/s / 101.33 kPa を保持 | 破棄予定 (エビデンス) |
| `run_0048_fix2_single_col` | (A)(B) 後, 単一種 N₂, 柱 IC, cfl 0.5, 3000 step | **健全**: P 101.0–102.0 kPa, T 305–1045 K, 軸追従 | ref (単一種基準) |
| `run_0049`–`run_0057` | (A)(B) 後の多成分: coupling 0/1/2, 一様組成 9 種, 平面, relax 1, control 0 | 異組成ジェットは**軸上ジェット核で P 148 kPa・T 1150 K に暴走** (1000 step)。一様組成 (`run_0053`) は健全、coupling 依存なし、relax 1 で悪化 → 多成分固有の第 3 の問題 | 破棄予定 (エビデンス) |
| `run_0058`–`run_0061` | dual-time (`unsteady 1, dualTime 1`, dt 2e-7 / 1e-6) と precond 0 定常 | **dual-time (`run_0059`/`run_0061`, 2 ms) と precond 0 (`run_0060`) は健全** → precond × 分離化学種の擬似時間不整合 (C) を特定 | ref (dual-time 代替の実証) |
| `run_0062_spdt1` / `run_0063_spdt2` / `run_0064_spdt0_bitcheck` | (C) `speciesPrecondDt` 1 / 2 / 0 の A/B (多成分, 柱 IC, cfl 0.5, 1500 step) | **mode 1・2 は P 100.2–102.6 kPa, T≤1045 K で健全**; mode 0 は 91–148 kPa / 1144 K で暴走 → mode 1 を既定に | ref (C の検証) |
| `run_0065_mix_col_cfl1_6k` | 全修正後の混合 (化学 OFF): 2 次, cfl_pseudo 1, relax 0.5, 柱 IC, 6000 step | 有界 (P 99.3–102.2 kPa, T 304–1045 K, NaN なし)、ジェットは混合中 (軸 Y_H2 x/d 20–50 で 0.0234→0.0143, T 655 K)。`check_convergence`: NOT CONVERGED (発達途上: rms_roUx/roe falling, rms_ro/roK/roOmega plateau) | ref (混合の起点) |
| `run_0066_mix_col_cfl2_cont` | run_0065 res_6000 から cfl_pseudo 2 で継続 | **step 932 で NaN 発散** (cfl_pseudo 2 は不可; 位置は `forge_run.log`) | 破棄予定 |
| `run_0067_mix_col_cfl1_cont` | run_0065 res_6000 から cfl_pseudo 1 で 20000 step 継続 (計 26000 step) | 有界 (P 101.2–102.1 kPa, T 304–1045 K, NaN なし)。軸 Y_H2: x/d 10 0.0178 / 20 0.0076 / 30 0.0046 / 50 0.0044 (下流はまだ発達中)。`check_convergence`: NOT CONVERGED (rms_ro/roe/roK/roOmega falling, rms_roUx/roUy plateau 4–6e-6)。次: さらに継続 → 反応 ON | active (混合場の最新) |
| `run_0068_react_cfl1` | **反応 ON** (jac 2, ji 5, coupling 2, tci 0), run_0067 res_20000 から cfl_pseudo 1, 10000 step | 有界。step ~6000 で x/d≈20 に自着火 (OH 3e-3)、火炎基部が x/d 10.7 (8000) → 7.7 (10000) と上流へ移動中、Tmax 1930 K。`check_convergence` NOT CONVERGED (過渡) | ref (着火の起点) |
| `run_0069_react_cont` | run_0068 res_10000 から 20000 step 継続 (tci 0 = 層流化学) | 火炎基部が x/d 4.6 (4000) → 3.4 (12000) → **0 (16000 以降, リップ付着火炎)**。Tmax 1800–2020 K, 軸 T max 1593 K (z/d≈30)。`check_convergence` NOT CONVERGED (全列 plateau)。TCI なし RANS の典型 (平均温度で層流反応率を評価し着火過大) → PaSR へ。図 `compare_cabra_20000.png` | ref (tci 0 の結論) |
| `run_0070_react_pasr` | 同条件で **PaSR (`tci 1`, tauChem 1, C_mix 1)**, run_0067 混合場から 30000 step | 着火は同じく x/d≈21 (5000)、火炎基部の上流伝播が遅くなるだけで 7.8 (10000) → 2.9 (20000) → 1.2 (30000) と**ほぼ付着**。軸 T は tci 0 と同様 (T>600 K が z/d 12.7, 実験 14)。`check_convergence` NOT CONVERGED。図 `compare_cabra_30000.png` | ref (PaSR C_mix 1) |
| `run_0071_react_pasr_cmix4` | PaSR C_mix 4 (τ_mix 4 倍 → κ 小) の感度, 30000 step | 推移は C_mix 1 とほぼ同一 (21.2 → 7.9 → 4.4 → 1.9 → 1.4 → 1.2)。**PaSR の κ は基部の上流伝播を止めない** (近傍せん断層の平均場が既に可燃・高温)。NOT CONVERGED。図 `compare_cabra_30000.png` | ref (PaSR 感度の結論) |
| `run_0072_react_li` | **機構 A/B (plan §5.1 P0-3)**: run_0068/0069 と同一条件・同一 IC (run_0067 `res_20000`)、機構だけ Li 2004 (`tools/mechanisms/h2co_li2004_cantera.yaml`) に交換, tci 0, 30000 step | **付着推移は Jachimowski とほぼ同一**: 着火 x/d (Y_OH>2e-4) 20.8 (6k) → 10.5 (8k) → 7.2 (10k) → 3.3 (16k) → **0.0 (20k 以降 `ignx` が閾値上 0 で維持)** (Jach は累積 26k で 0)。**Li への交換だけでは付着は解消せず、平均場反応率 (TCI 欠如) が主因という仮説を強く支持** (機構 1 本の A/B・擬似時間非収束なので「確定」ではない)。NOT CONVERGED、出口流束・Tmax は DRIFTING、NaN なし | ref (**機構 A/B の結論**) |
| `run_0073_react_li_cont` | run_0072 `res_30000` (付着済) から Li 機構のまま +30000 step 継続。目的: 付着後の統計定常確認 | `ignx` 0・T_max 1569 K・出口 Y_H2O は STEADY だが、**出口質量流束は 0.028–0.044 kg/s (平均 0.036 = 理論値 0.037) を周期 15–20k step で振動し DRIFTING** (環状逆流は step 5000 以降消失)。z/d=9 の T は snapshot 間 ≤34 K で定常、**z/d=26 は最大 183 K 変動**で実験比較に使えない。全残差プラトー = 擬似時間リミットサイクル。`check_convergence` NOT CONVERGED | ref (付着後の準定常性の記録) |
| `run_0074_react_li_faropen` | run_0073 `res_30000` から、**外周 (physID 5) を slip → `outlet_statPress`** (Ps=Pt=101325, Tt 1045; codex レビュー 2 の切り分け順序 ①), Li, +30000 step | **改善せず・不採用**: 外周で Uy −18〜+9 m/s の出入りが暴れ (coflow 3.5 m/s に対し非物理)、出口正味流束 −0.004〜0.072 kg/s (末尾 35 ノード逆流) で DRIFTING。静圧固定の側面境界は低マッハでは数十 Pa の圧力差で大流速を作り物理的 farfield にならない。内部は z/d 9/26 とも末尾 ΔT 10–13 K まで静まり残差 falling。`ignx` 0・T_max 1564 K STEADY | ref (外周圧力境界は不採用) |
| `run_0075_react_li_precond0` | run_0073 `res_30000` から **`lowMachPrecond: 0`・cfl_pseudo 0.5** (外周 slip; 切り分け順序 ③ リミットサイクルが前処理起因かを見る), Li, +30000 step | **前処理起因ではない**: 出口質量流束 0.034–0.042 kg/s (±11 %, precond 2 の ±20 % より小さいが消えない) で DRIFTING、逆流なし、全残差プラトー。z/d 9 は ΔT ≤8 K で定常、z/d 26 は ΔT 63–116 K。`ignx` 0・T_max 1569 K・出口 Y_H2O STEADY | ref (precond 無罪の記録) |
| `run_0076_react_li_extdomain` | **領域拡張** (`mesh/cabra_ext.msh`: Lx 0.25→0.40 m, R 0.10→0.20 m, nx 480, ny_co2 85, 110k ノード; `setup_cabra_case.py --mesh`), run_0073 `res_30000` から `interp_field.py` で cross-mesh 補間 (拡張域は最近傍 = coflow 値), 外周 slip, precond 2, Li, 30000 step (切り分け順序 ②) | 品質: `check_mesh_quality` は node 変換 h5 を読めない (CONNE 形式) ため cell 変換で判定 → **PASS** (AR max 104, skew 0)。**step 4617 で NaN 発散** (`res_nan_4618.h5`): 発散ノード 151 個は x −4.5〜0.7 mm・r 2.2〜4.6 mm = **ジェット管リップの付着火炎根元**、拡張域 (r>100 mm, x>250 mm) は coflow 値で健全。nx 300→480 で出口近傍の軸方向格子幅が 0.20→0.07 mm と 3 倍細かくなりリップ近傍の格子条件が変わったのが疑わしい (領域拡張の A/B としては交絡)。rms_roe は step 3000 から緩やかに上昇し 4615 で急変 | 破棄予定 (発散記録) |
| `run_0078_react_li_extdomain2` | 同じ領域拡張 (Lx 0.40, R 0.20) だが**近傍解像度を旧メッシュに一致** (`--nx 354 --ny_co2 82`: 出口の軸方向初期幅 0.20 mm・r=30 mm 外側の初期幅を旧値に合わせ、外側だけ粗く延長), run_0073 から補間, cfl 1.0, Li, 30000 step | **拡張しても振動は消えない**: 出口質量流束 (R=0.2 の理論値 0.143) が 0.174 (5k) → 0.144 (10k) → 0.112 (15k) → 0.174 (20k) → 0.146 → 0.114 と**周期 15000 step で ±22 %、内側 r≤0.1 (0.036–0.050) と外側が同位相**、逆流なし、x=0.25 面と出口面で同値 (柱全体が一様に脈動)。z/d 9 ΔT 27–45 K、z/d 26 56–139 K。`ignx` 0・T_max 1572 K STEADY。→ 領域の閉じ込めではなく、**速度入口+静圧出口の低速 coflow 柱に対する擬似時間反復 (局所 dt) のピストンモード**と解釈 | ref (領域拡張の結論) |
| `run_0079_mixfrac_smoke` | **CMC Phase A smoke**: run_0067 混合場から化学 OFF・`mixfrac` ON (Bilger ξ・分散輸送・χ̃), 300 step (`setup_cabra_case.py --mixfrac 1`) | NaN なし。`xi` はカーネル vs numpy 3e-8 一致、Y_H2 由来 ξ との差 ≤0.068 は出口リップ層流層の差分拡散 (物理)。分散 0–0.042 (実現可能性違反 9e-13)、χ̃ ≤3.8e3 1/s、`check_mixfrac.py` | ref (Phase A smoke) |
| `run_0080_mixfrac_dev` | 同条件で 5000 step (分散の発達場) | せん断層で sqrt(ξ''²)/ξ̃ 0.36–0.49、χ̃ が z/d 5 の 531 → z/d 30 の 4.4 1/s (文献の桁)。実現可能性違反 3e-11。**CMC の初期場** | ref (混合場+分散) |
| `run_0081_cmc_frozen` | **CMC Phase B 凍結検証**: run_0080 から化学 OFF・CMC 受動 (`--cmc 1 --couple 0 --cmcchem 0`), 300 step (分散復元バグで ξ''²=0 = δ-PDF だが混合線保存の判定には無関係) | 初版は条件付き T 998–1087 K (保存形輸送のドリフト) → 非保存補正後 **全 node 1045.00 K**、`cmc_dY` ≤0.017 (差分拡散分)。コスト 403 ms/step (混合分率のみ 61) | ref (**Phase B 検証**) |
| `run_0082_cmc_react` | **CMC 反応 ON・結合 ON 初回 (couple 1 = source coupling)**: run_0080 混合場から Li 2004, jac 2, ji 5, `cmc {nEta 41, couple 1, chem 1}`, cfl 1.0, 6000 step (**注: リスタートの分散復元バグで ξ''²=0 から開始**) | 条件付き空間は着火 (T_Q max 1615 K, 前線 x/d 18.8→6.4) するが**平均場が燃えない** (T_max 1123 K, 軸 z/d 26 840 K vs 実験 1410, `cmc_dY` 0.07 定常) = source coupling の構造的欠陥 (条件付き状態が燃え切ると ω→0)。652 ms/step。NaN なし | ref (couple 1 不採用の根拠) |
| `run_0083_cmc_react_couple2` (削除済) | couple 2 (緩和ソース) の 3 版: 初期 Ỹ_pdf=0 → step 2 NaN; エンタルピー緩和の非整合 → step 2 NaN; 反応熱規約修正後も差分拡散差 (0.017) を反応熱換算し Q̇ ±1e11 で step 2 NaN | **couple 2 不採用** (plan §9 (5)) | 破棄済 |
| `run_0083_cmc_react_couple3` | **couple 3 (平均 Y,T を PDF 積分状態へ α=0.05/step でブレンド)**, 分散復元修正後, 他は run_0082 と同一 | 全置換版 step 14 NaN → α 版 step 19 NaN → 往復バグ修正後も rms_roe が 25 step で 300 倍成長 (リップの差分拡散非整合, `run_0085`) → 停止 | 破棄予定 (結合デバッグの記録) |
| `run_0084_cmc_ab_a0` / `_p0` | couple 3 α=0 (往復のみ) と couple 0 (受動) の 60 step A/B | 修正前: α=0 でも step 40 NaN (更新前原始量と更新後 ρ の混用) → 修正後: 残差が受動と 3 桁一致 (2.611e-6 vs 2.609e-6) | ref (往復整合の証拠) |
| `run_0085_cmc_diag` | couple 3 α=0.05 の 40 step 診断 (5 step 出力) | リップ直後 x=0.2 mm r=2.4–2.6 mm で T 884→510 K (ξ 0.34–0.46, T_Q 1045 = 未着火)、P max 124 kPa: **差分拡散の平均場が混合線から 130 K 以上ずれ、ブレンドと綱引き** | ref (差分拡散非整合の証拠) |
| `run_0086_mix_sc0` | 混合場を **`speciesDiffusionMethod: 0` (定数 Sc)** で作り直し: run_0080 から化学 OFF・mixfrac ON, 3000 step | 完走。run_0098 の同一 step と全量 ~1 K 一致 → **fp32 化学と Q restart 経路を長期で検証** | ref |
| `run_0100_cmc_c5_fp32_cont2` | run_0099 の続き (couple 5, さらに 16000 step): 平均場の加熱が飽和するかの確認 | step 10000 で中断 (couple 6 採用で不要)。軸 z/d 26 は 1153 K 止まり | 破棄予定 |
| `run_0101_cmc_c6` | **couple 6** (Y ブレンド + 燃焼領域ゲート付きの h̃→h̃_pdf 緩和, 陰的経路) を run_0099 の 16000 step 状態から 8000 step。平均場 T が PDF 診断 T に追従するか | 完走。**平均 T = PDF 診断 T に数 K で一致** (z/d 9 半径 +2 K)、火炎基部 x/d 9.4〜10.9 (実験 ≈10)。ただし軸 ξ が z/d 20→30 で 0.51→0.02 に崩落 (実験 0.62→0.41): `xi_flux_check.py` で ∫ρuξdA が z/d 5 で +60 %・z/d 30 で −70 % と**非保存 = 定常反復の柱モードの過渡** (dual-time run_0077 は ±8 % 保存) | ref |
| `run_0102_cmc_c6_dualtime` | run_0101 の状態 (場 + Q) から **dual-time (dt 1e-5 s, 4000 step = 40 ms)** で couple 6。柱モードを消して下流 (z/d ≥ 20) の ξ・T を実験と比較する本命 | **失敗**: 10 ms で条件付き場が消炎、平均 ξ が z/d 5 以降で ~0 に崩落 (ξ フラックス 1.6→0.4 g/s)、ξ=0 の coflow が 770 K まで過冷却 (couple 6 の熱持ち越し qdebt が原因 → 持ち越し廃止)。原因切り分けは run_0103–0106 | ref (失敗記録) |
| `run_0103_dt_mix_sdm0` / `run_0104_dt_mix_sdm1` | dual-time 純混合 (chem 0, cmc 0) を run_0077 から 600 step: sdm 0/1 の A/B | 両者ビット同一・ξ フラックスは消炎過渡を除き保存 → 輸送は無罪 | 破棄予定 |
| `run_0105_dt_cmc_passive` | dual-time + CMC 受動 (couple 0, 平均場は laminar 化学) を run_0077 から 600 step | ξ フラックス保存 (1.51〜1.60 g/s), 質量流量 37 g/s 一定 → CMC 輸送は無罪 | 破棄予定 |
| `run_0106_dt_c6_from0077` | dual-time + couple 6 (持ち越し廃止版) を run_0077 から (Q は混合線初期化) | step 200 で質量流量 172〜233 g/s・出口逆流・P ±5 % = α ブレンドの音響過渡 → 中断。couple 6 は dual-time 不可 | 破棄予定 |
| `run_0107_dt_c7_from0077` | **couple 7** (物理緩和 τ_c=1e-4 s のソース) dual-time を run_0077 から 2000 step (20 ms) | **成功**: P 101〜102 kPa, ξ フラックス保存 (1.5〜1.8 g/s), 条件付き火炎基部 x/d 22→8.4 (T>1200 K は x/d 9.7; 実験 H/d≈10), 半径 T 平均差 +27/+32/+20/+42 K (z/d 9/11/14/26), 軸 T は実験と重なる (`cmc_pdf_mean_2000.png`)。凍結域 (T<Tfreeze) に初期場の OH が残る欠陥 → 修正 (couple 7 は Tfreeze を通す) | ref (**CMC 検証の基準**) |
| `run_0108_dt_c7_cont` | run_0107 の続き (修正後バイナリ) 4000 step = 40 ms: T_c=1045 K の基準場、火炎基部の定常確認 | 完走。基部 x/d 8.1〜8.4 で定常、半径 T 差 +24〜+42 K。ただし ξ_Ω の離散偏り (希薄側 −0.5〜−2 %) が緩和ソース経由で燃料シンクになっており、ξ フラックスが z/d 50 で 1.30 g/s に減少 (修正は 71daa48e) | ref |
| `run_0110_cmc_tc1030` | T_c = 1030 K, 5000 step (50 ms) dual-time (修正前バイナリ) | 火炎が x/d 20→26 に後退して 15 ms で吹き飛び、その後 ξ が z/d 10 以降で 0.008 に崩落 (m_ξ 0.4 g/s) = **離散 β-PDF の 1 次モーメント偏り × 緩和ソース (1/τ_c) の燃料シンク** (plan §9 (17))。T_c 応答としては無効 | ref (失敗記録) |
| `run_0111_dt_c7_xifix` (AWS) | 修正 71daa48e (目標組成を混合線に沿って ξ̃ に整合) の検証: run_0108 の状態から 1500 step dual-time, AWS g5 (`~/forge-chem/case/48.cabra_h2n2/`) | 完走 (1500 step)。ξ フラックス z/d 50 が 1.30→1.44 g/s と回復 (シンク解消)、基部 x/d 9.3〜9.4、半径 T 差 −15/−6/+1/+46 K (z/d 9/11/14/26) | ref (**修正後の 1045 K 基準**) |
| `run_0112_cmc_tc1045` … `run_0116_cmc_tc1075` (AWS) | **T_c 応答曲線** (修正後): run_0111 の状態から T_c = 1045/1030/1060/1015/1075 K で各 3000 step (30 ms) dual-time, 直列 (`run_tc_sweep.sh`, ログ `tc_sweep_aws_2026-09-12.log`) | (実行中) | active |
| `run_0099_cmc_c5_fp32_cont` | run_0098 の step 4000 (場 + `cmc_Q_4000.bin`) から fp32 化学 (`cmc.fp32: 1`, 最適化後バイナリ) で 16000 step 継続。Q(η) restart 経路の初使用と fp32 の長期 A/B | (実行中) | active |
| `run_0087_cmc_sc0` | run_0086 から CMC couple 3 (α 0.05), sdm 0, Li, 6000 step | **step 322 NaN**: 未着火のままリップ (x 0.4 mm r 2.7 mm, ξ 0.31) の T が 1034→400 K, P 95 MPa。管壁伝導・Le≠1 で T は混合線 T(ξ) と一致しない → **h のブレンドは非物理 (couple 3 不採用)** | 破棄予定 (根拠) |
| `run_0088_cmc_c4_diag` | **couple 4** (組成のみ α=0.05 ブレンド + その反応熱 −Σc_sΔρY_s をエネルギーへ), sdm 0, 400 step (100 毎出力) | NaN なし、P 100–102 kPa、rms_roe 横ばい (~1.2)、T_Q max 1356 K (条件付き空間で着火開始)。リップ付近で max|ΔT| 349 K (受動対照 run_0090 と比較) | ref (couple 4 診断) |
| `run_0089_cmc_c4` | couple 4 本番 (旧バイナリ) | 偽発熱 (混合線に沿った H2O 変化の生成エンタルピー) + doChem 逆転バグで無効 → 停止 | 破棄予定 |
| `run_0090_cmc_passive_ctrl` | couple 0 (受動) 400 step: ブレンドの影響を分離する対照 | 基底流の過渡は最大 157 K (x 12 mm r 3.5 mm)。couple 4 (旧) との差 +309 K は偽発熱と判明 | ref (対照) |
| `run_0091_cmc_c4_diag2` | couple 4 (反応熱=混合線離脱分のみ, doChem 修正後) 200 step, `cmc_xiOm/xipdf` 診断付き | ξ_pdf=ξ (差 0.003), P 101–102 kPa, NaN なし, リップ T は受動比 −89 K (減速中) | ref (couple 4 修正版の診断) |
| `run_0092_cmc_probe` | `CMC_PROBE_NODE=18815` (リップ節点) で Q(η) をダンプ, 3 step | **初期化直後に内点 Q が 0/NaN** → doChem 引数逆転を特定。修正後は内点が混合線のまま保存 | ref (バグの証拠) |
| `run_0093_cmc_c4` | couple 4 (熱=混合線離脱分) 本番 | step 1487 NaN: 平均 T 3477 K (T_Q 1392) — 平均場の数値的ずれ補正まで発熱計上 | 破棄予定 (根拠) |
| `run_0091_cmc_c4_diag2` (再利用) | couple 4 の熱源 A/B 600 step: D_pdf 時間差分 → ±2e9 W/m³ の偽吸発熱; **条件付き反応熱の PDF 平均 (最終形)** → 非負、P 101–110 kPa、条件付き着火 x/d 39–41、平均 OH>2e-4 x/d 38 | ref (熱源の決定) |
| `run_0094_cmc_c4` | couple 4 最終形 (熱源=条件付き反応の PDF 平均) 本番 | step 1000: T_max 1506 K, 平均 OH>2e-4 x/d 15.9, 条件付き着火 x/d 24.8, **P 96–186 kPa** → step 1291 NaN (擬似時間の発熱が速すぎ流れが追従不能) | 破棄予定 (根拠) |
| `run_0095_cmc_c4_cap` | + 発熱リミッタ ΔT ≤ 10 K/step (持ち越し), dtScale 0.5 | step 2000 で P 83–147 kPa → step 2654 NaN | 破棄予定 (根拠) |
| `run_0096_cmc_c4_cap2` | couple 4, 上限 2 K/step, dtScale 0.25, 6000 step | step 3000 まで P 100–112 kPa, 平均 OH x/d 36.6 → step 4000 で T 6000 K / P 9–127 kPa → step 4300 NaN。**直接加算は不採用** | 破棄予定 (根拠) |
| `run_0097_cmc_c5` | **couple 5** (発熱を陰的化学ソース経路 `res_roe` へ, 上限 10 K/step), dtScale 0.5, 6000 step | **完走・NaN なし・P 100–102 kPa**。火炎基部 (OH>2e-4) x/d 22.7 (2k) → 18.4 → 16.6 → 14.7 → **13.8 (6k)** と減速しつつ上流へ、条件付き着火前線 14.1 が追従。平均 T_max 1063 K・軸 z/d 26 883 K (加熱は遅い)。`check_quasisteady`: T_max STEADY、着火位置は上流へ移動中 (DRIFTING)。半径 T は z/d 9 で平均差 +10 K、z/d 14 で −198 K、z/d 26 で −353 K (平均場の加熱が未完) | ref (**結合方式の決着**) |
| `run_0098_cmc_c5_long` | couple 5 の長期 run: run_0086 から 20000 step (2000 毎出力), Q 永続化 (`cmc_Q.bin`) 付き, double 化学 | 完走。火炎基部 (Y_OH>2e-4) は x/d 7.8〜14.5 を柱モードと共に往復 (DRIFTING)。**PDF 診断 T (`cmc_pdf_mean.py`) は実験と整合 (z/d 11 半径平均差 +5 K, 軸 z/d 26 1328 vs 実験 1410) だが平均場 T は −120〜−220 K** = couple 5 の熱受け渡しの構造欠陥 (plan §9 (13)) | ref |
| (注) `run_0079_react_li_dualtime_cont` | 別セッションが作成した run (本セッションの `run_0079_mixfrac_smoke` と番号衝突)。表の管理は当該セッションに委ねる | — | 他セッション |
| `run_0077_react_li_dualtime` | run_0073 `res_30000` から **dual-time** (`unsteady 1, dualTime 1`, dt 1e-5 s 固定 = ジェット部 CFL 6.5, nStepInner 5, cfl_pseudo 0.5, precond 2, 外周 slip), Li, 10000 step = 物理 0.1 s (coflow 通過時間 71 ms の 1.4 倍), 出力 10 ms 毎 (切り分け順序 ④) | **物理時間では定常化する**: z/d 9・26 の T は 30 ms 以降 10 ms あたり 0–4 K、最終場 (100 ms) と ≤3 K で一致 (擬似時間では同点が 183 K 動いた)。出口質量流束は 0.032→0.041→0.033→0.038 と周期 ~50 ms で減衰する緩やかな柱モード (±12→±7 %、平均 0.036 = 理論値 0.037)、逆流なし。`check_quasisteady`: ignx/tmax/exit_yout_H2O STEADY、exit_massflux DRIFTING (減衰途中)。**擬似時間反復では減衰しない柱の調整モードが物理時間では減衰** → Cabra の下流比較は dual-time 場で行う。`compare_cabra_10000.png`: 軸 T>600 K は z/d 13.3 (実験 ≈14)、OH は z/d 0 から (付着) | ref (**④ の結論: dual-time で定常化**) |
| `run_0079_react_li_dualtime_cont` | run_0077 `res_10000` から dual-time をさらに +100 ms (同設定, 物理 100–200 ms) — 出口流束の減衰を最後まで見て定常参照場を作る | **内部場は完全に定常** (z/d 26 の T は 100 ms 間で最終場と ≤1 K)。出口質量流束は減衰が止まり **0.0367 ± 0.002 kg/s (±5.5 %) の小振幅物理振動** (周期 ~20 ms, 平均は理論値 0.037): `check_quasisteady` exit_massflux **OSCILLATING → 平均±振幅で報告**、ignx/tmax/exit_yout_H2O STEADY。**dual-time 参照場 = `res_10000` (t=200 ms)** | ref (**Cabra 反応 ON の定常参照場**) |

## 機構着火遅れ比較 (2026-09-05, plan §5.1 P0-3)

`ign_delay_mech_compare.py` (Cantera, Cabra 混合線・断熱定圧 1 atm・τ = max dT/dt) による Jachimowski 9sp20r の 1000–1100 K 妥当性検証。成果物: `ign_delay_mech_compare.{csv,png,log}`。

| 機構 | τ_min @T_c 1015 K | @1045 K | @1075 K |
| --- | --- | --- | --- |
| Jachimowski 9sp20r (forge 現行) | 6.50 ms | **1.46 ms** | 0.73 ms |
| Li 2004 | 11.32 ms | 1.69 ms | 0.85 ms |
| Burke 2012 | 35.55 ms | 3.49 ms | 1.27 ms |
| GRI3.0 H₂ subset (既知不良) | 103.5 ms | 18.9 ms | 3.04 ms |

**2026-09-05 訂正**: 初版は coflow のモル分率が合計 1.1 (N₂ 0.8537, codex レビュー 2 で発覚) で Cantera が正規化していた (旧値は `ign_delay_mech_compare_WRONGCOFLOW.csv` に保存)。上表は N₂ 0.7537 (balance) の再計算値。

ξ_MR は 0.025–0.045 (文献 ≈0.05 と整合)。**Jachimowski は Cabra 検証実績のある Li 2004 より 14–43 % 早く、Burke 2012 の 1/2.4 (1045 K)〜1/5.5 (1015 K)** — 検証済み機構に対し系統的に早着火側 (0-D 判定は max dT/dt、CFD の火炎基部は Y_OH>2e-4、実験は OH 発光で、定義の対応は未検証)。GRI subset の 10 倍以上の遅れは Benim 2020 の「使用不可」報告と整合。Burke/Li の Cantera YAML は `solver_density_cuda/tools/mechanisms/h2_burke2012_cantera.yaml` / `h2co_li2004_cantera.yaml` で、forge は Li をそのまま読める (`run_0072`)。Cantera 側は各機構付属の熱力学、forge は CEA `species_db.yaml` を使うため逆反応・反応熱は同一ではない (ホスト 0-D 照合は plan §9 (7))。
