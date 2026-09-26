# case/40.nozzle_design_tool — ノズル設計ツール (forge_design) の検証ケース

`design/forge_design` (ノズル設計ツール — [親計画](../../plans/active/tooling-nozzle-design-tool.md) /
[Phase 0 子 plan](../../plans/active/tooling-nozzle-phase0-foundation.md)) の評価パイプライン検証用ケース。
問題定義 YAML (`problem_bell_smoke.yaml`) から「区分構成壁 → 構造化 quad msh4.1 →
`convertGmshToForge` → 品質ゲート → 準 1D 等エントロピー IC → forge (SLAU+SST+陰解法) →
収束判定 → 推力メトリクス」を 1 コマンドで実行する:

```bash
PYTHONPATH=design python3 -m forge_design.evaluate.runner \
    case/40.nozzle_design_tool/problem_bell_smoke.yaml \
    case/40.nozzle_design_tool/run_NNNN_<slug> [--steps N] [--prepare-only]
```

条件はいずれも ③ベル (Pt 4 MPa / Tt 1500 K / 背圧 20 kPa / rt 10 mm / ε 9 /
θa 20° / L/rt 7)、メッシュ 221×65 (14,080 cells)、品質 VERDICT PASS (AR max 19 / skew max 0.42)。

## 計算 run 一覧

> **連番衝突の注記 (2026-08-13)**: 並行セッションが同日に番号を採ったため **0071/0072 が二重発番**。
> 壁関数系列 = `run_0071_E3` / `run_0072_repdiag`、設計ツール系列 = `run_0071_node_e2e_regr` /
> `run_0072_bell_top_smoke`。スラッグで区別する (以後の新規 run は 0073 から採番)。

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_bell_smoke` | E2E 初回スモーク。SST の roK/roOmega を IC 未設定 (既定 0) | step 2–6 で rms_roK/roOmega に NaN → EOS 床洗浄後プラトー (既知指紋の再現)。パイプライン自体は完走し η=0.981 | 削除済み(2026-08-31) (NaN 指紋の記録) |
| `run_0002_bell_smoke_sstic` | IC に roK/roOmega (k=1, ω=18000) を追加、壁第一セル y+~70 (壁関数) | **NaN なし**。残差 1.1–1.7 桁低下後プラトー (`CONVERGENCE_VERDICT.txt`: NOT CONVERGED/stalled — cell モード atomicAdd ノイズ床水準)。**推力は 4000→12000 step で 0.002% ドリフト** (F=1961 N, η=0.9790, ṁ=1.299 kg/s, `metrics.json`)。`residual_history.png` あり | active (Phase 0 E2E 基準) |
| `run_0003_bell_L6` | dv 応答確認: L/rt=7→6 (`problem_bell_L6.yaml`) | η=0.9708, F=1944.7 N (短縮 → 発散損失増で η 低下 ✓)。品質 PASS・NaN なし・プラトー同傾向 | active (dv 感度の記録) |
| `run_0004_bell_L9` | dv 応答確認: L/rt=7→9 (`problem_bell_L9.yaml`) | η=0.9837, F=1970.5 N (延長 → η 単調増 ✓)。ṁ は 3 run で 0.03% 一定 (スロート同一の整合) | active (dv 感度の記録) |
| `run_0005_bell_smoke_node` | node 再計算の初回 (冷間 IC + implicit + 細壁 1e-3) | step 2 から全列 NaN。後の切り分けで真因 = 壁第一セル過細 + レシピ不一致と確定 → [調査 plan](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md) §2.5 | 削除済み(2026-08-31) (発散の記録) |
| `run_0006_bell_node_warm` | **node 初収束**: 実証レシピ (壁 5e-3・warm start from run_0002・explicit RK3 cfl0.1・1次・katoLaunder) | **VERDICT: PASS / ALL PASS** (全列 3.0–4.6 桁低下)。η=0.9809, F=1964.8 N (cell run_0002: η=0.9790 と 0.2% 差) | active (node 基準) |
| `run_0007_bell_node_2nd_imp` → `_2nd_expl` | 2 次化・陰解法化の試行 (run_0006 収束場から) | **いずれも step 4–10 で発散** — node×SST の陰解法と 2 次精度は未解決のソルバ課題 (case/29 の実績も 1 次 explicit のみ)。plan §2.5 参照 | 削除済み(2026-08-31) (発散の記録) |
| `run_0008_bell_node_runner` | runner 実装後の E2E 一気通貫確認 (`--warm-from`, 24000 step) | **VERDICT: PASS / ALL PASS** (3.4–5.1 桁・全列 falling)。η=0.9817, F=1966.4 N — node レシピが runner 既定で完結 | active (旧 node 1次 explicit 基準) |
| `run_0009_node_sst_imp_repro` | 課題1再現試行: 1次+陰解法 cfl_pseudo 2 (bndFirstOrder あり), 300 step | 「step 4 発散」は再現せず (300 step 健全)。旧発散記録は 2 次との交絡と判明 | 削除済み(2026-08-31) (切り分け記録) |
| `run_0010_node_2nd_repro` | 課題2再現試行: 2次+explicit cfl0.15 (bndFirstOrder あり), 300 step | 300 step 健全。反実仮想 (bndFirstOrder 除去, scratchpad) は step 8–10 で発散 = **課題2の因子は bndFirstOrder 欠落と確定** | 削除済み(2026-08-31) (切り分け記録) |
| `run_0011_node_2nd_imp` | 2次+陰解法 cfl2 の短尺確認 (300 step) | 300 step 健全 (全列 1.4–2.0 桁低下) — ただし長尺は下記 0012 で発散 | 削除済み(2026-08-31) (切り分け記録) |
| `run_0012_node_2nd_imp_full` | 2次+陰解法 cfl2 の 12000 step (軸修正**前**) | **step 10612 で発散** (roK がベル部近軸で e-fold ~2000 step の指数成長 → roOmega inf)。課題1の実態 = 遅発性 | 削除済み(2026-08-31) (発散の記録) |
| `run_0013_node_imp_long` | 1次+陰解法 cfl2 の 12000 step (軸修正**前**) | **step 7875 で発散** (同モード)。= 課題1は convMethod 非依存、軸行真空化×SST k シートが真因 ([plan §2.6](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | 削除済み(2026-08-31) (発散の記録) |
| `run_0014_node_2nd_expl_long` | 2次+explicit cfl0.15 の 12000 step (bndFirstOrder あり) | **完走・全列 2.5–2.8 桁 falling** = 課題2は bndFirstOrder で解決済みの確証 | active (課題2解消の証拠) |
| `run_0015_node_imp_axisdir` | **nodeAxisDirichlet=1** + 1次+陰解法 cfl2, 12000 step | 完走。rms_ro 1.4e-9 / roK 2.5e-7 (floor)、roK 成長モード消滅。軸床ピン 0 ノード。η=0.9816 (1次 explicit 0.9817 と一致) | active (軸修正の効き) |
| `run_0016_node_2nd_imp_axisdir` | **最終目標: 2次+陰解法 cfl2 + 軸 Dirichlet** (runner E2E, 12000 step) | **VERDICT: ALL PASS**・quasisteady ALL STEADY。η=0.9907, F=1984.4 N, ṁ=1.2923。**12000 step ≈ 20 秒** (explicit レシピ ~7 分の ~20 倍) | active (**node 生産基準**) |
| `run_0017_cell_kato` | cell 対照: run_0002 + katoLaunder (warm from run_0002) | η=0.9790 = run_0002 と同一 → **node−cell の η 差 +1.2% は katoLaunder では説明されない** (2次化の node/cell 離散差) | active (η 帰属の対照) |
| `run_0018_bell_L6_node` | L スイープ node 取直し: L/rt=6 (壁 5e-3 修正版 yaml, warm from run_0016) | 完走 (4.0–6.1 桁低下, roOmega は床で振動 2.9e-3↔5.6e-3)。**η=0.9835** | active (dv 感度 node) |
| `run_0019_bell_L9_node` | 同 L/rt=9。cross-geometry warm は cfl4 直投入で発散 → **段階起動** (soft explicit 3000 step → 陰解法 12000 step) | **VERDICT: PASS**。**η=0.9957** | active (dv 感度 node) |
| `run_0020_node_runner_cfl4` | 2次+陰解法 **cfl_pseudo 4** + bndFirstOrder の E2E (旧既定) | **VERDICT: PASS / ALL PASS** (3.2–4.6 桁)。η=0.9907 (cfl2 の run_0016 と一致) | active (bndFirstOrder あり世代の記録) |
| `run_0021_node_2nd_imp_nobfo` | **bndFirstOrder 撤去の主検証**: 全域 2次+陰解法 cfl4+軸 Dirichlet | 12000 step NaN なし。rms_ro ~1e-7 プラトー (入口/収縮部の微小リミットサイクル, ro 相対 ~1e-5)・SST 4.1–4.2 桁。quasisteady **ALL STEADY**。η=0.9896 | active (全域 2 次の主検証) |
| `run_0022_node_runner_nobfo` | **現 runner 既定** (全域 2次+陰解法 cfl4, bndFirstOrder なし) の E2E 基準 | run_0021 と同等 (η=0.9896, ALL STEADY, プラトー同性状) | active (**runner 既定基準**) |
| `run_0023_bell_L6_nobfo` | L/rt=6 を現既定で取り直し (warm from run_0021) | **VERDICT: PASS (converged)**。η=0.9822 | active (dv 感度 node 現行) |
| `run_0024_bell_L9_nobfo` | L/rt=9 を現既定で取り直し (段階起動: soft explicit 3000→陰解法 12000) | rms_ro 9.2e-9 (実質床, flat 判定)・ALL STEADY。η=0.9942 | active (dv 感度 node 現行) |
| `run_0025_node_wallsmooth` | 壁列エントロピー市松の種切り分け (warm start 壁 T 平滑化) | res_0 は市松 5K に清浄化されるが **12000 step 後に run_0022 とビット同一の 206K 市松が再生成** = 市松は定常解固有 ([plan §2.8](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | active (壁市松の証拠) |
| `run_0026_node_su2axis` | **SU2 流軸対称 (`axisymMethod: 1`)** の一次検証 (全域2次+陰解法 cfl4, nodeAxisDirichlet **なし**, warm from run_0021 収束場) | **軸行が真に健全** (床0・kシート消滅・軸線滑らか)。12000 step 完走・η=0.9904・ALL STEADY。ただし implicit は喉部近軸 limit cycle で rms_ro ~9e-5 プラトー (explicit は 3e-7 到達 = 空間健全) → [SU2化 plan](../../plans/active/axisymmetric-su2-source-formulation.md) | active (axisymMethod 1 基準) |
| `run_0029_node_adiabfix` | **断熱壁熱流束リーク修正の反実仮想** (現既定構成, 修正のみ, warm from run_0022) | **壁 T 市松 mean 206→43.7K (−79%)・壁温 1631→1195K**。η=0.9884 (−0.12%)。残差は `NOT CONVERGED` プラトーだが、4000/8000/12000 の η・ṁ・壁温は手計算で不変。以後の node 断熱壁 run はこの物理が既定 | active (**断熱修正基準**) |
| `run_0030_cell_5e3` | **cell を現行 5e-3 メッシュで取得** (壁温三者比較用, cold start, `problem_bell_smoke_cell.yaml`) | η=0.9813, ṁ=1.3035。壁温 mean 1185.8K = node (1195.9K) と同水準 → **壁温 −240K は node 固有でなく forge 共通** | active (cell 同一メッシュ基準) |
| `run_0028_su2_sst` | **SU2 v8.5 クロスチェック** (同一メッシュを gmsh で .su2 化, 同一 BC, RANS-SST/ROE/MUSCL/Euler implicit+FGMRES, 手順書テンプレ) | 15000 iter + 固定 CFL10 継続 5000 iter。全残差は未収束プラトー (`rms[RhoE]≈10^1.06`, `rms[ω]≈10^0.30`)。軸健全・壁 T 市松なし、η_CF=0.9573。入口コア乱流量は SU2 `k=0.975/ω=17160` vs forge `1/18000` でほぼ一致する一方、SU2 low-Re と forge automatic 壁処理、Roe/SLAU、SST 生産補正、軸離散が未整合の粗比較 | active (SU2 対照) |
| `run_0030_su2_sst_wallfunc` | run_0028 固定 CFL 最終場から継続し、SU2 `STANDARD_WALL_FUNCTION` だけを追加 | 第一内点が forge に接近 (`U=925 vs 843m/s`, `k=9.93e3 vs 8.54e3`, `μt/μ=8.3 vs 5.6`)。ṁ=1.2934もforge 1.2928に一致、診断 η=0.9725 (forge 0.9884)。それでも壁 T は面積重み **1422K**、forge は1193K。SU2 は Crocco–Busemann 断熱回復温度を明示設定するが forge SST 壁関数には熱閉包がないことを確証。`wall_temperature_compare.png`; 残差は `NOT CONVERGED` の診断run | active (熱壁関数の対照) |
| `run_0031_su2_sst_constk` | run_0028 固定 CFL 最終場から継続し、分子熱伝導率を forge と同じ 0.0257 固定へ変更 | SU2 壁 T は1426→**1441K** (+15K) で、forge の低温差を説明しない。forge の Sutherland μ + 固定 λ はベル第一内点で分子 Pr≈1.50–1.87となり air (`Pr≈0.72`) と非整合。残差は `NOT CONVERGED`; `residual_history.png` | active (熱物性の対照) |
| `run_0027_node_rfloor` | **r 床 (`axisRFloor: 3.0e-4`, ユーザ提案)** の検証 (method 0, nodeAxisDirichlet **なし**, 全域2次+陰解法 cfl4) | 軸健全 (床0)・η=0.9906・rms_ro ~1e-7 プラトー。床帯縁に有界の k 帯 (~230; scratchpad soak +12000 step でビット不変の平衡)。素朴な面床は閉性破れで発散 → 離散閉性面積 (A_closure) で解決 | active (axisRFloor 基準) |
| `run_0032_cell_yp1_lowre` | y+≈1 メッシュ調整 1 発目 (`wall_first_frac 2e-4`, 221×81, AR 191 PASS)。cell wf=0+constPr, warm from run_0030 | bell y+ mean 1.84 (狙い未達)・T_w 1378K・η 0.9780・ALL STEADY | ref (frac 感度点) |
| `run_0033_cell_yp1_lowre` | **y+≈1 確定メッシュ** (`frac 1e-4`, AR 381 PASS)。cell `wallTreatmentSST:0`+`thermCondMethod:1, prandtlLam:0.72` (Prt は既定 0.85), warm from run_0032 | bell y+ mean 0.94・T_w 1368.2K・η 0.9779・ṁ 1.2993・ALL STEADY (残差はプラトー) | active (Prt 感度対照) |
| `run_0034_node_yp1_lowre` | node 同条件 (wf=0, Prt 0.85)。**発散切り分けで 2 バグ修正** (双対 CV 重心/体積の float32 桁落ち→wall_dist ガベージ / interp_field の node 照会点) — [plan §2.9](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)。段階起動 (陰 cfl0.5 1次 3000→陰 cfl4 2次, stageA/ 同梱)。warm は node 場 (run_0029) 必須 | 全列 2.8–3.5 桁低下 (rms_ro 3.4e-9)・ALL STEADY。T_w 1365.9K = cell と 2.3K 差・τ_w ~1% 差・ṁ 1.2959 (5e-3 世代の wf=0 ṁ −4.5% 異常は消滅) | active (Prt 感度対照) |
| `run_0035_su2_yp1_lowre` | **SU2 v8.5 low-Re を同一 y+1 メッシュで** (run_0028 cfg 流用: adaptive 15000 + 固定 CFL10 5000 iter, PRANDTL_LAM 0.72 / PRANDTL_TURB 0.9 既定) | rms[Rho] 10^-5.95 (固定 CFL 相でも全列低下継続)。bell y+ mean 0.93・**T_w 1414.2K**・τ_w forge±1–3%・ṁ 1.3017・η 0.9796。**チェンバで T_w 1576K = Tt+76K の非物理超過あり** (断熱壁は Tt 超え不可; `wall_temperature_compare_yp1.png`) | active (**SU2 y+1 対照**) |
| `run_0036_cell_yp1_prt09` | **cell y+1 正式基準**: run_0033 + `turbulentPrandtl: 0.9` (SU2 既定に整合; Prt 0.85→0.9 で T_w +20.5K, scratchpad D5) | **T_w 1388.9K・η 0.9779・ṁ 1.2993**・ALL STEADY | active (**Step1 cell 基準**) |
| `run_0037_node_yp1_prt09` | **node y+1 正式基準**: run_0034 + `turbulentPrandtl: 0.9` | **T_w 1387.2K (cell と 1.7K 差)**・ṁ 1.2959・ALL STEADY。η は出口積分 0.9905 だが**出口列アーティファクトで過大** — 内部列積分 0.970–0.977 ([plan §2.10](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | active (**Step1 node 基準**) |
| `run_0038_node_5e3_tawwf` | **Step2 熱的閉包の主検証 (旧「弱閉包」設計, superseded)**: node 5e-3 wf=1 + `sstThermalWallFunction: 1` (+constPr/Prt0.9), warm from run_0029 | **訂正 (2026-08-11)**: 報告していた「bell 壁温 1417.9K = SU2 と 4.1K 一致」は**壁ノードの実状態ではなく出力配列 (bvar `Ts`) に上書きした $T_{aw}$ 診断値**だった — 実状態 $T[W]$ (bell 平均) は **1195.3 K で OFF 基準 1196K とほぼ不変** (旧弱閉包は場にほぼ効いていなかった)。η/ṁ の一致は有効 (保存量不介入のため)。旧設計は superseded、後継は [`turbulence-sst-adiabatic-taw-fluxmodel.md`](../../plans/accepted/turbulence-sst-adiabatic-taw-fluxmodel.md) (SU2 式流束モデル置換, 再検証待ち) | ref (**旧設計の記録・要再検証**) |
| `run_0039_cell_5e3_tawwf` | 同 cell 5e-3 (ghost 閉包→弱閉包に改訂後) | 壁温 1489.5K = **+70–90K 過大** (cell 代表点=第一セルの T が node/SU2 同高さ比 ~100–160K 高い — cell wf=1 BL 熱監査 follow-up)。初版 ghost Dirichlet は Tt 飽和で却下 | active (熱的閉包 cell の限界記録) |
| `run_0040_node_yp30_tawwf` | **Step3 壁関数系列 node (旧「弱閉包」設計, superseded)**: (`problem_bell_yp30.yaml` frac 1e-2, bell y+ mean 98, AR 3.8 PASS) wf=1+熱的閉包+constPr/Prt0.9, warm from run_0038 | **訂正 (2026-08-11)**: 報告値 1405.2K も run_0038 と同じ理由で出力 $T_{aw}$ 診断値であり実状態ではない (実状態は bell 平均 1134.2K、run_0044 の OFF 相当と近い)。η 出口積分 0.9835 も §2.11 の出口列欠陥修正前の値で二重に無効 (η 出口=内部列は run_0051 参照)。ṁ 1.2864 は有効。後継: `turbulence-sst-adiabatic-taw-fluxmodel.md` | ref (**旧設計の記録・要再検証**) |
| `run_0041_cell_yp30_tawwf` | 同 cell (`problem_bell_yp30_cell.yaml`, bell y+ mean 71), warm from run_0039 | 壁温 1387.2K (±20K 市松)・η 0.9802・ṁ 1.3029・ALL STEADY。**cell 代表点バイアスは y+~70 では消滅** (5e-3 のバッファ層代表点固有) | active (Step3 cell) |
| `run_0042_su2_yp30_wallfunc` | **SU2 STANDARD_WALL_FUNCTION を同一 y+30 メッシュで** (low-Re 15000 → wf 継続 5000+10000 iter) | 壁温 **1418.9K**・η 0.9673・ṁ 1.2868。+10000 iter で全量不変 (準定常; rms[Rho] −4.5 プラトー, `residual_history.png`)。`wall_temperature_compare_yp30.png` | active (**Step3 SU2 対照**) |
| `run_0043_node_yp1_outletfix` | **node 出口列欠陥の根治検証** (run_0037 と同一入力の A/B): `outlet_statPress` の bvar `Psb` 動的化修正 ([plan §2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | 出口列 sag **19.4%→0.39%** (超音速行 0.40%)。**η 出口積分 0.9754 = 内部列 0.9747–0.9755** (アーティファクト解消)。ṁ 1.2959 列間ばらつき消滅・ALL STEADY。壁温は run_0037 と ≤0.85K 差 (出口隣接除く) | active (**node y+1 正 (η 込み)**) |
| `run_0044_node_yp30_outletfix` | 同修正の y+30 wf 系 A/B (run_0040 と同一入力) | sag 13.6%→0.17%。**η 出口積分 0.9835 (偽) → 0.9687 = 内部列 0.9679–0.9687**、SU2 wf (run_0042: 0.9673) と +0.14%。全列 2.8–5.1 桁低下・ALL STEADY。ṁ 1.2864 不変 | active (**Step3 node 正**) |
| `run_0045_node_yp1_outletfix_cont` | run_0043 の res_12000 から +12000 step 継続 (準定常確認) | 場は run_0043 と不変 (出口 P 最大 6e-4)。残差は低レベルリミットサイクル (rms_ro 1.1e-8, rms_roUx 3.5e-5) で ALL STEADY | active (継続確認) |
| `run_0046_node_yp30_offregr` | `sstEnergyWallFunction` 実装後の OFF 回帰 (run_0044 と同一入力・新バイナリ) | 場の差 rel ~1e-5 = run-to-run atomicAdd 床と同水準 (run_0047 で確認) → **OFF ビット等価** | 削除済み(2026-08-31) (回帰確認のみ) |
| `run_0047_node_yp30_repro` | 同一バイナリ再走で run-to-run 床を測定 | rel 2e-6〜3e-5 (roK 最大) — run_0046 の差と同桁 | 削除済み(2026-08-31) (床測定のみ) |
| `run_0048_node_yp1_isoT_ref` | **等温壁 Ts=1000K の y+1 low-Re 真値** (run_0043 入力 + wall_isothermal, warm from run_0043) | q_w: スロート −1.7 MW/m²・ベル −0.15〜−0.6 MW/m²。残差 2 桁低下プラトー・q_w 静定 | active (**等温壁真値**) |
| `run_0049_node_yp30_isoT_qwwf` | **`sstEnergyWallFunction:1` (Kader 原式) node y+30 wf** (run_0044 入力 + wall_isothermal Ts=1000K, warm from run_0044) | vs run_0048: **チャンバ/スロート −9%** (実用域)、**ベル部 (M≈4 冷却壁) +37〜+265% (積分 +87%) 過大 = 非圧縮 Kader T⁺ の圧縮性限界を実測** ([plan §8](../../plans/active/turbulence-sst-thermal-flux-model.md) フォローアップ) | active (**圧縮性限界の記録**) |
| `run_0050_node_yp1_outletic` | **出口修正の最終形 (ic 参照) A/B** (run_0037 と同一入力): Psb 入力専用化 + node outlet の p_tilde/勾配閉包を内部値参照に改訂 ([plan §2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md) 改訂) | sag 0.40% / η 出口 0.9754 = 内部列 0.9747–0.9756・ṁ 1.2959 — **前実装 (run_0043) と同一の根治**。Psb 分布指定 (inletProfile CSV) が保持される | active (**node y+1 正 (最終形)**) |
| `run_0051_node_yp30_outletic` | 同最終形の y+30 wf A/B (run_0040 と同一入力) | sag 0.18% / η 出口 0.9687 = 内部列・SU2 wf +0.14%・ṁ 1.2864 — run_0044 と同一 | active (**Step3 node 正 (最終形)**) |
| `run_0052_node_yp1_dofonly_regr` | **境界勾配 DOF-only 化の一般化検証** ([plan](../../plans/accepted/architecture-node-boundary-gradient-dof-only.md)): run_0050 と同一入力を全 primitive owner-state 化後の binary で再実行 | sag 0.40% / η 0.9754 / ṁ 1.2959 — **run_0050 と完全一致** (outlet 特例の一般化は無害と確認) | active (**DOF-only 化 正 (y+1)**) |
| `run_0053_node_yp30_dofonly_regr` | 同一般化 + `sstThermalWallFunction:1` (新 Taw 流束モデル置換, run_0051 と同一入力) | **step 1686 で roOmega NaN 発散** (壁非接触の第一内部行 11 点)。診断用、削除済み(2026-08-31) | 削除済み(2026-08-31) (発散記録) |
| `run_0054_cell_bitinv_regr` | **cell 不変性検証**: run_0036 と同一入力を新 binary で再実行 | run_0036 比相対差 0.6–2.3% だが同一 binary 再実行 (run_0056) でも 0.7–3.7% (=[[cell-atomicadd-nondeterminism]] のカオス的ばらつき範囲) — cell 経路は本セッションの変更で一切呼ばれないことをコード上も確認済み | active (cell 不変性の根拠) |
| `run_0055_node_yp30_dofonly_tawoff` | 診断: run_0053 と同一入力で `sstThermalWallFunction:0` に戻し一般化 GG/LSQ 単独の健全性を切り分け | **ALL PASS** (全列 3.5–5.1 桁低下) — 発散は Taw 流束置換機構に起因すると確定 | active (切り分け証拠) |
| `run_0056_cell_bitinv_repro` | run_0054 の同一 binary 再実行 (run-to-run 床測定) | 相対差 0.7–3.7% — run_0036 対比とほぼ同水準、cell 不変性の根拠 | 削除済み(2026-08-31) (床測定のみ) |
| `run_0057_node_yp30_taw_repedge` | 診断: Taw 流束置換を代表点 (Normal_Neighbor) エッジのみに制限する修正版で再試行 | **step 1480 で同型発散** — 後にユーザレビューで root cause 確定: 置換流束が T[W] 非依存になり壁半 CV の復元項を喪失、補償なき外部吸熱として壁ノードが EOS 床まで異常冷却 (step1000 で T[W] min 1100.7→27.1 K [全辺 run_0053] / 15.5 K [本 run])。roOmega NaN は下流症状 ([plan §4.3 訂正・§9](../../plans/accepted/turbulence-sst-adiabatic-taw-fluxmodel.md)) | 削除済み(2026-08-31) (発散記録) |
| `run_0058_node_yp30_taw_outputonly` | **最終形の検証**: Taw 注入を全面撤回し出力専用化 (`sstThermalWallFunction: 1`, run_0055 と同一入力) | **ALL PASS** (3.5–5.1 桁)。場は OFF (run_0055) と atomicAdd 床内で一致 (rel ≤5e-5) = W–I 熱流束は Taw 不変。T[W] min 1100.6 K (健全)。差分は壁面出力 Tsb のみ: 1412.7 K (Taw モデル値, SU2 wf 1418.9 K と 6 K 差) | active (**Taw 出力専用の生産 baseline 正**) |
| `run_0059_node_yp30_taw_su2coupled_lowcfl` | **mode 2 (SU2 coupled, experimental) 短尺試験**: warm restart (run_0055 res_12000) + cfl_pseudo 0.5, 200 step | NaN なし完走だが T[W] −0.7 K/step 単調冷却。壁 k_wf ピンは機能 (k[6564]=9142) も μt[W] 不変 9.5e-4 | 削除済み(2026-08-31) ([plan §10b](../../plans/accepted/turbulence-sst-su2-taw-coupling.md)) |
| `run_0061_node_yp30_mode0_regr` | mode 0/1 経路の回帰 (新バイナリ, run_0055 と同一 IC 1000 step) | run_0055 res_1000 と node 再現性帯 (~2e-6 rel) で一致 = 既定経路不変 | 削除済み(2026-08-31) |
| `run_0062_node_yp30_taw_su2c_lowcfl_long` | mode 2 の平衡確認 (run_0059 と同構成 5000 step) | T[W] min 1125→113.6 K 単調降下 (旧 run_0053 と同軌道)。T_I は加熱されず BL 加熱で釣り合う経路なし → mode 2 現形は不採用 | 削除済み(2026-08-31) |
| `run_0063_node_yp30_taw_su2c_mutw02` | + 実験 env `FORGE_TAW_WALL_MUT_SCALE=0.2` (壁側 μt を SU2 実効値相当へ) | 冷却速度半減も単調降下継続 (漸近 ~500 K 見込み=非物理)。壁 μt 縮小だけでは不十分と動的確証 | 削除済み(2026-08-31) |
| `run_0064_node_yp30_mode0_regr3` | mode 3 実装後の mode 0 回帰 (run_0055 同一 IC 1000 step) | run_0055 res_1000 と再現性帯 (~2e-6 rel) で一致 = 既定経路不変 | 削除済み(2026-08-31) |
| `run_0065_node_yp30_taw_defectflux_lowcfl` | **mode 3 (defect-flux 閉包) 短尺**: warm (run_0055) + cfl_pseudo 0.5, 200 step | T[W] が単調**上昇** 1125→1246 K (幾何減衰, 外挿≈1400 K)、内点 986 K 不動、NaN なし = 負帰還設計どおり | 削除済み(2026-08-31) |
| `run_0066_node_yp30_taw_defectflux_prod` | **mode 3 生産条件** (cfl_pseudo 4.0, 12000 step, warm) | T[W] が step1000 で **1418.2 K にプラトー** (Taw_diag 1416.5+ε, 12000 まで不変)。SU2 primitive 壁温 1428 K と **10 K 一致 (実状態同士)**。非壁場は run_0055 と rel mean ~1e-4。`check_quasisteady` **ALL STEADY**、出口 massflux/推力積分は baseline 比 ~1e-5 で不変。残差床 6–12 倍は出口リップ 16 ノードの既存 ω サイクル増幅 (完全局在, 非ブロッカー) — [plan §10c](../../plans/accepted/turbulence-sst-su2-taw-coupling.md) | active (**mode 3 正式採用の生産基準** — 以後の断熱壁 node run は `sstThermalWallFunction: 3` を既定とする) |
| `run_0067_node_yp30_taw_defectflux_kwfwall` | mode 3 + 壁ノード roK_wf Dirichlet 拡張の A/B (リップ ω 振動の壁 k 浮上仮説の検証) | 壁温・場は run_0066 と同一。rms_roOmega 0.141→0.165 で**改善せず → 仮説棄却**、mode 3 は最小構成 (壁 k ピンなし) に確定 | 削除済み(2026-08-31) (A/B 記録) |
| `run_0068_node_yp30_taw_defectflux_omgwf` | mode 3 + `nodeOmegaWfDirichlet: 1` (irep ω ピン, SU2 SetTurbVars 両変数化) | リップ ω サイクル完全消滅 (rms_roOmega 9.4e-4) だが **τ_w=ρ_rep·u_τ² が y+1 真値の 1.237 倍に過大** (ピンなし 0.945 / SU2 0.985)。第一内点 μt が SU2 の 3.3 倍 (リミッタ不発) → BL 過混合 → 摩擦過大 = 推力 −0.18% の正体。**正式化見送り** — [plan §10c](../../plans/accepted/turbulence-sst-su2-taw-coupling.md) | active (ω ピン棄却の根拠) |
| `run_0069_node_yp30_mode0_omgwf` | mode 0 + ω ピン (推力シフトの帰属分離) | 推力 309.623 N ≈ run_0068 → **−0.18% は ω ピン閉包そのものの効果、mode 3 と独立**。ただし冷壁との組合せは残差悪化 (rms_roOmega 0.23) | 削除済み(2026-08-31) (帰属記録) |
| `run_0089_closure_spalding` | §5.2 ③ **整合実験**: `FORGE_WF_INV_LAW=spalding` + **`FORGE_WF_CLOSURE_LAW=spalding`** ($g$→`wf_pk`,`roK_wf`)。`run_0088` 場から 12000 step | $\tau_w$=**1457.93 Pa** = y+1 内部基準比 **1.0027** (② の 1.0063 より 1.0 に近い)、壁温 1411.39 K、$q_w$=0 (断熱)。`wall_model_tau` **STEADY**。②/③ 比: `wf_pk` **+7.3%** / `roK_wf` **−7.8%** ($y^+$~92 は $g$ 反転域。交差点 $y^+$≈41.8) — [plan §8](../../plans/active/turbulence-reichardt-gap-residual.md) | active (§5.2 ③ の正) |
| `run_0087_invlaw_reichardt` | §5.2 連成 A/B の**対照**: `run_0066` 場から 12000 step、env OFF | $\tau_w$=**1388.33 Pa** = y+1 内部基準比 **0.9548**、壁温 1410.82 K。`wall_model_tau` **STEADY** (正式評価関数へ整合済) | active (A/B 対照) |
| `run_0088_invlaw_spalding` | §5.2 連成 A/B: 同 restart + `FORGE_WF_INV_LAW=spalding` (逆解き則のみ、$g$ は Reichardt 固定、E3 OFF) | $\tau_w$=**1463.11 Pa** = y+1 内部基準比 **1.0063** (+0.6%)、壁温 1411.75 K (+0.9 K)、$q_w$=0 (断熱)。`wall_model_tau` **STEADY** → **$\omega$ pin 1.2352・E3 1.1414 と違い過大化が +0.6% に収まる** — [plan §7](../../plans/active/turbulence-reichardt-gap-residual.md) | active (§5.2 の正) |
| `run_0072_repdiag` | **代表点幾何診断** (`FORGE_WF_REP_DIAG=1`、`run_0066` 場から 100 step) | `best_cos`=0.956・**接線オフセット/y=0.306**・`wall_dist[irep]/y`=1.046。**要素ベース**凍結場カウンターファクトで法線化すると $\tau_w$ が **0.610 倍** (旧 IDW の 0.679 は近傍数依存の産物で撤回)。y+1 内部基準比 0.9549→**0.5825** → 候補 A は単独では採らない — [代表点 plan](../../plans/accepted/turbulence-node-wf-representative-point.md) | active (§3 の正) |
| `run_0071_E3` | **E3 A/B ($\omega$ 生産を壁法則整合 $S_{\rm wf}^2$ へ)**: `FORGE_WF_OMEGA_SOURCE=1`、`run_0066` 場から 12000 step、他は生産構成と同一 | NaN なし・`check_quasisteady --quantity wall_tau` **STEADY**。**y+1 内部基準比 = 1.1414 に過大化** (単一評価関数 `wall_tau_eval.py`。従来の「≈1.061」は $\tau_w$ 比を 0.945 に掛けた推定で**過小だった**) → **採用不可** (ただし $\omega$ ピンの 1.237 より遥かに軽い)。壁温 1414.0→1421.1 K。case/26 側は 0.7615→0.8445 で改善 — [plan §5](../../plans/accepted/turbulence-node-wf-omega-source.md) | active (**E3 不採用の根拠**) |
| `run_0070_node_yp30_mode3_lsq` | **LSQ 勾配の本番検証** (`mesh.gradLSQ: 2`, 生産 mode 3 構成の A/B, 対照 run_0066) | 完走・NaN なし。**解は実質同一** (τ_w/真値 0.945 = GG と同値・壁温 1450.9K・ṁ/推力とも rel ~1e-5)。ただし**残差床は GG の 1.5-2 倍高い** → 構造化・低歪みメッシュでは LSQ の利得なし、既定 `gradLSQ: 0` 据置が妥当 — [LSQ plan](../../plans/active/discretization-lsq-gradient.md) 検証3 | active (LSQ A/B 記録) |
| `run_0071_node_e2e_regr` | **設計ツール E2E 回帰 (2026-08-13, 現バイナリ)**: runner 既定 (5e-3 smoke yaml, warm from run_0029) — [Phase 0 plan](../../plans/active/tooling-nozzle-phase0-foundation.md) クローズ用 | 完走。warm start から全列さらに 2.6–4.5 桁低下後プラトー (rms_ro 2.4e-9 = ノイズ床)・quasisteady **ALL STEADY**。**η=0.9734・ṁ=1.2928** (=run_0029)。旧 runner 基準 η 0.9896 (run_0022) との差は node 出口列欠陥の根治 (run_0043 系) による正常化 | active (**runner E2E 現行基準**) |
| `run_0072_bell_top_smoke` | **TOP 放物線壁の E2E スモーク** (`problem_bell_top.yaml`: bell_type top, θn20°/θe8°/L7 — [Phase 2 plan](../../plans/active/tooling-nozzle-moo-loop.md) ステップ 2, warm from run_0071) | 品質 PASS (AR 7.7/skew 0.42)・**VERDICT: PASS/ALL PASS (全列 3.4–6.0 桁)**・**ALL STEADY**。**η=0.9731・ṁ=1.2928** (同条件 Hermite run_0071 と −0.03%/一致 — 族の差は小、dv 感度は Phase 2 DOE で取る) | active (**TOP 壁基準**) |
| `run_0073_moo_smoke` | **MOO ドライバのスモーク** (DOE4+infill1, 各評価 = soft 段 3000→本段 12000, [Phase 2 plan](../../plans/active/tooling-nozzle-moo-loop.md) ステップ 4)。旧版で 2 バグ検出 (不始動流の偽アトラクタ / stage 最終場未ダンプ) → 修正 | 修正後 **5/5 PASS** (~33s/評価, 全評価 conv PASS・ALL STEADY・ṁ 比 0.9858±0.0003)。η=0.952 (L5.1)〜0.978 (L10.9)。`ledger.jsonl`/`pareto.json` | 削除済み(2026-08-31) (ドライバ検証の記録) |
| `run_0074_moo_top_r1` | **③ベル TOP dv の本番 MOO キャンペーン** (DOE24 + EHVI infill 2×13 = 50 評価, seed=run_0072, ref=(−0.90,13), 評価は eval タグ別サブディレクトリ) | **49/50 PASS・28 分** (conv PASS 37 / STALLED 12・全 PASS 点 ALL STEADY)。唯一の FAIL = `inf_11_0` (θn30/θe4/L5 角) は **ṁ/ṁ_1D=1.034 の物理ゲートで排除** (残差・準定常は健全な非物理解)。**パレート 19 点: η 0.953@L5 → 0.980@L10、前線ほぼ全点で θe→下限 4°**。基準 run_0072 は同 L で前線比 −0.2%。`ledger.jsonl` / `pareto.json` / `pareto_r1.png`。**前線ポリッシュ (2026-08-13)**: `opt/polish.py` で前線 19 点を +12000 step チャンク継続 → **19/19 settled (|Δη|<3e-4)・前線不変** (`<tag>_p1`, `polish.jsonl`, `pareto_polished.json`, `pareto_final_rao.png`) | active (**Phase 2 本番 r1・確定前線**) |
| `run_0075_rao_check` | **Rao 照合** ([Phase 2 plan](../../plans/accepted/tooling-nozzle-moo-loop.md) §6, 同一評価系+ポリッシュ, `problem_bell_rao.yaml`): rao80 = case/29 digitize チャートの ε=9 点 (θn33.8°/θe12.1°/L5.97)、rao100 = θn流用の意図的非最適対照 (θe8°/L7.46) | **合格**: rao80 η=0.96608 (settled) = 前線補間 0.96661 と **+0.053% 差で前線上**。rao100 η=0.97312 は前線が +0.28% 支配 (誤 θn は正しく支配される)。図 = `../run_0074_moo_top_r1/pareto_final_rao.png` | active (**Rao 照合の判定**) |
| `run_0076_fixedL597` | **固定 L=5.9713 の (θn,θe) 直接最適化** (Rao 判定の精密化 — 前線ガタつき/Rao超え疑いへの回答。サロゲート支援 exploit+explore 各3ラウンド 6 評価, warm from rao80) | 最適は **θn24.5–28.5×θe8–11 の平坦台地** (max η=0.96741 @θn25.8/θe10, settled)。**rao80 はその +0.133% 下 = 判定帯 ±0.15% 内 → 「Rao は前線上」確定** (Rao 超えではない)。θe 下限張り付きも固定 L では最適 θe=10° = θe(L) 減少は Rao 傾向の再現。infill +20 (計70評価) との確定前線 23 点 = `../run_0074_moo_top_r1/pareto_merged_final.{json,png}` | active (**Rao 判定精密化**) |
| (追加可視化) | **形状+マッハコンタ比較図**: Rao点・同L最適・前線各L(6.49/7.56/8.64/10.05) の壁重ね書き + 6パネル同一カラースケール M コンタ (`_p1` ポリッシュ済み場, griddata 構造格子補間で三角形分割縁アーティファクト回避) | `run_0074_moo_top_r1/shape_compare.png` / `mach_contour_compare.png`。全ケース M=1 白線がスロートから立ち上がり壁まで滑らかに発達、異常構造なし | active (**Rao vs 前線 可視化**) |
| `run_0077_rao80_hh` | rao80 を **Rd=0.382 (Huzel–Huang 正規値)** で再計算 (0.4 との比較) | 内部衝撃波・軸 M 過大とも**同一に再現** (axisM_exit 4.73, ΔM −0.49)。η=0.9661 も不変 → 円弧半径の 0.4/0.382 差は無関係 | active ([調査ノート](../../notes/investigations/nozzle-top-internal-shock-diagnosis.md)) |
| `run_0078_su2_rao80` | **SU2 v8.5 RANS-SST を rao80 と同一メッシュ (msh→su2 直変換)・同一 BC で** — 内部衝撃波が forge 固有かの判定 | **forge と定量一致**: axisM_exit 4.63 (forge 4.70, +1.5%)・衝撃波トレース同一軌跡 (差 0.3–0.7mm)・強度成長同傾向。15000 iter, rms[Rho]≈−5.3 プラトー (既知性状)。**成果図 `cross_compare_2x2.png` / `cross_compare_profiles.png`** | active (**forge 無罪の確証**) |
| `run_0079`〜`run_0085` (forge 切り分け 7 本) | **夜間バッチ**: cell A/B (0079) / axisymMethod=1 (0080) / Rd=0.8/1.5 (0081/0082) / Ru=0.75/3.0 (0083/0084) / **厳密 Euler cell+slip** (0085)。各 run に `diag.{json,txt}` | cell 4.72・Euler 4.72 = **離散化/粘性/壁関数とも無関係**。Rd↑で過膨張緩和 (4.70→4.35) = 円弧曲率起因と整合、Ru は無関係 (4.70/4.69)。0080 は陰解法深収束未達 (既知) で証拠除外 | active (切り分けマトリクス) |
| `run_0086_su2_euler` | **SU2 Euler** 同一メッシュ (非粘性 2×2 の完成) | **8.4k iter で rms[Rho] 10⁻¹¹ 深収束**。axisM_exit 4.62・ΔM −0.42 = forge Euler (4.72/−0.42) と一致。**結論: 内部衝撃波 (円弧→放物線接合起因) と軸 M 過大 (kernel 過膨張) は物理。forge 無罪** | active (**非粘性対照**) |
| `run_0950_sglsq_s2_gg` / `run_0951_sglsq_s2_lsq` → 延長 `run_0952_sglsq_s2_gg_ext` / `run_0953_sglsq_s2_lsq_ext` (2026-09-26、**AWS** `~/forge-pgrad-new/case/40.nozzle_design_tool/` のみ) | **S2 物理 A/B (スカラー勾配 gg / lsq)** [plan](../../plans/accepted/gradient-scalar-lsq-unification.md) §6 S2。起点 main ワークツリー `run_0045_node_yp1_outletfix_cont/res_12000` (sha256 `15e47af0…`、nozzle.h5 `849660af…`) から `restart_field.py` で 12000 + 延長 12000 step。バイナリ `f99f236d` (AWS `~/sglsq/forge_f99f236d`、sha256 `98703ca1…`)、`FORGE_CUDA_BLOCKSIZE=128`/`_SMALL=128`、双子の差は `mesh.scalarGradient` の 1 行のみ (`s2_setup.py`、`IC_FROM.txt` に起点 sha256 と設定差)。**起点からの設定差** (両双子に共通): 廃止キー `mesh.nodeAxisDirichlet: 1` を削除 (2026-08-16 廃止、書くと起動エラー)、`slauWallNormalChi` 実効 1 | **物理ゲート (延長末尾、`s2_eval_case40.py` = `thrust_metrics`、起点 run_0045 で η 0.97542・ṁ 1.29594 を再現)**: η_CF 0.97794/0.97795 (∈ 0.978 ± 0.003)、ṁ 1.29947 (∈ 1.29–1.30)、η・ṁ とも STEADY。**lsq−gg**: Δη_CF +0.000 %、Δṁ 0.000 %、壁温 L∞ 0.62 K (上限 0.1 %/0.2 %/15 K)。起点 0.9754 → 0.9779 の移動は両双子共通 (設定差)。**収束**: 両方 plateau (RISING なし)。`check_floor_ratio --start run_0045` は両方 NOT BACK ON FLOOR (rms_roOmega 56×)、現行 gg 初段の床基準では両方 PASS。解釈は F 判断待ち。**現行判定 (2026-09-26 ユーザ決定で基準を現行 gg 双子に差し替え)**: lsq/gg 末尾残差比 ≤ 1.5 で PASS (`case/09.Taylor-Green/_g0_lsq_seam/S2_GGREF.txt`・`s2_evidence/JUDGEMENTS.txt`)。以下の起点床の判定は旧基準の記録 | active (S2) |

**壁半 CV 凍結場収支診断スクリプト** (plan [turbulence-sst-su2-taw-coupling §10](../../plans/accepted/turbulence-sst-su2-taw-coupling.md) の §11 診断の実体):
`wall_budget_frozen.py` (forge 収束場の壁ノード収支を A/B/C/D+E 方式で再構成)、
`su2_wall_budget.py` (SU2 run_0042 場での同再構成)。パスはファイル内定数 (run_0055/run_0042 前提)。

**軸処理 3 方式の比較 (全域 2 次+陰解法 cfl4, bndFirstOrder なし)**:

| 方式 | 軸ノード | row1 k 帯 | 収束 (rms_ro) | η_CF |
| --- | --- | --- | --- | --- |
| `nodeAxisDirichlet: 1` (現 runner 既定) | 解かずミラー | 3.4 | ~1e-7 プラトー (run_0022; 1 次 or bndFirstOrder 世代なら 1e-9 PASS = run_0015/0016) | 0.9896 |
| `axisymMethod: 1` (SU2 流) | 真に解く | **1.35** | implicit ~1e-4 プラトー (explicit 3e-7) | 0.9904 |
| `axisRFloor: 3e-4` (r 床+閉性補正) | 解く (環状近似) | 230 (有界, 帯定義: y∈[0.5,1.5]mm × x∈[30,60]mm の max k) | ~1e-7 プラトー | 0.9906 |

※ 収束値は check_convergence の実測 (全域 2 次では 3 方式とも NOT CONVERGED/プラトー判定であり、
計量は quasisteady STEADY で担保)。η は同一 warm 系統・同一メッシュだが cell (run_0002) は壁 1e-3
世代メッシュのため単純比較不可。**node−cell 差 +1.1% は出口列アーティファクトと判明 (2026-08-11,
[plan §2.10](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md))**。
2026-08-05 の断熱壁修正 (run_0029) 以後は node の η が ~0.12% 下がる (漏れ加熱の除去)。


**dv 感度まとめ (η–L トレードオフの応答確認)**:
L/rt = 6 / 7 / 9 → η_CF は cell (2026-08-03 前半, 壁 1e-3 メッシュ) 0.9708 / 0.9790 / 0.9837、
**node 全域 2次+陰解法 (現既定, run_0023 / 0021 / 0024, 壁 5e-3)** 0.9822 / 0.9896 / 0.9942
(bndFirstOrder あり世代 run_0018/0016/0019 は 0.9835 / 0.9907 / 0.9957)。いずれも物理的に
正しい単調応答で、node−cell はほぼ一定の +1.1〜1.2% オフセット (トレンド保存)。オフセットの
帰属: katoLaunder ではない (run_0017)・軸修正でもない (修正前 run_0012 も 0.991)。
**2026-08-11 訂正: この +1.1〜1.3% は node 出口列の系統的不整合による出口積分アーティファクト**
([plan §2.10](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) — 内部列で再積分すると
node η は cell と整合する。node の η 出口積分値は修正まで生産値に使わない。

**node モードの状態 (2026-08-04 更新)**: 残っていた 2 課題は解決した
([plan §2.6–2.7](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md))。
2 次発散・陰解法の遅発性発散は**いずれも node 軸対称の軸行真空化が種** (ベル部で軸 CV が
radial 圧力平衡からデカップルし EOS 床まで過膨張 → 偽せん断 → SST k シート) — ソルバ新機能
**`nodeAxisDirichlet: 1`** (軸ノードを radial 隣接からの対称 Dirichlet に置換) で根治。
当初レシピに入れていた `bndFirstOrder: 1` (境界隣接 1 次化) は軸病変のマスキングであり、
壁境界層の精度を損なうため**撤去** (ユーザ指示 2026-08-04)。**現 runner 既定 =
全域 2 次 (convMethod 1) + 陰解法 (blockDPLUR, cfl_pseudo 4) + nodeAxisDirichlet**。
1 評価 12000 step ≈ 20 秒 (旧 explicit レシピ比 ~20 倍)。残差は形状により入口/収縮部の
微小リミットサイクルで ~1e-7 プラトーになるため (cell 全域 2 次と同格)、**計量は
`check_quasisteady.py` STEADY を確認して使う**。warm start は引き続き必須で、
**同メッシュ node 場からが最良**。異形状/cell 場からの warm は cfl4 直投入で初手発散する
ことがあり、その場合は段階起動 (soft explicit 3000 step → 陰解法、run_0019/0024 の手順)。
run_0001–0004 は壁 1e-3/2e-3 世代のメッシュ (現 problem yaml は全て 5e-3 に統一済み)。

注記: 残差プラトーの深掘り (真の定常収束品質) は Phase 2 (Rao 照合) で扱う。η の妥当性
(0.98 前後はベルノズルとして自然なオーダー) も Phase 2 の照合が正式判定。



**壁温・η_CF の生産値確定 (2026-08-11, 壁温真値 3 段階キャンペーン完了。※壁温の wf+閉包系列は
2026-08-11 (2) の訂正で要再検証、下記注意④参照)**:
y+≈1 low-Re 三者基準 (run_0033–0037)・熱的壁関数 `sstThermalWallFunction` (run_0038–0039,
旧「弱閉包」設計 [archived plan](../../plans/archived/turbulence-sst-thermal-wall-function.md),
後継 [`turbulence-sst-adiabatic-taw-fluxmodel.md`](../../plans/accepted/turbulence-sst-adiabatic-taw-fluxmodel.md))・
y+≈65–300 壁関数系列 (run_0040–0042) による挟み込みの結論:

| 量 (③ベル Pt4MPa/Tt1500K/ε9) | 生産値 | 根拠 |
| --- | --- | --- |
| **ベル部壁温 T_w (断熱)** | **1400 ± 15 K** (理論 T_aw ≈1400、y+1 low-Re と SU2 wf のみで支持) | y+1 low-Re (無汚染): forge 1387–1389 / SU2 1414。SU2 wf (無汚染, forge 側バグと無関係): 1419–1422。**forge node wf+閉包 (旧弱閉包設計) の 1405–1418K は出力 T_aw 診断値であり実状態ではないため生産値の根拠から除外** (④参照)。旧 wf 値 1193K は熱閉包欠落、旧 cell 1767K は熱物性誤り |
| **η_CF** | **0.978 ± 0.003** | 最良解像 (y+1) で cell 0.9779 / SU2 0.9796 / node 0.9754 (出口欠陥修正後 run_0043, 出口積分=内部列)。壁関数系列は 0.967 (SU2)〜0.981 (cell) に散り、メッシュ/壁処理依存が残る |
| ṁ | 1.29–1.30 kg/s (±0.5%) | 全系列。y+30 wf 系 (SU2/node) は −1% 側 |

**注意**: ① node 出口列欠陥は **2026-08-11 に根治済み**
([plan §2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md), run_0043–0045) —
node の η_CF は出口積分値をそのまま使ってよい (修正後の正値: y+1 low-Re 0.975 / y+30 wf 0.969)。
修正前 run (run_0042 以前) の node η 出口積分値 (0.9835–0.9905 系) は +1.3% 過大のため使わない。
② 壁温は `wallTreatmentSST: 1` +
`sstThermalWallFunction: 1` + `thermCondMethod: 1, prandtlLam: 0.72` +
`turbulentPrandtl: 0.9` を生産構成とする (node 一次)。
④ **`sstThermalWallFunction` の位置付け確定 (2026-08-11 (2)〜(3) の経緯)**: 旧「弱閉包」設計
(run_0038/0040) で報告していた壁温は壁ノードの実状態ではなく出力配列の $T_{aw}$ 診断値であり、
実状態は OFF 基準とほぼ不変だった (`run_0038` bell 平均 1195.3K vs OFF 1196K)。その後試した
W–I 流束注入は発散で棄却され (run_0053/0057)、**最終形 = $T_{aw}$ は「壁関数が再構成した物理壁面
温度のモデル値」として境界出力 (Tsb) 専用** (run_0058 で検証済,
[plan §0](../../plans/accepted/turbulence-sst-adiabatic-taw-fluxmodel.md))。**壁温の報告は
「モデル壁面温度 Tsb」(forge 1412.7K vs SU2 wf 1418.9K = モデル値同士で 6K 差) と「壁ノード実状態
T[W]」(~1270K, 壁半 CV 平均温度) を明確に区別して行うこと** — SU2 の壁関数壁温も同種のモデル値で
あり、モデル値同士の比較は正当。η/ṁ は全経緯を通じ無汚染。
③ 格子系列は壁処理レジームが
異なるため Richardson/GCI は非適用 (単調収束系列でない)。

**y+≈1 low-Re 三者基準解 (2026-08-11, Step 1 完了)**: 専用メッシュ 221×81 `wall_first_frac 1e-4`
(AR 381 / skew 0.42 PASS, bell y+ mean ≈0.93) で forge cell / forge node / SU2 の低 Re 三者比較を実施
(`problem_bell_yp1*.yaml`, run_0033–0037)。熱物性は SU2 と整合 (`thermCondMethod: 1, prandtlLam: 0.72`,
`turbulentPrandtl: 0.9` — **Prt 0.85→0.9 で壁温 +20.5K の強感度**を発見し 0.9 を正式条件とする)。

| 量 (bell 部) | forge cell (run_0036) | forge node (run_0037) | SU2 (run_0035) | 理論 T_aw |
| --- | --- | --- | --- | --- |
| T_w 面積平均 [K] | 1388.9 | 1387.2 | 1414.2 | ~1400 |
| T_w @x=30mm [K] | 1388.5 | 1391.1 | 1416.6 | 1406 |
| τ_w @x=30mm [Pa] | 1707 | 1754 | 1756 | — |
| ṁ [kg/s] | 1.2993 | 1.2959 | 1.3017 | — |
| η_CF | 0.9779 | 0.9905→**内部列 0.970–0.977** | 0.9796 | — |

結論: ① cell/node は壁温 1.7K 差・τ_w ~1–3% 差で相互整合し、理論 T_aw 1400K の −1% 以内。
② SU2 は bell で +1% 高めだが**チェンバで Tt+76K の非物理超過**があり厳密な正解ではない —
壁温の真値帯は **~1390–1415K (理論 1400K ± 1%)** とみなす (wf=1 の 1193K に対し熱的閉包欠落
−210K が本質、Step 2 で対処)。③ η_CF の低 Re y+1 値は **cell 0.978 / SU2 0.980 / node 内部列
0.975** に収束 — 従来の壁処理依存幅 0.955〜0.988 の下限 0.955 は「未解像 low-Re」、上限 0.988 は
「node 出口アーティファクト」でありいずれも真値でない。④ 発見・修正したソルバ/ツールバグ 3 件
(双対 CV 幾何 float32 桁落ち・interp_field node 照会点・node 出口列不整合 [根治済み 2026-08-11,
run_0043–0045]) は
[plan §2.9–2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md) 参照。
比較図: `run_0035_su2_yp1_lowre/wall_temperature_compare_yp1.png`。

**壁温の三者比較と原因 (2026-08-11 更新)**: SU2 low-Re は解析回復温度に近く (1431K)、forge は
cell/node とも約−240K (1186/1196K)。判別実験 D2 で forge `wallTreatmentSST=0` は同じ SU2 low-Re と
壁温1K差・η 0.25%差へ整列する。ただし両者とも `y+≈19–40` でlow-Reとして未解像なので、これだけで
`wf_pk` を欠陥と断定してはならない。追加対照 `run_0030_su2_sst_wallfunc` では、x=30mm第一内点の
`μt/μ=6.88` (forge wf=1 cellは3.80)・T=1015K (forge 997K) と近いまま、SU2はCrocco–Busemann
熱閉包により壁温1427Kを得た。したがって**壁面出力の約240K差の直接原因はforge SST壁関数に
断熱回復温度の熱閉包がないこと**であり、「BL全体のμt不足だけが原因」という帰属は棄却する。
一方 η 0.955〜0.988 の真値と壁応力側の妥当性は、`y+≈1` low-Re と `y+>30` wall-function の
格子系列で未決。元のD1–D3記録と限界は
[調査メモ](../../notes/sessions/wall-temperature-three-way-analysis.md)を参照。

**レビュー対応の実装 (2026-08-11)**: ① `check_quasisteady.py` の step パースバグ修正
(res_outlet_*/拡張子数字の混入 — 既報 ALL STEADY は修正後ツールで再確認済み・結論不変)。
② **`Taw_diag` 診断出力**を SST 壁関数に追加 (Crocco 型 T_aw = T_rep + Pr^{1/3}·U_t²/(2cp)):
現既定 node 構成で **ベル部 mean 1420.2K = SU2 壁関数の壁温 1422K と 2K 差** — 熱的閉包を
適用すれば SU2 と一致することをフォワード確認。壁状態への適用は保存整合の plan 化後。
③ **`thermCondMethod: 1`** (constant-Pr: k=μ(T)cp/prandtlLam, 既定 0=従来) を追加 —
スモーク (scratchpad D4) で壁温 −25K シフト (SU2 の逆向き +14.7K と符号整合、主因でないことも一致)。

