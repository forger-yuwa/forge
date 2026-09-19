# 05. sod shock tube

## 3D periodic 化検証 (median-dual M4 §4.5)

2D/1D sod を 3D box ([0,1]×0.1×0.1, 200×8×8 hex) に押し出し、transverse (y,z) を Cartesian periodic にして node periodic を transient 衝撃で検証。**物性・IC は run_0002_slau_cpg と同一** (cpg 物理単位: 左 P1e6/T2000, 右 P1e5/T400, roe=P/(γ-1), cp1038.8/γ1.4)。**教訓**: 非次元 `initial:"sod"` (P=1) は forge thermo と非互換 (T≈0.003K で roe 破壊→P 均一) なので使わない。cfl は run_0002 同様 0.2。

| run_* | 目的・設定 | 結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_node3d_periodic_sod` | **node + y,z periodic**, run_0002 同一 config (cpg 物理単位 PL1e6/TL2000・PR1e5/TR400, cfl0.2, SLAU, RK3) | **完走 (200 step NaN なし)・物理 sod 構造 (ro 0.84..1.92, P 1e5..1e6, 衝撃右伝播)・spanwise std=0〜1e-6 完全均一** = node periodic が transient 衝撃で機能。node は cell より前面がクリーン (決定的) | **検証完結** |
| `run_sod3d_cell_ref` | cell + y,z periodic (同 config, 参照) | 完走・物理 sod。滑らか域 spanwise std=2e-7 均一、衝撃前面のみ std=0.5 波打ち (cell atomicAdd ノイズ、periodic とは別) | ref |
| `run_node3d_periodic_regr3dgeo` | 3D median-dual 幾何堅牢化 (plan architecture-median-dual-3d-double-geometry) の新 converter で sod3d.msh を再変換し再実行 | 変換診断は健全 (体積 relErr 1.2e-7・閉性 3.6e-6)・1000 step NaN なし。ただし copy した config が非次元 sod IC (下記 IC バグ) で場は非物理 — **幾何検証の本体は thinwall 旧新比較 (plan 参照) であり本 run は補助** | 削除済み(2026-08-31) |
| `run_node3d_slip_lowcfl` | **IC バグ調査の記録**: 当初 `initial:"sod"` (非次元 P=1) は forge CPG thermo (cp1038.8,R296.8) と非互換で T≈0.003K→roe 破壊→**P 均一の異常状態**。これに node checkerboard が乗り発散していた | 削除済み(2026-08-31) (IC バグ記録) |
| `run_0009_aws_node_smoke` (EC2) | **AWS P1 スモーク** ([plans/active/tooling-cloud-gpu-env.md](../../plans/active/tooling-cloud-gpu-env.md)): g5.xlarge (A10G) コンテナで gmsh 生成→convert→node 実行。IC は `gen_sod3d_ic.py` で物理単位 cpg を焼き込み (組み込み `initial:"sod"` は現行も非次元で使用不可) | 完走・NaN なし・mesh 品質 PASS。**run_0011 (同一コード・同一入力のローカル RTX3060) と最大相対差 1.6e-6 / L2 1e-7 で一致** = AWS 環境妥当 | ref (EC2 上) |
| `run_0010_aws_cell_smoke` (EC2) | 同上の cell 版 | 完走・NaN なし。run_0012 と最大相対差 2.2e-6 / L2 1.2e-7 で一致 | ref (EC2 上) |
| `run_0011_aws_repro_node` | run_0009 の**ローカル再現** (EC2 生成の sod3d.h5 を bit 同一入力に使用、現行 HEAD native ビルド) | AWS↔ローカル一致の判定根拠 (上記) | ref |
| `run_0012_aws_repro_cell` | run_0010 のローカル再現 | 同上 | ref |
| `run_0013_monitor_explicit` | run_0012 の複製。**console モニタ行** ([plans/accepted/architecture-runtime-monitor-line.md](../../plans/accepted/architecture-runtime-monitor-line.md)) の非定常陽解法 (RK3, 固定 dt) 検証: `t`/`dt`/`maxCFL` + ms/step + 残差要約。`_baseline/` は旧バイナリ比較用 | 完走・NaN 0。`t/dt/maxCFL 0.07` 表示、log 1512→312 行、CSV 旧×2/新×2 列別比較 (`csv_compare.txt`) 合格 | 破棄予定 (ログ形式検証) |
| `run_0014_monitor_explicit_adaptdt` | run_0012 + `deltaT.control: 1` (CFL 0.2 適応) + `monitorInterval: 5`。モニタ行の `dt` (使用)/`maxCFL`/`dt_next` の刻み対応と `totalTime` 加算順修正の検証 (codex result Major 1)。`_baseline/ _baseline2/ _new2/` は CSV 比較用 | 完走・NaN 0。step 1: `t 1.0969e-06 dt 1.10e-06 maxCFL 0.20 dt_next 1.09e-06` (t = 使用 dt)。`csv_compare.txt` 合格 | 破棄予定 (ログ形式検証) |
| `run_0020_limA_prod` / `run_0021_limC` / `run_0022_limD` | **リミッタ無次元化の回帰** (`run_node3d_periodic_sod` と同一入力、`limiterMatchRecon`/`limiterScaled` のみ変更: A=0/0 生産既定, C=1/1 無次元化, D=1/2 比の形)。旧キー `LESorRANS`/`LESmodel` は `model: "none"` に置換 | **C/D が A より厳密 Riemann に近い**: L1 誤差 ro 7.048e-3→6.908e-3, P 3.212e-3→3.062e-3, u 6.048e-3→5.909e-3。衝撃後 ρ 1.9913→**2.0517** (厳密 2.2240) で鈍りが減る。厳密解上限超え 0 個・P 範囲外 0 個 | ref ([../../plans/active/convection-node-wall-reconstruction.md](../../plans/active/convection-node-wall-reconstruction.md) §4.22) |
| `run_0023_v2_bitcheck` / `run_0024_v2_bitcheck2` | **W2 V2: 既定パス不変の確認** (`run_0020_limA_prod` と同一入力を新バイナリで 2 回) | 旧 vs 新 **1.19e-4** < 同一バイナリ 2 回走の**ノイズ床 1.51e-4**。node 3D periodic は同一バイナリでもビット再現しないので判定はノイズ床比較で行う | ref (plan §4.24 V2) |
| `run_0025_perA` / `run_0026_perC` / `run_0027_perD` | **周期経路に `limiterScaled` を配線した後の再実施** (codex plan-3 Critical 1: 周期 node は `limiterPeriodic_d.cuh:115` が旧 Venkat を呼び `limiterScaled` を読んでいなかった = `run_0021/0022` は W1 のみの結果だった) | **C/D が A より厳密 Riemann に近い**: L1 ro 7.048e-3→6.914e-3 (C) / 6.905e-3 (D), P 3.212e-3→3.072e-3 / 3.062e-3, u 6.048e-3→5.927e-3 / 5.908e-3。衝撃後 ρ 1.9913→2.0524 / 2.0517 (厳密 2.2240)。配線前後で ρ 9.9e-4 / P 677 Pa 動いており配線が効いていることを確認 | ref (plan §4.22 の正) |
| `run_0028`/`run_0029_g1per{A,C}`, `run_0034`/`run_0035_g1per{A,B}200`, `run_0030`〜`run_0033_ksweep_*` | **周期経路の G1 検査 (W1h) と `venkatK` スイープ** (`limiter_g1_check_periodic_d` 追加。許容幅に絶対床 `1e-6·q_ref`) | **① と ② の優劣が逆転**: ① ρ 49712 (最大逸脱 22.7×)・他 0 / ② K=1.0 は全変数悪化 (ρ 86674, 速度 153157, 最大 308×) / **K を下げると単調改善し K=0.05 (SU2 既定) で ρ 逸脱 0**・残る速度 P は最大 3.4×/1.9×。L1 精度は K にほぼ非依存 (6.914e-3→6.905e-3)。→ **`venkatK` 既定 1.0 は不適、0.05 にすべき** | ref (plan §4.26) |

**結論**: node sod 3D periodic が **安定・物理的・spanwise 完全均一**で動作 (検証完結)。当初「node-explicit-shock 限界」とした発散は**誤診**で、真因は壊れた IC (非次元→P均一) + 過大 CFL だった (user 指摘で判明)。periodic 実装は transient 衝撃でも健全。

### 厳密 Riemann 解との照合 (t=5e-5s, 3パネル density/pressure/velocity)
図: `sod3d_density_compare.png`, `sod3d_compare_3panel.png`。**node が厳密解と一致、cell が横方向不安定で破綻**:
- **node**: 密度(膨張波・接触 ro0.94・衝撃後~1.98(理論2.224をやや平滑化)・右0.842)・圧力(p*4.4e5→1e5)・流速(u*~500)とも厳密 Riemann に一致。理論最大密度2.224 を超えるノード=0、P≤IC範囲。横方向は衝撃面で微小波打ち(x=0.54 std0.43)だが有界・物理的。
- **cell**: 同一 IC・config なのに rms_roUy=10.13 の横方向不安定で右状態が汚染 (ro 1.59≠0.842, P 6.6e5≠1e5, 局所 ro 4.213 で理論超え1256ノード, P 1.43e6 非物理)。cell atomicAdd 非決定性 / SLAU carbuncle が grid-aligned 衝撃で暴走。
**結論**: median-dual **node** が 3D Sod を厳密解通り正確に解く (cell より優秀)。periodic も滑らか域で完全均一。残課題: 衝撃面の微小波打ち (SLAU 衝撃感度、roe 等で改善余地)。
