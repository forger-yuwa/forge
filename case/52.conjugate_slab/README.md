# case/52 共役スラブ (CHT の V1 検証)

**1 次元純伝導の共役解**で、forge (流体) + `solid_shell.py` (固体) の**外部弱連成ループ**を解析解に照合する。
計画は [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md) §6 V1、
仕様は [`methods/boundary.md`](../../methods/boundary.md)「共役熱伝達 (CHT)」。

## 問題

静止した流体層 (厚さ $H$=0.01 m、$k_f$=0.0241 W/mK 一定) の上面を 350 K に固定し、
下面を固体 (全抵抗 $R_s=t/k_s+R_{back}$=0.2 m²K/W、背面 300 K) と連成させる。直列抵抗なので

$$R_f=\frac{H}{k_f}=0.414938,\quad q=\frac{T_\infty-T_b}{R_f+R_s}=81.3090\ \mathrm{W/m^2},\quad T_w=T_b+qR_s=316.2618\ \mathrm{K}$$

**流れは無い** (静止・重力なし) ので、解析解は圧力に依らない。

## 実行

```bash
python3 gen_mesh.py                      # 一様 20x16 quad -> template/mesh.h5 (静止 IC, p0=1013.25 Pa)
python3 ../../solver_density_cuda/tools/cht_loop.py run_0001_v1_slab \
    --template template --forge <forge> --solid solid.json \
    --phys-id 3 --phys-name wall_bot --max-iter 15 --tol-K 1e-3 --tol-rel 1e-3
python3 verify_v1.py run_0001_v1_slab    # VERDICT
```

## この case の設計上の注意 (実測に基づく)

- **壁近傍クラスタメッシュは使わない**。純伝導の擬似時間緩和は最小セル幅で律速され、
  case/24 由来の $d_1$=8.5e-5 m メッシュでは **40000 step でも $q$ が定常値の 6 倍**だった。
  一様 $h$=6.25e-4 m にすると 10000 step で定常。
- **低圧 (p0=1013.25 Pa) にして熱拡散率 $k/(\rho c_p)$ を上げる**。局所 $\Delta t$ は音速律速で変わらないので、
  $p_0$ を 1/100 にすると緩和に要る step 数がほぼ 1/100 になる (101325 Pa では 20000 step でも $q$=14 W/m²、
  1013.25 Pa では 10000 step で 120.3 W/m² = 厳密解 120.5 の −0.16 %)。解析解は圧力に依らない。
- **`cfl_pseudo` は 5**。case/24 の旧設定 (`cfl_pseudo: 100`) は**現在の HEAD では step 100 以内に NaN**
  (P 床・T 床に張り付いてから `ro` が非有限になる EOS 床洗浄)。クラスタメッシュでは 20 まで、
  本 case の一様メッシュでは 5 で安定。1 次/リミッタ無しでも NaN になるのでスキーム依存ではない。

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_v1_slab` | **V1 本体**: 外部弱連成ループ (11 反復, 各 10000 step)。`--flux q_compact`, $R_s$=0.2, 背面 300 K, 上壁 350 K | **VERDICT PASS**: $T_w$=316.2659 K (解析解 316.2618、温度上昇 16.27 K に対し **0.025 %**)、**両側 $q$ の不一致 0.0053 %** (81.3252 vs 81.3296 W/m²)、$q$ の解析解差 0.0199 %。壁面の $T_w$ ばらつき 6.1e-4 K。履歴 `cht_history.csv`、最終壁温 `Tw_final.csv`、各反復 `it_NNN/` | active (**V1 の根拠**) |
| `run_0002_v1_insolver` | **ソルバ内連成 (Phase 2a, `conjugate: {mode: local1d, ...}`)**: 同じ固体 ($R_s$=0.2)、20000 step、`interval: 50`、`warmup: 200` | $T_w$ 316.234 K まで来て更新量 3e-4 K/更新 (まだ緩和中) → `run_0003` へ継続 | 中継 |
| `run_0003_v1_insolver_cont` | 上の継続 (40000 step)。**再開は `conjugate_Tw_3.csv` → `wall_profile_3.csv`** で壁温を引き継ぎ (`ints: {conjugate: 1, wallProfile: 1}`) | **$T_w$=316.2371 K (解析解比 −0.152 % of rise)、両側 $q$ 不一致 0.0002 %** (81.1856 / 81.1855 W/m²)、更新量 1.5e-5 K で静定。外部ループ (run_0001, 316.2659 K) との差 **0.029 K = 0.18 % of rise** = 同じ壁温での流体側離散解の差 (±0.2 %) の範囲 | active (**Phase 2a の根拠 / V4**) |
| `run_0004_pinfix` / `run_0005_qeff` | `run_0003` の設定で 8000 step。`0004` = 等温壁エネルギーピンの修正後、`0005` = `iface_q_eff` 診断の追加後 (2026-09-20、本表への記載漏れを 2026-09-21 に補完) | `run_0005`: 下壁 `q_eff` 81.18368 / `q_compact` 81.18324 W/m² (差 −5.2e−6 相当)。`res_wall_{bot_3,top_4}_8000.h5` | ref |
| `run_0006_ctrl_newbin` | `run_0005` と同一入力を **2026-09-21 のバイナリ** (`iface_q_eff` から蓄積項 $H_wR_\rho$ を引く版) で再計算。`fx` A/B の対照 | 下壁 `q_eff` 81.18349 / 上壁 −81.51316 (`q_compact` −81.51314 に一致。旧 −81.50878) | active |
| `run_0007_fxhalf` | 上と同じ + `FORGE_NODE_FX_HALF=1` (node 内部双対面の `fx=0.5`)。[`discretization-node-face-weight-midpoint.md`](../../plans/active/discretization-node-face-weight-midpoint.md) §6 (c) | 下壁 `q_eff` 81.18345 (対照比 −5e−7) = **1 次元では不変**。`check_convergence`: `NOT CONVERGED (stalled/plateau)` = 対照と同じ (純伝導の丸め床) | active |

### 収束の扱い (AGENTS.md「収束確認」)

**`check_convergence.py` の VERDICT (必ず参照する)**:

| run | VERDICT |
| --- | --- |
| `run_0001_v1_slab/it_010` | `NOT CONVERGED (stalled/plateau)` |
| `run_0002_v1_insolver` | **`PASS (converged)`** |
| `run_0003_v1_insolver_cont` | `NOT CONVERGED (stalled/plateau)` (`rms_ro` init 1.23e-10 → fin 1.05e-10, drop 0.3 dec flat) |

**この case では残差ベースの VERDICT を「解が悪い」根拠として読まない**。静止・純伝導なので
残差は**初期から機械ゼロ近傍** (`rms_ro` 1e-10, `rms_roe` 4e-5) で、下げ幅が取れずツールは
`stalled/plateau` を返す。同じ設定の連続 run で `PASS` (run_0002) と `NOT CONVERGED` (run_0003) に
割れるのもこのため (run_0003 は既に静定した場から始まるので低下桁数が出ない)。

**したがってこの case の合否は次の 3 つで判定する** (いずれも上の表に数値を記載):

1. **解析解との照合** — `verify_v1.py` の VERDICT (`PASS`: $T_w$ 0.025 %、両側 $q$ 0.0053 %)。
2. **界面の静定** — 外部ループは `dTw_max` 7.2e-5 K・`res_rel` 7.2e-6 が 2 反復連続で許容以下。
   ソルバ内は更新量 max|d$T_w$| が 1.5e-5 K で一定。
3. **流体側の緩和** — 壁温固定の緩和試験で 10000 step 後に $q$ が 4 桁一定 (本 README 上部)。
