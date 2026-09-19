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

**収束の根拠**: 外部ループは `cht_history.csv` の `dTw_max` 7.2e-5 K・`res_rel` 7.2e-6 が 2 反復連続で
許容以下 (`--tol-K 1e-3 --tol-rel 1e-3`)。各反復の CFD は**壁温固定での緩和試験**で 10000 step 後に
$q$ が 4 桁一定になることを確認済み (上記)。**`check_convergence.py` はこの case では意味を持たない**:
静止・純伝導なので `rms_ro` が初期から機械ゼロ近傍で、ツールは `stalled/plateau` を返す
(判定は解析解との照合と界面残差で行う)。
