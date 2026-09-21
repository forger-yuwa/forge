# case/57 遷移平板 (ERCOFTAC T3A / T3B)

層流から乱流への遷移を forge が再現できるかを見る検証ケース。遷移モデルは Langtry–Menter 2009 の $\gamma$–$\tilde{Re}_{\theta t}$
(`turbulence.transition: lm2009`)。設計判断と合否基準は [`plans/active/turbulence-transition-lm2009.md`](../../plans/active/turbulence-transition-lm2009.md) §6、
式は [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §11。

## 条件

平面 2D の平板。前縁の 0.04 m 上流に入口を置き、助走は slip、$x\ge0$ が断熱 no-slip 壁。上面 slip、出口は静圧。
流速は $M$ 0.2 (69.44 m/s, 300 K, 101325 Pa) に固定し、**単位長さあたりのレイノルズ数を粘性で実験に合わせる**。

| ケース | $Re/m$ | 粘性 [Pa·s] | 入口 $Tu$ | 入口 $\mu_t/\mu$ | 入口 $k$ [m²/s²] / $\omega$ [1/s] |
| --- | --- | --- | --- | --- | --- |
| T3A | 3.6e5 | 2.269e−4 | 3.3 % | 12 | 7.87 / 3401 |
| T3B | 6.27e5 | 1.3033e−4 | 6.5 % | 100 | 30.56 / 2759 |

実験値: `ref/t3a_exp.dat` ($x$ [mm], $C_f$, $Tu$ [%]。OpenFOAM の T3A チュートリアル同梱の ERCOFTAC データ)。T3B の実験値は未入手
(ERCOFTAC のサーバに接続できなかった。2026-09-22)。

## メッシュ

`gen_mesh.py` が gmsh の構造格子から `.msh` / `.su2` / node 用 `.h5` を作り、`check_mesh_quality.py` を通す。格子系列は両方向 $\sqrt2$ 倍ずつ。

| タグ | 節点数 (壁法線方向の層数) | 第一層 [µm] | 品質 VERDICT |
| --- | --- | --- | --- |
| `plate_coarse` | 12880 (56) | 28.3 | PASS (AR ≤ 1000, skew ≤ 0.90) |
| `plate_base` | 25350 (78) | 20.0 | PASS (AR 最大 537) |
| `plate_fine` | 50490 (110) | 14.1 | PASS |

## 手順

```bash
python3 gen_runs.py --run run_NNNN_t3a_sst --mesh plate_base --model sst                 # 一様場から段階起動 (完全乱流の SST)
python3 gen_runs.py --run run_NNNN_t3a_lm  --mesh plate_base --model lm --ic-from run_NNNN_t3a_sst --reset-turb --main-steps 100000
python3 tools/cf_plate.py --forge run_A:ラベル --su2 su2_t3a_lm:SU2 --exp ref/t3a_exp.dat --out cf.png   # forge と SU2 を同じ式で Cf 化
python3 tools/cf_plate.py --series run_A run_B --series-out run_B/cf_series.csv                           # 報告量の時系列
python3 ../../solver_density_cuda/tools/check_quasisteady.py --series-csv run_B/cf_series.csv --series-cols x_onset,cf_min,x_end,cf_max,cf_x0.3,cf_x0.9,cf_x1.3
python3 ../../solver_density_cuda/tools/check_lm_kernel.py run_A      # output.level 2 の LM run: カーネルと numpy 参照の突き合わせ
```

- **出発場に依らない**: 完全乱流の SST 場からそのまま遷移モデルを入れても (`run_0011`)、$k$/$\omega$ を入口値に戻してから入れても (`--reset-turb`, `run_0004`→`run_0005`)、
  同じ解に落ちる (遷移開始・終了位置と $C_f$ 3 点が 4 桁一致)。「乱流境界層の中では $\gamma$ を減らす項が働かないので層流域が戻らないのでは」という事前の懸念は、平板では当たらなかった。
- $C_f$ は forge・SU2 とも**壁から第一内部節点の速度**で $\tau_w=\mu u_1/y_1$ として出す ($y_1^+<1$)。基準動圧は入口の一様流。
  遷移開始 $x_{onset}$ は $C_f$ が最小になる位置、遷移終了 $x_{end}$ はその下流で $C_f$ が最大になる位置。
- **全 run が `NOT CONVERGED (stalled/plateau)`**。床は単精度の丸め (`rms_ro`/⟨ρ⟩ が float32 eps の 0.8〜1.6 倍。SST 単体・層流でも同じ) で、残差は全域に一様に分布する
  (`run_0015`、`tools/residual_map.py`: 上位 10 節点の寄与 5 %、前縁 10 mm 以内 0.0〜0.1 %)。当初「前縁の特異点」と書いたのは誤り。
  したがって結果は**未収束・報告量は準定常**として扱い、量は `check_quasisteady.py --series-csv` の VERDICT で判定する。
- 局所 $y_1^+$ は前縁を含む全壁で 粗 1.67 / 基準 1.31 / 細 1.02 が最大 (いずれも前縁の節点)。1 を超える壁長は 0.34 / 0.12 / 0.08 %。run 一覧の「$y_1^+$ 最大」は $x>$ 0.02 m の値。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `su2_t3a_lm` | **SU2 8.5 SST (V2003m) + LM の参照解**。`plate_base` と同一メッシュ・同一条件、30000 反復、CFL 10 | rms 低下: ρ −4.0 桁、k −6.2 桁、ω −3.1 桁、γ −3.5 桁、$\tilde{Re}_{\theta t}$ −2.3 桁。$x_{onset}$ 0.3675 m、$x_{end}$ 0.8635 m、$C_f$(0.3/0.9/1.3 m) = 0.00227 / 0.00441 / 0.00418。`restart_flow.csv` | ref |
| `su2_t3b_lm` | 同 T3B 条件 | 実行中 | active |
| `run_0001_regr_none_head` / `run_0002_regr_none_new` | **`transition: none` の非退行**。変更前 (HEAD) と変更後のバイナリで SST 1 次 300 step、各 3 反復 (反復分は scratchpad)。変更前バイナリ自体が反復でビット一致しないのでノイズ床で判定 | `check_field_regress.py --boundary`: **VERDICT PASS** (全量 ノイズ床 × 2 以内、最大比 1.83) | 破棄済み |
| `run_0003_t3a_sst` | T3A・`plate_base`・遷移モデルなしの SST。一様場から soft → mid → ramp (cfl 1, 2) → main (2 次, cfl 5, 20000 step) | 前縁から完全乱流。$C_f$(0.3/0.9/1.3) = 0.00489 / 0.00412 / 0.00390。`NOT CONVERGED (stalled/plateau)` (平均流 1.3–1.5 桁で床、k 3.0 桁) | ref (LM の出発場) |
| `run_0004_t3a_lm` | **最初の LM run**。`run_0003` の流れ場 + $k$/$\omega$ を入口値に戻して ramp (cfl 1) → main 20000 step。`output.level 2` (診断場つき) | 20000 step では遷移後半が未発達 ($x_{end}$ 1.03 m)。`check_lm_kernel.py`: **VERDICT PASS** (カーネル vs numpy 参照、全 8 量で差 ≤ 6e−6、相関の反復 平均 3.3 / 最大 6 回) | ref (`run_0005` の出発場) |
| `run_0005_t3a_lm_cont` | 上の継続 80000 step (通算 100000)。`cf_series.csv`・`QUASISTEADY_VERDICT.txt` | **$x_{onset}$ 0.3595 m (SU2 比 −2.2 %)、$x_{end}$ 0.8496 m (−1.6 %)、$C_f$ 最小 0.002240 (+1.8 %)、$C_f$(0.3/0.9/1.3) = 0.002289 / 0.004413 / 0.004172 (SU2 比 +0.8 / +0.1 / −0.2 %)**。実験の $C_f$ 最小は測定点 0.395 m (測定間隔 0.1 m) で、それに対し −9 %。全 7 量 **STEADY** (末尾窓の drift ≤ 1.9 %、最も動く $x_{end}$ は 50000 → 100000 step に 0.843 → 0.850 m)。`NOT CONVERGED (stalled/plateau)` = 丸めの床 | active |
| `run_0012_t3b_sst` / `run_0013_t3b_lm` | **T3B** ($Tu$ 6.5 %、$Re/m$ 6.27e5)。SST 段階起動 → LM 100000 step (k/ω は戻さない) | 遷移開始 0.090 m ($Re_x$ 5.6e4)・終了 0.189 m、$C_f$($x$ 0.1/0.2/0.6/1.3) = 0.00520 / 0.00533 / 0.00450 / 0.00390、$y_1^+$ 最大 0.75 ($x>$ 0.02 m)。全 8 量 **STEADY**。`NOT CONVERGED (stalled/plateau)`。SU2 (`su2_t3b_lm`) との照合は SU2 完走後 | active |
| `run_0014_t3a_lm_unitcheck` | **カーネル単体検査 (強化版)** と初期出力・restart の確認。`run_0011` の準定常場から cfl 0.01・50 step・level 2 | `check_lm_kernel.py`: **VERDICT PASS** (11 量・全節点、最大で許容の 0.21 倍、相関反復 平均 3.3 / 最大 6 回、下限の作動 0 %)。`res_0` の遷移 5 変数が有限、`roGamma`/`roReth` は継続元とビット一致 (入口ピン 2〜5 節点のみ 1e−7)。`LM_KERNEL_VERDICT.txt` | active |
| `run_0015_t3a_lm_resmap` | **残差の床の所在**。`run_0011` から 2000 step・level 2 (平均流の残差場つき) | `RESIDUAL_MAP.txt`: 平均流の残差は全域一様 (上位 10 節点 5 %、前縁 10 mm 以内 0.0〜0.1 %)。`res_roOmega` だけ壁第一層に局在 (上位 10 節点 38 %) | active |
| `run_0006_t3a_sst_coarse` / `run_0007_t3a_lm_coarse` | 格子系列 (粗)。SST 段階起動 → LM 100000 step | $x_{onset}$ 0.3415 / $x_{end}$ 0.8209 m、$C_f$ = 0.002282 / 0.004400 / 0.004153、$y_1^+$ 最大 0.63。全 7 量 **STEADY**。`NOT CONVERGED (stalled/plateau)` | active |
| `run_0008_t3a_sst_fine` / `run_0009_t3a_lm_fine` | 格子系列 (細)。SST 段階起動 → LM 140000 step | $x_{onset}$ 0.3581 / $x_{end}$ 0.8694 m、$C_f$ = 0.002288 / 0.004431 / 0.004179、$y_1^+$ 最大 0.32。**基準 → 細の変化: 遷移開始 −0.4 %、$C_f$ +0.0 / +0.4 / +0.2 %**。全 7 量 **STEADY**。`NOT CONVERGED (stalled/plateau)` | active |
| `run_0010_t3a_lam` | 遷移・乱流モデルなしの層流 (`run_0003` の場から 60000 step)。層流ソルバの検査 | $C_f$ の Blasius 比 −1.8 % ($x$ 0.05) 〜 +2.4 % ($x$ 1.4 m)。残差の床は LM・SST と同じ (`rms_roe` 5e−2) → 床は乱流・遷移モデルと無関係。`NOT CONVERGED (stalled/plateau)` | active |
| `run_0011_t3a_lm_noreset` | **出発場の感度**。`run_0003` (完全乱流) の $k$/$\omega$ を戻さずに LM を入れて 100000 step | `run_0005` と遷移開始・終了・$C_f$ 3 点が 4 桁一致。全量 **STEADY** | active |
