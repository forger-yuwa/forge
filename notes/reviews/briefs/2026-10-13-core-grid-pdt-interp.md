# 諮問ブリーフ: 前処理の刻みの拡大だけを止めた A/B (§4.87、判定不能) の解釈と次の一手

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.82〜§4.87)
- **date**: 2026-10-13
- **諮問の条件**: AGENTS.md「エスカレーション条件」3 (判定不能)・4 (登録に無い run の前)・7 (解釈の確定の前)
- **ユーザの決定**: 「1 かな」(時間項の前処理でなぜ縮流部の残差が増えるのかを調べる)。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.f64`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
grep -n '^### 4\.8[6-7]' plans/active/tooling-nozzle-core-grid.md     (§4.86・§4.87 を sed で 30 行ずつ)
grep -n 'implicit_defect_correction_block_precond_d\|implicit_defect_correction_block_d\b\|__global__ void' solver_density_cuda/cuda_forge/timeIntegration_d.cu | head -30
sed -n '1276,1300p;1345,1480p' solver_density_cuda/cuda_forge/timeIntegration_d.cu     (前処理の block カーネル)
sed -n '1790,1830p' solver_density_cuda/cuda_forge/timeIntegration_d.cu                 (カーネルの選択と引数)
sed -n '15,60p' solver_density_cuda/cuda_forge/lowMachPrecond_d.cuh
sed -n '1,30p' case/45.isobutane_m6_d155/_band_ab/core_grid/pdt_judge.json
```

標準の block カーネル (`lowMachPrecond` < 2 のとき) と前処理の block カーネルの違いは、上の grep で両方の位置を確かめ、必要な範囲を 60 行以内ずつ追加で読んでよい。

## 1. 結果 (§4.87)

- 2 本とも `run_0608` res_60000 から point・cfl_pseudo 1・`lowMachPrecond: 3`・各 2 万 step、バイナリは同じ (sha 4e6e6c14、診断のスイッチ入り)。S+ (`run_0692`) は拡大あり、S0 (`run_0693`) は `FORGE_DIAG_PRECOND_NO_DT_SCALE=1` で拡大なし。資格 5 項とも満たした (刻みの比 6.03 で期待どおり、S+ は §4.83 の P3 を 0.9997〜0.9998 で再現)。
- 減衰率 d = q̄(16k〜20k)/q̄(10k〜14k): throatwall (x/r_t −2〜0・壁から 0〜10 節点目) の ρ・ρE は S+ 3.12・2.98、S0 0.95・1.02。縮流部 (x/r_t −8〜−2) の ρ・ρE は S+ 1.46・1.46、S0 1.05・1.03。縮流部の ρω は S+ 1.41、S0 0.65。
- 判定: 判定不能 (中間、S0 の d が A の ≤ 0.8 にも B の ≥ 1.1 にも当たらない)。
- 縮流部の ρ の q: S+ 677 (2000) → 676 (10000) → 1300 (20000)、S0 477 → 679 → 817。
- `check_convergence`: S+ NOT CONVERGED (rms_ro・rms_roe RISING)、S0 NOT CONVERGED (全列 falling、rms_ro 0.9 桁)。

## 2. 登録外の記述

- `=3` は右辺を変えないので、S0 と §4.83 の C (`run_0688`、`=0`、バイナリ d58ed787 = 4e6e6c14 から診断のスイッチを除いたもの、同じ起点・設定) は同じ離散方程式を解く。step 20000 の縮流部の ρ の q は S0 817・C 49 (C は 217 → 73 → 49 と下がり続ける)。刻みを広げずに Γ の時間項を入れても、縮流部の残差は C の約 17 倍の水準で横ばい。
- ただし前処理の腕は、標準の block カーネルではなく `implicit_defect_correction_block_precond_d` (物理の FVS + Γ の時間項 + Sherman–Morrison) で解かれる。Γ 以外の違い (粘性・軸対称のソース・行列の精度など) は確かめていない。`precondEps: 1.0` なら β = 1 で Γ = I・Δτ′ = Δτ になるはず (`lowMachPrecond_d.cuh` の U_r = min(c, max(|u|, εc)))。M6 の移植の回帰 (float plan §6.57 の R2) は β = 1 で 1 step の境界の修復を確かめたが、反復の振る舞いは見ていない。

## 3. 判断してほしいこと

推奨を 1 つに絞り、事前の読み方つきで:
1. §4.87 の判定不能をどう扱うか (どこまで言えるか)。
2. 次の一手。候補 (縛られなくてよい): (a) `=3`・`precondEps: 1.0` (Γ = I・拡大なし) 対 `=0` で、前処理のカーネルそのものが反復を悪くするか (カーネルの違い) を分ける。(b) Γ の時間項の大きさを変える別の切り分け。(c) 低マッハの線をここで閉じ、何を結論として残すか。
3. ユーザの元の問い (流れ方向に細かいとライン陰解法が遅くなる仕組み) と今回の調査のつながりについて、今言えること 2〜3 文。
