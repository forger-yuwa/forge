---
name: sern-eval
description: ⑤ SERN (case/46) の評価 run を組む・回す・判定するときの手順。起動レシピ (MOC 初期値・梯子・本段 CFL)、収束判定の 5 条件、3D 固有の既知欠陥を現行の実測値で示す。「SERN を回す」「SERN の設定」「3D SST が落ちる」「CFL を上げたい」「MOO を回す」と言われたときに使う。
---

# /sern-eval — ⑤ SERN の評価 run

正本は [`plans/active/tooling-nozzle-sern-chain.md`](../../../plans/active/tooling-nozzle-sern-chain.md)
(§4.12 起動・§4.16 CFL・§4.16.1 収束判定・§5.1 残作業) と
[`case/46.sern_design/README.md`](../../../case/46.sern_design/README.md) の run 一覧。
**過去 run の YAML をそのままコピーしない** — 却下済みの設定 (断熱壁・旧トポロジ) が残っている。

## 1. 起動レシピ (2D)

> **検証済み (2026-09-19, R-b)。** 「梯子 1500×3 + 本段 500 = 5000 step」を 3 設計 (`L_cowl` 1.2 / 3.0 /
> θ_r0 21°) で長時間参照 (4000×3 + 本段 12000 = 24000 step) と直接比較し、**全比較量が事前宣言した許容内**:
> `C_T_with_shear` 相対差 2.400e-6 / 3.413e-6 / 2.991e-6 (A/B/C)、`C_M` 4.842e-5 / 5.956e-5 / 6.354e-5、
> 壁圧 L2 ≤ 1e-3、`vis_turb` L2 ≤ 9e-4、剥離なし、質量不均衡 ≤ 0.002 %。
> run は `case/46.sern_design/run_0194`–`run_0199`、詳細は
> [`plans/active/tooling-nozzle-sern-startup.md`](../../../plans/active/tooling-nozzle-sern-startup.md) §5.1 R-b。
> **注意**: `vis_turb` の「500 step で DRIFTING」は準定常判定の話で、参照との**差**は 1e-3 以下 (上表の比較)。
> 単一 run の `check_quasisteady` は別途必ず見ること。

| 項目 | 生産値 | 根拠 |
| --- | --- | --- |
| `gas.model` | `frozen_tp` | 排気 = CEA 凍結組成 EXH、外気 = AIR。熱力学は NASA-9 |
| `spec.wall_thermal` | 等温 1000 K | 断熱は M∞10 で回復温度 4600 K になり残差がプラトー (R4b) |
| `mesh.ic` | **`moc`** | 一様 IC は入口状態を全域に置き出口で 17 倍ずれる。MOC 場を外挿して埋める |
| 暖機 (層流, 1 次) | **1500 step, CFL ramp 0.1 → 0.3 → 1.0** | `opt.warm_lam_ramp`。runner が forge を 3 回起動する |
| soft (SST, 1 次) | **1500 step, CFL 1.0** | `opt.soft_cfl` |
| mid (SST, 2 次) | **1500 step, CFL 1.0** (`soft_cfl` と共有) | 次数と CFL を同時に上げない |
| 本段 (SST, 2 次) | **500 step, CFL 5.0** | 上限は 6 (12 で発散)。合計 5000 step |
| `evaluate.implicit_relax` | **使わない** | 0.7 は定常解を 0.05 % 動かす。CFL を上げる方で速くする |
| `evaluate.p_min` | 20.0 | 既定 1.0 Pa は M∞10 の膨張で割れる |
| `evaluate.outlet_kind` | **`outflow`** (3D) | `statPress` の Ps は実出口圧より桁違いに低く、node の亜音速壁列から unstart |

**CFL の上限は段で 50 倍違う**。暖機の入りは 0.1、本段は 5〜6。梯子は外せない (一様 IC からの直接起動は cfl 0.1 でも発散)。
**soft と mid は `soft_cfl` を共有していて独立に振れない** (段ごとの上限は未測定)。
`opt.warm_lam_cfl: 0.1` は ramp の初段として残る (ramp を書かないときの固定値)。
**過去の「梯子は全段 CFL 0.1」という記述は誤り** — R-b で実際に回した入力は soft/mid とも 1.0 (codex plan-2 M2)。

## 2. 収束の判定 (5 条件すべて)

`metrics/sern_gates.py` の `evaluate_gates` が返す `verdict` を**必ず貼る**。

1. 保存場に NaN/Inf 無し・正値
2. **rising な残差列が無い** (プラトーまたは低下中)。旧「全列プラトー」は実装バグ (3 桁低下列が素通り) を
   修正したうえで `require_residual_plateau` を**既定 OFF に降格**し、`gates["plateau"]` として診断のみ報告する。
   プラトーは収束の十分条件ではない (残差の大きさを問わないので float32 の更新消失やクランプで止まった状態と
   区別できない) ため、4. の床検査と 5. の量の定常性で補う。<br>
   **rising 判定はジッタとの比較込み (2026-09-19, R-h)**: 末尾窓の上昇桁数がその列自身の散らばり σ を超えた
   ときだけ rising。これが無いと、プラトーが数倍の幅で揺れる列 (SERN の `rms_roY1` が典型) で判定が
   **同一設定の run 同士でも PASS/FAIL に割れる**
3. 残差列の桁が揃う (1 列だけ中央値の 1e6 倍以上でない)
4. 床・下限への張り付き 0 ノード (`pMin`/`tMin` 50 K/`roMin`/`roOmega` 1e-20/`k`≤0)
5. `C_T_with_shear` ほか力係数が STEADY

残差 3 桁低下は課さない (1–2.5 桁でプラトーする性質)。**「まだ低下中」を許すと過渡途中の値を定常値として採る**。
CFL を上げるのは収束を速めるためで、最終的に収束していることが条件。

## 3. 3D の既知欠陥 (評価に使う前に確認)

- **側面バンド末端の壁** (R4e、未解決): `mesh_sern3d` が幅外流体の下流端を `vehicle_side` で塞ぎ、
  遠方境界との交線で最大 5.72 MPa (外気 2851 Pa)。**この欠陥がある run は評価に使わない**
- **`L_sw = L_cowl` は退化**: カウル後縁と側壁後縁が交線を共有し三重点で発散。`mesh3d.L_sw` を 0.8 などにずらす
- **`rms_roOmega` の張り付き**: `sym` 上 1 ノードで `k=0`・`roOmega` 下限 → 交差拡散が ω 下限で割られ RMS 1e15。
  条件 2/4 で落ちる
- `vehicle_top` は等温粘性壁が既定 (`evaluate.vehicle_top_kind: slip` で旧挙動)

## 4. 道具

- `solver_density_cuda/tools/mirror_symmetry.py` — 対称面で鏡像し全幅を復元。複数の壁面出力を
  `patch` id つきで束ねられる。ParaView は生成される `.xmf` を開く
- `python3 solver_density_cuda/tools/check_convergence.py <run_dir>` / `check_quasisteady.py` / `check_mesh_quality.py`

## 5. 報告

run パス (リポジトリルートからの相対)、ゲートの verdict、`case/46.sern_design/README.md` の run 一覧を更新。
**`forge_run.log` は段ごとに上書きされる**ので、そこの `detectNaN` 文字列で成否を判定しないこと (誤判定の常習)。
