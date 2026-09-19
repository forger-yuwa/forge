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

> **⚠ 2026-09-19 時点で確定していない。** codex plan レビュー (NO-GO, M6) で、下表の
> 「梯子 1500×3 + 本段 500」は**一度も回していない組み合わせ**と判明し撤回した。
> 500 step 時点では `vis_turb` が未発達 (局所 DRIFTING)。**受理基準 (`require_residual_plateau`) にも
> 実装バグ**があり、3 桁低下してなお下降中の列が素通りする。再検証まで下表は**暫定値**として扱い、
> 生産採用しない。改訂後の残作業は [`plans/active/tooling-nozzle-sern-startup.md`](../../../plans/active/tooling-nozzle-sern-startup.md) §5.1。

| 項目 | 生産値 | 根拠 |
| --- | --- | --- |
| `gas.model` | `frozen_tp` | 排気 = CEA 凍結組成 EXH、外気 = AIR。熱力学は NASA-9 |
| `spec.wall_thermal` | 等温 1000 K | 断熱は M∞10 で回復温度 4600 K になり残差がプラトー (R4b) |
| `mesh.ic` | **`moc`** | 一様 IC は入口状態を全域に置き出口で 17 倍ずれる。MOC 場を外挿して埋める |
| 梯子 (暖機/soft/mid) | **各 1500 step, CFL 0.1** | 1000 は設計によって `rms_roY1` が上昇。3 設計で確認 |
| 本段 | **CFL 5.0, 500 step** | 上限は 6 (12 で発散)。頭打ちは cfl 6 で 500 step (0.001 % 以内) |
| `evaluate.implicit_relax` | **使わない** | 0.7 は定常解を 0.05 % 動かす。CFL を上げる方で速くする |
| `evaluate.p_min` | 20.0 | 既定 1.0 Pa は M∞10 の膨張で割れる |
| `evaluate.outlet_kind` | **`outflow`** (3D) | `statPress` の Ps は実出口圧より桁違いに低く、node の亜音速壁列から unstart |

**CFL の上限は段で 60 倍違う**。立ち上がりは 0.1、発達後は 6。梯子は外せない (一様 IC からの直接起動は cfl 0.1 でも発散)。
`opt.warm_lam_ramp: [0.1, 0.3, 1.0]` で暖機の CFL を段階昇圧できる (runner が forge を複数回起動)。
**soft と mid は `soft_cfl` を共有していて独立に振れない** (段ごとの上限は未測定)。

## 2. 収束の判定 (5 条件すべて)

`metrics/sern_gates.py` の `evaluate_gates` が返す `verdict` を**必ず貼る**。

1. 保存場に NaN/Inf 無し・正値
2. ~~**全残差列がプラトー**~~ **実装バグあり (2026-09-19)**: `check_convergence.py:108` は 3 桁低下した列に
   状態文字列を付けないので、**3 桁落ちてなお下降中の列が素通りする**。プラトーは収束の十分条件でもない
   (残差の大きさを問わないので float32 の更新消失やクランプで止まった状態と区別できない)
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
