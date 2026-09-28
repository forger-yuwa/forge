# 63. 軸対称 CHT の流れあり検証 — 管内層流の Graetz 問題 (壁温一様)

計画: [`plans/active/boundary-cht-axisymmetric-graetz.md`](../../plans/active/boundary-cht-axisymmetric-graetz.md)。
**判定基準・パラメータは plan §6 の登録文そのまま**で、本 README は入力の作り方と run の所在だけを書く。

## 幾何と条件 (登録値、plan §4.1)

```
 r=R ┌─ wall_up (3, 断熱) ─┬──── wall_heat (4, 共役 / 等温) ────┬─ wall_down (5, 断熱) ─┐
     │ inlet (1)          │   固体殻 R..R+0.5 mm (k_s 100)       │                       │ outlet (2)
 r=0 └────────────────────┴──────────── axis (6) ───────────────┴───────────────────────┘
    x=−10 mm            x=0                                   x=172.8 mm           x=182.8 mm
```

| 量 | 値 |
| --- | --- |
| 流体 | 定数物性 (`viscMethod: 0`・`thermCondMethod: 0`)、μ 4.085818e-5、k_f 5.700284e-2、c_p 1004.5、γ 1.4 |
| 流れ | R 1 mm、M 0.05 (U_m 17.359 m/s)、Re_D 1000、Pr 0.72、Pe 720。入口 `inlet_uniformVelocity` + 放物 `inlet_profile_1.csv` |
| 熱 | T_in 300 K、加熱区間 (x⁺ 0 … 0.12) の壁: 共役 (固体殻外面 Robin h 1e8・T_c = T_in + ΔT) または非連成の等温 T_in + ΔT |
| IC | Poiseuille (放物速度・線形圧力・300 K)。`make_run.py` が VALUE にパッチ |
| 基準 | `graetz_ref.py` (古典 Graetz、forge と独立)。`--selftest` PASS |

## ツール

| ファイル | 役割 |
| --- | --- |
| `graetz_common.py` | 登録値 (物性・寸法・physID)。case/61 の `axcht.py` を流用 |
| `graetz_ref.py` | 基準解 (march = 古典 Graetz、ellip = 有限 Pe)。`--selftest` |
| `gen_mesh.py` | 流体メッシュ (N_r 16/32/64、軸方向も入れ子で細分)。`mesh/graetz_r<N>.h5` |
| `check_dry.py` | 乾式 1 step の壁ダンプで境界の帰属を検査 (壁が r=R にあるか) |
| `gen_solid.py` | 固体殻 (`--ns` 層数、`--dT`)。`mesh/solid_<乾式 run>_s<層>_dT<ΔT>.h5` |
| `make_run.py` | run ディレクトリ (`dry` / `iso` / `cht`、`--dT`、`--cfl`)。固体の ROBIN/TC と ΔT を照合 |
| `eval_graetz.py` | `snap` (V-g1・V-g2) / `series` (準定常・末尾変動・反復の不確かさ) |
| `compare_nu.py` | `pair` (V-g4・V-g5・固体層) / `grid` (V-g3 格子間差) |
| `test_eval_graetz.py` | 評価器の負例試験 (合成データ、10/10 PASS: `test_eval_graetz.out`) |

## 生成手順

```bash
cd case/63.graetz_cht
for n in 16 32 64; do python3 gen_mesh.py --nr $n; done
python3 make_run.py dry run_0001_dry_r16 --mesh mesh/graetz_r16.h5      # 1 step 回して check_dry.py
python3 gen_solid.py --wall run_0001_dry_r16/res_wall_heat_4_1.h5 --dT 10
python3 make_run.py cht run_NNNN_... --mesh mesh/graetz_r16.h5 --solid mesh/solid_run_0001_dry_r16_s8_dT10.h5 --dT 10
```

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry_r16` / `run_0007_dry_r32` / `run_0008_dry_r64` | 乾式 1 step (非連成・加熱壁 300 K)。壁ダンプ座標 = 固体界面節点の出所 (ローカル float) | `check_dry.py` PASS (3 本) | ref (固体の入力) |
| `run_0002_smoke_iso_r16` | 起動確認 (ローカル float、等温 ΔT 10、cfl 5、200 step) | step 200 で P 負・\|u_r\| 119 m/s に崩壊 | 破棄予定 (登録外) |
| `run_0003_smoke_iso0_r16` | 同 ΔT 0 | 同様に崩壊 (P −650 kPa) → 熱の段差は主因でない | 破棄予定 (登録外。起動 A/B の根拠) |
| `run_0004_smoke_walls_{adiab,iso}_r16` | 切り分け: 壁を全部断熱 / 全部等温 (cfl 5) | 断熱は安定 (\|u_r\| 1.3e-3)、等温は崩壊 | 破棄予定 (登録外。同上) |
| `run_0005_smoke_iso0_cfl{2,1}_r16` | 同 ΔT 0、cfl 2 / 1 (400 step) | どちらも安定 | 破棄予定 (登録外。同上) |
| `run_0006_smoke_iso10_cfl2_r16` | 評価器の通し確認 (ΔT 10、cfl 2、400 step、未収束) | 評価器が端から端まで動くこと | 破棄予定 (登録外) |
| `run_0009_ab_cfl2_r64_cht10` / `run_0010_ab_cfl1_r64_cht10` | **起動 A/B** (plan §6 起動設定): AWS FP64・N_r=64・共役 ΔT 10・`cfl_pseudo` 2 / 1、10000 step | 両方破綻せず (ρ・P・T 物理的、界面残差単調減少) → **CFL 2 を採用** (plan §5.1 #5b) | ref (起動設定の根拠) |
| `run_0011`〜`0013` (G1 等温 ΔT 0/5/10, N_r 32) / `run_0014`〜`0016` (G2 共役 8 層 ΔT 0/5/10, N_r 32) / `run_0017`・`0018` (共役 16 層 ΔT 0/10) / `run_0019`・`0020` (N_r 16 ΔT 0/10) / `run_0021`・`0022` (N_r 64 ΔT 0/10) | **本番** (plan §6)。AWS FP64・600000 step・`cfl_pseudo` 2。各 run に `CONVERGENCE_CHECK.txt`・`EVAL_SNAP.txt`・`EVAL_SERIES.txt`・`graetz_nu_600000.csv`・`CHT_*_VERDICT.txt`、比較は `CMP_Vg{3,4,5}.txt`・`CMP_solid.txt` | 収束 12 本 PASS、V-g2 PASS (0.77 %)、V-g4/g5/固体層 PASS、G-if・G-cons PASS、**V-g1 FAIL 2 件・V-g3 保留 (N_r 64 未定常)** — plan §5.1 #6 | active (h5 は AWS `~/forge-graetz/case/63.graetz_cht/`) |
| `run_0023_g2_dT0_r64_ext` / `run_0024_g2_dT10_r64_ext` | N_r 64 の延長 (plan §5.1 #6e): `run_0021`/`0022` の step 600000 から再開 (流体ビット一致・壁温・固体状態)、+600000 step | ALL STEADY (u_it 8.2e-8)、G-if/G-cons PASS、V-g3 PASS (`CMP_Vg3_ext.txt`) | active |

**AWS 上の保存 (2026-09-29 に削減)**: `~/forge-graetz/case/63.graetz_cht/run_*` は各 run の最終ステップの場・壁・固体ダンプと `conjugate_state_4.h5` だけを残し、中間スナップショットと `conjugate_iface_log_4.csv` は削除、`residual_history.csv` は gzip。判定の出力 (`*.txt`・`graetz_*.csv`・`momentum_balance_series.csv`) は手元に回収済み。
