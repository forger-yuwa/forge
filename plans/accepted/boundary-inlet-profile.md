# 入口分布プロファイル (inlet profile) — 非一様入口境界値の指定

## メタ

- **area**: `boundary`
- **status**: `done`
- **related_docs**: [`methods/boundary.md`](../../methods/boundary.md)
- **created**: `2026-06-28`
- **owner**: `CFD Dev`

## 1. 目的

入口境界 (`inlet_*`) の境界値を**一様でなく分布として指定**できるようにする。第一の用途は
**壁法則 (law of the wall) の境界層/チャネル速度プロファイル** $u(y)$ を入口に与えること。
さらに「時には x、時には y、時には xyz 全方向」で補間したいという要望に応え、補間方向を柔軟に選べるようにする。

## 2. 設計

inlet_* カーネルは **per-face の境界値 `bvar` (`Uxb[ib]`,`Uyb[ib]`,`Uzb[ib]`,`Ttb`,`Ptb`,...) をそのまま読む**。
これらは現状 `bcond::bcondInitVariables` で config の一様値 (`valueTypes==1`) に初期化される。
→ **`bvar` を face 重心座標に応じて非一様にセットすれば、kernel 無改修で分布入口が実現**できる。

### 2.1 有効化と入力
- `bcondConfig.yaml` の対象 inlet に `ints: {inletProfile: 1}` を付ける。
- run dir の `inlet_profile_<physID>.csv` を読む。**1 行目ヘッダで補間方向と量を指定**:
  - 先頭の連続する `x`/`y`/`z` 列 = **補間座標**。1 列 (例 `y`) なら **1D 線形補間**、3 列 (`x y z`) なら **3D 最近傍**。
  - 残り列 = **bvar 量名** (`Ux Uy Uz Tt Pt Ts Ps ro k omega` 等; その inlet が持つ bvar のみ反映)。
  - 例 (1D-y): `y Ux Uy Uz` / `1.0 0 0 0` / `1.02 30 0 0` ... 例 (3D): `x y z Ux Uy Uz` ...

### 2.2 適用箇所
`main.cpp` の `readBcondConfig` (bvar セット) の直後・最初の `applyBconds` より前に
[`applyInletProfiles(cfg, msh)`](../../solver_density_cuda/boundaryCond.cpp) を呼ぶ。`msh.planes[ip].centCoords`
(face 重心) と `bc.bvar`/`bc.bvar_d` が揃った後に適用し、補間値を host bvar にセット → device に再アップロード。
未指定 inlet は一様のまま (挙動不変)。`cfg.gpu==0` では host のみ。

### 2.3 補間
- **1D**: 補間軸でテーブルを昇順ソートし線形補間。範囲外は端値クランプ。
- **3D**: テーブル各行との 3D 距離最小 (最近傍)。

### 2.4 壁法則 helper
[`solver_density_cuda/tools/gen_inlet_walllaw.py`](../../solver_density_cuda/tools/gen_inlet_walllaw.py):
Reichardt 合成則 $u^+ = \frac{1}{\kappa}\ln(1+\kappa y^+) + 7.8[1-e^{-y^+/11}-\frac{y^+}{11}e^{-y^+/3}]$ で
チャネル (上下壁) / 片側 BL の $u(y)$ を生成し `inlet_profile_<physID>.csv` を書く。$u_\tau$ は中央/外縁速度 = `Uc`
になるよう二分法で求める。引数: `--physID --ylo --yhi --Uc --nu(=visc/ro) --mode channel|bl`。

## 3. スコープ

- **やる**: 任意方向 (x/y/z/xyz) のテーブル補間で inlet bvar を非一様化。壁法則 helper。
- **やらない**: 時間変化する入口 (turbulent inflow generator は別途 `inlet_fluctVelocity`)。
- ~~k/omega 等は bvar に存在すれば反映するが、RANS scalar 入口の per-face 化は別経路
  (`applyRansScalarBoundaries`) で未対応。~~ 決着 (2026-09-08, §8): k/ω・化学種 `Y{s}` とも per-face bvar を
  読む経路になっており、node の化学種だけ境界ノードにピンされていなかったのを修正した。

## 4. 検証

`case/18.backstep/run_inlet_profile_test`: backstep inlet (physID1, y∈[1,3]) を `inlet_uniformVelocity` +
`inletProfile:1`、`gen_inlet_walllaw.py --Uc 49 --nu 8.5e-4 --mode channel` で壁法則生成。
**ログ `[applyInletProfiles] ... set 3 quantities ... 297 faces`、入口 Ux(y) が壁法則分布** (壁 y=1,3 で 0、
中央 y=2 で増加; 30 step で node が profile に緩和: 中央 45.1 vs 目標 49, 形状一致)。一様なら全面 49 のところ
graded 分布になることを確認。

## 5. 影響範囲

- `solver_density_cuda/boundaryCond.cpp` (`applyInletProfiles` 追加), `boundaryCond.hpp` (宣言),
  `main.cpp` (呼び出し 1 行), `tools/gen_inlet_walllaw.py` (新規)。
- 既存挙動: `inletProfile` 未指定の inlet は一様のまま **不変**。

## 6. 完了条件

- [x] `applyInletProfiles` 実装 + main.cpp 呼び出し + helper
- [x] backstep で壁法則分布を検証 (run_inlet_profile_test)
- [x] methods/boundary.md に節追記、本 plan を `accepted/` へ

## 7. 残作業 (優先順)

| 優先 | 項目 | 状態 |
| --- | --- | --- |
| 低 | 2〜3 軸の分布は最近傍のみ (2D 線形/双線形補間は未実装)。必要になったら `applyInletProfiles` に追加 | 未着手 |
| 低 | node で入口ノードが周期境界の seam にも属する場合、`periodicNodeGather` (残差除外の後) が相手の残差を合算する。両側とも入口なら 0 同士で無害。入口∩周期のケースが出たら gather 後に再除外する | 未着手 (該当ケース無し) |
| ~~低~~ | ~~`speciesImplicitCoupling: 1` (scalar-DPLUR) では隣接の δ(ρY) が Jacobi sweep でピンノードに漏れる~~ | 決着 (2026-09-08, codex High 指摘): `species_dplur_sweep_d` でピン行を δ=0 に拘束 (方式 1/2 とも) |
| 低 | `inlet_Pressure_dir` の多成分対応 (`Yb` を受け取り混合則で状態を組む)。現状は単成分熱物性 | 未着手 (必要になったら `inlet_Pressure_d` と同形に) |
| 低 | 時間変化する入口分布 (CSV の時系列) は非対応。合成乱流は `inlet_fluctVelocity` | 未着手 |

## 8. 変更ログ

- 2026-06-28: 実装・backstep 検証・accepted。
- 2026-09-08: **全温 / 組成 (H2O) 分布の運用整備**。(1) 亜音速入口 `inlet_Pressure` の `Tt`/`Pt`/`Y{s}`/`k`/`omega`
  列、超音速入口 `inlet_uniformVelocity` の `ro/Ux/Uy/Uz/Ps` 列が per-face で効くことを case/16 (node/cell, TP
  MIXDRY+H2O) と case/46 (node CPG) で確認。(2) **node の化学種入口 Dirichlet が境界ノードに効いていない
  バグ**を発見・修正 (`species_dirichlet_boundary_d` node 分岐でノードをピン + `speciesPinResidual_d_wrapper`
  で残差除外; ghost 専用だったため一様値の変更も分布も場に入らなかった)。cell は不変 (run-to-run ノイズ内)。
  (3) 生成・照合ツール `tools/gen_inlet_profile.py` (式/測定表 → CSV、CPG/TP の Tt,M,Ps→ρ,U,Ps 換算、
  `verify` で h0 由来 T0・Y と照合)。(4) 手順書 `procedures/inlet-profile.md` と skill `forge-inlet-profile`。
  検証 run: `case/16.nozzle_wys/run_0317_inletprof_tt_h2o` (Tt +30 K, Y_H2O 山: T0 max 0.3 K 差, Y 1e-7),
  `run_0318_inletprof_h2o_cell`, `run_0319_speciespin_node_regression`, `case/46.sern_design/run_0089_inletprof_supersonic_tt`
  (ρ/U 1e-5 一致)。codex レビュー (read-only) の指摘 6 件に対応: 連成陰解法 sweep のピン行拘束 (High)、
  `inlet_Pressure_dir` 多成分非対応の明記、ツールの TP 既定種 (N2)・成分ごとの範囲検査・`verify` の Ps/Ts
  対応と全点 0 列のスケール、`applyInletProfiles` の applied/IGNORED 報告と Y 検査、手順書の誤記 (CSV 欠落は
  エラー終了・`set N` はヘッダ列数)。
