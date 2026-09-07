---
name: forge-inlet-profile
description: forge の入口に分布 (全温 Tt / 全圧 Pt / 組成 H2O 等 / k・ω / 超音速入口の ρ,U,Ps / 入口境界層) を与えるときの手順。inletProfile CSV を gen_inlet_profile.py で生成し、起動ログと verify で場に入ったことを確認する。「入口に分布を与えたい」「全温分布」「H2O 分布」「入口プロファイル」「入口境界層を与えたい」と言われたときに使う。
---

# /forge-inlet-profile — 入口分布を与える

正本は [`procedures/inlet-profile.md`](../../../procedures/inlet-profile.md) (仕様は `methods/boundary.md`
「入口分布プロファイル」)。手順:

1. **入口種別を確認する** (`bcondConfig.yaml` の `kind`)。超音速ノズルでも入口は亜音速なら `inlet_Pressure`。
   - `inlet_Pressure` (亜音速): 分布にできる列は `Tt Pt Y{s} k omega`。**Tt・H2O 分布はそのまま列に書ける**。
   - `inlet_uniformVelocity` (超音速・全量固定): 列は `ro Ux Uy Uz Ps Y{s} k omega`。**`Tt` 列は無効**なので
     `gen --Tt --M --Ps|--Pt` で ρ/U/Ps に換算する。
   - `inlet_Pressure_dir` は多成分非対応 (単成分熱物性)。凝縮モーメントは分布不可。時間変化する分布も不可。
2. **run を複製**して `run_NNNN_<slug>` を作り、入口に `ints: {inletProfile: 1}` を付ける。`floats` の一様値は
   残す (CSV に書かない量と換算のデフォルト)。段階起動で run を複製するたびに **CSV も一緒にコピー**する
   (restart でも毎回読む。`inletProfile: 1` のまま CSV が無いと起動時エラー終了)。
3. **CSV を作る**: `python3 solver_density_cuda/tools/gen_inlet_profile.py gen --run RUN --physID N --axis y --range lo hi
   --Tt "式" --Y H2O="式" [--plot]` (測定値は `--table file.csv`)。化学種は名前で指定し、残りは他種が受け持つ。
   gen が出す各列の min/max と bcond 一様値で桁を確認する。座標はメッシュの単位 (m)。
4. **起動ログを確認**: `forge_run.log` の `[applyInletProfiles] physID=N ... applied: Tt Y0 Y1` に意図した列が
   全部あるか (`IGNORED` に出た列はその種別に無い列名)。行自体が無ければ `inletProfile: 1` 漏れ。CSV 欠落や
   Y の範囲/和の不正は起動時エラー終了。
5. **結果を照合**: `gen_inlet_profile.py verify --run RUN --physID N --plot`。Tt は `VALUE/h0` 由来の T0 と比較
   (`output.level` 1 以上)。cell は第 1 セル値なので境界層内は過渡で差が残る。node は境界ノード値。
6. 通常どおり `check_convergence.py` / `check_quasisteady.py` の VERDICT、run パス、case README の run 一覧。
   全温を上げると同じ Pt で質量流量が変わるので設計比較では ṁ を併記する。

既知の落とし穴: 2026-09-08 より前のバイナリは node で化学種の入口 Dirichlet が境界ノードに効かない (Y が一様の
まま)。古いバイナリの run は再実行する。
