# 諮問ブリーフ: chem のブランチの node 低マッハ前処理の修正 (`876ec8b4`) をこのブランチへ移植する計画 (2026-10-12)

- 依頼者: 主セッション (エスカレーション条件 6: `cuda_forge/` の数値の振る舞いを変える編集の前、条件 4: plan に無い修正の前)
- plan: `plans/active/architecture-float-state-double-geometry.md` (§5.1 #25・#26、§6.57 と結果)
- ブランチ: このブランチ `feature/nozzle-wall-fit-and-pipeline` (HEAD c11a1cc3)。元の修正は `origin/feature/chemistry-finite-rate` の `876ec8b4` (2026-09-05、case/48 Cabra で発見・修正)。
- ユーザの方針: 他のブランチのソルバの修正は必ず取り込み、A/B で要否を決めない (2026-10-11、メモリ「always-incorporate-solver-fixes」)。ユーザの目標: case/45 の NS を速くしたい。低マッハで収束が悪いという見立てを `lowMachPrecond 3` (左辺だけの前処理) で確かめたい。

## 0. 問い
1. 下の移植の範囲・やり方・回帰の確かめ方に穴はないか。
2. TP の出口の修正 (下の (b)) は case/45 の生産の値を変えうる。取り込みの後に何を測って記録すべきか。

## 1. 観測事実
- §6.57 (`run_0756〜0759`、判定 `_band_ab/cold_pair/rlg_judge.json`): G0・point・FP64 (`~/forge-fgeom11-fp64`、commit b0303547)・`precondEps 1.0` (β = 1) で `lowMachPrecond` 0 対 3 の 1 step。P3 では入口・出口の境界の節点の更新が P0 の 0.02〜1 %、壁の ρ の更新が 20 % 違い、拘束の成分 (軸の ρv、壁の ρu・ρv・等温壁の ρE) が 0 でない (最大 1.4e-4・1.7e-3・2.2e-3・2.69)。内部も 0.03〜1.1 % 違う。結果 A。
- `876ec8b4` のソルバの部分 (git show で読める): `timeIntegration_d.cu` 57 行 (前処理の block カーネル `implicit_defect_correction_block_precond_d` に、標準の block カーネルと同じ node の行の処理 = 境界の半割面の粘性の対角を除く `isNode`、軸の ρv の行、壁 3 行・等温壁の ρE の行の切り離し。rank-1 項は g[row] = 0)、`boundaryCond_d.cu` 27 行 (`outlet_statPress_d` の亜音速の特性の構成で TP の γ_mix とコンフィグの γ の取り違えを直す。超音速・CPG はビット不変)、`setDT_d.cu` 35 行・`speciesTransport_d.cu` 4 行・`solverConfig.{cpp,hpp}` 10 行・`variables.hpp` 4 行 (`time.deltaT.speciesPrecondDt`、既定 1: 前処理のとき化学種は拡大前の Δτ で進める、`dt_local_sp`)。
- このブランチへの cherry-pick は `timeIntegration_d.cu`・`speciesTransport_d.cu`・`solverConfig.{cpp,hpp}` で衝突する (ブランチが分かれてから両方が大きく変わった)。case/48 Cabra のファイル・chem の plan も同じ commit に入っている (移植しない)。
- chem のブランチの記録 (メモリ): 修正後は cell と前処理 0・1 はビット不変。TP の亜音速の出口の修正で BK (case/47) は残差が 0.5 % 変わった (再現の揺れ 1e-5)。

## 2. 移植の案
- ソルバの部分だけを手で移植する (docs の該当節も)。前処理の block カーネルの引数に node の印 (`isNode`・壁/軸/等温壁の印) を足し、標準の block カーネルと同じ行の処理にそろえる。chem のブランチの実装をそのまま写し、このブランチで後から入った変更 (hoop の r 重み、`implicitSolvePrecision`、`qAccumulatorFP64`、診断のフックなど) と整合させる。
- 回帰:
  - R1: `lowMachPrecond 0` の経路のビット不変 (case/45 G0 の B0 とpoint、各 10 step、FP64 と float、移植の前後のバイナリで res の全量を比べる)。ただし TP の出口の修正で case/45 の出口の壁際の亜音速の節点は変わりうる → その節点の外はビット一致、変わった節点の数・場所・大きさを記録。
  - R2: §6.57 を移植後のバイナリで回し直して結果 B (境界・軸・壁・内部の d が P0 と 1 % 以内か再実行の幅の中、拘束の違反なし)。
  - R3: CPG のケース (例 case/08 の bump か case/48 の平板) でビット不変。
- その後に試験 P の本番 (`lowMachPrecond 3` で収束の速さ) を登録し直す。

## 3. 読んでよいファイル
- `git show 876ec8b4` (ソルバの部分)、このブランチの `solver_density_cuda/cuda_forge/timeIntegration_d.cu` (`implicit_defect_correction_block_precond_d` と標準の block カーネル)、`boundaryCond_d.cu` の `outlet_statPress_d`、`setDT_d.cu`
- float plan §5.1 #25・#26、§6.57 と結果
- `plans/accepted/time_integration-lowmach-preconditioning.md`
