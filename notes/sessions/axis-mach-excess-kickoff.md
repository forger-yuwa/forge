# 依頼: 軸対称計算で軸上の M が +0.1〜0.2 % 高く出る原因の追及 (別セッション用の開始プロンプト)

以下をそのまま新しいセッションに貼る。正本 (残作業・観測値) は plan の §5.1 #17 で、ここは写しとポインタだけ。

---

軸対称ノズル (case/45.isobutane_m6_d155, M6) の計算で、試験部の軸上だけ M が設計値より +0.1〜0.2 % 高く出る。
オーバーシュート (局所的な行き過ぎ) とは別の「軸対称計算の軸付近の誤差」と見ている。数値誤差なら原因を突き止めたい。

- 作業ツリー: `/home/sano/work/forge-integ-1005` (ブランチ `feature/nozzle-wall-fit-and-pipeline`)。元のツリー `/home/sano/work/forge` は別ブランチなので、そちらで探さないこと。
- 正本: `plans/active/verification-m6-axis-wave-mesh-su2.md` の **§5.1 #17** (観測値・既知の関連・切り分け候補)。
  観測の出所: `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` §9 の 2026-10-06 「軸付近の格子依存」。
- 使える run (同じ壁 run_0092 で格子だけ違う): `case/45.isobutane_m6_d155/run_0105_ns_coarse_cfl1` (粗)・`run_0103_ns_finemesh_pass_cfl1` (細)・
  `run_0106_ns_finemesh3_pass_cfl1` (第三水準)、設計壁の Euler `run_0114_euler_pin_G1_recal_ext6k`。**AWS の `~/forge-wallfit/case/45.isobutane_m6_d155/` にある**
  (最終場と res_0 だけ残っている。ローカルには run_0114 などの一部のみ)。AWS は skill `forge-aws-run` (起動は毎回ユーザに確認)。
  バイナリは `~/forge-wallfit-bin/solver_density_cuda/build/forge` (共有の `~/forge-integ` は他セッションが作り直すので使わない)。
- 進め方: #17 は担当 F。**まず切り分けの計画 (どの比較で何を棄却するか・合否の数値) を plan に書き、上位の判断役に諮ってから**計算を回す
  (AGENTS.md「モデル分担とエスカレーション」; 判断役は `~/.config/forge/diagnose-backend` を確認)。
  候補: (a) 軸側間隔だけを細分した格子系列 (mesh2d の `axis_gap_frac` を `wall_first_frac_throat` と併用できるようにする必要あり)、
  (b) 同一格子・同一 BC の SU2 比較 (`procedures/su2-cross-check.md`)、(c) 軸上の値の取り出し方 (直読と偶関数外挿 `axis_extract.axis_curve_node(mode="evenfit")`) の差、(d) cell 離散化との比較。
- 既知の関連 (先に読む): memory `node-axis-value-position-study`・`node-axisym-axis-dirichlet`・`axisym-pref-gauge-hoop-fix`・`axisym-rweight-closure-fp32`、
  `plans/accepted/architecture-axisym-axis-singularity.md`、`plans/active/axisymmetric-su2-source-formulation.md`。
- ソルバの数値カーネル (`solver_density_cuda/cuda_forge/`) を変える前は必ず上位に諮る。
