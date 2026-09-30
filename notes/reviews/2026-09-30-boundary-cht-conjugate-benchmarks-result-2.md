# codex レビュー: boundary-cht-conjugate-benchmarks (result)

- **plan**: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md)
- **stage**: `result` (diff base `feature/cht-axisym-graetz`)
- **date**: 2026-09-30
- **commit**: `64312448` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m2
- **extra**: `case/64.conjugate_pipe_wall/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

A の「FP64・node・流れ場固定」という限定結果は支持する。収束・準定常・界面収支と、代表 run の参照計算を再現できた。
ただし、A2 r16 の原因断定と安全停止の説明には根拠不足があり、判定条件・結論の文書同期も必要。現状のまま `accepted` へ移すことには反対する。

1. **Major — A2 r16 の FAIL を「観測差が許容を超えた」「参照側でなく forge の粗格子誤差」と断定している。**

   根拠：[plan §5.1 #7b](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:147)、[case README:17](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/README.md:17)、[現在仕様:350](/home/sano/work/forge-cht/methods/boundary.md:350)。

   実際の [A2 r16 の L3 再判定:4](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/ab_levels3/EVAL_CONJ_a2_r16_L3.txt:4) は、加熱区間の熱流束について、

   - 観測差：**1.8751 %**
   - 不確かさ：**0.22475 %**
   - 合計：**2.09985 % > 許容 2 %**

   である。登録規則による `FAIL` は正しいが、**観測差単独は許容内**。参照細分化後にも保守的な合格条件を満たさなかったことは、参照側の寄与を排除した証明にはならない。A1 r16 と同じ強さでは原因帰属できない。

   **対案：** A2 は「参照を追加細分しても不確かさ込みの合格条件を満たさない。forge の粗格子誤差が支配的であることを示唆する」と記す。`FAIL` は維持し、「差が許容を超える」「参照側でなく」という断定を削除する。この限定で閉じるなら追加 run は不要。

2. **Minor — A2 r16 の安全停止を「偽発火」とする根拠がない。**

   根拠：[現在仕様:350](/home/sano/work/forge-cht/methods/boundary.md:350) は偽発火と記す。一方、[停止条件](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:1029) は更新量の連続増加と累積増幅を検出するもので、報告された `6.1e-4 → 1.3e-3 K` はその条件に整合する。

   `case/64.conjugate_pipe_wall/run_0008_a2_r16/CONVERGENCE_CHECK.txt` も、step 32300 で **`NOT CONVERGED`**。別の `Df_scale` にした run の収束や、不動点の理論上の不変性だけでは、元の反復がその後自然に収束したとは示せない。

   **対案：** 「`Df_scale=5` で更新量増大を検知して安全停止。`Df_scale=20` の別 run で収束を確認。同一格子の `Df_scale=5` 収束対照はない」とする。

3. **Minor — 最新の結論・事後改訂・完了記録が文書間で一致していない。**

   根拠は次のとおり。

   - [plan §1](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:29) は A2 r64 を判定不能としているが、最新の `run_0010_a2_r64/EVAL_CONJ_L2.txt` は全項目 `PASS`。`Q_up/Q_tot` は差 `4.6411e-4 + U 2.4454e-4`。
   - [README:3](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/README.md:3) の「判定基準は plan の登録文そのまま」は、[§4.6 の事後改訂](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:111)と整合しない。数値許容幅は同じでも、`U` の構成は変更されている。
   - [§5.1 #3c](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:140) に完了済みの準定常評価が残作業として残り、[§6.1](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:159) には既実施の result レビューが未記録。[計画索引:28](/home/sano/work/forge-cht/plans/README.md:28) も分割前の説明である。

   **対案：** 最新結論を「事後改訂した同一有限領域の比較条件で、A1/A2 の r32・r64 が許容内」に統一する。旧判定は履歴として区別し、残作業表・レビュー記録・変更ログ・索引を同期する。C の未完了事項は後継 plan に残す。

検証で確認できた範囲は以下のとおり。run の所在はすべて `case/64.conjugate_pipe_wall/` 配下で、恒久索引は [README の計算 run 一覧](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/README.md:9)。

| run | 最新 `EVAL_CONJ_L2.txt` | 再確認した前提ゲート |
|---|---|---|
| `run_0005_a1_r16` | FAIL（一部判定不能） | 全 PASS |
| `run_0006_a1_r32` | PASS | 全 PASS |
| `run_0007_a1_r64` | PASS | 全 PASS |
| `run_0009_a2_r32` | PASS | 全 PASS |
| `run_0010_a2_r64` | PASS | 全 PASS |
| `run_0011_a2_r16_df20` | FAIL（一部判定不能） | 全 PASS |

前提ゲートは、単一区間の全保存量について `check_convergence.py` の **`PASS (converged)`**、保存された50時点の全系列について公式準定常判定関数の **`ALL STEADY`**、G-if・G-cons の **`VERDICT: PASS`** を再現した。圧縮 CSV の読み込みと結果保存の抑止だけをメモリ上で適応し、判定ロジックは変更していない。最終場の `VALUE/*` に NaN/Inf はなく、3格子の品質は **`VERDICT: PASS`**、最大 AR 31.3・最大 skewness 0。

A1 r32 は参照計算自体も再実行し、掲載値と `VERDICT: PASS` を再現した。参照解の自己検査も PASS。一次データ台帳は手元の **249件すべて一致**した。ただし、AWS のみにある中間場からの系列再生成は今回確認していない。

指定 diff にソルバ本体・`design` の変更はない。参照解の軸面積、界面熱流束の符号、`per rad` の積分に新たな明白な誤りは見つからなかった。共通ツールの `--abs-scale` 未指定時は、旧版との100系列比較で差がなく、不正尺度も拒否した。これは float32・cell・周期 seam の精度検証を追加したことにはならない。

C は保留が妥当である。例えば `case/65.conjugate_flat_plate/run_0007_c1_n64/` は、保存記録で流体収束・G-if とも **`NOT CONVERGED`**、界面残差最大 **677.8 W/m²**。後継 plan がこれを未解決として保持している点は適切。

**推奨は、上記3点を番号順に修正し、A の限定結果として `accepted` へ移すこと。** 事後改訂を事前登録の総合 PASS と扱わず、C は後継 plan で継続する。ファイルは変更しておらず、本レビューの指摘は plan 未反映。

指摘数: Critical 0 / Major 1 / Minor 2
