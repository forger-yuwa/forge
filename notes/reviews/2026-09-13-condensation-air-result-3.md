# codex レビュー: condensation-air (result)

- **plan**: [`plans/active/condensation-air.md`](../../plans/active/condensation-air.md)
- **stage**: `result` (diff base `feature/sern-design`)
- **date**: 2026-09-13
- **commit**: `79131b2d` (feature/condensation-air)
- **codex**: effort `high`, 4.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m2
- **focus**: 2026-09-13 result レビュー 2 回目 (NO-GO, M1 非有限入力の反転成功化 / M2 U_x の 2 倍床超過と判定方法 / M3 面 R_eff 床 / m4 文書 / m5 §5.1) の採用・反映 (§6.1, §9 2026-09-13 ②) が実装 diff・単体試験 (test_cond_air (i)(j))・3 反復ノイズ床 (run_0027/0030/0033) と最終バイナリ run_0035・読込試験 run_0034 で裏付けられているか。accepted にしてよいか
- **extra**: `case/34.arthur_n2_nozzle/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

前回 Major 3 件の修正と、3 反復ノイズ床・最終 run による再判定は確認できました。新たな Critical / Major は見つかりません。
ただし、回帰基準の適用範囲と文書の不整合を直してから `accepted` に移してください。

1. **Minor — 「3 反復以上」の新基準を、回帰全体には適用できていない。**

   **根拠:** [plan:178](/home/sano/work/forge-cond/plans/active/condensation-air.md:178) は同一バイナリ・同一設定の **3 run 以上**を要求します。しかし dry cell の床は依然 `run_0015o_dry_slip_cell_oldbin` と `run_0031_dry_slip_cell_oldbin_rep` の **2 run**、node の「ノイズ水準」は旧・新バイナリの比較です。[README:63](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/README.md:63)、[README:70](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/README.md:70)

   dry cell の `diff_res.py --factor 2` は再実行でも **PASS、exit 0**でした。したがって退行を示す指摘ではありませんが、「新基準で回帰完了」とはまだ書けません。

   **対案:** 移動前に dry cell/node の同一バイナリ反復を補い、変数別床と exit code を記録してください。空気の床を N2 比較にも使う点についても、適用根拠を §6 に明記してください。

2. **Minor — 文書同期は改善したが、現在仕様・結論に古い説明が残っている。**

   **根拠:**

   - [methods/condensation.md:757](/home/sano/work/forge-cond/methods/condensation.md:757) は、CPG 二相音速を適用しない理由を依然「潜熱フィットの整合修正が先」と説明しています。今回その修正は完了しており、残る理由は一般 EOS 固有系・Jacobian 側です。
   - [plan:234](/home/sano/work/forge-cond/plans/active/condensation-air.md:234) は Arthur 偏差を「N2 モデルの較正問題」と断定しています。原因未切り分けとした §5.1 #10 と矛盾します。
   - [methods/condensation.md:647](/home/sano/work/forge-cond/methods/condensation.md:647) の「原始量・ρe とも前ステップ値を保持」は広すぎます。[dependentVariables_d.cu:76](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/dependentVariables_d.cu:76) では速度を先に更新しており、失敗分岐は更新済み保存量を前ステップへ戻していません。保持するのは `T/P/sonic/Ht`、`roe` は反転による上書きをしない、という動作です。

   **対案:** 上記を実装・未確定事項に合わせて訂正してください。境界 kernel 試験・流束収支は [§5.1:146](/home/sano/work/forge-cond/plans/active/condensation-air.md:146) と [後続 plan:51](/home/sano/work/forge-cond/plans/active/condensation-followups.md:51) に戻っていることを確認しましたが、**未実施の後続課題**という扱いを維持し、「全件反映済み」を検証完了の意味で使わないでください。

前回指摘の裏付けは次のとおりです。

| 対象 | 今回の確認結果 |
|---|---|
| M1 非有限入力 | [反転入口:105](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:105) と成功条件に有限性検査あり。`test_cond_air` (i) の7例すべて `ok=false`・有限温度、exit 0 |
| M3 面 `R_eff` 床 | [面関数:77](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:77) に床なし。(j) の `R_eff<1` 6状態を含む掃引で最大相対誤差 **2.19e−16**、PASS |
| M2 ノイズ床 | `run_0027/0030/0033` の設定・初期場は同一。全ペア最大から JSON **全26変数を完全再現**。`Ux` 床 **1.3394476763948262e−4** |
| 場差の再判定 | `run_0024→0027`、`run_0017→0029`、`run_0027→0035` はすべて **PASS、exit 0** |
| 省略時既定 | [run_0034 のログ:86](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/run_0034_cfg_default_check/forge_run.log:86) で `condKantrowitz=0`、`condKantrowitzGammaMode=0` を確認 |

最終成果物は `case/34.arthur_n2_nozzle/run_0035_air_cpgcarrier_final/res_12000.h5`、索引は [case README の run 一覧](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/README.md:35) です。再解析は onset **2.207 in / 38.94 K**、壁圧比 **1.2135 / 1.3514 / 1.4547**。判定は以下でした。

- `check_convergence.py`: **NOT CONVERGED (stalled/plateau)**
- `onset_analysis.py --series`: **STEADY**、保存時刻の出口最小 `u_n/c` **5.498**
- `check_quasisteady.py --quantity machmax,pmax`: **STEADY**。これは onset・壁圧比の判定を代替しません。

対象24 run・300保存場の `VALUE/*` に NaN/Inf は0、使用メッシュは cell/node とも **VERDICT: PASS**。結論は**未収束の準定常比較**として支持します。

単体は既存実行ファイルを再実行し、`test_cond_air` と既存 host 3本が PASS。再ビルドは未実施です。`test_kwc` の device 部分は CUDA 初期化エラーで再確認できず、今回の PASS 実績には含めません。

**推奨は、上記1の回帰根拠補完 → 2の文書訂正 → §6.1への今回の記録・採否反映の順に完了し、`accepted` へ移すことです。** ファイルは変更しておらず、今回の指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 0 / Minor 2
