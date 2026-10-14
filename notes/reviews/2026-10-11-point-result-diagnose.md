# codex 諮問 (diagnose): point-result

- **brief**: [`notes/reviews/briefs/2026-10-11-point-result.md`](../../notes/reviews/briefs/2026-10-11-point-result.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `13e85cc6` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.1 min, rc=0
- **結論**: **次は C1――`time.deltaT.lineDtDirectionalCap` だけを 50→1 にした、同時期の対照付き・30,000 step の精度比較を行う。**
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/rpt_judge_s2.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

| 重大度 | 諮問事項 | 採否・根拠・対案 |
|---|---|---|
| Major | 「point では同じ解へ向かう」 | **要再検証**。`rpt_judge_s2.json` では4本とも `OVERALL: CHECK FAILURES ABOVE`、100,000〜120,000 step の θ_r は `DRIFTING`。正式判定は「判別不能」を維持し、「登録した窓で精度差が小さい」と記す。 |
| Major | 「B0 の遅いモードを丸めの揺れが押す」 | **仮説としてのみ採用**。point 化はライン・方向別 dt・ITJ を同時に変更している（[plan:1112](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1112)）。Gref の局所更新の自己相関と発生場所だけでは、G0 の精度差が同じモードであるとも、駆動がランダムな揺れであるとも確定しない。まず上限だけを変える。 |
| Major | C1 を旧 B0 の結果だけと比較する | **却下**。無介入の再実行 F0 が旧対照の許容範囲を外れた実績がある（[plan:1091](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1091)）。C1 は採用するが、`Cap=50` の陽性対照を同じ実行条件で取り直す。 |
| Minor | 「120,000 step まで約1/50」 | **表現を修正して採用**。保存・評価した窓に限定する。30,500〜99,500 step の場は削除対象なので全期間の上限ではない（[plan:1112](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1112)）。W4 の符号付き平均精度差は −0.030〜−0.031 ポイントで、B0 の正方向の差とは符号も異なる。 |

結論: **次は C1――`time.deltaT.lineDtDirectionalCap` だけを 50→1 にした、同時期の対照付き・30,000 step の精度比較を行う。**

第 1 仮説: B0 で観測した30,000 step の大きな精度差には、方向別 dt の拡大が関与し、上限1でその期間の差を小さくできる。 **確度: 中**
  
根拠:
- `case/45.isobutane_m6_d155/run_0487_tr_fp64/`・`run_0488_tr_f32/`・`run_0489_tr_f32b/` の記録では、B0 の W2 の平均精度差は θ_r(40/70/94) で +1.950 / +1.525 / +1.285 ポイント（plan §6.27）。
- point の `case/45.isobutane_m6_d155/run_0713_pt_d1/`・`run_0714_pt_d2/`・`run_0715_pt_f1/`・`run_0716_pt_f2/` は、W4 の `max|ΔF−ΔD|+N` が 0.03358 / 0.03282 / 0.03190 ポイント。ただし上記のとおり未収束・`DRIFTING`。
- [setDT_d.cu:158](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:158)〜194 は、ライン面を含む CFL 最大値から上限を課す。粘性 CFL 割引が無効の条件では、上限1で全面基準の dt に戻せる。ライン結合と ITJ は残るため、point 化より変更を絞れる。

反証条件: 実際に dt が縮小したことを確認し、`Cap=50` の大差も再現したうえで、`Cap=1` でも両窓・3断面に大差が残ること。この場合、「上限1だけでこの期間の差を抑えられる」を退ける。

第 2 仮説: ライン結合・ITJ・精度依存の演算／初期化と反復との相互作用が残り、dt の上限だけでは不十分。**未確認**。

第 3 仮説: dt の縮小は精度差の発達を遅らせるだけで、長期の差を除去しない。**未確認**。point でも W3→W4 で精度差が増えており、有限期間の小差だけでは除外できない。

判別 A/B:
- **設定の差は上限50対1だけ**。各上限で FP64×2・float×2、計8本を共通の Q32 入力から30,000 step。`lineImplicit=1`、`lineDtDirectional=1`、ITJ=5、CFL・緩和・バイナリ・環境変数を固定する。同時実行数と起動順も両条件で揃える。
- W1=20,000〜25,000、W2=25,000〜30,000。Δ・尺度 D・再実行幅 N は §6.48 と同じ定義とし、N は上限ごとに求める。
- **前提**: `Cap=50` で両窓・3断面とも `min sign(D)(ΔF−ΔD)−N ≥ 0.5|D|` を再現する。入力・実バイナリのハッシュ、各精度内の初期状態、有限性、縮流部と壁際の実効 dt も照合する。
- **結果 A**: `Cap=1` で全6条件が `max|ΔF−ΔD|+N ≤ 0.1|D|` → 第1仮説を支持し、「この期間でも上限変更だけでは不十分」を退ける。**長期の解の一致や遅いモードの同定とはしない**。
- **結果 B**: `Cap=1` でも全6条件が `min|ΔF−ΔD|−N ≥ 0.5|D|` → 第1仮説を退ける。符号も記録し、逆符号の大差を改善扱いしない。
- 中間・混在・陽性対照不成立は判別不能。全残差の `check_convergence` と θ_r・Q_w の `check_quasisteady` を添え、未収束なら有限期間の応答として扱う。

やらない方がよいこと: point の延長や C2 を同時に始めること、30,000 step の小差を恒久的な抑制と呼ぶこと、この結果だけで double 廃止・既定値変更・opt-in 削除を決めること。

呼び出し側の前提への異議:
- **共通の Q32 入力と、初期化後の状態の一致は別**。JSON には精度間の `res_0` の最大相対差として `roe=5.52e−4` が記録されている。現在の比較は初期化の精度依存も含み、step 内の丸めだけを分離していない。
- 「右辺を変えない」は、観測された停止状態が右辺ゼロの同じ解だという証明ではない。
- 「効かなかった介入」を一括して原因候補の除外根拠にしない。少なくとも pRef の正式判定は前提不成立による判別不能である（plan §6.47）。

不足情報: 対象の実 run ディレクトリはこのワークスペースになく、HDF5・残差全列・実バイナリを独立再検査できなかった。今回は判定 JSON・判定コード・plan・`setDT` を照合したレビューであり、forge-25 の一次データも未確認。**ファイル変更なし、plan 未反映**。呼び出し側で §5.1 #22・§6.50 に観測と仮説を分けて反映すること。
