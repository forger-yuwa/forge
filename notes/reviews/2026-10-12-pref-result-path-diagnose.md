# codex 諮問 (diagnose): pref-result-path

- **brief**: [`notes/reviews/briefs/2026-10-12-pref-result-path.md`](../../notes/reviews/briefs/2026-10-12-pref-result-path.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-12
- **commit**: `b1e4f000` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.7 min, rc=0
- **結論**: **次は T2 を「G0・上限 1・`pRef: 0` に固定した float／FP64 各 2 本の、登録停止規則までの比較」として実施する。**
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/rpc_judge.json`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/rcap_judge.json`, `plans/active/time_integration-line-implicit-speed.md`, `procedures/solver-settings.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical なし）

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | §6.55 の解釈 | **「判別不能」を採用。** JSON から再計算すると、結果 A の左辺／右辺は W1 が **6.405／9.138 K**、W2 が **13.801／13.632 K**。結果 B も成立しない。「登録窓で温度差を小さくした」は書けるが、「半減を実証」「単独では半減不可能」のどちらも不可。[判定値:116](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rpc_judge.json:116) |
| Major | S1「供給源＋増幅器」 | **作業仮説として採用、原因確定としては却下。** 末尾量子化の陰性は、初期化・途中の格納・右辺・左辺を分離しない。`pRef` 後の軌跡の差の減少率を、圧力演算が担う誤差の寄与率にはできない。これらを候補に残す。[float plan:1356](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1356) |
| Major | S2「FP64 の停滞と float の離れは同じ増幅器」 | **要再検証。** Gref・G0・G0c、停滞残差・精度間の温度差・下流の θ_r は別の観測である。近い場所で増えることは、同じモードの証明ではない。「方向別 `dt` の拡大との関連がある」までに限定する。[float plan:1322](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1322) |
| Major | D1 と T1〜T3 | **D1 を試験候補として採用し、停止規則を明示した T2 を優先。** T3 の 3,000 step では抑制と遅延を分けられない。上限は平均流と SST の両方の刻みを縮める。まず上限 1 の実用性を判定し、中間値の掃引は後にする。[設定:76](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:76)、[float plan:1366](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1366) |
| Major | D2「float の既定に `pRef`」 | **既定化は却下、個別の試験候補として保持。** 確認したのは一つの格子・起点・基準圧・短い窓での効果。`pRef` は運動量流束と hoop の圧力項を変更し、全演算の精度を改善するものではない。D1 の試験には混ぜない。[対流実装:618](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:618)、[軸対称実装:52](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:52) |
| Minor | 数字の表現 | **補正が必要。** W1 の「67〜77 % 減」は反復番号どうしを対応させた記述。全組合せなら **65.49〜77.42 %**。また上限 1 の θ_r の再実行幅は最大 **0.000183 ポイント**で、厳密なゼロではない。丸め表示と実測を分ける。[rpc の数値:118](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rpc_judge.json:118)、[rcap 判定記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rcap_judge.json) |

結論: **次は T2 を「G0・上限 1・`pRef: 0` に固定した float／FP64 各 2 本の、登録停止規則までの比較」として実施する。**

第 1 仮説: **上限 1 では、float が同じ上限の FP64 に対し、登録した停止時の値・全残差の許容範囲を満たせる。** 確度: 中  
根拠: `case/45.isobutane_m6_d155/run_0718_c1_d1`・`run_0722_c1_d2` と `run_0720_c1_f1`・`run_0724_c1_f2` の 20,000〜30,000 step では、θ_r の差と再実行幅を合わせた最大値は **0.000547 ポイント**。ただし保存済み判定は全本 **`OVERALL: CHECK FAILURES ABOVE`、θ_r・Q_w は `DRIFTING`**。支持するのは有限期間の小差までである。[float plan:1174](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1174)  
反証条件: 同じ停止規則で float が予算内に到達しない、または到達時の精度間差・全残差が事前登録した許容範囲を超える。

第 2 仮説: 上限 1 は離れを遅らせただけで、長く進めると許容差を超える。未確認。  
第 3 仮説: 平均流の改善と引き換えに SST の進行が律速となり、上限 1 が両精度とも実用的な停止条件へ届かない。Gref の ω の悪化が動機だが、G0 では未確認。

判別 A/B:

- **変える因子はビルド精度だけ。** G0・共通の Q32 起点・上限 1・`pRef: 0`・その他の実効設定を固定する。既存の上限 1 系列を継続利用する場合も、新しい `run_*` に再開し、入力・バイナリ・設定・再開直後の保存量を照合する。
- **停止規則は既存方針を尊重する。** 2,500 step ごとに出力し、`|欠損| ≤ 0.1 kg/s` と θ_r 各断面のドリフト `≤ 0.05 %／2 万 step` を 2 出力連続で満たしたところで停止。最大は共通起点から通算 20 万 step。「収束するまで」の無期限延長はしない。[既存の停止規則:490](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:490)
- **判定基準の提案**は同節の停止時比較を使う。各精度 2 本の全組合せで、θ_r 各断面の相対差 ≤ 0.05 %、符号つき Q_w の相対差 ≤ 0.1 %、末尾 2 万 step の各 `rms_*` 中央値と停止時の `2π·Σ|res_ro|` が FP64 の 2 倍以内。全残差、局所温度差、壁熱流束分布、δ_loc、到達 step も記録する。
- **A：両精度が到達し、全組合せでフラグなし** → 第 1 仮説をこの停止規則・格子・起点の範囲で支持する。**B：未到達またはフラグあり** → 「上限 1 だけでこの運用要件を満たせる」を退ける。欠損データは判定不能。
- `check_convergence` と各対象量の `check_quasisteady` の VERDICT はそのまま添える。A でも未収束なら、名称は **「停止規則での比較にフラグなし」**とし、「同じ定常解」とは呼ばない。

やらない方がよいこと: T1 の中間値掃引、T3 の短期無再発を根拠にした生産採用、根拠のない部分 double 化を先行させること。現行の上限 50・FP64 の停止値を真値に置くことも避ける。D1 が通った後の生産比較では、現行構成との**停止時の量・全残差・同じ停止規則までの総時間**を比べる。共有 GPU の経過時間だけで速度優位を決めない。

呼び出し側の前提への異議: **「`pRef` が供給源を減らし、上限 1 が増幅器を止めた」はまだ説明仮説である。** §6.55 では `Ep` が W1→W2 で **2.15〜2.20 倍**に増えている。8 本の記録は **`OVERALL: CHECK FAILURES ABOVE`**で、**`STEADY` は評価域の平均温度だけ**。局所差の定常化は示していない。現在書いてよい全体像は「圧力ゲージと反復設定の双方に精度差の発達が敏感で、末尾量子化だけでは登録窓の差を再現できなかった」までである。

不足情報: 対象 `run_*` の HDF5・CSV・個別 VERDICT はこの環境に無く、生データからの再検証はできなかった。今回は判定コード、保存済み JSON、計画記録を照合し、判定値を再計算した。**ファイル変更なし・plan 未反映**。上記は呼び出し側が事前登録するための提案である。
