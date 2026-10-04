# codex 諮問 (diagnose): d7275-height-ab-result

- **brief**: [`notes/reviews/briefs/2026-09-27-d7275-height-ab-result.md`](../../notes/reviews/briefs/2026-09-27-d7275-height-ab-result.md)
- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **date**: 2026-09-27
- **commit**: `c6370129` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.0 min, rc=0
- **結論**: **`run_0014_t26_h079` の実効設定とバイナリを確認し、`physProp.thermoFloat: 1 / 0` だけを変える継続 A/B で、残存する混合精度経路の寄与を調べる。**
- **extra**: `case/60.flatplate_d7275_m7/README.md`, `case/60.flatplate_d7275_m7/acceptance.json`, `case/60.flatplate_d7275_m7/tools/height_ab.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 判断 | 採否 | 理由・対案 |
|---|---|---|
| 高さ A/B は上端近接説を支持 | **限定して採用** | 記録された流れ残差の低下と圧力振幅の減少は支持する。比較域の残差増加は絶対量が小さく、これを覆さない。ただし局所化学種残差と指定された ρ・k・ω 時系列が不足し、事前登録の全項目を満たしたとは確定できない |
| (a) 一様初期値から段階起動し直す | **却下：現時点の次手として** | 初期ピークを大きくして低下桁数を増やしても、残る床の原因を識別できない |
| (b) 同じ設定で延長する | **却下：現時点の次手として** | 記録上は横ばい。先に、実装上残っている低精度経路を切り分ける |
| (c) 第一内部列の ω 残差を調べる | **採用** | 下記の `thermoFloat` A/B に絞る。G15 は未達のまま維持する |

結論: **`run_0014_t26_h079` の実効設定とバイナリを確認し、`physProp.thermoFloat: 1 / 0` だけを変える継続 A/B で、残存する混合精度経路の寄与を調べる。**

第 1 仮説: **「全域 FP64」ビルドでも残る `thermoFloat=1` の組成・熱力学評価の量子化が、自由流域の残差床と近壁 ω 残差を維持している。** 確度: **中。対象バイナリでの経路の有効性と、残差への因果関係は未確認。**

根拠:

- D-7275 の生成器が使う設定は `thermalMethod: 2`、`thermoHrefTemp: 298.15` で、`thermoFloat` を指定していない。現行コードの既定値は **1**。したがって、この生成設定ではハイブリッド経路が有効になる。[make_case.py:32](/home/sano/work/forge/case/56.gap_tp1187/tools/make_case.py:32)、[solverConfig.hpp:585](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:585)
- この経路は `flow_float` とは別に、**明示的な `float` で密度逆数・組成・組成和を計算し、その結果を double に昇格する**。さらに比熱比も float 除算で評価する。保存量配列を double にしても、この丸めは残る。[dependentVariables_d.cu:89](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:89)、[同:185](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:185)
- 読み取り専用の合成反例でも、真の組成を固定して密度だけ微小変化させると、この演算経路の組成に最大 **1.49e−8** の幅が生じた。double 経路は最大 **2.78e−17**。これは演算上の反例であり、対象 run の床を説明した実測ではない。
- 壁 ω は ν/y² に依存し、内部の消滅項は βρω²。小さな物性・状態変動が第一内部列の残差へ現れる経路は存在する。ただし寄与の大きさは未測定。[ransBoundary_d.cu:63](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransBoundary_d.cu:63)、[ransSource_d.cu:231](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:231)

反証条件: **実バイナリで既に `thermoFloat=0` だった場合、この仮説は対象外。1→0 が実際に切り替わり、切替過渡が落ち着いても対象残差と局所振幅が対照の 0.8–1.2 倍に残るなら、「この経路が床の主因」を棄却する。**

第 2 仮説: **近壁 SST の分離型 point-implicit 更新が反復振動を維持している。** 確度: 低。第一内部節点は実際に更新する DOF であり、壁ピンの残差として除外できない。以前の緩和 A/B は別の領域構成で判別保留だったため、この候補は除外されていない。[update_d.cu:380](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:380)

判別 A/B:

- **共通条件**: `case/60.flatplate_d7275_m7/run_0014_t26_h079` の最終保存場から、新規 run を二つ作る。両腕とも `restart_field.py --keep-src-dtype` で全保存量を同一にし、メッシュ・BC・物性 DB・SST・CFL・緩和・内反復・バイナリを固定する。
- **変更する一点**: A=`physProp.thermoFloat: 1`、B=`0`。既定値に頼らず両腕に明記する。
- **長さ**: まず各 5,000 step。3–4k と 4–5k の対象指標が 5% 超動くなら各 10,000 step まで延長し、それでも動けば判別保留。
- **見る量**: 全残差列、自由流域と第一内部列の成分別 √Σres²、固定節点の高精度な ρ・P・k・ω 時系列。両腕で `FORGE_OMEGA_BUDGET=1` と必要な `extraFields` を同じにし、`omg_trans + omg_prod − omg_dest + omg_cross` の収支も確認する。項別診断は既存実装にある。[variables.cpp:335](/home/sano/work/forge/solver_density_cuda/variables.cpp:335)、[ransSource_d.cu:238](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:238)
- **事前の読み**: B の対象残差が A の **1/10 以下**となり局所振幅も減少すれば、第1仮説を支持する。**0.8–1.2 倍に残れば主因説を棄却**する。流れだけ改善して ω が残れば、改善の帰属を流れに限定する。中間結果は保留。いずれも診断条件であり、G15 の代替合格条件にはしない。

やらない方がよいこと:

**Major — 低下桁数を稼ぐ目的の再初期化、別方程式の段の連結、未収束 `run_0012` の参照床への採用。** 判定器は「開始値」ではなく、**判定系列のピーク／末尾中央値**を使う。開始履歴を変えると同じ床でも判定が変わり得るため、床の解消とは区別する。対案は上記 A/B で残差そのものを比較する。[check_convergence.py:182](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:182)

呼び出し側の前提への異議:

- **Major — 「FP64 なので丸め床ではない」は受け入れない。** 明示的な float 経路が残る。対案は精度を配列型だけでなく実行経路で監査する。
- **Major — 高さ A/B の「支持」と「全項目確認済み」は分ける。** `height_ab.py` の領域別集計は最終場一枚で、共通して存在する残差フィールドだけを処理する。記録の旧上端域には化学種残差がなく、ρ・k・ω 時系列も未取得。対案は**「流れ残差・圧力振幅について支持、未取得項目あり」**と記録する。[height_ab.py:73](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/height_ab.py:73)、[acceptance.json:182](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:182)
- **Major — Cary の ω 残差 0.3 は、今回の 0.6 の許容根拠にならない。** ピーク基準なので PASS／NOT CONVERGED の違いは説明できるが、残差の絶対値は無次元化されておらず、格子・体積・物性が異なる。Cary の生成設定は `thermalMethod: 0` で、今回の TP 経路も通らない。対案は各ケースの局所項別収支と状態変動で床を評価する。[residualMonitor_d.cu:152](/home/sano/work/forge/solver_density_cuda/cuda_forge/residualMonitor_d.cu:152)、[Cary gen_runs.py:43](/home/sano/work/forge/case/59.flatplate_cary_m6/gen_runs.py:43)

不足情報: **対象 run の実効 YAML、実バイナリのソース差分・識別情報、残差履歴、場、VERDICT 原本がローカルにない。** 独立確認できたのはコードと生成設定、合成反例まで。記録上の高さ B の判定は **`NOT CONVERGED (stalled/plateau)`**、St は **`ALL STEADY`** であり、G15 未達を維持する。ファイル変更・`forge` 起動なし。**plan 未反映。**
