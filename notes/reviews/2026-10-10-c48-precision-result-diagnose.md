# codex 諮問 (diagnose): c48-precision-result

- **brief**: [`notes/reviews/briefs/2026-10-10-c48-precision-result.md`](../../notes/reviews/briefs/2026-10-10-c48-precision-result.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-10
- **commit**: `9da7430b` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.6 min, rc=0
- **結論**: **case/48 の追加 run は行わず、共通 double 座標に揃えた保存済み場で、動的 `je` と固定 `je` の後処理 A/B を一度だけ行う。**
- **extra**: `notes/reviews/2026-10-10-float-v6-redesign-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0／Major 4／Minor 1）。実測値と VERDICT は提供記録として扱い、コードで確認できた事実と区別する。

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 「共通 double 座標・同じ節点列で後処理した」 | **要再検証。実装が登録条件を満たしていない。** [抽出器:27](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/c48_prec_series.py:27) は各結果ファイルを読み、[読込処理:41](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/cooled_plate_eval.py:41) はその `/MESH/COORD` を使う。一方、[出力処理:141](/home/sano/work/forge-integ-1005/solver_density_cuda/output/output.cpp:141) は `geom_float` の座標を書いている。共通 double 座標を要求した [plan:648](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:648) と不整合。節点対応を確認し、共通入力の double 座標から抽出列・積分距離を固定して再集計する。影響の大きさは未確認。 |
| Major | 「準定常閾値を超えたことが判別不能の原因」 | **単独の説明として却下。** 提供値では δ*(0.6) の D＋E＝8.35e−4＝1.67τ。準定常条件だけを通しても、変動幅を含む比較で判別不能が残る。[判定器:78](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/c48_prec_judge.py:78) の分岐を合成入力でも確認した。**準定常条件と D＋E 条件の二つを満たさない**と記録する。 |
| Major | 「case/45 の % 級の偏りは case/48 では観測されない」 | **窓平均の記述に限って採用、原因の除外には却下。** 「この4本・登録窓・指定14量の抽出値では、窓平均に % 級の差は観測されなかった」まで。case/45 自体も両腕 `NOT CONVERGED / NOT ALL STEADY` で、比較は停止時の値である（[plan:572](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:572)）。SST、commit、残差集積、面演算のいずれも候補から外れない。#15 の候補順を今回だけで変更しない。 |
| Major | 判定器が前提条件を保証する | **要再検証。** [判定器:67](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/c48_prec_judge.py:67) は `check_convergence` の実行失敗・VERDICT 欠落を記録するだけで判定不能へ送らない。また、量の有限性検査は登録窓だけ（:48–50）、完走判定は最終 step だけ（:38–39）。対案は、実行成功、24点全体の有限性、正式ツールの正常終了を前提検査へ入れること。今回これらが実際に失敗した証拠はないが、判定 JSON だけでは否定できない。 |
| Minor | float の揺れを一律に「相対1e−4、FP64より3桁大きい」と要約 | **修正して採用。** `run_0057` の δ*(0.6) は幅8.7e−4で、FP64の1.6e−7に対し約5,400倍。[brief:29](/home/sano/work/forge-integ-1005/notes/reviews/briefs/2026-10-10-c48-precision-result.md:29) の例外を明記する。「時間の揺れ」は、ここでは定常計算の反復履歴の変動である。 |

結論: **case/48 の追加 run は行わず、共通 double 座標に揃えた保存済み場で、動的 `je` と固定 `je` の後処理 A/B を一度だけ行う。**

第 1 仮説: `run_0057` の δ*(0.6) に固有の大きな変動幅は、離散的な `je` の切替によって主に増幅されている。 確度: 中  
　根拠: 提供記録では、この断面だけ `je=74〜76` が変動する。[`station():96`](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/cooled_plate_eval.py:96) は `je` を整数で選び、その変更が積分上限だけでなく正規化に使う ρₑ・uₑ も変える。したがって、δ* の跳びを境界層全体の変動と直結できない。  
　反証条件: 共通座標・同じ保存場で `je=76` に固定しても、当該系列の相対変動幅が動的 `je` の場合の半分を超えて残るなら、「切替が変動幅の過半を説明する」を棄却する。

第 2 仮説: `je` の切替とは別に、float の反復状態に変動が残る。 確度: 中。提供記録では `je=78` で一定の x=0.3 でも θ・δ* が `OSCILLATING`、`je` に依存しない q_w にも変動がある。発生源を commit や原子加算まで特定する証拠はない。

判別 A/B: **ソルバ追加計算は0 step**。4本の24,000〜48,000 step、各13点を使い、両腕とも共通 double 座標・同じ節点列に揃える。

- A：登録どおり、各場から動的 `je` を選ぶ。
- B：x=0.6 の `je` だけを全4本・全時刻で76に固定する。ρₑ・uₑは各場の同じ節点から取り直す。
- 主判定は `run_0057` の δ*(0.6) の相対変動幅。B/A≤0.5なら端点選択による増幅を支持、B/A>0.5なら第1仮説を棄却する。併せて同じ閾値の `check_quasisteady` と D・E を記録する。

これは抽出感度の診断であり、Bの結果で旧登録判定を合格に変更しない。x=0.3 の不合格も別に残る。この確認後、case/48 の延長はせず、#15 は既定どおり `case/45.isobutane_m6_d155/run_0483_hp7_k1/`・`run_0484_fl1_f32_k1/` の結果に基づいて進める。

やらない方がよいこと: 合格を得るための窓延長・閾値緩和、固定 `je` による旧結果の上書き、「SSTは無関係」「軸対称が原因」という除外・確定、真因未確認での `qAccumulatorFP64` 対応。

呼び出し側の前提への異議: **観測された窓平均差の小ささは採用できるが、定常的な精度差の小ささや原因の除外には読み替えられない。** #15 に対しては「今回の平板構成で case/45 と同規模の窓平均差は再現しなかった」という探索上の情報に留まる。軸対称・TP・ライン陰解法・収縮部・格子条件との相互作用は絞れていない。

不足情報: 対象の `case/48.flat_plate_cooled_m4/run_0057_pv_f32a/`、`run_0058_pv_f32b/`、`run_0059_pv_f64a/`、`run_0060_pv_f64b/`、実YAML、HDF5、判定JSONはこの作業ツリーに存在せず、実測の再集計はできなかった。提供された VERDICT は、4本とも `CHECK FAILURES（NOT CONVERGED）`、float の厚さ群が `NOT ALL STEADY`、FP64 全量と全4本の壁量・積分量群が `ALL STEADY`、総合が「判別不能」である。

ファイル変更・forge 起動は行っていない。**plan 未反映**。反映先は `architecture-float-state-double-geometry.md` §6.22・§5.1 #10／#15。
