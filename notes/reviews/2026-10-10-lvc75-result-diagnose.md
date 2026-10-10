# codex 諮問 (diagnose): lvc75-result

- **brief**: [`notes/reviews/briefs/2026-10-10-lvc75-result.md`](../../notes/reviews/briefs/2026-10-10-lvc75-result.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **date**: 2026-10-10
- **commit**: `511cc614` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.0 min, rc=0
- **結論**: **次は候補(a)、マスク7を維持し、`time.deltaT.lineDtDirectional` の1／0だけを変えるA/Bを行う。**
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvc75_judge.json`, `notes/reviews/briefs/2026-10-10-lvc75-result-summary.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| Major | 今回の「支持」：**限定付きで採用** | 事前ゲートは I 75・R 9・S 3・V 5、本判定は59項目すべて合格。記録された分類を判定器の `pair()` に渡した再評価も「支持」だった。結論は登録どおり、**この条件・2000 step以内では、熱伝導の近傍Kを外す介入で非有限化を回避した**まで。[判定規則:425](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:425)。 |
| Major | 「マスク5は安定」「B1は定常値で止まった」：**却下** | B1・B2とも `check_convergence --segment: NOT CONVERGED`（区間0〜1999）。末尾500 stepの `rms_roK` の傾き×500は、それぞれ **+0.179、+0.228桁**。密度残差の低下だけでは評価できない。**「2000 step有限、両方未収束、長期挙動は未確認」**と書く。[結果:480](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/lvc75_judge.json:480)、[同:583](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/lvc75_judge.json:583)。 |
| Major | 「同じ場所の熱過渡を違う速さで動かした」：**要再検証** | 同じ壁際で温度変化が大きいことは、共通の原因や同じモードを示さない。Kを変えると補正の空間構造も変わる。さらに出発場は `NOT CONVERGED`。[親plan:411](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:411)。局在は観測点の選定に使い、物理的な熱過渡という解釈は保留する。 |
| Major | 「共通関数PASSだからK自体は正しく、組み合わせだけが問題」：**要再検証** | PASSは凍結物性・薄層モデルの微分についての結果。共通関数自身が実残差の厳密な微分ではないと明記している。[共通関数:94](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94)。今回の介入照合も、Kの変更内容の確認である。実残差・境界・更新との整合と、§6.8のC・D・Tは未解決のまま残す。 |
| Minor | faceh planの閉じ方：**後継への明示的な移管を採用** | [§5.1 #8・#9:66](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:66)を消化済みにしない。今回の結果と採否を追記し、#7を完了、#8・#9を後継のactive planへ移管したうえで、facehを「診断の区切り」として閉じる。#9の「マスク5でも壊れたら」という着手条件は外す。マスク5が有限でも、全部入りの整合確認は必要である。 |

根拠runは次のとおりです。数値は保存された判定JSONからの引用です。

| リポジトリルートからのrunパス | 記録された判定 |
|---|---|
| `case/45.isobutane_m6_d155/run_0570_lvc75_m7_a1/` | `DIVERGED`、最初の非有限：CSV 347／検知表示348 |
| `case/45.isobutane_m6_d155/run_0572_lvc75_m7_a2/` | `DIVERGED`、最初の非有限：CSV 279／検知表示280 |
| `case/45.isobutane_m6_d155/run_0571_lvc75_m5_b1/` | `FINITE`、2000 step、`NOT CONVERGED` |
| `case/45.isobutane_m6_d155/run_0573_lvc75_m5_b2/` | `FINITE`、2000 step、`NOT CONVERGED` |

主要な判定成果物は [lvc75_judge.json](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/lvc75_judge.json)、恒久索引は [case README:194](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/README.md:194)です。

結論: **次は候補(a)、マスク7を維持し、`time.deltaT.lineDtDirectional` の1／0だけを変えるA/Bを行う。**

第1仮説: **方向別dtによる大きな擬似時間刻みが、全部入りの更新で過大な補正を許しており、方向別dtを外すだけで今回の2000 step以内の非有限化を回避できる。** 確度: **中。ただしdtとの因果は未確認。**

根拠: 今回はDが不変で、Kの熱伝導項だけを除くと有限完走へ分かれた。熱伝導Dは常に入り、Kだけがマスクで切り替わる。[共通関数:132](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:132)。また、時間項はDに `V/Δτ` として入り、方向別dtは内部ライン面の音響制約を除外する。[時間項:868](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:868)、[設定仕様:75](/home/sano/work/forge-faceh/procedures/solver-settings.md:75)。

これは「Kが強すぎる」と確定した意味ではない。Dを残してKを外すと、隣接補正との相殺も失われるため、マスク5の改善をKの符号・式の誤りへ直結させられない。

反証条件: 介入が成立し、方向別dtありで破綻を再現した対照に対して、**方向別dtなしも各2本とも2000 step以内に非有限化する**こと。この場合、「方向別dtを外すだけで回避できる」を棄却する。

第2仮説: 実残差と近似LHS、壁拘束、ライン外結合、分離したSST更新の不整合が、Kを入れた補正で顕在化する。確度: **中・原因箇所は未確認**。第1仮説と併存し得る。

第3仮説: 製品の組立・演算精度の未解決項目が補正を変えている。確度: **低・未確認**。§6.8のC不一致・D判別不能3件・T保留は、今回のゲートでは解消していない。

判別 A/B:

- **変更点は1点**：A=`time.deltaT.lineDtDirectional: 1`、B=`0`。`lineImplicit: 1`、値3、マスク7、キー5、`cfl_pseudo: 4`、緩和0.7、5 sweep、ISP 0、段③の同一バイナリ、同一出発場を固定する。各2本、新規runで最大2000 step。既存runを新対照の代用にしない。
- **介入ゲート**：同じ初期状態で、残差・物性・拘束・Kが不変であり、実際の `dt_local` が変わったことを確認する。Dの差は、拘束前では時間項 `V/Δτ` の変更を丸め込みで再現し、拘束行ではその差が上書きされることを確認する。旧7/5ゲートの「D不変」をそのまま流用しない。
- **見る量**：全残差列、最初の非有限step、保存場の有限性とρ・P・Tの正値、壁隣接域の `dt_local` と補正量。有限側には区間0〜1999の `check_convergence --segment` のVERDICTを残す。
- **Aが2本ともDIVERGED、Bが2本ともFINITE** → 第1仮説の「この条件・期間での回避」を支持。平均流とSSTを含むdt変更の効果であり、熱伝導単独の機構や長期安定性は確定しない。
- **A・Bとも2本ともDIVERGED** → 第1仮説の十分性を棄却。第2・第3仮説は残る。
- Aが破綻を再現しない、Bの結果が分かれる、ゲート不成立 → 判別不能。

方向別dtを外すと、同じstep数でも擬似時間の進み方が変わる。したがって、Bの温度変化が小さいだけで「熱過渡が解消した」と判定しない。

やらない方がよいこと: マスク5の本番化、Kの列の追加削除、CFL・拘束・精度の同時変更、別の出発場への変更。候補(c)も今回は先行させない。薄層近似と実残差の微分が一致しないこと自体は設計上あり得るため、測る成分・凍結範囲・拘束方向を定めずに「不一致ならバグ」とする試験にはできない。親plan §6.16の残差全体の応答を、熱伝導Kだけの誤差とも扱わない。

呼び出し側の前提への異議: H1の局在と残差増大は観測だが、「同じ熱過渡」「外側反復の増幅率」は解釈である。H2の「K自体は正しい」も、薄層モデル内の照合という範囲を超えている。今回の支持は、これらを確定していない。

facehは、今回の結果・採否と未解決事項の移管先を記録し、resultレビューを完了した時点で閉じることを勧める。後継planには今回のdt比較と#8・#9を残す。値2・3と診断切替の残置、既定値0は、[親plan §6.17:419](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:419)の決定を維持する。

不足情報: 対象runのCSV・HDF5・`CONVERGENCE_SEGMENT.txt`・実行バイナリは確認したローカル保存先にない。独立に確認できたのはコード、plan、判定JSON、および記録分類からの判定分岐であり、生データの再測定や53ファイルのハッシュ再照合ではない。局所補正の大きさと実効dtの分布も未確認である。

ファイル変更・forge起動なし。**plan未反映**。採否・結果・後継への移管は、依頼どおり呼び出し側で反映すること。
