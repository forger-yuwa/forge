# codex 諮問 (diagnose): lvcgeom-result

- **brief**: [`notes/reviews/briefs/2026-10-10-lvcgeom-result.md`](../../notes/reviews/briefs/2026-10-10-lvcgeom-result.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **date**: 2026-10-10
- **commit**: `39a6e7a8` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 4.2 min, rc=0
- **結論**: ?
- **extra**: `case/45.isobutane_m6_d155/lvcgeom_pregate.py`, `case/45.isobutane_m6_d155/lvcgeom_judge.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcgeom_judge.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論は、**§6.9 の棄却を「この修復だけでは十分でない」に限定して採用し、次は同一バイナリでマスク 7/5 を比較し直す**ことです。

コード・判定 JSON を確認し、記録された分類から判定器の分岐を再評価して同じ「棄却」を得ました。事前ゲートは I=71・R=9・V=4、本判定は50項目すべて合格です。ただし、AWS の CSV・HDF5・帳簿はローカルにないため、生データの再測定ではありません。

採否表:

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 「精度の不整合は主因でない」 | **却下**。「座標差の修復だけでは非有限化を回避できない」を採用する。記録された最初の非有限 step は旧292/198、新356/234で、すべて `DIVERGED`。他の精度依存や複合原因を排除する試験ではない。[登録:335](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:335)。step 50 の −1.8% は観測値として残せるが、増幅率や機構の同一性の根拠にはしない。新の2本は旧最大292の2倍を超えず、登録した「遅らせる」の条件も満たさない。[判定器:256](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcgeom_judge.py:256)。 |
| Major | 「1331節点が壁法線の1列で非有限 → Thomas の枢軸破綻」 | **要再検証**。同じ起点の記録は **1ライン121節点**で、1331=121×11。[親plan:184](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:184)。複数列か、複数変数に重複する節点を数えたのか、抽出を確認する必要がある。対案は変数別の非有限節点集合、その和集合、実際の `line_cells/line_offsets` による列別集計。さらに LAYOUT2 は分解失敗を検知すると旧補正へ戻すため、「検知された枢軸破綻がそのまま列全体のNaNになる」という経路ではない。[実装:2508](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2508)、[同:2539](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2539)。非有限の伝播と発生源を分ける。 |
| Major | 候補(ii): step 0 のライン行列から外側反復の増幅を評価 | **次の一手としては却下**。D/K だけから外側反復の増幅行列は決まらない。実残差の微分、ライン外の結合、5 sweep、緩和、境界・SST の更新が必要。コードでもライン外の補正を後続 sweep の RHS に取り込む。[実装:958](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:958)。加えて起点 S0 は `NOT CONVERGED` である。[親plan:411](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:411)。対案は、まず単因子比較で対象を絞り、その後に実際の更新写像の応答を測る。行列の条件数や特異ベクトルを外側反復の不安定モードと呼ばない。 |
| Major | 候補(iv): 等温壁拘束×方向別dtの2×2 | **今回は却下**。値3では `lineViscCoupling >= 2` により拘束行が強制され、`implicitThermalJacobian` のビット2だけでは切り替わらない。[実装:1263](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1263)。拘束だけ外すと、隣が等温壁なら熱伝導Kを消す処理との整合も崩れる。[共通関数:141](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:141)。対案は拘束とdtを固定した7/5比較。既往の値0での壁拘束単独試験を出し直す必要もないが、値3との相互作用は未解決のまま残す。 |
| Major | §6.8・§6.9をもって「組立・精度を除外済み」として閉じる | **却下**。§6.8のC不一致、D判別不能3件、T保留は残る。[plan:276](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:276)。§6.9は係数変更の成立を示すが、これらを解消した試験ではない。対案は否定できた十分性と未解決事項を分けて処置する。 |

**結論:** 候補(i)、すなわち**段③の同一FP64バイナリ・同一出発場で、`FORGE_LVC_TERMS` の7/5だけを変える比較を新規登録して行う**。

**第1仮説:** この固定条件では、熱伝導の近傍Kを除くことが2000 step以内の非有限化を回避するのに十分である。  
**確度: 中。ただし探索的な根拠に限る。**

- **根拠:** 親planの `case/45.isobutane_m6_d155/run_0305_e1_m7/` は150 stepで非有限、`run_0306_e1_m5/` は2000 step有限。ただし登録上は比較無効で、中間場も欠けている。[親plan:279](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:279)。コードではマスクのビット2が熱伝導Kだけを切り替え、熱伝導Dは常に入る。[共通関数:107](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:107)。
- **反証条件:** 介入ゲートを通した新しい比較で、マスク7・5が各2本とも2000 step以内に非有限化すること。その場合、この回避仮説を棄却する。ただし同じ破綻機構とは断定しない。

**第2仮説:** 熱伝導Kを除いても、実残差と近似LHS、ライン外結合、分離したSST更新などの不整合によって非有限化する。確度: 低。共通関数は実残差の厳密な微分ではないが、原因箇所は未確認。[共通関数:94](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94)。

**判別 A/B:**

- **変更点:** A=`FORGE_LVC_TERMS=7`、B=`5`だけ。値3、キー5、方向別dt、上限なし、CFL4、緩和0.7、5 sweep、ISP0、面エンタルピー既定を固定する。出発場は `run_0183_ns_coldmesh_tw300_ext/res_100000.h5`。各腕2本、新規runでA1/B1/A2/B2、最大2000 step。
- **事前ゲート:** 同じ凍結状態の7・5・7再実行を採取する。入力・状態・物性・dt・拘束・接続を照合し、Dは再実行でビット一致するなら腕間もビット一致を要求する。残差・初期RHSは既存の再実行ノイズ規則を使う。Kの変更が自由なエネルギー行の熱伝導項だけであることを、独立計算した κ∂T_j/∂Q_j と照合する。許容と変更信号の解像条件は採取前に固定する。**旧試験のマスク0対値0の同値性は別の問いであり、この新規7/5比較の必須ゲートには含めない。旧INVALIDは保存する。**
- **見る量:** 全残差列、最初の非有限step・変数別節点集合、100 stepごとの保存場の有限性とρ/P/Tの正値。`check_convergence --segment` の区間付きVERDICTも保存する。有限完走を収束とは呼ばない。
- **分岐:**  
  **Aが各2本DIVERGED、Bが各2本FINITE** → この条件・期間で熱伝導K除去による回避を支持。Kの式の誤りやマスク5の長期安定性は証明しない。  
  **A・Bが各2本DIVERGED** → 第1仮説を棄却。熱伝導Kなしでも破綻する経路が残る。  
  Aが破綻を再現しない、腕内で結果が分かれる、ゲートが外れる → 帰属不能または判別不能。事後延長しない。

局在は、この比較に添える**観測点の選定**に使います。まず既存データの抽出を修正・確認し、両腕で同じ列と隣接列を見るべきです。121節点単位の番号付けなら、4960〜4945は列40の壁側、4961は列41の先頭に相当しますが、最終確認は接続情報で行ってください。step 100の温度上昇は既に大きな有限変化なので、微小な不安定モードや最初の発生源の証明には使えません。

**やらない方がよいこと:** マスク5を本番設定にする、Kの密度列をさらに削る、CFL・dt・拘束・精度を同時変更する、非有限化の遅延を安定化と呼ぶ、列全体のNaNだけからThomasの枢軸を原因と決めること。

**呼び出し側の前提への異議:** 「精度の問題ではない」「1331節点＝1列」「step 0のD/Kだけで外側反復を判定できる」は受け入れません。壁際のT・ρの変化も、エントロピーの不安定モードという解釈とは分けて記録してください。

本planは、面エンタルピーとLHS座標差の**各修復単独では回避できなかった、という診断の区切り**として閉じられます。ただし、閉じる前に§6.9の結果を記録し、残作業表へ①新規7/5比較、②§6.8のC・D・T未解決、③実残差・更新写像との整合確認を残すか、後継のactive planへ明示的に移管してください。「原因解明」「製品経路PASS」としては閉じられません。値2・3と診断切替は、親plan §6.17の残置決定を維持します。

**不足情報:** AWSの `case/45.isobutane_m6_d155/run_0560_lvcgeom_old_a1/`〜`run_0566_lvcgeom_old_dump2/` の生データ、局在抽出スクリプト、非有限節点の重複除去後の集合とライン対応です。現時点では、局在と最初の非有限発生源を独立に確定できません。

ファイル変更・forge起動なし。**plan未反映**。採否と事前登録の反映は、依頼どおり呼び出し側が行う前提です。
