# codex レビュー: time_integration-line-viscous-jacobian-faceh (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `0be79325` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m2
- **focus**: §4・§6 (諮問の採否を反映済み) と、台本 lvcfh.sh・判定 lvcfh_judge.py が §6 の文言どおりかを点検する。台本は AWS の ~/forge-wallfit/case/45.isobutane_m6_d155 で動く (cold_cfl.py・cold_pair.py は同じディレクトリ)。判定器は人工入力で支持・棄却・負の値・列の欠けを確かめた。コードは変えない。
- **extra**: `case/45.isobutane_m6_d155/lvcfh.sh`, `case/45.isobutane_m6_d155/lvcfh_judge.py`, `notes/reviews/2026-10-10-lvc-faceh-design-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
§4 の単因子 A/B と、§6 の結論を「2000 step 以内の非有限化回避」に限る方針は妥当です。  
ただし判定器は、証拠欠損や事前登録と異なる共通設定でも「支持」を返せます。以下を投入前に修正してください。

1. **Major — 証拠保存の不備をゲートが検出できない**

   根拠: [lvcfh_judge.py:175](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:175) は格子実体ではなく `MESH_SHA.txt` 同士を比較し、186 行は帳簿のサイズしか見ません。[99 行](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:99)では `DIVERGED` を返した時点で保存場の検査を飛ばします。

   **メモリ上で `main()` まで実行した結果**、全腕の `nozzle.h5`・`res_0`、A 側の `res_100`・`res_nan_*` がなく、帳簿がヘッダーだけ、B 側の収束判定記録もない入力が、**全ゲート合格・支持・終了コード 0** になりました。HDF5 読み込みは正常な模擬配列に置換した、判定制御の試験です。

   対案: 分類とは独立した証拠ゲートを設ける。格子の実ハッシュ、初期場、破綻までに保存されるべき場、`detectNaN` に対応する出力の可読性、帳簿の節点・記録範囲を検査する。`FINITE` 側には区間が確認できる収束 VERDICT を必須とする。ただし **`PASS` 自体は要求しない**。

2. **Major — 腕同士の一致だけでは、事前登録した条件を保証しない**

   根拠: [lvcfh_judge.py:178](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:178) は全設定を `a1` と比較するだけです。[cold_cfl.py:121](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/cold_cfl.py:121) は親の設定を継承し、185 行で保証する緩和値も「親と同じ」であって `0.7` ではありません。

   模擬試験では、全腕を共通して **CFL 0.5・緩和 0.1・上限 50** に変更しても、全ゲート合格・支持・終了コード 0 でした。台本は CFL 4 を指定しますが、継承される緩和・sweep・上限についてはこの穴が残ります。

   対案: `a1` の生成前に親設定・入力を固定し、生成後に §6 の期待値と照合する。少なくとも方向別 dt、キー 5、値 3、CFL 4、緩和 0.7、sweep 5、上限なし、ISP 未設定、2000 step・出力間隔 100 を検査する。メッシュ品質も「各腕で同じ記録」だけでなく、当該格子に対応する許容内の VERDICT を確認する。

3. **Major — 欠損・破損入力を理由付き `INVALID` にする処理が未完成**

   根拠: [lvcfh_judge.py:66](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:66) は `KeyError` を捕捉せず、[191 行](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:191)は `INVALID` の run にも無条件で `records()` を呼びます。

   実行確認では、`step` 欠損と `phase` 欠損はそれぞれ **`KeyError`**、CSV 不在時の `records()` は **`TypeError`** でした。また、列数が多い壊れた行も `broken_rows=0` で通りました。例外終了は判定結果の保存を保証しません。

   対案: 必須ヘッダー・重複ヘッダー・行の列数を検査し、ファイル欠損や解析失敗を run ごとの理由付き `INVALID` に変換する。記録不能な run は `records()` を省略し、異常系でも今回の判定 JSON を出して終了コード 1 にする。人工入力試験は関数単体に加えて `main()` まで通す。

4. **Minor — 最初の非有限 step を取り違える**

   根拠: [lvcfh_judge.py:100](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:100) は CSV よりログを優先します。**CSV の最初の非有限が step 20、ログが step 122** の模擬入力で、`diverged_step=122` になりました。§6 の「最初」と一致せず、遅延倍率の記録にも影響します。

   対案: CSV とログの検出 step を別々に保存し、代表値は最小値にする。両者が異なる場合も記録する。

5. **Minor — §1 のマスク比較の断定が、§3 の留保より強い**

   根拠: [plan:18](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:18) は「マスク 7 で発散し、5 で回る」と一般化していますが、§3 は `run_0306` の比較を無効・探索的観測としています。今回もマスクを比較しないため、熱伝導 K の必要性は確定できません。

   対案: §1 を「既往の探索ではその差が観測されたが、帰属は未確定」に修正し、問いを §6 と同じ「値 3・マスク 7 における面エンタルピー精度切替の効果」に統一する。

**推奨は、上記 1 → 2 → 3 → 4 → 5 を修正して事前登録を固定し、現在の主の 4 本だけを実行することです。**

親 plan は本線不採用で閉じていますが、今回は限定した診断であり、point 仕上げの `NOT_SUPPORT` とも別の問いです。面エンタルピーの切替は[実装:401](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:401)で確認でき、値 3 が厳密 Jacobian ではなく近似であることとも矛盾しません。solver 変更を伴わない今回は、cell・周期・他ケースへの回帰拡大や派生量の定常性判定を先に追加する必要はありません。

旧 run・初期 HDF5・AWS バイナリはローカルにないため、既往の破綻 step と方向微分値は独立再計算できていません。コード照合、シェル構文検査、メモリ上の判定器試験を実施しました。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 2
