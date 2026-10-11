# 諮問: faceh・dt の 2 つの plan の result 段レビューの Major の採否と、閉じる前の再確認の範囲 (2026-10-11)

ユーザ決定 (2026-10-11):「閉じていいし、後片付けしていいよ」。2 つの plan (`plans/active/time_integration-line-viscous-jacobian-faceh.md`・`plans/active/time_integration-line-viscous-jacobian-dt.md`) を `accepted/` へ移し、AWS の自分の run を手元へ移して消す。既定はこの plan では決めない (ユーザ決定、変えない)。

## 1. 観測事実 (レビューの指摘)

レビュー本文: `notes/reviews/2026-10-11-time_integration-line-viscous-jacobian-dt-result.md` (GO-with-changes、C0/M3/m2)、`notes/reviews/2026-10-11-time_integration-line-viscous-jacobian-faceh-result.md` (GO-with-changes、C0/M1/m1)。どちらも「追加の GPU 計算は不要、記述と判定器を直し、保存済みの証拠で確かめ直してから閉じる」。

- dt M1: 未定常の到達時点を「止まる状態」「一致」と呼んでいる (dt plan 261・347 行付近)。記録では L7・L7r の Q_w は DRIFTING、L0 の θ_r 70/94 と L7r の θ_r 40 は TRANSIENT-UNSETTLED、4 本とも NOT CONVERGED。
- dt M2: 残差の大小から「解から遠い」「θ_r の正しい値は生産の側が有力」と推論している (dt plan 263〜264 行、§6.6 の 350 行も同型)。
- dt M3: `lvcmp_judge.py` の `assess` は途中の見張りの行の量の有限性を見ず (到達の行だけ)、`qs()` (check_quasisteady) の実行の失敗 (OVERALL の行なし) でも REACHED を返す。`segverdict` は `CONVERGENCE_SEGMENT.txt` のテキストだけを読み、残差 CSV の有無・一致を見ない。codex は模擬入力で 3 つとも通ることを確かめた。`lvcmp_rerun_judge.py` も同じ処理を使う。今回のデータが不正だった証拠ではない。
- dt m4: `lvcmp_rerun_judge.py` の CENSORED の床の窓が step 150000 を要求するが、ループは 0〜149999 (今回は REACHED なので影響なし)。
- dt m5: 再開の条件 3 の「流れ方向の細分化」は逆 (forge-25 の Grefc は粗くした)。§7 の「2000 step に限る」が古い。
- faceh M1: `lvcaudit_judge.py` は sweep の配列 (`rhs_s*`・`dqnew_s*`・`dqold_s*`) の存在だけを見て有限性を見ず、T 判定は同じビット列の NaN を一致として通す (codex が模擬入力で確認)。実データ (`run_0550`〜`0552`) の対象配列に非有限は無く、結論 (FAIL・T 保留) は変わらない、と codex。plan §6.8 の「次から検査する」が未実装。
- faceh m1: `methods/time_integration/implementation.md:543` の単体試験の説明 (行列全体の最大値で割る) が古い。`plans/README.md` の faceh の行が古い。

## 2. 期待値と出典

- AGENTS.md: Critical・Major の採否は上位に諮る (条件 5)。採否は plan §6.1 に書く。
- 判定器の性質: 「判定不能を合格にしない」(AGENTS.md の収束確認、`test_gate_bad_input.py` の考え方)。

## 3. 再現条件・データの所在

- dt の判定の記録: `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcmp_judge.json`・`lvcmp_rerun_judge.json`・`lvcdt_judge.json` (commit 済み)。
- run の実体: AWS B (`~/forge-faceh-audit/case/45.isobutane_m6_d155/run_0553〜0559`・`run_0567`・`run_0585` の写し)、AWS A (`run_0580〜0587`・`run_0550`)。いま手元 `~/forge-evidence/2026-10-11-faceh-dt/{A,B}/` へ sha256 照合つきで取り込み中。消すのは採否と確かめ直しの後。
- faceh の run (`run_0550〜0552`・`0560〜0566`・`0570〜0579`) は手元 `~/forge-evidence/2026-10-10-faceh/` (照合済み)。`run_0550` の linedump の残りは AWS にもある (取り込み中)。

## 4. 実施済みの操作と結果

- codex result 段を 2 本 (commit `274359b7`、diff base: faceh は分岐点 `3de780a4`、dt は `eea98f4b^`)。
- 同型の欠陥の横の確認 (grep のみ、未修正): `segverdict` と同じ「区間付きの収束のテキストだけを読む」形は `lvcdt_judge.py:228`・`lvc75_judge.py:216`・`lvcgeom_judge.py:215` にもある (これらは 2000 step の有限性の判定で、VERDICT が無ければ証拠の欠けとして INVALID)。「NaN のビット一致を一致として通す」比較 (`view(np.int64)` の `array_equal`) は `lvc*_pregate.py` にも多数あるが、pregate は plan-6 の Major で配列の有限性の検査を足してある (全配列かは未確認)。

## 5. 仮説 (私の採否の案)

全件採用し、GPU は使わない。

- dt M1・M2: 採用。§5.1・§6.4・§6.6・case README・`plans/README.md`・自分のメモリの文を「同じ停止規則による到達の時点で観測された差」に揃える。L0/L5 は「到達の時点の値の差 0.03 % 以内」。「解から遠い」「θ_r の正しい値は生産の側が有力」「止まる状態の再現」は消し、「到達の窓の rms_ro の中央値が 3.30・3.36 倍」「構成間の差が許容の目安 0.1 % を超えた」に限る。解の精度の順位は保留と書く。
- dt M3: 採用。`lvcmp_judge.assess` で、DIVERGED 以外は全行の deficit・θ_r・Q_w の有限性を要求、`qs()` の OVERALL の行が無い・rc が 0/1 以外は INVALID (DRIFTING などの正常な判定は記録として残す)、`segverdict` は残差 CSV を必須にして `check_convergence.py --segment` をその場で回し直し、保存のテキストの VERDICT と一致を要求し、CSV の sha256 を記録する。模擬入力 (途中の NaN・Inf、qs の失敗、CSV なし) で INVALID になることを確かめ、AWS B 上で `lvcmp_judge.py`・`lvcmp_rerun_judge.py` を保存済みの証拠で回し直し、VERDICT が変わらないことを確かめてから消す (変わったら止めて諮り直す)。
- dt m4: 採用。CENSORED の窓は最後の残差の step (BUDGET − 1) で終える。plan に窓の定義の変更を書く。
- dt m5: 採用 (文の訂正)。
- faceh M1: 採用。`lvcaudit_judge.py` の S に、3 本の書き出しの sweep の配列の形 (N, 5) と有限性を足し、T の比較の前に非有限を INVALID にする。NaN・空の配列の模擬入力で確かめ、手元の保存済みの証拠 (`run_0550〜0552`) で回し直して結論が変わらないことを確かめる。
- faceh m1: 採用 (methods と索引の文)。
- 同型の横展開: `lvcdt_judge.py`・`lvc75_judge.py`・`lvcgeom_judge.py` の区間付きの収束も CSV から回し直す形に直し、lvcdt は A 上で、lvc75・lvcgeom は手元の証拠の写し (作業用ディレクトリ) で回し直す。pregate の NaN のビット一致は、比較の前に全配列の有限性を見ているかを確かめ、穴があれば同じく直す。

## 問い

1. 上の採否 (全件採用) でよいか。却下・変更すべきものは。
2. 確かめ直しの範囲: dt M3 の判定器 2 本は AWS 上で回し直す。faceh M1 は手元。横展開 (lvcdt・lvc75・lvcgeom の収束の照合、pregate の NaN) まで回し直すべきか、コードの修正と模擬入力の確認で足りるか。
3. 確かめ直しで VERDICT が変わらなければ、そのまま accepted へ移してよいか (result 段をもう一度回す必要はあるか)。
