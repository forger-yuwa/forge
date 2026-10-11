# 諮問: faceh・dt の result 段レビューの Major の採否と確かめ直しの範囲 (diagnostician、Fable)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-dt.md`](../../plans/active/time_integration-line-viscous-jacobian-dt.md)・[`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **諮問先**: `diagnostician` (ユーザがこのセッションの諮問先に選んだ。新しく起動して文脈なしで送った)
- **ブリーフ**: [`briefs/2026-10-11-faceh-dt-result-disposition.md`](briefs/2026-10-11-faceh-dt-result-disposition.md)
- **レビュー**: [dt result](2026-10-11-time_integration-line-viscous-jacobian-dt-result.md)・[faceh result](2026-10-11-time_integration-line-viscous-jacobian-faceh-result.md)
- **date**: 2026-10-11

結論・採否表・手順を転記する (一部の言い回しを短くした。採否は各 plan の §6.1)。

**結論**: 7 件とも採用 (dt M2・M3・m4 と faceh M1 は内容を 1 点ずつ変えて採用)。確かめ直しは「判定器を直す → 模擬の負例 → 保存済みの証拠で回し直して旧 JSON と鍵項目を diff」までで足り、GPU も result 段の再レビューも要らない。VERDICT・分類・到達 step・segment_verdict・quasisteady.overall・floors のどれかが変わったときだけ止めて諮り直す。

## 採否表

| 指摘 | 採否 | 理由 (諮問先が確かめた根拠) |
| --- | --- | --- |
| dt M1 | 採用 | `lvcmp_rerun_judge.json` の quasisteady: L7・L7r の Q_w DRIFTING、L0 の θ_r 70/94 TRANSIENT-UNSETTLED、4 本 NOT CONVERGED。「止まる状態」「一致」は言えない |
| dt M2 | 変更して採用 | 「残差が高い = 解から遠い」は根拠なし。加えて plan の「全部入りを使うと反復の不確かさが 0.1 % を超える」は数値に反する: 全部入り 2 本どうし (L7 vs L7r) の差は θ_r 0.02〜0.06 %・Q_w 0.08 % で 0.1 % 以内。超えているのは構成間の系統差 (0.5〜0.8 %) |
| dt M3 | 変更して採用 | `assess` は途中の行の有限性を見ず、`qs()` の失敗でも REACHED、`segverdict` はテキストだけ、を確認。ただし判定器の中で `check_convergence --segment` を run ディレクトリで回すと証拠を書き換える (下記)。実データは 4 本とも全行の 5 量が有限なので、厳格化しても分類は変わらない見込み |
| dt m4 | 変更して採用 | outer_begin は step 0〜N−1 (run_0571 で確認、main.cpp:3482)。CENSORED の窓は [BUDGET−1−20000, BUDGET−1] (REACHED と同じ 20001 行) にそろえる |
| dt m5 | 採用 | 「細分化」は逆 (Grefc は粗くした)。§7 の「2000 step に限る」は §6.3〜§6.6 と矛盾 |
| faceh M1 | 変更して採用 | sweep の配列に加え `scalar_visc_line` も S の有限性の対象外。T のループで CTRL_ARRAYS の各名前について x・y・z の 3 本とも有限を要求し、非有限なら INVALID (同型を一括で潰す)。S には監査の dump の sweep 15 配列と scalar_visc_line の形・有限性を足す |
| faceh m1 | 採用 | 文書の訂正のみ |

同型の横展開:
- `lvcdt_judge.py`・`lvc75_judge.py`・`lvcgeom_judge.py` は自分で CSV を読み 2000 行の有限性を見ているので「CSV なしで通る」穴は無い。残るのは「テキストがその CSV から作られた保証がない」だけ。テキストの `[last step N]` == NSTEP−1 == CSV の最後の step と、CSV の sha256 の記録で十分。回し直しは任意。
- pregate 4 本の NaN のビット一致: 比較する配列はすべて比較の前に有限性を見ている。穴なし (回し直し不要)。

## 手順 (GPU なし)

1. `lvcmp_judge.py` を直す。判定の規則 (閾値・分岐) は触らない。
   - `assess`: DIVERGED 以外は全行の deficit・θ_r 3 つ・Q_w が存在し有限 (派生列は要求しない)。
   - `qs()`: rc ∉ {0,1}、OVERALL の行なし、4 列ぶんの分類の行が無い → INVALID。DRIFTING 等は記録のまま。
   - `segverdict`: `residual_history.csv` と `stage_manifest.json` の実在、テキストの `[last step N]` == CSV の最後の outer_begin の step、REACHED では N ≥ reach_step、CSV の sha256 を記録。再計算は scratch に CSV と manifest を複製して回し、VERDICT のトークン・判定区間名・last step の 3 つだけを照合する。ツールの sha256 も記録する。食い違いは INVALID (合わせに行かない)。
2. 模擬の負例 (途中の行の NaN/Inf・qs の失敗・CSV なし・last step 不一致・CENSORED 風の末尾 CSV) で INVALID、正常な末尾 CSV で CENSORED の窓が通ることを確かめる。
3. AWS B で `lvcmp_judge.py`・`lvcmp_rerun_judge.py` を回し直す。旧 JSON を退避し、VERDICT・arms.*.class・reach_step・segment_verdict・quasisteady.overall・floors.*.rms_ro・dev_vs_L0_pct・floor_ratio_vs_L0 を diff して記録する (追加の記録は変更に数えない)。
4. faceh M1: `lvcaudit_judge.py` を直し、NaN・空配列の負例を確かめ、手元の `run_0550〜0552` で回し直す (出力は dump の中に書かれるので退避して diff)。
5. 横展開はコードの修正と模擬の負例まで必須。回し直しは複製 (symlink でなく) で任意。

## accepted へ移してよいか

鍵項目が不変なら、result 段を回し直さずに accepted へ移してよい。§6.1 の対応欄に diff の記録のパスと判定器の commit を書く。判定の規則 (閾値・分岐・窓の定義以外) を変えていないことが条件。回し直しが要るのは、鍵項目が変わった、閾値・分岐を触った、segverdict の再計算で保存テキストと VERDICT が食い違った、のどれか。そのときは進めずに諮る。

## 案の穴 (採用)

- dt plan の「全部入りを使うと θ_r の反復の不確かさが 0.1 % を超える」は消す。正しくは「全部入り 2 本と生産の到達時点の差 (θ_r +0.48〜0.74 %・Q_w +0.72〜0.80 %) は、全部入りどうしの差 (≤ 0.08 %) と生産・粘性入りの差 (0.03 %) より 1 桁大きい、再現する系統差。どちらが解に近いかは判定していない」。
- 「止まる状態」を両 plan・case README・`plans/README.md`・メモリから拾い、「到達の時点」「到達の窓の中央値」に揃える。forge-25 の引用 (推測・参考と明記済み) はそのまま。
- 「床」は §6.5 で定義した「到達の窓の rms_ro の中央値」の名前としてだけ使う (rms_roK/roOmega は RISING)。
- 判定器の出力先は commit 済みの記録を上書きするので、退避してから回す。
