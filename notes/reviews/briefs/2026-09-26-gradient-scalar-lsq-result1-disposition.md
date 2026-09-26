# 諮問ブリーフ: result 1 回目 (NO-GO, C0/M5/m2) の採否と、追加 A/B の事前規則 (2026-09-26)

plan: `plans/active/gradient-scalar-lsq-unification.md`。レビュー: `notes/reviews/2026-09-26-gradient-scalar-lsq-unification-result.md`。
前回諮問: `notes/reviews/2026-09-26-gradient-scalar-lsq-s2-result-diagnose.md` (S2 = 条件付き・未合格)。

## 決めてほしいこと

指摘 7 件それぞれの採否 (採用なら下の対案で良いか)、および M2・M3 の追加 A/B の判定規則をここで固定すること。

## 読んでよいファイル

レビュー本文、plan §5.1・§6、`solver_density_cuda/tools/{check_floor_ratio.py,stage_manifest.py}` (該当関数)、`case/09.Taylor-Green/_g0_lsq_seam/s3_perf.py`、`solver_density_cuda/tools/check_passive_budget.py` の先頭 docstring。巨大ファイル・h5 は読まない。

## 呼び出し側での再現 (観測事実)

- **M1 再現**: 必須 5 列の合成 CSV で、(a) 前半 80 点 1e-6・末尾 20 点 0 → `check_floor_ratio` PASS (`check_convergence` は NOT CONVERGED)、(b) 参照が全点 Inf・対象が定数 1e-6 → PASS。
- **M4 再現 (原因も特定)**: `s3_perf.py` が起動するバイナリのファイル名は `forge_f99f236d` なので、プロセス名は `forge` ではなく `pgrep -xc forge` に数えられない。case/48 の初回 S3 が自分の #5a・#5b (同じバイナリ) と重なって汚染され、それを検出できなかったのはこのため。取り直し (`S3_case48.txt`) の時点で他セッションの `forge` (case/58、`forge-cht` のバイナリ = 名前 `forge`) は 0 本だったが、自分の別 run が無かったことは記録上保証されない。
- **M5 確認**: `stage_manifest.py` の `HARD_PATTERNS`・`YAML_HARD_PATHS` に `scalarGradient` が無い。

## 採否の案 (親)

| # | 案 |
| --- | --- |
| M1 | **採用**。`check_floor_ratio.py` で `analyze` の列ごとの判定不能 (「末尾窓の代表値が 0」等) を不合格へ伝播、参照系列の NaN/Inf・必須列・床 > 0 を検査。`--selftest` に 2 例を追加し、`test_gate_bad_input.py` の流儀で負の試験。修正後に S2 の全床比判定を取り直す |
| M2 | **採用**。実装直前版 `36d8ba03` を AWS でビルドし、case/48 で `run_0950_sglsq_s2_gg` と同一入力・同一実効設定 (`scalarGradient` キーは旧版に無いので書かない、chi は既定のまま) で 48000 step (`run_0957_sglsq_s2_base36d8`)。事前規則案: 現行 gg (0950+0952 の連結 48000) と比べ、`check_floor_ratio --start run_0950+0952` の末尾平均比 ≤ 1.5 (全列) かつ物理量 (Cf・q_w・δ*・θ 3 station) 差 ≤ 0.1 % → **A**: 床移動は本 plan の commit 群より前から (本 plan の実装は gg 経路の床を変えていない)。**B**: 床比 > 1.5 または物理差 > 0.1 % → 本 plan の commit 群が gg 経路を変えた疑い → 実装を止めて調べる。床移動の真因の特定 (09-12 以降の 96 commit の bisect) は本 plan の範囲外とし、別途 F 項目に切り出す |
| M3 | **採用**。(a) `check_passive_budget.py --mode fct` を `run_0954/0955` に適用。(b) nSub 感度: 同じ FCT smoke を `nSubIterDualTime 30` で gg・lsq 各 1 本 (`run_0958/0959`)。事前規則案: nSub 15→30 で各双子自身の変化量 (gg15 vs gg30、lsq15 vs lsq30) が lsq−gg 差 (nSub 15) の 1/10 未満なら「lsq−gg 差は反復不足ではなく作用素差」と記録。そうでなければ反復不足が混ざるので差の審査は nSub 30 の双子で行う。旧上限超過の記録は残す |
| M4 | **採用**。監視を `nvidia-smi --query-compute-apps=pid` で GPU を使う PID 全体と測定対象 PID を照合する方式に替え、走行中に測定対象以外の GPU プロセスが 1 つでも現れたら取り直し。検出試験 (測定中にダミーの GPU プロセスを立てる) を通してから S3 を 2 case とも取り直す |
| M5 | **採用 (暫定修正 + #2g は既定化の前提のまま)**。`YAML_HARD_PATHS` に `mesh.scalarGradient` を追加し、明示 gg↔lsq の段を別区間にする (省略と明示 gg は別キーになる = 分けすぎる側に倒れるので安全)。起動順・バイナリ id・legacy の扱いは #2g の別 plan (未起票) で。#2g 起票は本 plan の Phase 2 前提として残す |
| m6 | **採用**。`methods/discretization.md`・`methods/gradient.md` を「既定 gg / opt-in lsq (実装済み・検証未完了) / cell は GG 固定」に更新、`procedures/solver-settings.md` に `mesh.scalarGradient` を追記 |
| m7 | **採用**。case/48 README に `run_0956_sglsq_s2_gg_chi0` の行、case/16 README に 0476 継続の結果、§5.1 に「床移動の真因」を独立項目 (F) として追加 |

## 仮説

- M2 は A になる見込み (S1 (1) で gg 経路は実装前後で 1 step の全配列が旧新ノイズ規則内、ただし長時間の床は未確認)。
- M3 は nSub の影響が小さく作用素差が主、という見込み (未検証)。
