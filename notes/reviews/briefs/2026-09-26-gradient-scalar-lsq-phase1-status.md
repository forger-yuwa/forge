# 諮問ブリーフ: スカラー LSQ Phase 1 の現状解釈 — result 2 回目に回す前 (2026-09-26)

plan: `plans/active/gradient-scalar-lsq-unification.md` (§5.1 #5, #5a–#5l, §6, §6.1)。
前回: result 1 回目 NO-GO (`notes/reviews/2026-09-26-gradient-scalar-lsq-unification-result.md`)、採否 (`notes/reviews/2026-09-26-gradient-scalar-lsq-result1-disposition-diagnose.md`)。

## 決めてほしいこと

1. 下の状態で、result 2 回目 (`codex_review.py --stage result`) に回してよいか。回すなら「Phase 1 の判定として何を主張するか」を 1 文で (例: 「opt-in lsq は gg と比べて物理量を変えず、性能劣化 ≤ 2 %、既定化は保留」)。
2. FCT smoke (#5h = B) と床移動の真因 (#5l、F) を、Phase 1 の完了条件から外して Phase 2 の前提 (または別 plan) に移してよいか。外すなら plan のどこにどう書くか。
3. Phase 1 を `done` にせず `in_progress` のまま Phase 2 の前提作業 (#2g 起票) に進むべきか。

## 読んでよいファイル

plan 全文、`case/09.Taylor-Green/_g0_lsq_seam/{M2_case48.txt,M3_fct.txt,FLOOR_RECHECK.txt,S3_case48.txt,S3_case39.txt}`、4 case README の `run_095x_sglsq_*` 行。h5 は読まない。

## 観測事実 (result 1 以降)

- #5f (M1): `check_floor_ratio.py` の誤合格 3 種を修正 (自己試験 9 例)。修正後ツールで S2 の全床比判定を取り直し、判定は不変 (case/39 のみ両方 PASS、判定不能なし)。
- #5g (M2) = **A**: 実装直前版 `36d8ba03` と現行 `f99f236d` gg を case/48 で同じ 24000 + restart + 24000 手順で比較。延長段の残差床比 0.989–1.092、Cf・q_w・δ*・θ (3 station、step 12000–24000 平均) の差 ≤ 0.023 %、旧版も STEADY・起点床を外れる (rms_ro 1.83×)。
- #5a: chi 0 でも case/48 の起点床に戻らない (B、chi 単独説明は棄却)。床移動の真因は未特定 (#5l、F)。
- #5h (M3) = **B**: FCT smoke の nSub 15→30 の自己変化が lsq−gg 差の約 2 倍。`check_passive_budget --mode fct` は gg/lsq × nSub 15/30 の 4 本とも同値の remainder 8.45e-2 で FAIL (閉合 3e-10・総量照合 4e-9 は通過)。
- #5i (M4): S3 の監視を solver PID と GPU PID の照合に変更。取り直し結果と競合検出試験は下に追記。
- #5j (M5): `stage_manifest.py` に `mesh.scalarGradient` の hard キー、YAML 解析不能時の非連結、回帰試験 4 例。
- #5k/#5l: 文書 (methods・procedures) と run 索引を更新。
- S2 の物理ゲートと lsq−gg 上限 (FCT smoke 以外): case/48 ≤ 0.016 %、case/40 Δη 0.000 %・壁温 0.62 K、0476 Y0 1.2e-5・Xi 3.1e-4 (継続後 1.6e-5/3.4e-4)、0482 onset 差 0・壁 p/p0 lsq−gg 0.003 %・vs 0482 0.16 %、case/39 全量 STEADY。

## S3 取り直しと競合検出試験 (追記)

- 競合検出試験: 測定 solver が GPU に乗った 5 s 後に別名バイナリ (`forge_36d8ba03`) の run を割り込ませた → `_s3_gg_a` で「solver 以外の GPU プロセス [33711]」を検出し取り直した (`S3_v2_driver.txt`・`S3_v2_inject.txt`)。
- 取り直し: case/48 比 lsq/gg **0.997** (中央値 gg 1.541・lsq 1.537 ms/step)、case/39 **1.018** (2.549・2.596)。どちらも上限 1.05 以内 (`S3_case48.txt`・`S3_case39.txt`)。旧監視版の結果 (1.010 / 1.017) と整合。ポーリング 0.5 s より短い競合は排除できない。

## 仮説

- Phase 1 の実装 (opt-in lsq) は、物理量・性能の面で gg と区別できる差を作っていない。未合格の収束ゲートは起点作成以降の共通の変化で、実装直前版でも再現する (case/48 で確認、他 case は未確認)。
- FCT smoke の保存収支 FAIL は dual-time 手順側 (起点の定常場からの切替) の問題で、スカラー勾配とは独立 (gg でも同値)。
