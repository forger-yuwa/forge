# 諮問ブリーフ: Graetz result レビュー (NO-GO C0/M4/m1) の採否と閉鎖文面の確認 (2026-09-29)

AGENTS.md エスカレーション条件 **5** (Major の採否)。plan `plans/active/boundary-cht-axisymmetric-graetz.md` §6.1 の `result` 行に採否、§5.1 #6g に訂正を書いた。
レビュー: `notes/reviews/2026-09-29-boundary-cht-axisymmetric-graetz-result.md`。

## 読んでよいファイル
plan の §6.1 result 行、§5.1 #6f・#6g の行、`case/63.graetz_cht/temp_reproduce.py` (判定部)、`temp_reproduce_x0_uncertainty.py`、
`MANIFEST_primary_data.sha256` (先頭 5 行)、`run_0014_g2_dT0_r32/TEMP_REPRODUCE.txt`

## 採否 (私の案、ユーザ決定を含む)
- M1: **ユーザ決定 (2026-09-29) で例外として受け入れて閉じる**。FAIL は保持、主張を限定 (§6.1 の文面)。status を done にして accepted へ移す。
- M2/M3: 採用。x=0 の不確かさ 2.84e-3 K > 1.67e-3 K (N_r 16/32/64/128 細分化、N_r=32 の全内部節点で比較) → 温度再現 A/B は判定不能。参考として A/B の幅の差 0.135 K は不確かさの 48 倍。
- M4: 採用。最終場・圧縮残差を手元に回収しハッシュ台帳 (181 ファイル)。手元で check_convergence 14 本・eval snap 3 対を再実行して一致。中間スナップショットと界面ログは復元不能と明記。
- m5: 採用。methods/boundary.md の検証範囲に「case/63: 流れのある軸対称 CHT、差し引き Nu 0.82 % 以内、前提検査 V-g1 FAIL (例外受け入れ)」を閉鎖時に追記。

## 諮りたいこと
1. 上の採否で閉鎖してよいか (ユーザ決定の例外受け入れを前提に)。閉鎖に足りない記録は何か。
2. 閉鎖時の主張文 (plan §1 の完了状態・methods の 1 段落) に含めるべき限定・落とすべき表現。
推奨を 1 つに絞ること。
