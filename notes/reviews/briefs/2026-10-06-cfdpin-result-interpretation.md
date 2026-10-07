# 諮問: CFD ピン plan の result 段の解釈 (accepted にしてよいか、何を限定として書くか)

日付 2026-10-06。codex (diagnose)。エスカレーション条件 7 (result 段の解釈を確定する前)。
plan: `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` (§1〜§6、§5.1 #8〜#11g、§6.1 の全諮問・レビュー、§9 全履歴)。前回 result 段レビュー (NO-GO): `notes/reviews/2026-10-05-tooling-nozzle-cfd-pinned-initial-line-result.md`。

## 最終設計と判定 (要約; 詳細は §9 2026-10-05〜10-06)

- 生産レシピ: CFD ピン初期線 (run_0062 凍結) + 当てはめ壁 (joint) + 細分格子 (ni 2000 × nj 97、壁解像 PASS、格子ゲート: 粗→細 +1.08 % で FAIL → 第三水準で細→第三 +0.21 % 合格) + 出口較正を生産格子の Euler でやり直し (Md_moc_offset −4.16e-4 → +3.770e-4、Euler G1 E2 run_0113/0114 で 5.999998) + C2 (k_f 1.054129・r_t 76.6539 mm)。
- result NO-GO の対応 #8〜#12: 実装・試験済み (δ_E 取り違えの訂正、入力契約、pw_ramp、報告ツールの正式壁解像)。
- 最終 NS `case/45.isobutane_m6_d155/run_0116_ns_recal_final` + 延長 `run_0117_ns_recal_final_ext` (cfl 1、各 60000 step): 出口コア M 5.998887 STEADY (登録 6.000 ± 0.02 % 合格)、δ_E/δ_C 0.9998 STEADY、オーバーシュート η0.1 0.0079 % STEADY、出口半径 0.7749995 m、壁解像 PASS (3.6 %)、**波 η0.1 0.0065 % (閾値 0.01 % 内) は DRIFTING のまま** (延長後も; 標本間隔 A/B #11g で探索の粗さは原因でない、ゆっくり減少) → ユーザ決定「A」: 未達のまま記録して ④ へ。check_convergence は全 run NOT CONVERGED (plateau)。
- 凝縮 ④ `run_0118_ns_recal_final_cond`: 凝縮 4 量 STEADY (開始 x 57.94・S_max 16.86・出口コア g 3.13e-4・出口コア M 5.9864)、凝縮残差 plateau (RISING なし)。
- 報告: run_0117・run_0118 の nozzle_report + pptx。
- 途中の環境変更: AWS 共有バイナリ消失 → 同じソースを自分の worktree でビルド (sha256 6b47811b…、δ_E 同等性 −0.003 %)。
- 本 plan 外に残すもの: 凍結線 (run_0062、粗い Euler 格子) の格子依存は未確認; Euler 較正の格子依存 (G0→G1 −0.00079) は分かったが G2 は未実施; 波の準定常未達。

## 問い

Q1. この状態で plan を `done`・accepted にしてよいか (波 η0.1 の DRIFTING 未達・NOT CONVERGED plateau・凍結線の格子依存未確認を「限定」として書く形で)。それとも accepted 前に必須の作業があるか。
Q2. accepted にする場合、§8 完了条件と変更ログに書くべき限定事項の文言 (過大主張を避ける)。
Q3. result 段レビュー (codex) に回す前に plan 本文で直すべき箇所 (§4・§6 が #11a〜#11g の追加で読みにくい、など)。
