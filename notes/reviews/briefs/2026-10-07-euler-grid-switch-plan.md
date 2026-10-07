# 諮問: Euler の格子を壁に寄せない配点に切り替え、出口較正と MOC の V5 をやり直す — 実装方法・手順・順序 (§4・§6 の新設)

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 1。ユーザから「時間がかかりすぎ」との声があるので、1 回で §4・§6 の骨子を決めたい (後で plan 段レビューは回す)。
plan: `plans/active/verification-case45-euler-total-enthalpy.md` (§9 の 2026-10-07 の E2 とユーザ決定)、`plans/active/discretization-moc-axis-limit-and-corrector.md` (§6 V5・V5b、§9)、`plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md` (§6 U4)、`plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md` (§5.1 #11f・#11h、出口較正の経緯)。HEAD 7fb7041e。

## 観測事実

- **E2** (`case/45.isobutane_m6_d155/_band_ab/euler_t0_e2_eval.json`): 同じ壁 (単調壁 + MOC の新しい軸処理)・2000 × 97・等エントロピー IC・soft 3000 + 本段 54000 で、半径方向の配点だけを変えた。B (全域で第 1 間隔の比 0.005、軸側 0.0187 r_w) は窓 42000〜54000 の全時点で全領域の \|T₀ − 1600\| ≤ 0.103 K、時間の幅 ≤ 0.1 K。A (G1: 1.3e-5・スロート 4.5e-6、軸側 0.0888〜0.1000 r_w) はスロートで > 182 K、時間変動が大きい。登録の判定は判別不能 (A が整定の前提を満たさない)。両腕とも本段 NOT CONVERGED stalled/plateau、NaN なし。
- **ユーザ決定** (2026-10-07): Euler の格子を壁に寄せない配点 (B と同じ) に切り替える。NS は今の G1 の壁際のまま。G1 の Euler で決めた出口較正 (`Md_moc_offset` +3.770e-4、run_0113+0114) と MOC の V5 をやり直す。単調壁の E′ の採用は取り消さない。
- **今のコード**: Euler の `prepare` と NS の `prepare_ns` は同じ `mesh` ブロックを読む (`runner_axismach.mesh_params`、2026-10-06 の諮問で「Euler 経路が鍵を黙って無視していた」のを直した経緯)。Euler の問題 YAML (`problem_d155_euler_pin_G1_recal*.yaml`) は NS と同じ G1 の鍵 (wall_first_frac 1.3e-5・wall_first_frac_throat 4.5e-6・ブレンド) を明示している。
- **出口較正の経緯** (accepted plan §5.1 #11f): 粗い Euler 格子 (1100 × 65) の較正 (−4.16e-4) を、生産 NS 格子と同じ格子パラメータの Euler (G1) でやり直して +3.770e-4 にした (出口コア M 5.999207 → 6)。#11h: 出口コアの軸付近の量は半径方向の配置に敏感。B の格子は軸側の間隔が G1 (0.0888) より細かい (0.0187)。
- **MOC の V5・V5b** (G1 の Euler): 7 本すべてが量の準定常の前提を満たさず保留。延長しても量が動き続けた。同じ壁で IC により P 傾きが 0.05 違った。
- **NS 側の予定**: 上流の多項式化の U4 (今の MOC に固定して `poly` を評価、cross-mesh 移送) と、MOC の V5′ (その後に MOC の変更だけ)。どちらも出口較正の値を使う。
- **好みと制約**: ユーザは「標準と決めた方式はコードの既定にする」方針 (過去の YAML の再現のために旧既定を残さない)。AWS の Euler G1 級は 2 本並列で本段 54000 が約 10〜15 分。NS は数時間級。

## 問い

1. **配点の切り替えの実装**: (a) Euler の問題 YAML の鍵を書き換える、(b) Euler 専用の鍵 (例 `mesh_euler:`) を設け、`prepare` は既定で全域 0.005 の等比にする (NS の `mesh` の壁際の鍵を Euler で黙って使わない・無視もしない形)、(c) 別案。標準を既定にする方針との整合。
2. **出口較正のやり直し** (E4 として登録): 格子 (ni 2000 × nj 97 の配点だけ変えたものでよいか、軸付近の配置の違いを NS とどう整合させるか)、IC (G1 の収束場は異常があるので使わない? 等エントロピー IC から段階起動?)、長さ、出口コア M の評価位置と準定常の前提 (E2 で得た絶対の幅の考え方をどう使うか)、較正の式 (#11f と同じ)。
3. **MOC の V5 のやり直し** (V5d として登録): 腕 B (今の単調壁・legacy MOC) と 腕 M (単調壁・analytic + converge) を新しい格子で。IC の扱い、本数、長さ、準定常の前提 (V5 の追加登録の条件のままでよいか)、出口較正の値 (今の値か、E4 の新しい値か)。
4. **順序**: E4・V5d と NS (U4・V5′) の順序。出口較正が変わると設計壁が変わるので、NS の評価をどこで回すと最少の NS の本数で済むか (U4 と V5′ を分けて評価するという前回の諮問の方針は維持したい)。
5. 各 §6 の判定の骨子 (事前登録の文言の要点)。

## 読んでよいもの

上記 4 つの plan、`case/45.isobutane_m6_d155/{euler_t0_e2.py, euler_t0_e2_eval.py, moc_v5_euler.py, moc_v5_euler_eval.py, README.md}`、`_band_ab/{euler_t0_e2_eval.json, moc_v5_euler_eval.json, moc_v5b_ext_eval.json}`、`design/forge_design/evaluate/runner_axismach.py` (`mesh_params`・`prepare`・`prepare_ns`)、`design/forge_design/meshing/mesh2d.py`、`procedures/nozzle-design-workflow.md`。編集は禁止。
