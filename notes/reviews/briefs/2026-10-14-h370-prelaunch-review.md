# 諮問ブリーフ: V-h370 の投入前のレビュー (判定器の修正 #36g と出口の列の被覆の許容) (2026-10-14)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1・3・4: 事前登録の判定の具体化の確認、新しい run の前)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §6 V-h370・§5.1 #36e〜#36g
- 前の諮問: `notes/reviews/2026-10-14-h370-judge-criteria-diagnose.md` (Major 4、全件採用 → #36g)
- 実装: commit `fdf3b728` (このレビューの作業ツリーはその commit のまま)。手元の試験: `case/45.isobutane_m6_d155/test_h370_launch.py` 111 項目・`design/tests/run_h370_prelaunch_tests.py` 77 項目・`run_stop_rule_tests.py`・`run_euler_h_ref_tests.py` ALL PASS。forge はまだ回していない。

## 観測事実

1. **dry-run (AWS B、b5b53a1b)**: Euler・NS とも品質 PASS、壁の保存・復元・壁の出力の設定も通った。変換後の座標で δ_in 最大 7.964e−4 mm (許容 1e−3 mm 以内) だが、NS の最後の列 x = 95.29740274407894 が Euler 参照の端 95.29740274239427 を 1.68e−9 r_t 越えて被覆 False (n_uncovered 1)。**原因**: Euler の格子の msh を既定の 10 桁で書くため (出口 7.306010214129155 m → 7.306010214 m、÷ r_t で観測値と一致)。NS は 17 桁。
2. **被覆の許容 (#36g)**: `deltastar.euler_cover_range`/`euler_covers` に 1 本化。`remap_edge="missing"` の経路だけ、端の外側に tol = 1e−4 × 端の列の間隔 (生産の Euler 2000 × 97・throat_refine 5 の出口で 6.46e−6 r_t ≈ 5e−4 mm)。許容の内側で端を越える列は端の列の値で評価。10 桁の丸めは全列で最大 6.5e−9 r_t。旧い run_0174 の 0.07 r_t の不足は欠測のまま。既定の clamp の経路は今のまま。適用: 抽出・`delta_in_check`・停止規則の帯と δ の記録・受け渡しの関門・判定器の A/B。
3. **判定器 (#36g、前の諮問の Major 4 への対応)**:
   - Euler の 6 点の M の準定常は「予備の資格」(`--euler-prequalify`、`h370_euler_prequalify.json`)。最終の資格は σ_ref に繋ぎ、判定不能・未整定なら `update_conditions_met` を偽にする。
   - σ_ref: 2 窓の全 17 出力 × 代わりの参照 4 つ (Euler の末尾の窓 [40000, 60000] の 40000・45000・50000・55000、固定の 60000 と合わせて 5)。代表の選び方は予備の資格の時点で固定。全 20 出力にすると抽出 361 回 ≈ 11.6 h、代表は 89 回 ≈ 2.9 h。代表の間の変化を見るため、末尾の窓の全 21 出力の sha と場の差 max|ρU_x − ρU_x,pin|/max を記録だけで残す。停止の出力 × 37000・49000 の検査は「スクリーニング」として残す。
   - 比べる列は座標と補間できる範囲から先に固定 (`fixed_columns`)。その中の非有限・欠測・held はその出力の比較を合格にしない。
   - 2 窓の 17 出力を列挙して照合し、欠けて再抽出で補えなければ判定不能 (`two_window_check` も行の欠けで不合格)。
   - 分類: 結果 A = max D + σ_ref ≤ 0.01 mm かつ全出力で合否が同じ、結果 B = min D − σ_ref > 0.01 mm、または合否の違いが代わりの参照でも保たれる、それ以外は判定不能。
   - 前の諮問の反例 3 つ: 中央 1 列の欠測 → 結果 B (合否の違いが 4 つの代わりの参照で保たれる。諮問の予想は判定不能)、停止の 1 行だけ → 判定不能、途中の 1 出力だけ 0.030 mm → σ_ref 0.030 で判定不能・最終の資格なし。

## 実装者が判断を求めた点 (呼び出し側は受け入れる案)
- (i) σ_ref の代表 4 出力 (全 20 にしない)。
- (ii) 腕 A の比べる列の中の欠測を「合否の違い → 結果 B」と読む (代わりの参照で保たれることが条件)。
- (iii) σ_ref が 0.01 mm を超えても A/B の分類を変えなければ資格ありとする (§6 の「判定を左右するなら資格未達」の文言どおり)。σ_ref ≤ 0.01 mm を別の条件として足すかは方針。
- (iv) A/B が分類に至らない (旧い参照の抽出の失敗など) ときも更新に進めない (安全側)。
- (v) 被覆の許容 1e−4 × 端の列の間隔。代わりに Euler の問題に `msh_digits: 17` を入れて丸めをなくす案もあるが、Euler の格子の sha が変わる。

## 問い
1. 上の修正で、前の諮問の Major 4 が解消したか。投入前に直すべき Critical / Major が残っていないか。
2. (i)〜(v) の採否。
3. 被覆の許容が「範囲の不足を隠す外挿」にならないか (旧い参照の 0.07 r_t の不足は欠測のまま、を含めて)。

## 読んでよいファイル
- `case/45.isobutane_m6_d155/h370_judge.py`・`h370_run.py`・`test_h370_launch.py`
- `design/forge_design/metrics/deltastar.py` (`euler_cover_range`・`euler_covers`・`_load_structured`)・`delta_in_check.py`、`design/forge_design/evaluate/stop_rule.py` (`required_window_steps`)・`design/forge_design/feedback/deltastar_loop.py` (`two_window_check`・`coverage_check`・`handoff_gate`)
- `design/tests/run_h370_prelaunch_tests.py`
- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §6 V-h370・§5.1 #36e〜#36g
- `notes/reviews/2026-10-14-h370-judge-criteria-diagnose.md`
