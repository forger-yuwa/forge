# codex 諮問 (diagnose): wall-temperature-profile-result

- **brief**: [`notes/reviews/briefs/2026-10-13-wall-temperature-profile-result.md`](../../notes/reviews/briefs/2026-10-13-wall-temperature-profile-result.md)
- **plan**: [`plans/active/tooling-nozzle-wall-temperature-profile.md`](../../plans/active/tooling-nozzle-wall-temperature-profile.md)
- **date**: 2026-10-13
- **commit**: `4ea3626d` (feature/nozzle-deliverables)
- **codex**: effort `high`, 1.2 min, rc=0
- **結論**: accepted へ移す前に、壁温分布だけを変えた prepare-only の A/B で、積分法の δ_r の差が生成された物理壁まで伝わることを確認する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: accepted へ移す前に、壁温分布だけを変えた prepare-only の A/B で、積分法の δ_r の差が生成された物理壁まで伝わることを確認する。

第 1 仮説: **[Major・要再検証] 現在の検証は、積分法が温度表を読むところまでで止まり、初期物理壁への反映漏れを検出できない。** 実装不具合があるとは未確認。確度: 高（検証の不足について）
  根拠: `design/tests/run_wall_temperature_tests.py:320` は `res["Tw"]`、同 `:370` は `thermal_bc` の記録を照合している。一方、`design/forge_design/evaluate/runner_axismach.py:1183` には δ_r の平滑化、同 `:1320` には物理壁の構築という後続処理がある。提示された V1b の 5.7e-14 K と V2 の 1.5e-5 K は、この後続処理で補正が失われても成立し得る。
  反証条件: 同じ設計壁・実効 `scale` で、温度分布による δ_r の差と、生成物理壁の半径差が、平滑化・壁形状近似の許容誤差内で対応する既存の試験結果が示されること。

第 2・第 3 仮説: 無し。

判別 A/B: A は既存の温度表、B はスロート温度を保ち下流側の表の温度だけを +100 K。設計壁・`scale`・補正設定を固定し、寸法逆算は介さず `--init-integral --prepare-only` を各1回、NS は0 step。積分法の δ_r から独立に求めた期待半径差 Δr_ref と、保存された物理壁の半径差 Δr_out を、同じ実寸 x で比較する。比較区間は径方向補正をそのまま適用する拡大部とし、ε = max（1e-8 m、両ケースの壁形状近似の誤差上限の和）を事前登録する。
  → **判定A**: max|Δr_ref| > 10ε、かつ max|Δr_out − Δr_ref| ≤ ε なら、反映漏れの疑いを除外できる。
  → **判定B**: 期待差が十分あるのに不一致なら、完了判断を保留し、δ_r から物理壁までの経路を調べる。期待差自体が小さければ判別不能であり、合格にしない。

やらない方がよいこと: この穴を埋める目的で V2 の NS を延長すること。壁ノードの温度拘束を再確認しても、初期物理壁の生成経路は検証できない。

呼び出し側の前提への異議: 「V1〜V4 合格」と「積分法の初期物理壁まで配管の失敗を塞いだ」は分けるべき。V2 を温度拘束の検査に限定する解釈は採用する。`case/70.m5_air_wt_demo/README.md:25` の main **NOT CONVERGED** は、その限定された合否を否定する根拠にはしない。

不足情報: 上記の δ_r 差と物理壁差の照合値、および壁形状近似の誤差上限。run の生データと VERDICT 原本は指定に従い未読で、数値は提示値・README 記載として扱った。plan 未反映（読み取り専用の諮問）。
