# codex 諮問 (diagnose): deltastar-ns-rebuild-route

- **brief**: [`notes/reviews/briefs/2026-10-13-deltastar-ns-rebuild-route.md`](../../notes/reviews/briefs/2026-10-13-deltastar-ns-rebuild-route.md)
- **plan**: [`plans/active/tooling-nozzle-deltastar-band-root-zone.md`](../../plans/active/tooling-nozzle-deltastar-band-root-zone.md)
- **date**: 2026-10-13
- **commit**: `cf35b7c1` (feature/nozzle-deliverables)
- **codex**: effort `high`, 1.7 min, rc=0
- **結論**: **Euler 参照の受入れを「完走・NaN なし」から「収束 PASS・対象量 STEADY」へ変更し、まず `run_0012` の反復数だけを変える A/B で、δ_E が参照場の過渡に依存していないか確認する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **Euler 参照の受入れを「完走・NaN なし」から「収束 PASS・対象量 STEADY」へ変更し、まず `run_0012` の反復数だけを変える A/B で、δ_E が参照場の過渡に依存していないか確認する。**

第 1 仮説: **【Major・要再検証】未収束の Euler 参照から作った δ_E が壁に転写され、出口 M の誤差を残す、または改善理由を取り違える。** 確度: **中**。受入れ条件の欠落は確認できるが、新しい Euler が実際に未収束になるかは未確認。

  根拠:
  
  - `plans/active/tooling-nozzle-deltastar-band-root-zone.md:121` は、補正の入力となる Euler の合格条件を「完走・NaN なし」とし、収束 VERDICT は記録だけとしている。
  - 同 plan `:98–100` では、その Euler から抽出した δ_E を寸法・壁の入力にする。`design/forge_design/feedback/deltastar_loop.py:343–349` の処理は δ_E の尺度換算であり、入力に含まれる参照場の誤差を除去しない。
  - 同 plan `:147–150` の最終 NS の判定だけでは、この問題を防げない。**未収束の Euler の流量を定数として分母に置いても、NS が定常化すれば流量比は STEADY になり得る。**
  - 前回レビュー `notes/reviews/2026-10-13-deltastar-band-root-zone-diagnose.md:41` も、旧参照 `run_0003_euler_36k` の `NOT CONVERGED` を報告している。ただし、今回その判定原本は閲覧していない。

  対案: **参照場の受入れゲートを追加する変更を採用**する。旧設計の `run_0012` と新設計の `run_0013` の双方について、同一の方程式・BC・離散化の区間で `check_convergence --segment: PASS`、出口コア M と流量について `check_quasisteady: STEADY` を要求する。

  反証条件: 下記 A/B の双方が判定を通り、参照の反復数を増やしても出口コア M・流量・抽出 δ_E が事前の識別閾値以内であること。この場合、参照場の過渡を今回の主因候補から外す。

第 2・第 3 仮説: 挙げない。現資料では、補正 1 回・Re⁻⁰·² 換算・cross-mesh IC のどれが支配的かを順位づける証拠がない。

判別 A/B: **変えるのは旧設計 Euler の本段反復数だけ**。A＝計画の 36000 step、B＝同じ計算を別 run に複製して同一設定で追加 12000 step。NS 入力は `run_0010` に固定し、両 Euler から同じ修正済み抽出を行う。

  見る量は両判定ツールの VERDICT、出口コア M、流量、共通座標での δ_E 全分布。識別閾値は事前に **|ΔM| ≤ 0.00025、流量の相対差 ≤ 1e−4、max|Δδ_E| ≤ 1e−4 r_t** と置く。これは診断用の閾値であり、最終 NS の M 精度を保証する換算ではない。

  → **A・B とも PASS／STEADY かつ差が閾値内なら**、参照の過渡が支配的という第 1 仮説を棄却する。  
  → **A が判定不合格、または延長で閾値を超えて変わるなら**、36000 step の参照をそのまま採用する方針を棄却する。B も不合格なら判別未完了であり、過渡同士の差を補正の物理的効果とは解釈しない。

やらない方がよいこと: 参照場の受入れを確定する前に、出口 M のずれを「補正 1 回では不足」と断定して、補正反復や係数調整を始めること。

呼び出し側の前提への異議: **カードから実寸を保って作り直す判断は採用**する。`case/70.m5_air_wt_demo/_sizing_ab/sizing_ab.json:30,57,63` は、A が旧寸法を再現し、拡大部長の誤差が −0.485 mm → +4.422 mm になることを確認できる。ただし、この寸法 A/B は流れの収束や補正の物理的妥当性を証明しない。また、最終 NS の V5 だけでは、入力側の参照場に対する前回の不足情報は解消しない。

不足情報: `run_0010` の対象量の準定常 VERDICT、および `run_0012`・`run_0013` の収束／準定常 VERDICT と判定区間。閲覧制限に従い独立確認していない。**plan 未反映**。変更先は §4.3・§5.1 #5・§6 V3／V5 を想定し、反映は呼び出し側が行う。
