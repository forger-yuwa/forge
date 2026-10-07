# codex 諮問 (diagnose): cfdpin-result-interpretation

- **brief**: [`notes/reviews/briefs/2026-10-06-cfdpin-result-interpretation.md`](../../notes/reviews/briefs/2026-10-06-cfdpin-result-interpretation.md)
- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **date**: 2026-10-06
- **commit**: `e2210e0e` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.7 min, rc=0
- **結論**: **`active` を維持し、実効入力検査・報告・完了条件を是正した版を result 再レビューへ回す。追加 CFD は、その前提作業として要求しない。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 論点・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **現状のまま `done`・accepted：却下** | [plan:137](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:137) は「§6 を満たす」を要求するが、V3′ は保留、最終 dry の波は `DRIFTING`。ユーザ決定はトレードオフ採用と「未達のまま④へ」であり、判定の合格化ではない。**実装完了・採用判断・性能検証の結果を分けて完了条件を書き直す。** |
| **Major** | **前回 M2「入力契約」の対応完了：要再検証** | [cfd_initial_line.py:100](/home/sano/work/forge-integ-1005/design/forge_design/feedback/cfd_initial_line.py:100)、[同:137](/home/sano/work/forge-integ-1005/design/forge_design/feedback/cfd_initial_line.py:137) はガス条件を `prepare_info` から照合し、入口 BC の実効全温を照合しない。読み取り内容だけを変更した試験で、入口 `Tt=1600→1500 K` でも provider は受理し、`gas` を照合済みと記録、熱力学ハッシュも不変だった。**実効 BC・組成・熱力学条件とメタデータ／利用側を照合し、不整合を拒否する。** 今回の実 run が誤設定だったという指摘ではない。 |
| **Major** | **⑤報告の完了扱い：却下** | [nozzle_report.py:208](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:208) の波は `100(M/Md−1)` の平滑残差だが、[build_pptx.py:125](/home/sano/work/forge-integ-1005/design/forge_design/report/build_pptx.py:125) は「圧力波」と表示する。また [同:149](/home/sano/work/forge-integ-1005/design/forge_design/report/build_pptx.py:149) は末尾幅だけを掲載する。実際の両 PPTX に `DRIFTING`・`STEADY`・「未達」の記載はない。**Mach 波へ改称し、量別 VERDICT・判定区間・未達継続の決定を掲載して再生成する。** |
| **Major** | **格子差からの原因確定・精度保証：却下** | [plan:182](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:182) の「出口 M もほぼ格子収束」は最終場比較、[同:185](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:185) の「Euler 分約2割、残りは境界層側」は異なる方程式・壁条件の差を加算的に分解した解釈。**観測した格子感度までに限定し、残差分の帰属は未分離とする。** δ_E の旧形状上の格子ゲート合格を、最終形状・凍結線・波の格子独立性に拡張しない。 |
| **Minor** | **§4・§6 の整理：採用** | [plan:61](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:61) と [methods:355](/home/sano/work/forge-integ-1005/methods/design/overview.md:355) は旧 offset のまま。**現行値 `+3.770e-4`・固定 Euler 参照 `run_0114`・生産格子を本文に集約**し、旧条件は履歴と明示する。§1 の「非劣化を確認」も実際の採用理由へ訂正し、§5.1 は完了／未達／移管を更新する。 |

結論: **`active` を維持し、実効入力検査・報告・完了条件を是正した版を result 再レビューへ回す。追加 CFD は、その前提作業として要求しない。**

第 1 仮説: 完了を妨げている主因は、限定付きの採用判断と検証合格の混同、および入力検査・報告への反映不足である。  
確度: **高**

根拠: 実データは `/home/sano/work/forge/` 側で確認した。以下はリポジトリ相対パス。

- `case/45.isobutane_m6_d155/run_0117_ns_recal_final_ext/`  
  `quantities_series.csv` を正式ツールで再判定。末尾5枚、**step 40000〜60000**：
  - 出口コア M **5.99888699 — STEADY**。登録範囲内。
  - δ_E/δ_C **0.99983458 — STEADY**。
  - オーバーシュート **0.0079397 % — STEADY**。
  - Mach 波 **0.0065105 % — DRIFTING**、drift **5.4 %/tail**。値の範囲は **0.0061323〜0.0066641 %**。
- `case/45.isobutane_m6_d155/run_0118_ns_recal_final_cond/`  
  `cond_series.csv` を登録の末尾5枚、**step 14000〜18000**で再判定：onset・S_max・出口コア g・出口コア M は **すべて STEADY**。出口コア M の末尾平均は **5.9864071**。これは凝縮評価の登録合格であり、dry の M6 許容範囲への合格ではない。
- 両 run の保存済み収束判定は **`NOT CONVERGED (stalled/plateau)`**。凝縮関連列に `RISING` はない。
- 全保存場の `VALUE/*` を確認し、NaN/Inf は0。保存済みメッシュ品質は **`VERDICT: PASS (AR<=5000, skew<=0.90)`**。
- 最終報告の正式壁解像は **`VERDICT: PASS (壁解像)`**、超過面積 **3.5 %**、最大 y₁⁺ 約 **1.505**。これは登録した超過許容5 %での合格であり、全壁で y₁⁺≤1 の意味ではない。

反証条件: 入力不整合を拒否する検査が成立し、提出報告と完了条件が量別 VERDICT・例外決定・検証範囲に整合していることが確認できれば、この完了阻害は解消する。

第 2・第 3 仮説: 追加しない。波の物理的原因は、今回の証拠から確定できない。

判別 A/B: **入口 BC の `Tt` だけを変える provider の拒否試験、CFD 0 step。** A＝1600 K、B＝1500 K。他のメタデータ・snapshot・利用側条件は固定する。  
→ **A受理・B拒否なら**実効入力契約は成立。**両方受理なら**メタデータの一致を実効条件の一致と誤認している。今回の読み取り差し替え試験は**両方受理**だった。修正後は同じ試験を回帰条件とする。

やらない方がよいこと: 波の無条件延長、`--drift` の事後緩和、較正の追加、凍結線の取り直し、未達を合格へ書き換えること。凍結線の格子依存・Euler G2 は、限定付き実装完了のために直ちに追加する必要はない。

呼び出し側の前提への異議:

- 「探索の粗さは除外済み」は強すぎる。#11g は登録上**中間・保留**。細標本でも `DRIFTING` が残り、主因説は支持されなかった、まで。
- 「凝縮④合格」は採用する。ただし **登録4量の準定常条件と残差非上昇条件の合格**に限定する。
- **限定付きの実装完了として閉じる方針は採用可能**。上記修正と result 再レビュー後、§8・§9 には次の文言を推奨する。

> 本 plan は CFD ピン機能と case/45 の指定レシピの実装・評価を完了した。Hall 対比の全量非劣化、流れ場の残差収束、最終設計の格子独立性を証明したものではない。V1 点上ゲート不合格、V2 許容超過、V3 不採用、V3′ 保留の履歴を保持し、採用根拠は2026-10-05のユーザによるトレードオフ判断とする。最終 dry の Mach 波は閾値内だが DRIFTING であり、2026-10-06のユーザ決定により未達のまま凝縮評価へ進んだ。dry・凝縮とも残差は NOT CONVERGED。凝縮の登録4量は指定末尾区間で STEADY。凍結線の格子依存、Euler 較正の G2、波の準定常達成は未確認として残す。

不足情報: ローカル配置には最終2 run の `residual_history.csv` がなく、残差判定は保存済み VERDICT の確認に限る。result 再レビューには元 CSV と実行時の判定条件を揃えること。未達・未確認事項を別 plan へ移管するなら、その節と担当も明記すること。

**ファイル変更・forge 実行は行っていない。提案は plan 未反映。**
