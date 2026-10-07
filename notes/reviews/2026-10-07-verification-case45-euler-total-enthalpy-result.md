# codex レビュー: verification-case45-euler-total-enthalpy (result)

- **plan**: [`plans/active/verification-case45-euler-total-enthalpy.md`](../../plans/active/verification-case45-euler-total-enthalpy.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-07
- **commit**: `aca4edb1` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m2
- **focus**: E1〜E6 と E4V の結果 (§9) が事前登録を満たすか、§9「§8 の整理」(諮問で訂正済み) が証拠の範囲に収まっているか、accepted に移してよいか
- **extra**: `case/45.isobutane_m6_d155/euler_ref_de_sensitivity.py`, `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

E2「判別不能」、E4 段1「判別不能」、E4V「実務較正として合格」という区別は、保存済みデータと整合します。  
訂正後の §9「§8 の整理」は、原因未確定・適用範囲限定という解釈で妥当です。移行前に以下の記録を修正してください。

1. **Minor — 未確定事項の扱いを §5.1 に残す必要がある。**  
   根拠: [残作業表](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:55)では result レビューだけが未完了ですが、[§9](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:212)には原因候補の未分離、起動過渡との因果関係、初期線への離散化誤差などが未確定として残っています。  
   **対案:** §5.1 に「今回の完了条件外として保留する事項」と再調査条件を記載し、MOC 採否は対応 plan に引き継ぐ旨を明示してください。追加 CFD を今回の完了条件にする必要はありません。

2. **Minor — 評価窓全体の範囲と最終時点の値が混在している。**  
   根拠: [§9①(b)](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:211)の「本段42000〜54000で −234.9〜+292.6 K」は、実際には `run_0161_euler_t0cluster_g1` の **54000時点**の値です。`_band_ab/euler_t0_e2_eval.json` の同じスロート壁域では、53000時点に **+294.462559 K** があります。  
   **対案:** 「最終54000時点で −234.9〜+292.6 K」と限定するか、窓全体なら「−234.9〜+294.5 K」に訂正してください。E2 の判別不能という結論は変わりません。

検証結果の照合は次のとおりです。以下の run はすべて `case/45.isobutane_m6_d155/` 配下です。

| 項目 | 確認結果 |
|---|---|
| E1 | 保存記録の IC 偏差約0.011 K、soft 後の数百 K の偏差、B→C 最大1 ULP は記述と整合。soft 内の初期化と反復途中は未分離という限定が適切。 |
| E2 | `run_0161`・`run_0162` の記録はともに **NOT CONVERGED (stalled/plateau)**。不合格列は停滞のみ。A の時間幅条件未達による「判別不能」は正しい。 |
| E3 | `mesh_euler` の独立読込・既定値・移行エラーは §4 と整合。読み取りだけで実行できる既存テスト部分では、NS の設定・生成座標・接続の維持を確認。既存失敗3件の免除も明記済み。 |
| E4／E4V | `run_0163` の全温幅 **0.1065316824 K** は不合格。`run_0164` は最大幅 **0.0689320056 K**、出口平均 **5.9999996762**、最大偏差 **3.0574e−6**。保存系列から判定を再現。出口の記録は両窓 **STEADY**、残差は **NOT CONVERGED**であり、実務較正限定の採用は登録どおり。 |
| E5 | 保存された13点から幅を再計算し、N0 **3.7518e−5**、N1 **1.5851e−5**、N2 **4.0780e−5**。全件 **PASS**。これは参照時点への感度の合格であり、収束の証明ではない。 |
| E6 | 別ワークツリーの `run_0062_euler_wallfit_fit_r1_ext6k/res_6000.h5` を直接再計算。偏差 **−0.058304〜+0.072583 K**、絶対偏差0.1 K超・1 K超とも0節点を再現。 |

`run_0174` は全温の13枚窓の幅 **0.126446 K** が未達です。「E5 の参照感度は合格、較正場としては採用しない」という訂正を支持します。

`main...HEAD` の差分から本 plan 関連の実装を確認し、`methods/design/overview.md` と手順書の格子分離・現行較正値も照合しました。新規 run の残差 CSV・HDF5 は手元にないため、その一次出力の再検査と全テスト再実行までは行っていません。E1/E2/E4/E5 は保存済み評価記録による照合です。

**推奨は、上記2点を修正し、result レビュー記録・`status`・`plans/README.md` を同期して accepted に移すことです。** 原因同定や MOC 生産採用まで完了したとは扱わないでください。ファイルは変更していません。

指摘数: Critical 0 / Major 0 / Minor 2
