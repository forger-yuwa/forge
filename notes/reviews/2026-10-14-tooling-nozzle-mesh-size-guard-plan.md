# codex レビュー: tooling-nozzle-mesh-size-guard (plan)

- **plan**: [`plans/active/tooling-nozzle-mesh-size-guard.md`](../../plans/active/tooling-nozzle-mesh-size-guard.md)
- **stage**: `plan`
- **date**: 2026-10-14
- **commit**: `88675245` (feature/nozzle-cooled-wall-settling)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

課題の同定と、実効節点数による事前ガードは妥当です。既存の accepted plan は Euler 格子の分離までで、今回の規模制限とは重複しません。  
ただし、V1b の一般化と判定精度、要求カード生成中のガードに穴があり、現状のまま上限 2000 を既定化することは推奨しません。

1. **Major — V1b の「比 3 が合格 → 全ケースで上限 2000」は成立しない。**

   **根拠:** [対象 plan:93](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mesh-size-guard.md:93) は、case/70 の 1800→600 を NS 6000→Euler 2000 の根拠にしています。しかし軸方向の配点は、領域長・スロート細分・局所細分から決まります（[mesh2d.py:53](/home/sano/work/forge-deliv/design/forge_design/meshing/mesh2d.py:53)）。NS/Euler の点数比が同じでも、Euler の局所 Δx や形状の解像度は同じになりません。accepted plan も「2000×97 は全ケースの汎用既定の根拠ではない」と明記しています（[verification-case45-euler-total-enthalpy.md:47](/home/sano/work/forge-deliv/plans/accepted/verification-case45-euler-total-enthalpy.md:47)）。

   **対案:** V1b の結論を case/70 の条件内に限定してください。既定変更の検証には、実際に `mesh.ni > 2000` となる対象形状で、従来の Euler 格子と 2000 の直接比較を追加する。適用範囲を領域長・形状・配点条件で記録し、2000 を精度保証付きの普遍的上限として扱わないこと。解そのものの格子依存が原因なら、§4.1 の高次補間では解消しません。

2. **Major — 準定常ゲートが比較許容差に対応しておらず、時間変動を格子差と誤認し得る。**

   **根拠:** [対象 plan:91](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mesh-size-guard.md:91) は δ_E・流量比に 1e-4 の許容を設定しますが、準定常判定の窓・閾値を固定していません。参照する [qs_series.py:6](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/qs_series.py:6) の既定は drift 5%、osc 10% です。

   保存済み CSV に対して、`check_quasisteady.py --tail 0.4 --drift 5e-5 --osc 5e-5` を再実行した結果は以下でした。

   | run（`case/70.m5_air_wt_demo/` 配下） | 流量の VERDICT | 末尾40%の流量の相対最大最小差 |
   |---|---|---:|
   | `run_0013_euler_ns2` | `OSCILLATING` | 2.61e-4 |
   | `run_0010b_ns_c2b_cont60k` | `DRIFTING` | 6.40e-5 |

   `check_convergence.py` の再判定は両者とも **`NOT CONVERGED (stalled/plateau)`**。最終 HDF5 の `VALUE` 19 データセットには非有限値はありませんでした。これは新しい V1b の不合格を意味しませんが、「同じ段取りなら床に達する」「既定の STEADY なら比較できる」という前提を支持しません。

   また、出口 M・流量が落ち着いても局所の δ_E は保証されません。accepted plan には、Euler の延長だけで δ_E が最大 1.5e-4 r_t 変化し、帯位置も切り替わった記録があります（[tooling-nozzle-deltastar-band-root-zone.md:114](/home/sano/work/forge-deliv/plans/accepted/tooling-nozzle-deltastar-band-root-zone.md:114)）。

   **対案:** 各量について評価窓と時間変動の許容を事前登録し、格子差の許容より十分小さくしてください。δ_E も固定した NS 入力に対する複数の Euler 時刻から評価する。判定には既存の `check_quasisteady --series-csv` を使い、変動が許容を覆う場合は「判定不能」とする。未収束場の利用は条件付き感度評価として明記し、収束解の精度確認とは分けるべきです。

3. **Major — `prepare` のガードだけでは、要求カード生成中の大規模格子確保を防げない。**

   **根拠:** `estimate_mesh` は [wt_card.py:871](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:871) で `_mesh_ar_full` を呼び、必要なら繰り返します。その内部は [wt_card.py:716](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:716) で `generate_axisym_mesh` を直接実行します。NS が 6000×97 なら、`prepare` に到達する前に **582,000 節点**を確保します。§4.3 の「生成器の出力で事前提示」では、この確保より後になります。

   さらに、点数計算を単に `len(_radial_fracs_capped(...))` として共通化すると、[mesh2d.py:103](/home/sano/work/forge-deliv/design/forge_design/meshing/mesh2d.py:103) のリスト生成で、判定前に導出 nj 分のメモリを使います。実行確認では指定例の nj は **201、節点数は402,000**でした。

   **対案:** 共通ガードを `prepare`・`prepare_ns` に加え、要求カードの全格子品質見積もりの直前にも置いてください。超過時は点数の見積もりを返して停止し、全格子は生成しない。点数計算自身も巨大な径方向配列を作らず、超過を早期判定できる実装にする。V3 に「要求カード経由でも生成関数を呼ばない」と「極小 `axis_gap_frac` でも巨大確保しない」を追加してください。

4. **Minor — 明示した Euler `ni` の優先順位が、設計と試験に抜けている。**

   **根拠:** 現行コードは `overrides.mesh_euler.ni` があれば自動設定を行いません（[wt_card.py:1497](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:1497)）。一方、§4.1 の代入式と V2 にはこの例外がなく、実装次第でユーザ指定を黙って上書きします。

   **対案:** 新しい式は「明示指定がない場合だけ」と固定し、明示値は保持したまま節点数ガードに渡してください。V2 に上限以下・上限超過の明示指定を追加する。実装ステップには、現在仕様である `methods/design/overview.md` の更新も明記してください。

推奨は、**規模ガードを先に完成させ、Euler の既定変更を限定された検証範囲で導入する**方針です。実装前の修正順は **③ガードの対象経路 → ②判定精度 → ①既定値の適用範囲 → ④明示指定の保持**とします。

本変更は設計ツール側で閉じるため、周期境界や block-DPLUR のソルバ回帰を広く追加する必要はありません。node の対象ケースで、既存のメッシュ品質ゲート・初期化・段階起動を維持し、各検証の VERDICT を記録する方が適切です。ファイル変更・新規 CFD 実行は行っていません。

指摘数: Critical 0 / Major 3 / Minor 1
