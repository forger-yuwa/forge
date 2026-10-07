# codex レビュー: tooling-nozzle-cfd-pinned-initial-line (result)

- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-06
- **commit**: `24e5f475` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

主要な検証値と報告書の修正は確認できました。未達を残す採用判断と性能認定を分ける方針は妥当です。
ただし、CLI の実行不具合、凍結源の入力検査、残作業の移管記録を直すまでは `accepted` に移せません。

1. **Major — NS pass の CLI が報告生成前に `NameError` で終了する。**  
   [`deltastar_loop.py:298`](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:298) の `sys.exit(main())` より後に `_design_report` が定義されています。CLI ではその定義に到達せず、`run_pass`／`run_pass0_integral` の末尾で失敗します。計算とファイル書き込みを mock に置き換えた実行で、`NameError("name '_design_report' is not defined")` を再現しました。  
   **対案:** `__main__` ブロックを全関数定義の後へ移動し、CLI の両経路について報告呼び出しまでの回帰試験を追加する。「報告失敗でも chain を続行する」という手順書の契約も確認してください。

2. **Major — `inletProfile` を無視し、異なる実効入口条件を「照合済み」とする。**  
   [`cfd_initial_line.py:83`](/home/sano/work/forge-integ-1005/design/forge_design/feedback/cfd_initial_line.py:83) は `floats` の一様値だけを読み、`ints.inletProfile` と CSV を検査しません。[入口分布の仕様](/home/sano/work/forge-integ-1005/procedures/inlet-profile.md:11)では CSV が実際の Pt・Tt・組成を上書きします。  
   `run_0062` の読み取り内容だけに `inletProfile: 1` を追加すると、provider は受理し、熱力学ハッシュも元と同じ `cb08f7398369bad5` でした。#13a の実効条件照合には穴が残っています。  
   **対案:** 今回の一様入口契約では `inletProfile != 0` を明示的に拒否する。受理・拒否の試験と現在仕様の記述を同期してください。

3. **Major — 未確認事項の移管が「追記予定」のままで、完了記録が閉じていない。**  
   [plan §8](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:163) は評価格子・`ni` 感度を移管するとしていますが、[移管先の残作業表](/home/sano/work/forge-integ-1005/plans/active/tooling-design-problem-campaign-recipe.md:171)に該当項目がありません。凍結線の格子依存・Euler G2・波の準定常達成も、§8 の未確認事項から追跡先を確定できません。  
   また、[case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:94) は出口 M を単に「合格」とし、`run_0124/0125` は「走行中」、`run_0126/0127` は未掲載です。ユーザ決定 B と最終監査結果を反映していません。  
   **対案:** 各未確認事項を移管先 §5.1 に優先順付きで登録するか、延期・対象外の決定を明記する。run 索引には延長 run、判定、成果物、出口 M の「下限境界上」を反映してください。

4. **Minor — 非有限の入口条件が照合を通過する。**  
   [`cfd_initial_line.py:125`](/home/sano/work/forge-integ-1005/design/forge_design/feedback/cfd_initial_line.py:125) の比較は NaN に対して不一致を検出しません。読み取り時の入口 Tt を NaN にすると、provider が受理することを確認しました。  
   **対案:** 比較前に Pt・Tt・組成・熱力学パラメータの有限性と必要な物理範囲を検査し、非有限値の拒否試験を追加する。

5. **Minor — 採用根拠の `+0.009 %pt` は Mach 波ではなく圧力波。**  
   [plan §1](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:26) の名称が誤っています。[保存済み比較結果](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/wallfit_euler_ab_fit_vs_pincal_diag_fixedcoef.json:543)では、軸の差は `P_wave` が **+0.009196 %pt**、`M_wave` が **+0.001495 %pt** です。後者は `tail_ok: false`、総合判定は保留です。  
   **対案:** 採用根拠を「軸の圧力波」に訂正する。今回の報告書で Mach 由来の量を改称した処置を、過去の圧力由来の量に適用しないでください。

実測との照合結果は次のとおりです。run の一次索引は上記の case README です。

| 対象 | 確認結果 |
|---|---|
| `case/45.isobutane_m6_d155/run_0117_ns_recal_final_ext/` | `quantities_series.csv` の末尾5枚で出口 M **5.998886995 — STEADY**、波 **0.00651046 % — DRIFTING**、δ_E/δ_C **0.9998346 — STEADY**。計画の記述を再現 |
| `case/45.isobutane_m6_d155/run_0118_ns_recal_final_cond/` | `cond_series.csv` の末尾5枚で開始位置 **57.9382**、S_max **16.8618**、出口 g **3.12749e−4**、出口 M **5.9864071**。再実行した `check_quasisteady.py` は **全 STEADY** |

両 run の保存済み残差判定は **`NOT CONVERGED (stalled/plateau)`**。最終 HDF5 の `VALUE/*` に非有限値はなく、再生成済み PPTX（14枚／15枚）には Mach 波表記と量別 VERDICT が反映されています。抽出器の V0・主要入力契約の読み取り専用試験、物理壁試験、格子パラメータ試験も合格しました。

ただし、参照した `/home/sano/work/forge/` のローカルコピーには元の `residual_history.csv` と `run_0126/0127` の実体がなく、全残差と軸細分監査の独立再計算は未実施です。

**推奨は、上記番号順に修正・確認し、未達事項の移管を完了してから、§8 の「実装・評価完了、性能認定は保留」という範囲で `accepted` に移すことです。** 本レビューではファイルを変更していません。

指摘数: Critical 0 / Major 3 / Minor 2
