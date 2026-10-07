# codex レビュー: tooling-nozzle-cfd-pinned-initial-line (result)

- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **stage**: `result` (diff base `integrate/main-2026-10-05`)
- **date**: 2026-10-05
- **commit**: `4215a6af` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 11.8 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **extra**: `case/45.isobutane_m6_d155/README.md`, `case/45.isobutane_m6_d155/_band_ab/c2pin_solve_pass2.json`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

抽出器・物理壁の追加試験と既定経路の互換性は確認できました。  
ただし、δ_E の取り違え、凍結入力の検査不足、凝縮結果の準定常証拠不足があり、現状の「P5 完了」は承認できません。

指定 diff を取得して確認しました。CUDA ソルバ本体の変更はありません。追加試験2本はともに FAIL 0。Hall＋`interp` は比較元に対し22万点で r・r′・r″ の差がすべて0でした。V3′ の非劣化判定が保留だったことも JSON と整合しています。ユーザによるトレードオフ採用自体は指摘対象にしていません。

run は別配置の `/home/sano/work/forge/case/45.isobutane_m6_d155/` を参照しました。以下の run 名はすべて `case/45.isobutane_m6_d155/` 配下です。最終 run の残差 CSV はこの配置に無く、収束については保存済み `CONVERGENCE_VERDICT.txt` を確認しました。

1. **Major — δ_E と緩和後の壁入力を混同し、閉包誤差と原因分析を約半分にしています。**

   [c2pin_hand0.py:18](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/c2pin_hand0.py:18) と [c2pin_solve.py:24](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/c2pin_solve.py:24) が読む `delta_r_next.csv` の第2列は、抽出 δ_E ではありません。[deltastar_loop.py:63](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:63) で ω=0.5 の緩和と再平滑化を受けた δ_next です。未緩和の平滑化抽出値は第4列です。

   保存 CSV から再計算すると：

   | run・出口 | 報告に使った δ_next/δ_C | 抽出 δ_E/δ_C |
   |---|---:|---:|
   | `run_0090_ns_c2pin_final` | 1.003582 | **1.007207** |
   | `run_0092_ns_c2pin_pass2` | 1.001861 | **1.003768** |

   生抽出値でも後者は **1.003614**。最終許容 ±1% には入りますが、「2 pass 目の予測 1±0.003 内」は誤りです。また手0の ΔM/ΔM_pred は、x=20/40/60/80/94.7 で **2.79/4.05/3.18/4.22/5.42 → 1.39/2.03/1.59/2.10/2.72**。変更ログの「1D 換算が3〜5倍過小評価」は、この取り違えを含みます。

   **対案:** 現行の緩和更新は δ_target と明記して保存し、δ_E を別列・別名称にする。閉包ゲートと手0を未緩和の抽出値で再判定し、§9・README・較正記録を訂正してください。

2. **Major — 「同じ形・同じガスの Hall 基準場」という凍結入力契約を実装が保証していません。**

   [cfd_initial_line.py:34](/home/sano/work/forge-integ-1005/design/forge_design/feedback/cfd_initial_line.py:34) は文字列で node／軸対称を判定し、Euler 判定は欠落を許す `prepare_info.viscous` だけです。実効 `viscMethod`、乱流、壁 BC、形状、ガス、Hall 基準場かどうかを照合しません。

   実測では `run_0062_euler_wallfit_fit_r1_ext6k` を **R=3、γ=1.4** に渡しても受理され、元の CFD 線に別条件の Hall M″=0.060605 が混ざりました。本来の値は0.024553です。読み取り内容だけを差し替えた検査でも、実効 config が `viscMethod: 2`・SST の入力を loader が受理しました。

   さらに snapshot 省略時は毎回最新場を選び、§5.1 #6c が要求する抽出値のハッシュも [runner_axismach.py:504](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:504) にありません。

   **対案:** YAML を構造として読み、実効設定・縮流部形状・熱力学条件・基準場の出所を照合する。生産入力は snapshot を明示し、snapshot／抽出線／熱力学条件のハッシュを保存してください。

3. **Major — 凝縮 `run_0093` の値は再現しますが、定常結果としての合格証拠がありません。**

   `run_0093_ns_c2pin_pass2_cond/res_12000.h5` から出口コア M=**5.985476**、onset=**59.39497**、S_max=**16.50177** を再現しました。しかし計算後の保存場は4000・8000・12000の3枚です。これらから再抽出した5量に正式ツールの `classify_series` を適用した結果は、すべて：

   `TRANSIENT-UNSETTLED — only 3 snapshot(s) (<4)`

   `res_0` は g=0 の dry 初期場で onset が存在せず、onset 系列へ加えると `NONFINITE` です。保存済み収束判定も `NOT CONVERGED` で、凝縮残差 `rms_rog_0`・`rms_roQ2_0` は `RISING`。これを単に plateau とまとめるのも不正確です。

   **対案:** 凝縮状態について十分な保存点を追加し、onset・出口 g・出口 M・S_max の同一区間の VERDICT を保存する。それまでは [plan:139](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:139) の⑤を「値は範囲内、準定常判定未完了」に戻してください。

4. **Major — 壁解像 FAIL の設計を完了扱いにし、提出報告ではさらに狭い範囲の別指標を表示しています。**

   `run_0094_ns_c2pin_pass2_ext6k` に `check_wall_resolution.py --groups wall` を再実行した結果：

   `VERDICT: FAIL` — y₁⁺最大 **13.392**、目標超過面積 **29.4%**。

   一方、同 run の `report/report.json` は最大 **8.671**・超過割合 **18.4%**。[nozzle_report.py:381](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:381) は速度差から壁応力を近似し、[同:399](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:399) で x>0 に限定して節点割合を返します。正式ゲートの全壁・面積割合とは別物です。

   **対案:** PPTX に正式な壁解像 VERDICT と全壁の数値を載せる。現状は暫定設計として扱い、壁解像または格子感度による裏付けを残作業の合格条件にしてください。旧 run も FAIL だったことは免除根拠になりません。

5. **Major — `joint` の物理壁が case/45 の縮流長に固定され、標準 geometry で使えません。**

   [wall_axismach.py:257](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:257) のランプは固定値 `(-11,-6)`。標準の L_U=3.5 では [同:370](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:370) により必ず拒否されます。δ_r=0 の構成でも再現しました。`prepare_ns` はこのランプを指定する経路を持ちません。

   **対案:** ランプを problem の明示設定として渡し、case/45 に `[-11,-6]` を登録する。他形状には適合する設定を要求し、対応範囲と短い縮流部の試験を追加してください。

6. **Minor — CFD ピンなのにアンカー出所が `hall` と記録されます。**

   [runner_axismach.py:496](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:496) は旧 `x_reach_cfd` だけで判定します。実際の `run_0094/prepare_info.json` も `initial_line.source="cfd"` と `anchor_source="hall"` が併存しています。

   **対案:** ピン経路を区別し、M・M′ は CFD、M″ は Hall という成分別の出所を記録してください。

7. **Minor — 現行仕様と残作業表が変更ログに追従していません。**

   [plan:80](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:80) の P2 本文は旧窓 `[-9,-6]`、[同:83](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:83) は撤回済みの出口 M 許容差と「x≤−1.5で縮流部不変」を残しています。#3〜#5 の状態も未更新で、本レビューで必要になった検証を置く行がありません。`check_plans.py` は PASS ですが、内容の整合性は検査していません。

   **対案:** §4・§5.1・§6 に現在有効な仕様と未完了条件を集約し、旧条件は変更履歴として明示する。README と `plans/README.md` も同期してください。

**推奨は、active を維持し、①δ の台帳と評価の訂正 → ②入力契約・汎用設定の修正 → ③凝縮と設計適格性の不足検証 → ④文書同期、の順で result レビューをやり直すことです。** dry の保存済み5評価量については正式判定の `STEADY` を再現できましたが、それだけで今回の完了条件全体は満たせません。ファイルは変更していません。

指摘数: Critical 0 / Major 5 / Minor 2
