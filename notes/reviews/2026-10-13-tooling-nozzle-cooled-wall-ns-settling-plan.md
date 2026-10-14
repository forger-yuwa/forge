# codex レビュー: tooling-nozzle-cooled-wall-ns-settling (plan)

- **plan**: [`plans/active/tooling-nozzle-cooled-wall-ns-settling.md`](../../plans/active/tooling-nozzle-cooled-wall-ns-settling.md)
- **stage**: `plan`
- **date**: 2026-10-13
- **commit**: `4e8b3672` (feature/nozzle-cooled-wall-settling)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **extra**: `case/70.m5_air_wt_demo/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

課題の同定は妥当で、前身 plan の accepted 範囲にも整定問題の解決は含まれていません。  
ただし、停止判定には誤合格する反例があり、CFL 掃引・帯固定・Euler 参照の受入れにも不足があります。以下を実装前に修正してください。

1. **Major — 5 枚の単一等比外挿では、遅いモードを見逃して合格できる**

   **根拠:** [plan:58](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-cooled-wall-ns-settling.md:58) は、同符号で減衰する増分から残りを推定し、最初の合格窓で停止します。しかし、複数モードに対する成立条件がありません。

   計画の式を、合成系列 `Mᵢ = 4.99 + 10⁻⁴·0.5ⁱ + 10⁻²·0.999ⁱ`、`i = 0…4` に適用して確認しました。

   | 判定量 | 値 |
   |---|---:|
   | 窓内変化 | 1.337e−4 |
   | 末尾の増分比 ρ | 0.7215 |
   | 推定残変化の絶対値 | 4.203e−5 |
   | 推定漸近 M | 4.999924 |
   | **真の漸近 M** | **4.990000** |

   増分はすべて同符号で減衰し、計画の M 条件と設計条件を通りますが、真値は不合格です。これは実 run の予測ではなく、**判定仕様の反例**です。逆に、完全な定数列では比が 0/0、小さな丸め変動では符号反転となり、現仕様では整定していても止まれません。

   **対案:** 等比外挿を単独の停止根拠にしないこと。窓長・開始位置を変えた漸近推定の安定性、独立した後続窓での予測確認、連続合格を要求してください。定数・微小変動には別の有界変動判定を定義し、S0 に定数列、微小振動、二重指数、減衰比が 1 に近づく系列を追加する。合否境界を推定不確かさがまたぐ場合は保留です。

2. **Major — S1 は残差が悪化した腕を「速い」として選び得る**

   **根拠:** [plan:87](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-cooled-wall-ns-settling.md:87) の採用条件は、4 量の整定と腕間差だけです。

   既存の `case/70.m5_air_wt_demo/run_0014d_ns_ns2_cfl5/` を、通常判定で PASS の `run_0014_ns_ns2/` を参照として `check_convergence.py --from-floor` で再評価すると、

   > `NOT CONVERGED (residual left the reference floor, or column mismatch)`

   でした。末尾平均は参照床に対して、`rms_roUx` **32.72 倍**、`rms_roUy` **15.05 倍**、`rms_roe` **5.32 倍**です。開始場は S1 と異なるため将来の腕の失敗を証明しませんが、**CFL 5 の加速実績を品質の裏付けにはできません**。

   また、float32 の commit は [update_d.cu:301](/home/sano/work/forge-deliv/solver_density_cuda/cuda_forge/update_d.cu:301) の加算で、小さい更新が消えます。場が動かないことと離散方程式の収束は別です。整定した腕が 1 本だけの場合、現仕様では腕間一致も検証できません。

   **対案:** S1 の採用前に全保存量の収束 VERDICT を要求し、継続区間・参照床の扱いを事前登録してください。少なくとも異なる CFL の有効な 2 腕で確認し、1 本しか成立しなければ判別未了とする。残差と場の停止が食い違う場合は更新量／ULP を調べる。なお、既存 `qAccumulatorFP64` は [main.cpp:3371](/home/sano/work/forge-deliv/solver_density_cuda/main.cpp:3371) で軸対称を拒否するため、そのまま代替に使えません。

3. **Major — 帯位置の窓内固定は、窓間の枝切り替わりを解決しない**

   **根拠:** [plan:68](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-cooled-wall-ns-settling.md:68)・[plan:90](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-cooled-wall-ns-settling.md:90) は、選択した一つの窓の内部だけを評価しています。実装は [deltastar.py:357](/home/sano/work/forge-deliv/design/forge_design/metrics/deltastar.py:357) で最初の閾値交差を選ぶため、梯子を平均しても、別の窓で別枝を選ぶ可能性が残ります。

   `case/70.m5_air_wt_demo/_band_c_ab/` の CSV を再集計すると、`run_0010` の 50k→55k の最大変化は c=1.6 で **9.50e−5 r_t** と S2 の閾値内ですが、50k の出口 δ_E は c=1.25 に対して **−9.96e−4 r_t** 変わっています。窓内変動の小ささは、抽出値の頑健性を保証しません。

   **対案:** 連続する非重複窓で帯をそれぞれ独立に選び、平均 δ_E の差にも使用域全体で **≤1e−4 r_t** を要求してください。共通帯での場の変化と、帯再選択による変化を分けて記録する。各スナップショットの hard/held 検査と、前身 V3e の残変化判定を維持し、不合格を平均で消さないこと。

4. **Major — Euler 参照の厳密な受入れが、前提として置かれたまま未計画**

   **根拠:** [plan:81](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-cooled-wall-ns-settling.md:81) は Euler の PASS を要求しますが、S0〜S3 にそれを達成・判別する工程がありません。再実行した VERDICT は次のとおりです。

   - `case/70.m5_air_wt_demo/run_0012_euler_c2b/`: `NOT CONVERGED (still converging — run more steps)`
   - `case/70.m5_air_wt_demo/run_0012b_euler_c2b_cont12k/`: `NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`

   さらに、[check_convergence.py:359](/home/sano/work/forge-deliv/solver_density_cuda/tools/check_convergence.py:359) は通常判定で不合格の参照を `--from-floor` に使うことを拒否します。単なる延長や床の指定で受入れ問題を解消できる計画にはなっていません。

   **対案:** S3 より前に独立した Euler 受入れ工程を追加してください。前後両設計について、残差停滞の切り分け、出口 M・流量の準定常、同じ NS を使った参照変更時の δ_E 感度を確認する。[既存の Euler 調査 plan](/home/sano/work/forge-deliv/plans/accepted/verification-case45-euler-total-enthalpy.md:43) の専用格子・参照適格性検査を出発点にし、case/70 に同じ原因があるとは決めつけないこと。

5. **Minor — 引継ぎ用の δ_E 再生スクリプトが現行コードでは動かない**

   **根拠:** [deltaE_series.py:30](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/deltaE_series.py:30) は旧実装の文字列が 1 回存在することを assert します。現行 `deltastar.py` 内の該当数を実測すると **0** でした。`one()` は [同ファイル:80](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/deltaE_series.py:80) で必ずこの処理を呼ぶため、前身の再生経路をそのまま使えません。

   **対案:** S2 の前提作業として実装前の文字列差し替えを除去し、現行抽出器を直接呼ぶ経路へ更新する。既存保存 CSV との比較で再現性を確認してください。

6. **Minor — 残作業表と影響範囲が、既に決めた方針に追従していない**

   **根拠:** [plan:77](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-cooled-wall-ns-settling.md:77) は原因を「温度境界層の発達」として完了扱いですが、§3.1 は連動の観測に限定しています。#2・#3 は方針未決の記述、§7 は未定のままです。共通抽出器を変更するなら M6 にも影響します。

   **対案:** #1 は「連動する量・領域の測定完了、原因未確定」と訂正する。#2・#3 を具体的な実装・検証作業に改め、case/70 専用処理と共通 API の変更範囲、既存単一スナップショット経路の非劣化確認を明示してください。

**推奨は、既存ソルバを使う段階的な検証計画として進めることです。** 実装前の優先順は、①停止判定の修正、②Euler 受入れ工程の追加、③S1 の残差・複数腕ゲート、④S2 の窓間再現性ゲートです。

case/70 を主検証にする選択は妥当です。CUDA 共通部を変更しない範囲で cell・周期ケースまで増やす必要はありません。S3 では前身 V5 のメッシュ品質・壁解像・寸法ゲートを維持し、同一メッシュ restart と再設計後の cross-mesh restart を明記してください。

今回、ファイル変更・forge の新規実行は行っていません。既存系列への判定ツール再実行、CSV 再集計、停止式の合成反例検証を実施しました。上記提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
