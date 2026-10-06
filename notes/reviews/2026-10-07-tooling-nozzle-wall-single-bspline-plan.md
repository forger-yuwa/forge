# codex レビュー: tooling-nozzle-wall-single-bspline (plan)

- **plan**: [`plans/active/tooling-nozzle-wall-single-bspline.md`](../../plans/active/tooling-nozzle-wall-single-bspline.md)
- **stage**: `plan`
- **date**: 2026-10-07
- **commit**: `56ce09ce` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **extra**: `notes/reviews/2026-10-07-wall-single-bspline-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

全域 1 本の x の 5 次 B-spline にする方針は妥当です。ただし、改訂後の W1 を満たす λ は今回の走査でも見つからず、入力仕様・下流連携・IC 移送にも実装前の修正が必要です。推奨は、方針を維持して以下の条件を先に確定することです。

目的は既存 plan と重複していません。`monotone-r2` は設計壁の形状拘束、本件は δ_r の平滑化と物理壁の表現統合です。ノットの和集合による加算、C² 接続での重複度 3、3 階微分の積に対する区間ごとの 3 点 Gauss 求積は数学的に成立します。スカラー関数と平面曲線の係数を区別した CAD 方針も適切です。[SciPy の B-spline 定義](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html)

今回の変更で CUDA の流束・block-DPLUR・周期処理を直接変更する必要はありません。検証は幾何試験と case/45 の node NS・凝縮を中心にするのが妥当で、cell・周期箱の追加計算は不要です。

1. **Major — 改訂後も W1 の実行可能性が確認できない。λ 選択の試算を実装より前に置くべきです。**

   根拠は [plan:60](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:60)、[W1:113](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:113) と、保存入力による独立試算です。

   別ワークツリー `/home/sano/work/forge` の `case/45.isobutane_m6_d155/run_0147_ns_mono_final/` にある `prepare_info.json`・`delta_r_initial.csv` から入力を復元しました。旧平滑化と保存 δ_r の最大差は **6.66e−16 r_t** です。

   | 平滑化 | 相対 RMS | 出口半径の旧方式との差 |
   |---|---:|---:|
   | 現行 | 1.90186% | — |
   | x 基底、λ＝3 | 2.00942%：残差条件内 | +5.937 µm |
   | x 基底、EDF 同等 λ＝17.58294 | 2.66090%：残差条件外 | +4.302 µm |

   さらに **λ＝0 と 10⁻⁸〜10⁴、0.05 decade 刻みの計242候補**を調べました。残差条件を満たす173候補のうち、ランプの高周波条件も満たすものは **0** でした。

   この試算では「0.5 r_t 移動平均」を中心付き連続平均とし、r″ の平均を両端の r′ から積分で評価しました。ランプの評価点を10,001点から40,001点に増やしても、最良候補の高周波最大値は旧方式の **1.00099218倍**、最大位置は x＝−11 のままでした。有限走査なので、全 λ に対する不可能性の証明ではありません。

   **対案:** §5 の最初に、フィルタの中心・端部処理・評価領域・数値許容差を固定した実行可能性試験を追加してください。合格 λ が存在することを確認してから統合実装へ進む。見つからなければ §4.1 の停止条件を適用し、閾値を結果に合わせて緩めないことです。

2. **Major — CSV の意味と有効範囲が未決定のため、二重平滑化と C² 接続の破綻を防げません。**

   [plan:64](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:64) は「明記する」としていますが、仕様自体はまだありません。現行の [deltastar_loop.py:88](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:88) は抽出値を平滑化し、緩和後にも再平滑化します。一方、[runner_axismach.py:910](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:910) の CSV は既に平滑化済みです。

   平滑化は冪等ではありません。同じ保存入力・x 基底・λ＝3で平滑化をもう一度掛けると、**出口半径がさらに −0.583 µm** 変化しました。「同じ関数を通す」だけでは W3 の同一壁を保証できません。

   また現行 CSV 読み取りは範囲外を端値で延長し、導関数を0にします（[同:923](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:923)）。許容されている部分範囲入力 `δ_r＝0.01＋0.001x、x∈[0,10]` を実行すると、x＝0の左右で δ_r′ が **0 → 0.001** と跳びます。これを厳密保持しながら C² の壁にすることはできません。

   **対案:** §4.1で入力を「生値」「旧平滑化済み CSV」「新しい係数・ノット」に区別し、平滑化を行う場所を確定してください。新形式は再平滑化せず復元し、旧 CSV は明示的な移行処理を経る仕様にする。新経路では δ_r が実際に適用される全域を覆う入力を要求し、不足範囲の黙示的な端値延長を拒否してください。W3 は**同じ入力の異なる受け渡し経路**を比較する試験に限定すべきです。

3. **Major — 「同じ呼び出しで使える」というインターフェースと実装範囲が不足しています。**

   [plan:76](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:76) の `r`・`x_in`・`x_e` だけでは、現行の NS 準備は動きません。

   - [runner_axismach.py:1059](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1059) は `validate()`、`:1142`以降は `x_throat`・`r_throat`・`kappa_throat`・`_dstar_hist`・`offset_mode` を使います。
   - [ic.py:68](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/ic.py:68) はスロート属性が無いと **黙って `(0,1)` に戻り**、面積比と亜音速／超音速の枝を決めます。
   - [nozzle_report.py:485](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:485) は壁クラスを受け取らず、CSV を読み直して補間します。§7の「報告は壁クラス経由」は現状と異なります。
   - §4.1で変更するとした `deltastar_loop.py`・`c2pin_solve.py` も、§5・§7の実装一覧から抜けています。

   **対案:** 壁クラスの必要属性と診断情報を列挙し、スロート量を統合後の spline から求める仕様を追加してください。保存係数を読む共通ローダーを報告にも接続し、上記ファイルを実装表へ追加する。対応する設計壁・オフセット方式も明示し、未対応の組合せは新キー指定時に拒否してください。

4. **Major — W4 が継承する IC 移送手順は、今回の壁移動にそのまま適用できません。**

   [plan:118](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:118) が参照する [monotone plan:252](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-throat-monotone-r2.md:252) は番号写像・段階起動なしです。しかし [ic_index_map.py:26](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ic_index_map.py:26) の移動上限は **1 µm**。今回の λ＝3では出口だけで **5.937 µm** 動きます。これは残差条件内の候補でも起きます。

   MOC plan には既に、変換後メッシュの移動・接続・反転の検査と別 IC による確認が追加されています（[MOC plan:130](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:130)）。本件単独で採用する場合にも必要です。

   **対案:** W4の前提として、**最終形状の変換後メッシュ**で移動量・境界対応・反転・壁距離・品質を確認し、IC 移送法を決める工程を追加してください。番号写像の上限だけを引き上げず、採用時は別 IC・段階起動との依存性確認を含める。`geom_float` も float32 なので（[flowFormat.hpp:7](/home/sano/work/forge-integ-1005/solver_density_cuda/flowFormat.hpp:7)）、倍精度 spline の誤差保証と変換後メッシュの健全性は別々に検査すべきです。

5. **Minor — W5 の許容差は量ごとの次元と零曲率を扱えていません。**

   [plan:119](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:119) は値・接線・継ぎ目に一律 `1e-12 r_t`、曲率に相対 `1e-9` を指定しています。r′ は無次元で、直管の曲率は0です。相対誤差だけでは直管や変曲点を判定できません。

   **対案:** 無次元座標で半径・r′・r″・曲率を別々に定義し、各量に絶対許容差と相対許容差を設ける。重複ノットでは左右極限を評価し、物理単位への変換試験と分離してください。

6. **Minor — 機能実装だけで完了する場合のレビュー工程が抜けています。**

   [残作業表:101](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:101) は W0 を列挙せず、result レビューも生産反映の行にだけ置いています。一方、[完了条件:136](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:136) は機能実装完了でも W0 と2回のレビューを要求しています。

   **対案:** `W0〜W3・W5 → 機能実装の result レビュー` を必須工程にし、生産採用時はさらに `W4 → 生産結果のレビュー → ユーザ判断` と明記してください。

**推奨は、全域 1 本の方針を維持し、①W1の実行可能性、②入力仕様、③連携インターフェース、④IC移送条件の順に plan を修正してから実装することです。** λ探索が不成立なら、その時点で止めるべきです。

保存 run の `CONVERGENCE_VERDICT_segment.txt` は **main 区間・最終 step 59999、NOT CONVERGED（stalled/plateau）** でした。今回の数値は幾何・平滑化の試算であり、流れの収束・一致を示しません。CFD は起動せず、ファイルも変更していません。指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
