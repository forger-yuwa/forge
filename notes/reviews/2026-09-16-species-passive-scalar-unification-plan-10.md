# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `a0424677` (feature/sern-design)
- **codex**: effort `high`, 5.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m1
- **focus**: 10 回目 (plan-9 NO-GO M1–M4/m1 の反映; §4.7 v9, 実装・単体済, commit a0424677)。変更: (1) 独立照合の始点 = 計算開始前の root-only 総量 ([passive] initial total, %.12e)、終点 = 全後処理後の確定状態の再集計; (2) 収支ゲート: FCT/クランプ記録を最終 step に結合 (不一致・欠落は FAIL)、モーメント成分はクランプ記録必須、全生値の有限性 + solver 側 nonfinite フラグ、独立照合 1e-6、失敗系試験 (実書式 parse→判定) を追加; (3) 次数ゲート: --nsub 必須、全物理 step の inner_iter 行と必要列の有限性、ゼロ始点→非ゼロは拒否、config の bdfOrder / nSub 倍 / passiveFct・SFR 一致、totalTime 有限; (4) 実現可能性ゲートは solver と同じ ρ_l(T) (H2O/N2 式) で全凝縮種、特異不整合・負値・非有限も違反; (5) r_H/a。問い: 収支・次数ゲートに残る偽 PASS/誤判定はあるか、§6-2/§6-6 の合否をこれらのツールの VERDICT に委ねてよいか。設計 (§4.1–4.8) の再指摘は不要 — ゲートの完備性だけ。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
対象は `a0424677`。前回の総量集計時点の修正は確認しましたが、収支・次数ゲートには偽 PASS と誤 FAIL が残っています。  
現状では §6-2／§6-6 の合否を、この 2 ツールの VERDICT に委ねられません。指摘はゲートの完備性に限定します。

1. **Major — 収支ゲートが、検査対象・記録の欠落を許している**

   **根拠:** [check_passive_budget.py:76](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:76) は「読み取れた成分」だけを検査します。設定から必須成分集合を確定せず、FCT 非作動時には閉合・総量差の検査も行いません。CSV が無ければ最終 step 照合を省略します（[同:115](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:115)）。

   実書式を使った `parse_lines → evaluate` で、以下を再現しました。

   | 入力 | 現行判定 |
   |---|---|
   | FCT 非作動、初期総量 `1` → 最終総量 `2`、補正ゼロ | PASS |
   | FCT 宣言は 5 成分、収支記録は `roXi` だけ | PASS |
   | 最終 step 不明、step 1 の記録だけ | PASS |
   | 累積 floor の `lo=NaN`、絶対量・相対量はゼロ | PASS |
   | solver の `nonfinite` フラグが欠落 | PASS |

   `lo/hi` は正規表現で取得しても検査辞書へ保存せず、`nonfinite` 欠落はゼロ扱いです（[同:44](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:44)、[同:53](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:53)）。「全生値の有限性・記録の完備性」は未達です。

   **対案:** config と実行時宣言から必須成分・FCT 作動条件・終了 step を確定し、欠落や解析失敗は拒否してください。全数値トークンと必須フラグを検査し、非 FCT の保存試験にも収支検査を用意するか、明示的に「判定対象外」としてください。

2. **Major — sub-iter ゲートが全化学種・全反復・全物理 step を保証しない**

   **根拠:** [analyze_moment_order.py:121](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:121) の必須列に `rms_roY*`、`rms_roXi`、凝縮種 `_1` 以降がありません。低下量不足を FAIL にするのは必須列だけです（[同:106](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:106)）。

   また、各 step の最初・最後以外の値を検査せず、期待 step 集合も CSV 内の `outer_*` 行から作ります（[同:73](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:73)）。

   ファイルを作らない模擬入力で、実際の `main()`／CSV 判定を通して確認しました。

   - 化学種残差が初回から最終へ **100 倍増大**しても `VERDICT: PASS`。
   - 中間 `inner_iter` に **NaN** があっても `VERDICT: PASS`。
   - ある物理 step の outer／inner 行が丸ごと欠落しても `VERDICT: PASS`。
   - 必須列が全反復で厳密ゼロの場合は、逆に `no usable sub-iter data` で FAIL。

   実 run でも、`run_0258_passiveD_order_bdf2_dt4e-6_nsub40/` の `rms_roY1` 最小低下は **0.99 桁**と表示されますが、この列の低下不足は FAIL 理由に入りません。

   **対案:** 必須列を config の全輸送成分から生成し、全 inner 行の有限性、期待する物理 step・最終 sub-iter の存在を検査してください。全反復ゼロは成立として扱い、途中だけ非ゼロになる列とは区別する必要があります。

3. **Major — 次数試験の設定照合が、不正な比較を受け入れる**

   **根拠:** [analyze_moment_order.py:146](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:146) は、双方の `nsub` が真値の場合だけ倍増を確認します。欠落・ゼロでは検査を省略します。刻みも有限・正であることを確認せず比を比較するため、NaN は拒否されません（[同:139](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:139)）。

   模擬入力の `main()` 実行で、**`nSub` 欠落、全 checkpoint の `dt=NaN` がそれぞれ `VERDICT: PASS`、exit 0** になることを確認しました。

   さらに、FCT／SFR は「4 run 間で同じ」ことしか確認しません。全 run が SFR 0 でも PASS でき、`solver`・`dualTime`・`passiveScalarScheme` は照合対象外です。これでは §6-6 の **FCT が実際に作動する試験**の証明になりません。

   **対案:** 必須設定の存在・型・範囲、config と checkpoint の刻み一致、3 水準の同一条件、倍増 run の差分限定を検査してください。FCT 有効試験には期待条件を明示する検証モードを設け、実行ログの作動記録まで照合すべきです。

4. **Major — 確定場の有界性・実現可能性を、solver と同じ条件で検査していない**

   **根拠:** [analyze_moment_order.py:40](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:40) が読むのは `rog/roQ*` ではなく **原始量 `g/Q*`** です。float32 の原始量変換を挟むため、許容境界付近で判定が変わります。

   同じ温度依存物性式・同じ不等式で実測すると、

   | 最終場 | 原始量からの違反数 | 保存量からの違反数 |
   |---|---:|---:|
   | `case/44.vitiated_air_wt/run_0289_passiveG_order3_bdf2_dt1.6e-5_nsub40/res_100.h5` | 26 | 18 |
   | `case/44.vitiated_air_wt/run_0281_passiveG_order1_bdf2_dt8e-6_nsub40/res_200.h5` | 76 | 70 |

   また、solver は **無次元 `x` または `y ≤1e-30`** を退化修復します（[condensationRealizability_d.cuh:21](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:21)）。ゲートは生の `Q1/Q2 ≤0` しか特異判定しません（[analyze_moment_order.py:175](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:175)）。`x=1e-31, y=1e-16` は solver の修復対象ですが、ゲートの `VERDICT: PASS` を再現しました。

   トレーサについても、収支ゲートは場を読まないため、**総量が正しくても局所的に `roXi/ro` が範囲外である状態を拒否できません**。

   **対案:** 必須全成分の確定保存量と正の密度・有限な温度から検査してください。solver と同じ退化条件・物性経路を使い、成分欠落、`roXi/ro` の上下限、モーメントの非負性・実現可能性を独立した必須ゲートにしてください。

5. **Minor — 同一終了時刻の検査が、float32 の時刻積算誤差で誤 FAIL する**

   **根拠:** `totalTime` は `flow_float` で毎 step 加算されます（[solverConfig.hpp:42](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:42)、[main.cpp:1870](/home/sano/work/forge/solver_density_cuda/main.cpp:1870)）。一方、ゲートは相対 `1e-6` で一致を要求します。

   同じ名目終了時刻 `0.0016 s` の `run_0256/0257/0258` は、checkpoint 上ではそれぞれ

   `0.001600000425 / 0.001599996001 / 0.001600001473`

   となり、`totalTime differs across runs` が発生しました。これらの run は他の条件でも FAIL ですが、この時刻判定は別の誤拒否要因です。

   **対案:** 時刻を double または開始時刻＋整数 step×刻みで管理し、ゲートもその情報と照合してください。単純に許容幅を広げる対応は推奨しません。

**推奨は一つです。設計を維持し、上記 1→4 の順にゲートを欠落時に必ず拒否する形へ修正し、5 の時刻整合も直してから正式な統合検証に進んでください。** 今回再現した偽 PASS／誤 FAIL を、判定まで通す回帰試験として追加する必要があります。

既存の収支単体試験は `ALL PASS` でした。一方、実 run の収支ゲートは `run_0140_passiveG_fct_step_seam/` が **FAIL**（補正合計比 `7.31e-6`）、`run_0289` も **FAIL** でした。旧 run の結果を修正版 solver の検証実績として扱ってはいません。

ファイル変更・新規 CFD 実行はしていません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
