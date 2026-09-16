# codex レビュー: species-passive-scalar-unification (plan)

- **plan**: [`plans/active/species-passive-scalar-unification.md`](../../plans/active/species-passive-scalar-unification.md)
- **stage**: `plan`
- **date**: 2026-09-16
- **commit**: `b437a45a` (feature/sern-design)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m0
- **focus**: 11 回目 (plan-10 NO-GO M1–M4/m1 の反映; ゲートの完備性のみ, 実装済 commit 済)。変更: (1) passive_gate_common.load_config が solverConfig.yaml から必須成分 (tracer/nCondSpecies)・FCT 作動条件・終了 step・名目終了時刻 dt×nStepOuter を確定; check_passive_budget.py は mode auto (fct / conservative [非 FCT 非定常: 総量変化 ≤ tol] / steady)、必須成分の欠落・余分・不完了・全トークン非有限・nonfinite フラグ欠落・FCT 宣言なのに作動記録なし を FAIL、確定場ゲートも実行; (2) check_passive_field.py: 保存量 (roXi/ro の [0,1], rog/roQ* の非負, solver と同じ ρ_l(T)・無次元・退化条件 x,y≤1e-30 の実現可能性), ρ>0, T 有限; (3) analyze_moment_order.py: config の存在/型/範囲、3 水準同一・倍増 run の差分限定、名目時刻一致と checkpoint 整合 (1e-9)、config 同一 (solver/dualTime/scheme/SFR/passiveFct/condModel/nCond/tracer)、--expect-fct (config + 各 log の [passiveFct] active)、全輸送列・全 inner 行の有限性・期待 step の outer/inner 行・ゼロ始点→非ゼロ拒否・全列の低下 ≥2 桁、確定場ゲート; (4) 物理時刻は totalTimeD (double) で積算し checkpoint に書く; (5) 射影の ρ_l は double の式 (float 表との差がゲート許容を超えないよう)。問い: ゲートに残る偽 PASS / 誤 FAIL はあるか、§6-2/§6-6 の合否をこの 3 ツールの VERDICT に委ねてよいか。設計の再指摘は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
対象は `b437a45a`。ゲートに偽 PASS と誤 FAIL が残っており、§6-2／§6-6 の合否はまだ委ねられません。  
設計の再指摘はせず、現行コード・既存 run・ファイルを作らない模擬入力で確認した問題に限定します。

1. **Major — sub-iter の必須列と反復履歴を保証していない**

   **根拠:** [analyze_moment_order.py:48](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:48) は、CSV に存在する列だけを検査します。config から必須列を生成していません。また、ゼロ始点の検査は最終値との比較だけです（[同:74](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:74)）。

   メモリ内の入力を実際の `main()` に通し、次の両方で **`VERDICT: PASS`、exit 0** を再現しました。

   - 化学種・モーメント・運動量・エネルギーの残差列を削除し、`rms_ro` だけ残す。
   - 各 step の残差を **`0 → 1 → 0`** にする。

   `inner` 番号も検査しないため、所定の最終 sub-iter まで記録された保証がありません。

   **対案:** config から流れ・化学種・受動種の必須列を生成し、欠落を拒否してください。各 step の反復番号と終了記録を照合し、ゼロ始点は**全 inner 行がゼロ**の場合だけ受理する必要があります。

2. **Major — 「nSub 以外は同一」という比較条件が未実装**

   **根拠:** 設定比較は [analyze_moment_order.py:123](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:123) の限定されたキーだけです。`convMethod`、`limiter`、`speciesImplicitCoupling`、`passiveImplicitCoupling`、緩和、メッシュ・IC・BC は比較されません。

   倍増 run だけ **`convMethod 1→0`、`limiter 2→0`、化学種 coupling `1→0`** に変更した模擬入力でも、`main()` は **`VERDICT: PASS`、exit 0** でした。これでは差を sub-iter 誤差として解釈できません。

   設定読込にも不一致があります。[passive_gate_common.py:20](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:20) は `bdfOrder` 欠落を **1** にしますが、solver の既定は **2** です（[solverConfig.cpp:450](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:450)）。「必須キーの存在確認」にもなっていません。

   **対案:** 実効設定を solver と同じ既定値で正規化し、許可する差分を `dt/nStepOuter/nSub` と出力先に限定してください。メッシュ・IC・BC の内容も照合し、整数項目は変換前に型・範囲を検査してください。

3. **Major — 収支ログの「全トークン・全期間の有限性」を検査していない**

   **根拠:** [check_passive_budget.py:24](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:24) は per-step 部分を読み飛ばし、[同:45](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:45) は過去の記録を同一成分の最新値で上書きします。有限性検査は残った辞書に対してだけです。

   既存単体試験の実ログ書式を使った `parse_lines → evaluate` で、次がともに **True＝PASS** でした。

   - 最終行の `per-step ... lo nan`。
   - step 50 の累積絶対補正が `nan`、step 100 の記録は有限。

   **対案:** 読込時点で全数値トークンを検査し、一度でも異常を検出したら最後まで失敗を保持してください。対象ログの解析失敗も黙って読み飛ばさず拒否し、今回の反例を回帰試験に追加してください。

4. **Major — 検査する保存場が終了 step の場とは限らない**

   **根拠:** [passive_gate_common.py:33](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:33) は存在する最大番号のファイルを選ぶだけで、`nStepOuter` や checkpoint 時刻と照合しません。収支側も CSV が無ければ完了照合を省略します（[check_passive_budget.py:178](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:178)）。

   config・収支ログは **100 step 完了**、保存場は **`res_0.h5`／時刻 0 だけ**という模擬入力で、収支ツール全体が **`VERDICT: PASS`、exit 0** になりました。さらに CSV を欠落させても PASS でした。

   **対案:** config・ログ・CSV・HDF5 を同一の終了 step に結び付け、終了場と完了記録を必須にしてください。`check_passive_field.py` 単独でも、古い場を最終場として受理してはいけません。

5. **Major — double 積算への変更だけでは時刻照合の誤 FAIL が解消しない**

   **根拠:** `dt` は依然 `flow_float` です（[solverConfig.hpp:44](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:44)）。[main.cpp:1870](/home/sano/work/forge/solver_density_cuda/main.cpp:1870) は丸め済みの `dt` を double に変換して加算します。一方、ゲートは YAML の double 値に対し、時刻を相対 `1e-9`、刻みを `1e-12` で照合します（[analyze_moment_order.py:137](/home/sano/work/forge/case/44.vitiated_air_wt/analyze_moment_order.py:137)）。

   実データと同じ変換で次を確認しました。

   | YAML の `dt` | float32 の実効 `dt` | 相対差 |
   |---|---:|---:|
   | `8e-6` | `7.999999979801942e-6` | `2.525e-9` |
   | `0.002` | `0.0020000000949949026` | `4.750e-8` |

   後者は `case/09.Taylor-Green/run_0148_passiveG_fct_step_seam/res_160.h5` の checkpoint にも記録されています。**新しい double 積算でも、この刻みの丸め差は残り、両検査の許容を超えます。**

   **対案:** 実行精度で丸めた実効刻みと開始時刻を正本にし、checkpoint の刻み・終了時刻を照合してください。名目刻みとの区別を明示し、実際の float32 刻みを使った正常系試験を追加してください。

6. **Major — 非 FCT の非定常 run を一律「総量不変」で判定するのは誤り**

   **根拠:** [check_passive_budget.py:172](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:172) は非 FCT の dual-time を一律 `conservative` にし、[同:124](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:124) は境界流束・物理ソースを考慮せず初期総量との差を検査します。

   これは閉じた無ソースの移流試験には適切ですが、§6-6 の凝縮ソース作動試験には不適切です。例えば初期総量 `1`、ソース積分 `0.001`、最終総量 `1.001` の保存的な解も、総量変化約 `1e-3` として拒否されます。

   **対案:** 総量不変の判定を、閉境界・無ソースと確認できる試験に限定してください。一般の非 FCT 非定常 run には境界流束・ソース込みの収支が必要です。記録が不足する場合は「判定不能」とし、物理的な総量変化を保存違反と判定しないでください。

**推奨は、設計を維持してゲート修正を先に完了することです。** 上記 1→4→2→5→6→3 の順で直し、正常系と今回の反例を `main()`・終了コードまで通す回帰試験にしてから正式検証へ進めてください。

既存の収支単体試験は `ALL PASS` でしたが、上記の抜けは検出しません。実 run の `case/09.Taylor-Green/run_0148_passiveG_fct_step_seam/` は現行ゲートで **`VERDICT: FAIL`**、補正合計比 **`7.30e-6`** でした。これは旧成果物の確認であり、修正版 solver の検証実績とは扱っていません。

ファイル変更・新規 CFD 実行はしていません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 0
