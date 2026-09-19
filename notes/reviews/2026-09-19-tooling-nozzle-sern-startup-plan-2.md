# codex レビュー: tooling-nozzle-sern-startup (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-startup.md`](../../plans/active/tooling-nozzle-sern-startup.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `fde8f3fc` (feature/sern-design)
- **codex**: effort `high`, 6.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m2
- **focus**: R-h (rising 判定にジッタとの比較を足した変更, check_convergence.py) と R-b の目的量比較結果

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
R-b は短縮レシピの有望さを示していますが、記載された設定・誤差上限には不一致があります。  
R-h は実際の残差爆発を `flat` に変えるため、「変更はすべてジッタの誤検出修正」という結論を支持できません。

1. **Major — R-h は本物の発散傾向を消す。回帰結果の解釈が誤っている。**

   根拠は [`check_convergence.py:120`](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:120)。長い窓全体の回帰勾配がジッタ閾値を超えないと、末端で急増していても `rising` を禁止します。

   実例は `case/37.pintle_nozzle/run_0009_slau_outflow_full/residual_history.csv`。`rms_roUx` は step 1103 の **2.77e-4 → 1106 の 7.02e5 → 1108 の 7.58e22** と爆発していますが、旧判定の `rising` が新判定では **`flat`** になります。全履歴の VERDICT は別列の Inf により `DIVERGED (NaN/Inf)` を維持しますが、**Inf 前の step 1106 までで再判定すると、全列有限・rising 列なし**でした。

   また、メモリ上の反例として「末尾400点で基調が0.3桁上昇し、±0.4桁の交互振動を重ねる」系列では、`D=0.306`、`σ=0.400`、2窓平均比 **1.414** に対して `flat`。初期ピークを含めると列判定まで `PASS` になります。**散らばりが大きいことは、持続的な上昇がない証明ではありません。**

   **対案:** R-h を未完了に戻す。末端の急増検出をジッタ判定から独立させ、緩やかな上昇はブロック集約した傾きとその不確かさで判定する。判定不能を自動的に `flat` にしない。上記実 run の有限な途中履歴、上昇＋振動、定常振動を回帰試験に追加してください。

2. **Major — R-b が検証したレシピと、plan に掲げるレシピが違う。**

   [`startup.md:141`](/home/sano/work/forge/plans/active/tooling-nozzle-sern-startup.md:141) は `soft/mid` を **CFL 0.1** としています。しかし実入力は [`problem_rb_prop_B.yaml:123`](/home/sano/work/forge/case/46.sern_design/problem_rb_prop_B.yaml:123) の **`soft_cfl: 1.0`**。A/C も同じで、保存された `residual_soft_cfl1.csv`・`residual_mid_cfl1.csv` と整合します。

   長時間参照も [`problem_rb_ref_B.yaml:117`](/home/sano/work/forge/case/46.sern_design/problem_rb_ref_B.yaml:117) の本段 **CFL 5.0**、暖機 ramp、`soft_cfl: 1.0` です。したがって今回の比較は、**新しい起動条件での短時間／長時間比較**です。§4.16.2 の旧レシピとの直接比較にはなっていません。

   **対案:** 実測済みの **暖機 ramp → soft/mid CFL 1.0 → 本段 CFL 5.0** を検証対象として明記し、参照も同条件の長時間版と定義する。撤回済みの表を現行仕様として残さず、レシピ・入力・run の対応を一本化してください。

3. **Major — R-f の化学種検査には、NaN を正常扱いする穴が残っている。**

   [`sern_gates.py:30`](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:30) の `FIELD_VARS` に `roY*` がありません。追加された床ゲートも [`sern_gates.py:195`](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:195) で非有限値を負値検査から除外し、組成和でも除外します。

   `run_0195_rb_ref_B/res_12000.h5` の読取り結果を**メモリ上だけで**変更し、`roY1[100]=NaN` とした試験では、`field_health.ok=True`、`floor_gate.ok=True` の両方を再現しました。実ファイルは変更していません。

   **対案:** 設定された全化学種について、存在・形状・有限性を先に検査し、不成立なら即不合格にする。その後に非負性・組成和を検査する。全点判定不能を表す `-1` も不合格にしてください。

4. **Major — R-e の「真因は IC 不整合で確定」は、R-h 後には根拠が成立しない。**

   [`startup.md:173`](/home/sano/work/forge/plans/active/tooling-nozzle-sern-startup.md:173) は旧 IC の `run_0182_sm1_B` と新 IC の `run_0191_mainB_5p0` の FAIL/PASS 差、初期／最終残差比を因果の証拠にしています。

   現判定器で再評価すると、**両方とも `residual_health.ok=True`、VERDICT は `NOT CONVERGED (stalled/plateau)`**。`rms_roY1` の回帰上昇量 `D` は旧 IC が **−0.129桁**、新 IC が **+0.140桁**でした。端点比の「5.4倍対1.5倍」から、旧 IC だけが継続的に上昇していたとは言えません。

   NASA-9 整合化自体は妥当ですが、**それが旧 FAIL の原因を解消したことは未証明**です。

   **対案:** 「熱力学的不整合は修正済み、残差判定差との因果は未確定」に戻す。同じ判定器・実行条件で旧／新 IC を反復比較し、末尾の分布・局所残差・保存収支から判断する。種輸送の切り分けを不要とする結論は保留してください。

5. **Minor — R-b の比較値は概ね再現できるが、記載された精度上限は満たしていない。**

   各 `force_history.csv` の最終値から、`|提案−参照|/|参照|` を再計算しました。run 索引は [case README](/home/sano/work/forge/case/46.sern_design/README.md:244) です。

   | 設計 | 比較した run（`case/46.sern_design/` 配下） | `C_T_with_shear` 相対差 | `C_M` 相対差 |
   |---|---|---:|---:|
   | A | `run_0198_rb_prop_A` / `run_0199_rb_ref_A` | 2.400e-6 | 4.842e-5 |
   | B | `run_0194_rb_prop_B` / `run_0195_rb_ref_B` | **3.413e-6** | **5.956e-5** |
   | C | `run_0196_rb_prop_C` / `run_0197_rb_ref_C` | 2.991e-6 | **6.354e-5** |

   [`startup.md:170`](/home/sano/work/forge/plans/active/tooling-nozzle-sern-startup.md:170) の「`C_M` 相対5e-5以内」と、変更ログの「`C_T_with_shear` 相対3e-6以内」は不正確です。

   再実行した判定は、6 run とも **残差: `NOT CONVERGED (stalled/plateau)`、力係数: `STEADY`**。保存されたメッシュ品質判定はすべて **`VERDICT: PASS`** です。`vis_turb` の平均・ピークも `classify_series` で `STEADY`、提案／参照の最終場の相対L2差は最大約 **9.0e-4** でした。短縮の有望さは支持しますが、これを残差収束の証明にはできません。

   **対案:** 「6桁再現」を実際の誤差値に置き換え、最終値比較か末尾平均比較かを固定する。R-b の残りの場・保存収支・連続窓の検証を完了するまで、生産レシピ確定とはしないでください。

6. **Minor — 正本・写し・現行仕様が食い違っている。**

   [`startup.md:10`](/home/sano/work/forge/plans/active/tooling-nozzle-sern-startup.md:10) は本体を正本としますが、本体の [`chain.md:597`](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:597) は依然 `require_residual_plateau` を既定 True としています。切り出し側でも §4 の旧条件と §5.1 の修正結果が併存しています。

   目的と node・平面・`m6_on` への限定は妥当で、既存の準定常ツールを利用する方針にも問題はありません。ただし、共通の `check_convergence.py` 変更は、この限定範囲を超えて影響します。

   **対案:** 起動・受理条件の正本をこの小計画へ移し、本体はリンク参照にする。旧条件は明示的に失効扱いとし、共通判定器の変更範囲と回帰対象を別途明記してください。

**推奨は、短縮レシピの昇格を止め、受理判定を修復してから既存 run を再評価することです。** 優先順は **R-h の発散見逃し修正 → R-f の非有限値検査 → レシピと正本の統一 → R-e の因果再検証・R-b の残検証**。現状の「R-h 完了」を前提に実装を進めることは承認できません。

ファイルは変更していません。本レビューの提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
