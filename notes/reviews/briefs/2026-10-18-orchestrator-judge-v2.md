# 諮問: `judge_run` の停止規則 v2 対応と、風洞ノズルの判定 (`wt_qualify`) での「解析の成立」の読み方 (2026-10-18)

plan `plans/active/tooling-design-orchestrator.md` §5.1 #5 (§4.2)。関連: `plans/active/tooling-convergence-and-wall-resolution-gates.md` §4.6・§5.1 #12
(`judge_run.py` の正本)、`plans/accepted/tooling-settle-convergence-and-monitor.md` (停止規則 v2・`check_settled.py`)。

聞きたいこと: 下の設計案 A の (1)〜(6) の可否と穴。特に (Q1) v2 の run の `convergence` を `check_settled` に置き換えてよいか、
(Q2) `quasisteady` の量を v2 の目的の量と 1 回の計算で兼ねる形でよいか、(Q3) `wt_qualify` の「解析の成立」に `RUN_VERDICT` を
どの粒度で入れるか (壁解像の FAIL の扱いが `wt_qualify` の補正の分岐と衝突する)。

## 1. 観測事実

作業ツリー `forge-cardmesh` (ブランチ `feature/card-local-mesh`)。2026-10-18 に `origin/feature/generator-defaults-launch-gate`
(forge-gate、`judge_run.py`・投入の関門・判定の仕様 `JUDGE_SPEC.json`) を取り込んだ (マージの commit は本ブリーフと同じ日)。
取り込みの衝突は解いた。forge-gate 側の試験 (`solver_density_cuda/tests/unit/test_judge_run.py`・`test_gate_bad_input.py` ほか) は、
試験の設定にリミッタの基準値を明示して ALL PASS (こちら側の 10-16 の変更で、基準値の記録も明示も無い段は区間を連結しない
`stage_manifest.limiter_ref_entry` → status unknown)。

**実測 (合成の run)**: `design/tests/run_settle_migration_tests.py` の (c) の run `rc_main` (停止規則 v2・REACHED、`run_case.sh` の
`CONVERGENCE_VERDICT.txt` は `TOOL: check_settled.py` で `VERDICT: SETTLED`・CONVERGENCE_EXIT 0) を取り込んだ後の `run_case.sh` で回すと、
`RUN_VERDICT.json` は `state final`・総合 FAIL (`convergence` FAIL = `check_convergence --segment` の低下の桁数; `physical`・`quasisteady`・
`wall_resolution` は合成の run なので判定不能)。基準の無い run `rc_legacy` も同じ形。下の 1. の食い違いの再現である。

取り込んだ後のコードで、停止規則 v2 の run (`settle_criteria.yaml` のある axis-Mach NS、case/71 の run_0011・run_0014 の形) を
`judge_run.py` で判定すると次のようになる (1. は上の実測、2.〜6. はコードを読んで確かめた事実):

1. `convergence` は全 run で `check_convergence.py` (低下の桁数、manifest があれば `--segment`) を回す (`judge_run.py:856` `check_convergence`)。
   v2 の run では `run_case.sh` (`solver_density_cuda/tools/run_case.sh:80-99`) が `check_settled.py --select` で道具を選び、本段は
   `check_settled.py --bind` で判定する。`check_settled` の契約は「低下の桁数は見ない」(`check_settled.py` docstring、plan settle §4.2)。
   したがって v2 の run の `RUN_VERDICT` の `convergence` は `SETTLE_VERDICT` と別の基準で決まり、整定 (SETTLED) した run が
   低下の桁数の不足で FAIL になりうる。
2. 停止規則 v2 の状態は `REACHED`・`CENSORED`・`FAILED` (`design/forge_design/evaluate/stop_rule.py:1081` `V2_STATES`)。`judge_run.py:151`
   `STOP_STATUS` に `FAILED` が無いので、v2 の `FAILED` は「既知でない」で INDETERMINATE になる。
3. 生成器の仕様: `runner_axismach.run_staged_ns` (`design/forge_design/evaluate/runner_axismach.py:3154`) は
   `J.spec_axismach(..., stop_rule=(sr if sr_use else None))` を書く。v2 の run では `sr_use = False` なので、仕様の量は
   出口コア M (criterion relative = `check_quasisteady` の既定) と質量の収支 r (criterion absolute) で、系列 `judge_series.csv` は
   本段の `res_<N>.h5` から作る (`judge.write_axismach_series`)。v2 の目的の量 (`settle_criteria.yaml` の `exit_core_mach`・
   `mass_balance`・`line_mach`・`wall_heat_nodes`・残差の見張り) とは定義・窓・閾値が別 (例: 窓が末尾 4 割 vs W 20,000 step の格子)。
4. 鮮度の照合 `collect_inputs` (`judge_run.py:1341`) は `settle_criteria.yaml`・`STOPPED_AT`・停止の時の判定の専用ファイル
   (`stop_rule_v2/SETTLE_VERDICT_stop_<step>.json`)・流量の出力を入力に入れていない。
5. `wt_qualify.py` の「計算の成立」(`design/forge_design/feedback/wt_qualify.py:288` 付近) は、停止規則 v2 の `REACHED` と停止の step の
   再評価 `SETTLED` を、停止の時の判定の記録との結び付き (`deltastar_loop._v2_binding`) で確かめる。承認のモード (#2c) では
   判定専用の checkout で `check_settled --end-step <停止> --no-write` を回し直す (`approval_rechecks`、同 `:853`)。
   `RUN_VERDICT.json` (非有限・物理性・solver の終了コードなど) は読んでいない。
6. `wt_qualify` の分岐 (`:946` 以降): 判定不能・計算の成立の FAIL → `incomplete` (escalate)。公差・壁解像・寸法の FAIL →
   `back_to_correction` (壁の形だけなら較正し直し → recalibrate、壁解像を含めば `mesh_first_layer` の補正を挙げて escalate)。
   `judge_run` の総合も壁解像を含む (`check_wall_resolution --step N --groups <全壁> --weight area`、真理値表 FAIL > INDETERMINATE > PASS)。

## 2. 期待値と出典

- plan orchestrator §5.1 #5 の合格の条件 (事前に書いたもの): 「run_0011 で RUN_VERDICT が SETTLE_VERDICT と整合、`SETTLED` でも壁解像 FAIL なら
  総合 FAIL、公差未達でも整定している例・必須の CAD の欠けの例の試験」。
- plan orchestrator §4.2: 判定を 要求の一致 / 解析の成立 (`judge_run`・停止規則 v2・壁解像・メッシュ品質・`launch_gate`・受け渡しの関門) /
  要求の達成 (qualification) に分ける。公差は `judge_run` の外。
- plan judge §4.6: `judge_run` は既存の判定ツールの共通の集約層で、数値の判定基準は各ツールのまま、設計の受理の条件は各ゲートに残す。
  §5.1 #14b (2) では停止規則 v1 の量を `stop_rule.judge_rows` そのもので判定する形 (criterion `stop_rule`) にした (定義を二重に持たない)。

## 3. 再現条件

- コード: `forge-cardmesh` の `feature/card-local-mesh` の HEAD (取り込みの後)。
- 対象の run: `case/71.m55_combustion_wt/run_0011_ns_v2` (停止規則 v2、REACHED・SETTLED)、`run_0014_chain_v2` (V6、REACHED・SETTLED、壁解像 PASS)。
  どちらも forge-gate より前に回したので `RUN_VERDICT.json`・`JUDGE_SPEC.json` は無い。

## 4. 実施済みの操作と結果

- 取り込み・衝突の解消・試験 (上の 1.)。`check_wall_resolution.py` は両側が同じ欠陥 (`--step`・`--groups` の欠けを最新・見つかった壁で
  PASS) を別々に直していたので、こちら側 (10-17: 対象の壁と step を評価の前にそろえる、`--groups` 無しは bcondConfig の no-slip 壁を全部必須・
  kind が wall* でない出力は対象外) を残した。forge-gate 側の試験「`--groups` なしは出口の出力を拾って判定不能」は、こちら側の仕様に合わせて
  「出口は対象外・壁 2 枚で PASS」に直した。
- `judge_run` の v2 対応はまだ書いていない (この諮問の後)。

## 5. 設計案 A (仮説)

(1) **道具の選択**: `judge_run` の `convergence` は `settle.select_judge(run)` (= `check_settled.py --select`、`run_case.sh` と同じ規則) が
`check_settled.py` を選んだ run では `check_settled.py <run> --no-write --json` (照合は既定で有効) を回し、終了コード 0 SETTLED → PASS、
1 NOT_SETTLED → FAIL、2 DIVERGED → FAIL、3 INDETERMINATE → INDETERMINATE。`check_convergence` へは落とさない。`check_convergence.py` を
選んだ run (基準の無い run・前の段と確かめられた run) は今のまま。終端の step は `check_settled` の既定 (STOPPED_AT、無ければデータの最後)。

(2) **停止の時の判定との整合**: 停止規則 v2 の run (STOP_RULE_STATE.json の method v2、reeval あり) では、(1) の判定し直しを
`stop_rule.stop_verdict_record` で結び付けた停止の時の判定と比べる (総合と各量の verdict)。食い違えば `convergence` を INDETERMINATE
(「停止の時の再評価と判定し直しが食い違う」; 判定器の版の違いなど)。結び付かない (専用ファイルの sha256 の不一致など) も INDETERMINATE。

(3) **準定常の量**: JUDGE_SPEC に criterion `settle` を足す。生成器 (`spec_axismach`) は v2 の run では `settle_criteria.yaml` の目的の量
(残差の見張りを除く) の名前を仕様の量として書き、`judge_run` は (1) と同じ 1 回の `check_settled` の結果の各量の verdict を写す
(2 回目の計算をしない; v1 の criterion `stop_rule` と同じ「定義を二重に持たない」)。仕様の量と基準の量の名前の完全一致を照合 (欠け・余りは
INDETERMINATE)。`exit_core_M`・`r` の今の量は v2 の run では書かない。

(4) **停止規則の状態**: `STOP_STATUS` に `FAILED` → FAIL を足す (v2 の状態)。

(5) **鮮度**: `collect_inputs` に `settle_criteria.yaml`・`STOPPED_AT`・`MAIN_STAGE_STARTED`・停止の時の判定の専用ファイルと、
`check_settled --json` が返す「使ったファイルの一覧と sha256」(流量の出力 CSV など) を足す。

(6) **wt_qualify の解析の成立**: 今の v2 の結び付き (停止の時の REACHED・SETTLED) に加えて、`RUN_VERDICT.json` の `state final`・
`check_fresh` = FRESH と、項目 `solver`・`nonfinite`・`physical`・`convergence`・`stop_rule` の PASS (または NOT_APPLICABLE) を要求し、
どれかが FAIL・INDETERMINATE・無い → `incomplete`。**`wall_resolution` と `quasisteady` は `RUN_VERDICT` から読まない** —
壁解像は `wt_qualify` が今どおり自分で回して補正の分岐 (`mesh_first_layer`) に使い、準定常の量は (3) で `convergence` と同じ判定なので
二重に数えない。承認のモードでは `RUN_VERDICT.json` は作業側のコードの自己申告なので、判定専用の checkout で `judge_run.py` を
別の出力先へ回し直す (`--out` を足す) か、項目ごとに既存の再検査 (`approval_rechecks`) で足りるかを決めたい。
`RUN_VERDICT` の総合 (壁解像を含む) は「`SETTLED` でも壁解像 FAIL なら総合 FAIL」を満たすが、wt_qualify の結果は壁解像の FAIL を
`back_to_correction` に送る (今どおり)。

迷っている点:
- (6) で `RUN_VERDICT` の総合をそのまま使わず項目を選ぶのは、`judge_run` の真理値表を迂回することにならないか。総合 FAIL でも
  `wt_qualify` が壁解像の補正へ送るのは §4.2 の「解析の成立が通ったときだけ要求の達成を評価する」と矛盾しないか
  (壁解像は §4.2 の表で「解析の成立」の側にも「要求の達成」の `accepted` の条件にも出てくる)。
- (1) で SETTLED を `convergence` の PASS と呼んでよいか (`check_settled` は残差の低下を見ない; 項目名と中身のずれ)。
- (3) で `line_mach` は wt_qualify の V5c (線の一様性の公差) と同じ線を監視するが、`judge_run` の準定常は公差を見ない (整定だけ) ので
  役割は分かれていると考えている。
