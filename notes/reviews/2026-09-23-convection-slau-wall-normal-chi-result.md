# codex レビュー: convection-slau-wall-normal-chi (result)

- **plan**: [`plans/accepted/convection-slau-wall-normal-chi.md`](../../plans/accepted/convection-slau-wall-normal-chi.md)
- **stage**: `result` (diff base `39526328`)
- **date**: 2026-09-23
- **commit**: `cec9c787` (feature/sern-design)
- **codex**: effort `high`, 6.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **focus**: 実装 diff (9 ファイル 106 行) と §6.2 の実測を突き合わせること。特に (1) 効果の帰属: 本計画の chi_n と、別 plan で同時に直した出口 BC (outlet_statPress -> outflow, tooling-sern-mesh-blocking B1f) を分離できているか。当方の切り分けは「step を揃えた比較 (対策なし step 800 / 対策あり step 1000) で出口の亜音速率が 7.1 % 対 9.0 % とほぼ同じ」の 1 点に依っている。(2) C_L 0.433 % / C_M 0.316 % がノイズ床の 5-6 倍で実効果と判明した件を、判定対象外として扱ってよいか。既定化の条件として何を要求すべきか。(3) V1-d を run_0439 の check_quasisteady ALL STEADY で合格としてよいか (check_convergence は同区間で stalled/plateau だが rms_ro 2.40e-10 = 通算 4.8 桁低下)。(4) 一次資料 Shima-Kitamura 2011 §III.K が入手できないまま実装を先行させた判断 (既定 0・ビット同一の opt-in) の是非。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
`chi_n` の CUDA 実装は概ね設計どおりですが、V1-c の診断コードが変更後の流束を計算していません。  
出口 BC との効果分離と V3 の受入条件も未完了であり、現状の「検証完了」は支持できません。

指定 diff を HEAD `cec9c787` まで確認しました。現在は **10 ファイル、98 行追加・7 行削除**です。質量流束だけの変更、再構成後の速度、壁 mask、境界半割面の除外、既定 0、config の拒否条件は設計と整合しています。追加部分に明白な符号・単位・絶対閾値の誤りは見つかりませんでした。周期・軸対称の実測確認はありません。

case/48 は別 worktree の実データで再検証でき、V2 は **`VERDICT: PASS`**。V3 は両側とも **`NOT CONVERGED (stalled/plateau)`**でした。11 スナップショットの各 5 地点の `Cf`・`qw` は、既存評価関数と `check_quasisteady.py` の判定関数により、指定閾値で全て **`STEADY`**。全スナップショットの `VALUE/*` に非有限値はありませんでした。case/46 の対象データはローカルに見つからず、SSH は sandbox に拒否されたため、その数値は独立に追認できていません。

1. **Major — V1-c は別の流束で判定されており、合格根拠が無効です。**

   **根拠:** [diag_wall_cv_budget.py:65](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_wall_cv_budget.py:65) は全速度から従来の `chi` を計算しています。`--equilibrium` も [同ファイル:182](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_wall_cv_budget.py:182) でその関数を呼びます。対して実装は [convectiveFlux_slau_d.inc.cuh:552](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:552) で `chi_n` を使います。

   診断関数を直接実行した反例では、壁側速度 0、内点接線速度 1000 m/s、両側音速 700 m/s、面法線速度 0、圧力差 900 Pa、面積 1 m² に対して、診断は `mdot=0`、変更後の式は `mdot=−0.642857 kg/s` です。検証したい補充項そのものが診断から抜けています。

   さらに、[plan:218](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:218) は「STEADY 後の最後の 3 dump」を要求しますが、結果は上昇途中の `run_0437/res_6000` だけです。V1-b の正式結果もありません。密度が床より上であることだけでは、`pMin`・温度処理まで不活性とは証明できません。

   **対案:** V1-c の PASS を撤回し、診断を実効 config・`chi_n`・前処理に合わせる。実カーネルの面流束との照合後、起点直後の V1-b と、定常化後 3 dump の V1-b/c を再判定してください。

2. **Major — 出口 BC と `chi_n` の独立性は、現在の比較では証明できません。**

   **根拠:** [plan:332](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:332) の step 800 対 1000 は同一時点ではありません。亜音速率 7.1% 対 9.0% の一点比較から「加速していない」とは判断できません。

   `run_0435` 対 `run_0437` は **flag 1 における BC 変更効果**を調べています。旧 BC での flag 0/1 比較と合わせても、**`outflow`＋flag 0** が欠けており、修正後 BC でも `chi_n` が必要か、両変更の相互作用があるかは未確定です。

   また、[runner_sern.py:279](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:279) と 3D runner は既定 BC を変更しています。両方 `outflow` の V3 が完走したことは、[plan:445](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:445) の「出口既定変更が生産レシピを壊していない」の比較根拠にはなりません。変更自体は別 plan B1f と procedures に明記されていますが、回帰の実証は別問題です。

   **対案:** 同一バイナリ・同一 IC で **BC 2 種 × flag 2 値**を揃え、共通 step／累積 CFL 区間で壁密度・収支・出口量を比較する。まず不足する `outflow`＋flag 0 を実施し、出口既定変更の回帰は flag を固定して別に判定してください。

3. **Major — V3 の一部だけを満たして、受入全体を PASS にしています。**

   **根拠:** [plan:238](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:238) 以降で必須にした次の結果が、§6.2 にありません。

   - case/16 の V2・V3。
   - SERN の衝撃足壁圧分布と、事前に固定した比較ノルム。
   - SERN の力係数時系列に対する指定閾値の準定常判定。
   - 格子感度と累積 CFL あたりの収束比較。
   - `cfl_pseudo: 0.2 → 0.4` による固定点確認。

   SERN の 3 本ずつの反復は再現性を調べますが、時間方向の未収束バイアスは除去しません。`C_T` の終端差が許容内でも、これらのゲートの代わりにはなりません。V0 の結果欄も config 拒否試験が中心で、面反転・等状態・非対象面の試験記録を追えません。

   §5.1 は case/16 を低優先度へ下げ、V1-b、衝撃足、格子感度、CFL 試験を独立した未完了項目として保持していません。周期・軸対称試験も「result 段まで」の約束が残っています。

   **対案:** 不足項目を §5.1 に戻し、V3 総合 PASS を保留する。指定したゲートの全結果、入力、抽出コマンド、VERDICT を揃えて再レビューしてください。

4. **Major — V1-d の準定常判定を、方程式全体の収束証明に置き換えています。**

   **根拠:** [plan:382](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:382) は `NOT CONVERGED` に対して、`rms_ro` の絶対値だけから「収束しきってノイズ床」と断定しています。記載値 `6.02e−6 → 2.40e−10` の低下は **4.40 桁**であり、4.8 桁でもありません。

   case/48 も「収束済みの場から」と説明していますが、起点
   `/home/sano/work/forge/case/48.flat_plate_cooled_m4/run_0025_B_tw300_y3_fx05/`
   の再判定自体が **`NOT CONVERGED (stalled/plateau)`**でした。未収束の参照をノイズ床と呼ぶ根拠が不足しています。

   **対案:** **V1-d 単独は、要求した密度系列が STEADY かつ全壁の床到達 0 を確認できれば合格にして構いません。** ただし「全方程式が収束した」とは分けてください。`run_0437–0439` の実効設定同一性を確認し、BC 変更後の同一区間を全残差列で判定する。ノイズ床の主張には別の裏付けが必要です。`run_0439` の省略していない密度系列・VERDICT も提示してください。

5. **Major — 「剥離縁だけに効く」という説明は実装から成立せず、力係数変化の扱いを誤らせています。**

   **根拠:** mask は [convectiveFlux_slau_d.inc.cuh:549](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:549) の壁フラグだけで、剥離を検出していません。付着境界層でも、壁法線方向の面に接線速度差があれば `chi_n` と `chi` は異なります。流束差が小さい理由としては、[plan:298](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:298) の **圧力差が小さい**という説明が適切です。

   `run_0437` の面数比 0.53% は位置を証明せず、その 3D ケースの比率を SERN 2D に転用もできません。[plan:439](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:439) の「そこは元々 `chi=0`」も、単に `chi` が変わった面数からは導けません。

   **対案:** 適用範囲を「全壁隣接内部面」と正確に記述し、ケースごとに `Δchi`・圧力差・流束差・壁圧寄与の位置を確認する。`C_L` 0.433%／`C_M` 0.316% は事前の `C_T` ゲート違反にはなりませんが、**生産利用上の精度評価から外してはいけません**。既定化には、これらを含む許容誤差の事前設定、準定常確認、少なくとも 3 水準の格子検証、独立した検証基準との比較を要求すべきです。

6. **Minor — 結果台帳と数値の出所が同期していません。**

   **根拠:** [plans/README.md:62](/home/sano/work/forge-sern-design/plans/README.md:62) はまだ「本命は未確定」。[case/46 README:342](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:342) は `run_0439` が実行中、次行は反復 b/c が未完了の記述です。case/48 の `_v2`・`_v3` も、この checkout の run 索引から追えません。

   case/48 の保存済み `cooled_plate_eval.json` における 5 地点の最大差は `Cf` **0.001804%**、`qw` **0.002388%**で、[plan:295](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:295) の 0.000%／0.011% と異なります。許容内という判断は変わりませんが、使用量・抽出方法を特定できません。

   **対案:** 台帳と計画一覧を同期し、評価定義・成果物の実パス・抽出コマンドを残す。`methods/index.md` の既存リンクは有効でしたが、methods の仕様説明には非周期・非軸対称という検証範囲も追記してください。

7. **Minor — 原典未入手の扱いを更新し、文献と今回の変更を式単位で区別してください。**

   **根拠:** [plan:59](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:59) は原典未読ですが、今回、Shima–Kitamura の**著者公開版**で §III.K と Fig.20 の説明を閲覧できました。面法線 Mach の比較は記載されていますが、Fig.20 の説明は運動量式の圧力項を指しており、今回の「壁隣接面・質量流束のみ」と同一の実験ではありません。[著者公開本文](https://www.researchgate.net/publication/258474939_Parameter-Free_simple_Low-Dissipation_AUSM-Family_scheme_for_all_speeds)

   **対案:** §3.1 と §5.1 #2 を更新し、変更する式・適用面・試験条件の相違を明示する。**既定 0 の実験機能として実装を先行した判断自体は許容できます。** ただし、文献の懸念を承知したうえで設定した格子感度・衝撃足の検証を省略する理由にはなりません。2023 年論文の抄録も収束悪化と格子形状への感度を報告しています。[Furusawa–Kitamura 2023](https://onlinelibrary.wiley.com/doi/10.1002/fld.5183)

**推奨は、既定 0 のまま `active` に留め、診断修正 → BC×flag の効果分離 → V1/V3 の不足ゲート完了、の順で再レビューすることです。** 現段階で `accepted` へ移すべきではありません。本回答はレビュー提案であり、ファイル変更・plan 反映はしていません。

指摘数: Critical 0 / Major 5 / Minor 2
