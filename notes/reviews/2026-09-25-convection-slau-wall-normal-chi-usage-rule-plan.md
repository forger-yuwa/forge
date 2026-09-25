# codex レビュー: convection-slau-wall-normal-chi-usage-rule (plan)

- **plan**: [`plans/accepted/convection-slau-wall-normal-chi-usage-rule.md`](../../plans/accepted/convection-slau-wall-normal-chi-usage-rule.md)
- **stage**: `plan`
- **date**: 2026-09-25
- **commit**: `615eafda` (feature/sern-design)
- **codex**: effort `high`, 7.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m1
- **extra**: `plans/accepted/convection-slau-wall-normal-chi.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

既定 0 の維持と、未解決の 2D 感度評価を独立させる目的は妥当です。
ただし、Roe 比較の物性条件と格子細分化の手順が成立せず、過渡を合格扱いできる判定も残っています。
現案のまま計算を投入しても、目的とする「正否」の根拠にはなりません。

1. **Major — Q3 の Roe は、対象の二成分 TP と同じ物理条件を計算しない**

   対象は `EXH/AIR` の二成分ですが、Roe は左右の面状態を **`sp[0]` の単成分物性**で再構成します。Roe 平均状態も `thermo_T_from_h(sp, 1, &Yone, ...)` であり、混合組成を使いません。根拠: [Roe の左面状態](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_roe_d.inc.cuh:159)、[Roe 平均状態](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_roe_d.inc.cuh:242)。

   したがって、[plan §4.4](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:70) の比較には流束方式以外の差が混入します。壁 CV への流束が初手で反転することも、力係数の精度を保証しません。また、同一物理の比較ができたとしても、最細格子 1 点で Roe に近いだけでは正しさを決められません。

   **対案:** 現行 Roe による Q3 を外し、独立精度評価は保留する。同じ混合物性・BC を扱う参照手法と、その格子・反復誤差を確認できるまで、Q1 の結果を「正否」に読み替えない。Roe の多成分化を行うなら、本件の「コード変更なし」とは別の計画にする。

2. **Major — 2D メッシャの `scale` は格子密度ではなく物理寸法を変える**

   [plan §4.2](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:58) は `scale` で細分化するとしています。しかし runner は **`scale=H`** を渡し、メッシャは接続を作った後に **`coords *= prm.scale`** を実行します。根拠: [runner](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:495)、[メッシャ](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern.py:325)。

   この値を半分にしても節点数は増えません。物理寸法を縮めれば、今回の粘性計算では Reynolds 数も変わり、固定形状の格子比較になりません。

   **対案:** `H_m`・形状・計算領域を固定し、`ni_*`・`nj_*` と第一層幅などの**実際の解像度キー**で 3 水準を定義する。生成コマンドだけでなく、節点数・局所接線間隔・第一内部ノード距離・形状不変の確認を記録する。新格子への場の移植元と `interp_field.py` による restart 手順も、§5.1 に明記する。

3. **Major — Q4 は未定常な過渡を「限界サイクル」に変更して通せる**

   [plan §4.5・§6](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:75) の「一度倍にしても STEADY でなければ限界サイクル」は成立しません。判定器は `DRIFTING`、`OSCILLATING`、`TRANSIENT-UNSETTLED` を明確に分けています。根拠: [判定実装](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_quasisteady.py:293)。単調増加する人工系列を実際に判定すると、長さを倍にしても `DRIFTING` のままでした。

   また、§4.1 の「2D 生産は flag 0 で収束」は既存記録と矛盾します。保存された本段 CSV を再判定した結果は、両側とも次のとおりです。

   `case/46.sern_design/_r3_m1m3/sern2d_conv/flag{0,1}/residual_history.csv`  
   **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**

   **対案:** `DRIFTING`／`TRANSIENT-UNSETTLED` は延長または判定保留とし、予算上限に達してもラベルを変更しない。`OSCILLATING` の場合だけ、平均・振幅の安定性を確認して比較する。準定常確認は Q1 の `C_L/C_M` だけでなく、Q2・Q3 を含む比較対象の全 run と、`C_T`・衝撃足の量にも適用する。「残差未収束」と「目的量が比較可能」を区別して記録する。

4. **Major — Q2 の `D_inter` だけでは CFL 独立性を示せない**

   [plan §6 Q2](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:106) は、累積 CFL を揃えた run 間差だけを判定します。両者が同じ過渡を追えば、この差は小さくなります。前計画はそのために **各 run 自身の `D_intra`** を要求していました。根拠: [前計画 V7](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:469)。

   `D_inter(1) ≤ 2 D_inter(0)` も、flag 0 の差が丸め床付近なら不安定な判定です。両側が ε を十分下回っていても落ち得ます。

   測定対象にも注意が必要です。既存 [v7_dist.py](/home/sano/work/forge-sern-design/case/46.sern_design/cad/v7_dist.py:13) は **3D 専用の壁 physID を固定**しており、2D の ID 1–3 は入口・出口です。組成・乱流列も読みません。前計画自身が組成の遅いモードを未解決として記録しています。

   **対案:** 各 CFL の独立した延長窓で `D_intra` と比較量の準定常を確認した後、`D_inter` を評価する。比は検出限界付きの補助情報に留める。壁集合は実際の BC から作り、第一内部ノード・組成・乱流量の監視条件を測定前に固定する。

5. **Major — Q1 の分岐は、格子極限について誤った結論を出し、未定義の場合もある**

   [plan §6 Q1](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:105) の比率 `>0.7` は、異なる極限への収束を意味しません。例えば

   \[
   \Delta(h)=0.04h^{1/4}
   \]

   は細分化でゼロになりますが、指定された 3 水準では `0.0400 → 0.03668 → 0.03364`、fine/coarse は **0.841**。`C_M` 許容 0.02 を超えるので、現規則はこれを「flag が収束解を変える」と判定します。

   逆に `0.10 → 0.06 → 0.03` は、単調減少・許容超え・比 0.3 で **(a)/(b)/(c) のどれにも入りません**。両フラグに共通する格子誤差は、フラグ間差だけでは検出できません。また、衝撃足 L2 は測定量にあるものの、分岐条件に反映されていません。

   **対案:** 各フラグの格子間変化とフラグ間差を別々に評価する。有限の 3 水準から言える結論は「検証範囲で許容内」「差が残る」「判定不能」に限定し、全ケースを覆う分岐と複数量の集約規則を定義する。衝撃足の壁圧・位置にも合否と未達時の処置を明記する。

6. **Major — 適用規則の条件が壁 CV 排出を識別していない**

   [plan §4.1](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:53) の条件は「床到達 > 0 **または** 高速な壁隣接内点」です。しかし、[`--summary` の実装](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_wall_cv_budget.py:141) は**全域**の床到達数を数えており、壁 CV に限定していません。

   速度条件も排出の十分条件ではありません。接線速度が大きくても、法線速度と圧力差がゼロなら、その面の質量流束はゼロです。実装でフラグが変える項は **`(χ_n−χ)ΔP`** であり、再構成後の面状態を使います。根拠: [SLAU 実装](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:541)。

   **対案:** 現在の二条件は「診断を始める兆候」に下げる。適用判断には、対象壁 CV の密度時系列・全接続面の正味流出・実効面状態での補充項の消失を使う。既知の 3D 接続構成以外、とくに周期・軸対称の生産構成へ一般化しない。

7. **Major — 引き継ぐ振幅評価は、両側の差の不確かさを過小評価する**

   [plan §4.2](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:61) が参照する #10b 式は、不確かさを **両 run の半幅の大きい方**と定義しています。根拠: [前計画 #10b](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:942)。

   位相が独立なら、差の変動には両側が寄与します。例えば平均差 0.015、半幅 0.003／0.004 では、現式は **0.019 で合格**ですが、対称な変動帯どうしの差の上限は **0.022** となり、`C_M` 許容 0.02 を超えます。

   **対案:** 瞬時変動帯を比較するなら両側の区間差を使う。時間平均の差を比較するなら、窓間の平均の安定性と平均値の不確かさを評価する。振幅を平均値の誤差と同一視せず、採用する比較式を測定前に固定する。

8. **Minor — `stage_manifest` は SLAU→Roe の切替を区間として分離しない**

   [plan §6](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi-usage-rule.md:109) は全 run を `--segment` で判定しますが、[`stage_key`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:184) の hard キーには `solver` がありません。同一 config の `solver: SLAU` を `ROE` に変えた呼出しで、**`stage_key(SLAU) == stage_key(ROE)`** を再現しました。

   **対案:** 将来 Q3 を実施するときは、Roe 側だけの独立した履歴・manifest を作り、SLAU の履歴を連結しない。共通ツールで保証するなら `solver` を hard キーに追加して検証する。

**推奨は、既定 0 と既知構成での限定的な opt-in を維持し、上記を反映した計画を再レビューしてから Q1・Q2 を実施することです。** Q3 による正否判定は、同一物理を扱える参照手法が成立するまで保留してください。3D の固定点確認も、2D Q2 で解決した扱いにせず、委譲先の未解決事項として残すべきです。

指定どおりファイルは変更していません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
