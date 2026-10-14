# codex 諮問 (diagnose): speedup-tests

- **brief**: [`notes/reviews/briefs/2026-10-12-speedup-tests.md`](../../notes/reviews/briefs/2026-10-12-speedup-tests.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-12
- **commit**: `cb199505` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.1 min, rc=0
- **結論**: **試験Pの最初の一手として、G0・共通起点・point・FP64で `precondEps=1` を両腕に固定し、`lowMachPrecond=0/3` の1 step比較により、前処理経路への切替が境界更新を壊さないか確認する。**
- **extra**: `plans/accepted/time_integration-lowmach-preconditioning.md`, `procedures/solver-settings.md`, `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | P を現状のまま20,000 step実行：**要再検証** | 通常の point 経路は node 境界半割面の粘性対角を除外するが、前処理経路は除外しない。前者には「退化距離による巨大対角→境界ノード凍結」の説明もある。[timeIntegration_d.cu:957](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:957)、[同:1342](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1342)。使用予定の元 commit `0b4dff4e` にも同じ差を確認した。**まず前処理を数学的に無効にした β=1 の比較で、経路切替自体の影響を測る。** |
| **Major** | 「FP64 同士なので前処理だけの比較」：**却下** | P0 の既定 `implicitSolvePrecision=0` は `ST=float`、P3 の行列・解法は `flow_float` なので FP64 ビルドでは double。さらに P3 には通常経路の軸・壁の拘束行処理が渡されない。[同:1764](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1764)、[同:1313](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1313)、[同:1122](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1122)。**設定一式の効果の比較と、低マッハ前処理への因果帰属を分ける。** |
| **Major** | P の「停止値への残りが半分なら効く」：**却下** | 保存済み JSON を再計算すると、`case/45.isobutane_m6_d155/run_0713_pt_d1/` の θ_r(40) は22,500 stepで0.077222607、25,000 stepで0.077078560となり、T2参照値0.077158535を通過する。その後30,000 stepでは参照値より−0.421%。**近接は過渡の通過でも起こる。** [rpt_judge_s1.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rpt_judge_s1.json)。距離は参考値に下げ、最終的な加速は同じ収束・準定常要件までの時間で判定する。20,000 stepは有限期間の応答試験に限る。 |
| **Major** | G の「差≤0.1%なら答え不変・AR受入条件成立」：**要再検証** | 参照の `run_0752_t2_d1/`・`run_0754_t2_d2/` 自身が未収束・DRIFTING。停止規則は残存反復誤差を上限0.1%に抑える保証ではない。[t2_judge.json:38](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/t2_judge.json:38)。現案で言えるのは「登録停止時の量の差」まで。**答え不変には両格子の収束・準定常確認が必要**で、壁熱流束分布・δ_loc・出口量にも採否基準を事前設定する。AR緩和の結論は今回のG0c・条件に限定する。 |
| **Major** | 同一機械で並行実行：**数値応答の収集は条件付き採用、速度判定は却下** | plan §6.11 も共有GPUの経過時間を参考扱いとし、採用判断には専有時の費用を使う。[float plan:491](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:491)。PとGを同じGPUで走らせた時間から加速率を出さない。**今回の短い確認は専有・逐次を推奨**。後の並行収集には資源余裕と同時負荷の記録が必要。速度は native・GPU専有で測る。 |

結論: **試験Pの最初の一手として、G0・共通起点・point・FP64で `precondEps=1` を両腕に固定し、`lowMachPrecond=0/3` の1 step比較により、前処理経路への切替が境界更新を壊さないか確認する。**

第 1 仮説: **P3では、低マッハ前処理そのものとは別に、node境界処理の欠落が更新を抑制・変形する。** 確度: **中**。処理の欠落はコードで確認済みだが、今回の入力での影響量は未確認。

  根拠: 上表の境界半割面処理と拘束行処理の差。P3の関数には `isNode`・壁／軸フラグの引数自体がなく、通常経路と同じ処理を実施できない。

  反証条件: β=1で、境界・軸・内部それぞれの拘束条件と自由成分の補正が、再実行の幅を考慮して実用上同等であり、境界に系統的な補正抑制が認められないこと。

第 2 仮説: 境界処理の問題を除けば、低マッハの擬似時間剛性が主要な律速で、P3が改善する。確度: **低**。case/23の `=2` はRHSも変えるため、case/45の `=3` を支持する直接証拠にはならない。

判別 A/B:

- 両腕とも同じ実バイナリ・G0・初期場・point設定・`precondEps: 1.0`。変更は **`time.deltaT.lowMachPrecond: 0` 対 `3` だけ**。各1 step、再実行を含め各2本。`nStepInner` は元設定を保つ。
- `precondEps=1` なら β=1、Γ=I、擬似時間倍率=1になる。[lowMachPrecond_d.cuh:18](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/lowMachPrecond_d.cuh:18)、[setDT_d.cu:268](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:268)。
- 初期化後の9保存量、更新前RHS、`dt_local` を照合し、**BCの再適用で隠れる前の最終補正 `dq`** を境界・軸・内部別、成分別に比較する。
- 判別用の暫定閾値を事前登録する。**A：**境界の自由成分の補正がP0の1/100以下に抑制されるなど、再実行幅の10倍を超える系統差がある → 「切替は純粋な前処理操作」という前提を退け、境界経路を先に修正する。**B：**全評価領域で補正の相対L2差が1%以下か再実行幅内で、拘束違反もない → 今回の入力で重大な境界更新障害が出るという第1仮説を退ける。中間・欠測は判別不能。Aでも境界欠落への単独帰属には、精度差などの影響を分ける必要がある。

やらない方がよいこと: **P3の発散を直ちに「低マッハ前処理は効かない」と扱うこと、停止値の通過を加速と呼ぶこと、PとGの共有GPU時間から速度順位を決めること。** Gを継続場から始めるなら、先行30,000 stepとその費用を総到達費用に含める。

呼び出し側の前提への異議: 「RHSを変更しない」はコード上支持されるが、「停止した場が同じ定常解」「同じstepの差は進み方だけ」とは言えない。参照FP64の保存済み判定は、再開後122,500〜142,500 stepについて次のとおり。

```text
check_convergence: OVERALL: CHECK FAILURES ABOVE
check_quasisteady: θ_r(40/70/94), Q_w は各量とも DRIFTING
                  OVERALL: NOT ALL STEADY
```

不足情報: 対象runの生CSV・HDF5と実バイナリはこの作業環境にないため、数値は保存済み判定JSONから検算した範囲である。実バイナリの境界挙動、Gの局所量の許容値、実行機械のGPU・メモリ空き状況は未確認。ファイル変更・forge実行は行っていない。**本提案はplan未反映**であり、呼び出し側が `architecture-float-state-double-geometry.md` §5.1 #25と検証節へ反映する。
