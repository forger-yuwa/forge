# codex 諮問 (diagnose): precond-port

- **brief**: [`notes/reviews/briefs/2026-10-12-precond-port.md`](../../notes/reviews/briefs/2026-10-12-precond-port.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-12
- **commit**: `c11a1cc3` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.2 min, rc=0
- **結論**: **出口修正を取り込み済みとして移植範囲を訂正し、未取り込みの node 境界処理を現行仕様に合わせて移植したうえで、その修正だけの前後を P3・β = 1 の 1 step で比較する。**
- **extra**: `plans/accepted/time_integration-lowmach-preconditioning.md`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/rlg_judge.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 5 / Minor 1）

| 重大度 | 対象・判定 | 根拠と対案 |
|---|---|---|
| **Major** | TP 出口の「未移植」判断：**却下** | **既に `37e2e6265` で移植済み**。現行の [boundaryCond_d.cu:606](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/boundaryCond_d.cu:606) は `ga_loc` を内部・出口音速、等エントロピー関係、Riemann 不変量に使用し、`:701` で `gamma` を渡している。試験元とされる `b0303547` にも同じ実装があり、同ファイルの HEAD との差分は空。**出口修正は「取り込み済み」と記録し、今回の R1 から出口節点の差分許容を外す**。 |
| **Major** | node 境界処理の移植：**採用。ただし旧コードの逐語的コピーは却下** | 境界半割面の粘性対角除外と拘束行処理は必要。ただし現行の [timeIntegration_d.cu:1737](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1737) は `nodeIsothermalEnergyBC != 1` のときだけ等温壁エネルギー行を拘束する。元コミットの無条件の `iso_wall_flag` 渡しをコピーすると、弱形式のエネルギー更新を止める。**現行標準経路の適用条件を正本とし、弱形式ではエネルギー行を残す。既存の `ge_*` も保持する**。 |
| **Major** | R2「修正後は P0 と全域 1% 以内」：**要再検証** | 境界欠落以外にも差がある。標準経路の TP 軸ソース `[2][0]` は [timeIntegration_d.cu:1085](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1085) の式、前処理経路は [:1382](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1382) の `chi_hoop` を追加する式。さらに P0 の既定 solve は `float`、P3 は `flow_float`。**境界処理の修復と、P0/P3 の演算子・精度の同等性を別の判定項目にする**。β = 1 では rank-1 係数がゼロなので、`g[row]=0` の正しさも検証できない。β < 1 の行列レベルの拘束確認を追加する。 |
| **Major** | `speciesPrecondDt` の移植と回帰：**採用、現案の検証範囲は不足** | 現行では [speciesTransport_d.cu:609](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/speciesTransport_d.cu:609) の `dplurSweepOnce` を化学種と受動種が共有し、周期・非周期の両経路がある。`:1559` は受動種からの呼び出し。ここを一律に `dt_local_sp` に替えると受動種まで変わる。**刻みポインタを呼び出し側から明示し、化学種の対象経路だけ切り替える**。`speciesImplicitCoupling=0` は別の `scalarTimeIntegration_d` 経路（`:1028`）なので、対応範囲を明記する。R1〜R3 に加えて、β < 1・多成分で刻みの実効値、`ρY` 更新、再正規化前後の収支を確認する必要がある。 |
| **Major** | §6.57 の「前提を満たした」判定：**要再検証** | [rlg_an.py:34](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/rlg_an.py:34) は `q != b+d` を検出しても記録だけで判定不能にしない。`:51` の RHS 比較も「再実行幅以内」ではなく、`max(10×P0再実行幅, 1e−6)` を許容し、P3 の再実行 RHS を検査していない。**前提違反は判定不能にし、両腕の再実行幅と実際の RHS 差を出力する**。今回の結果 A を否定する証拠ではないが、移植後の合格判定にはそのまま使えない。 |
| **Minor** | R3 の `case/08.bump`：**却下** | [verification/README.md:23](/home/sano/work/forge-integ-1005/procedures/verification/README.md:23) は入口 BC の問題と二次精度での残差停滞を理由に、回帰基準として未適格としている。**候補を `case/48.flat_plate_cooled_m4` に絞り、CPG・node・壁条件の実効設定を確認する**。 |

結論: **出口修正を取り込み済みとして移植範囲を訂正し、未取り込みの node 境界処理を現行仕様に合わせて移植したうえで、その修正だけの前後を P3・β = 1 の 1 step で比較する。**

第 1 仮説: node 境界処理の欠落が、入口・出口の更新抑制と壁・軸の拘束違反の主因である。確度: **高**。  
根拠: 標準経路は [timeIntegration_d.cu:962](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:962) で境界半割面の粘性対角を除外し、`:1123` 以降で拘束行を処理する。前処理経路の `:1342–1352` は除外せず、拘束フラグの引数もない。保存済み判定 JSON では、`case/45.isobutane_m6_d155/run_0756_lg_p0a`・`run_0757_lg_p0b` 対 `run_0758_lg_p3a`・`run_0759_lg_p3b` の入口・出口の更新比が、ρ で `2.09e−4`、ρE で `7.73e−4`。軸の最大 |dρv| は `1.40e−4`、壁の最大 |dρE| は `2.69`。ただし、この数値は **JSON の記録を確認したもの**で、生データからの再集計ではない。  
反証条件: 正しいフラグの配線・対象行・同一入力/RHSを確認しても、境界処理だけの修正で更新抑制が再実行幅を超えて改善しない場合、「更新抑制の主因」という説明は棄却する。

第 2 仮説: 修復後に残る P0/P3 の差には、TP 軸ソース Jacobian と solve 精度の違いが寄与する。確度: **中**。コード差は確認済みだが、寄与量は未測定。  
第 3 仮説: 提示しない。現在の証拠から「低マッハ剛性が case/45 の主な収束律速」とまでは判断できない。

判別 A/B: **唯一の変更は P3 の node 境界処理の有無**。同一ソース基準・同一ビルド条件から比較用バイナリを作り、G0・同一 restart・point・FP64・`lowMachPrecond: 3`・`precondEps: 1.0` を固定して各 1 step、各 2 回。化学種側の変更は両腕で同じ状態にする。新規 run に `Q0`、RHS、拘束前の要求更新 `d`、境界フラグを保存する。

- **結果 A**: 拘束対象の最大 |d| が 0 になり、入口・出口の更新抑制も再実行差の 10 倍を超えて改善する → 第 1 仮説を支持。残る P0/P3 差は第 2 仮説の対象とする。
- **結果 B**: 配線・適用を確認済みでも拘束違反が残る、または更新抑制が改善しない → この移植だけで症状を説明できるという仮説を棄却し、試験 P の長時間計算へ進まない。

これは**取り込みの要否ではなく、修復効果を判別する試験**である。1 step の結果から収束・加速は主張しない。

やらない方がよいこと: 出口修正を重ねて実装すること、R1 の差を「出口修正だから」と許容すること、R2 の 1% 超過を直ちに移植失敗とすること、β = 1 だけで実際の前処理と化学種刻みの検証を完了扱いすること。

呼び出し側の前提への異議: **「今回 TP 出口修正が生産値を変える」という前提は、提示されたソース基準では成立しない。** 今回記録すべきものは、取り込み済みコミット、試験バイナリのハッシュとソース対応、修正前後の P0 の全保存量差である。ビット差があれば、同一バイナリの再実行差と区別する。

仮に実際の生産バイナリだけが古いと判明した場合は、別の出口変更として扱う。その際は出口の分岐別面数・面積、γ、法線速度、逆流面積、質量・エネルギー流束、壁熱流束分布と積分値、既存評価量 θ_r を時系列で記録する。**出口変更の影響は内部へ伝播するため、「10 step 後も変更節点の外はビット一致」という条件は使えない。** 生産値の比較には両系列の `check_convergence.py` と `check_quasisteady.py` の VERDICT が必要。

不足情報: `run_0756〜0759` のディレクトリはこの環境に存在せず、`stage.h5`、実効設定、実バイナリとの対応を独立確認できない。確認できた成果物は [rlg_judge.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rlg_judge.json)、台本、判定器、[case README の run 索引](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:225)。**ファイル変更・forge 起動は行っていない。plan 未反映**。呼び出し側で `architecture-float-state-double-geometry.md` §5.1 #26・§6.57 に反映すること。
