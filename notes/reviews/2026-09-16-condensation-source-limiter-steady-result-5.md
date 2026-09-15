# codex レビュー: condensation-source-limiter-steady (result)

- **plan**: [`plans/active/condensation-source-limiter-steady.md`](../../plans/active/condensation-source-limiter-steady.md)
- **stage**: `result` (diff base `f1e40f41`)
- **date**: 2026-09-16
- **commit**: `a86cfde7` (feature/sern-design)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m3
- **focus**: 6 回目 (result-5 M1/m1 反映後)。前回までの再指摘は不要。集中点: (1) result-5 M1 の交差 restart (case/44 run_0177–0180, README「交差 restart」表, plan §9 先頭): cfl 2 の収束場を cfl 0.5 で restart しても動かない (cell g_0 max 6.2e-4 ≤ 床 2 倍 1.15e-3, node 反復床) こと、cfl 0.5 系の場を cfl 2 で restart すると 6000–8000 step で cfl 2 の場に収束することから『cfl 0.5 系の系統差は低 cfl 側の遅い過渡 (履歴) で固定点は Δτ 非依存』とした判断と、§6 (b) を『同一場からの交差 restart』で判定する手順への確定の妥当性。輸送/ソース分離の残差は計測していない (交差 restart で履歴寄与を直接示したため) — これで足りるか。(2) m1 の methods §4c 追記。(3) accepted に移してよいか。残件は followups F-cf10 (低 cfl 側の減衰の遅さの原因)。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **GO-with-changes**

交差 restart の実測は再現し、cell の **§6(b) ノイズ床ゲートは達成**しています。  
`methods` §4c の復旧規則追記も妥当です。  
移管前に、因果関係の断定・node の到達時刻・残作業表を修正してください。

## 指摘一覧

### 1. Minor — 「低CFL側の遅い過渡と確定」は交差 restart だけでは言い切れない

**根拠:** [plan §9](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:186) と [README](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:583) は、差を過渡と確定しています。

実測が示すのは、次の二点です。

- CFL 2 の保存場は、CFL 0.5 でも指定ノイズ床内に留まる。
- 低CFL系列の保存場は、CFL 2 に変更すると基準場の近傍へ移る。

これは**共通の場を指定精度内で維持できる証拠**です。しかし、低CFL側の丸めによる停滞まで排除していません。[更新確定処理](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationUpdateLimiter_d.cuh:69) は増分加算後に `flow_float` へ丸めるため、`θ_u>0` でも保存値が変わらない場合があります。4本とも正式な収束判定は `NOT CONVERGED` です。

**対案:** §6(b) の交差 restart 手順は採用し、結論を次に限定してください。

> ソース残差の明示的なΔτ依存を除去し、case/44 の交差 restart で共通の最終場を指定ノイズ床内で維持・再現した。低CFL系列の履歴差が残る機構は未同定。

**輸送・ソース分離残差の追加計測は、この実測ゲートの合格には不要**です。原因を確定する計測として、更新前増分と丸め後変化量の比較も含めて [F-cf10](/home/sano/work/forge/plans/active/condensation-followups.md:62) に残すのが妥当です。

### 2. Minor — node の「8000 stepで収束」は時系列と合わない

**根拠:** [README の `run_0180` 行](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:593)。`run_0170` に対する正規化最大差を再計算すると、

| step | `max|Δg_0|/max|g_ref|` |
|---:|---:|
| 8000 | `4.6242e-4` |
| 10000 | `6.9118e-4` |
| 18000 | `4.3214e-4` |
| 24000 | `1.2692e-4` |

node の反復床は `2.0325e-4` です。8000 stepでは床の2倍も超えており、その後に差が再増加しています。最終値が床内という結論は再現しました。

**対案:** 「8000 stepで収束」を「8000 stepまでに差の大部分が減少し、24000 stepの最終場は反復床内」に修正してください。plan §9 と F-cf10 の「6000–8000 stepで消える」も同期させてください。

### 3. Minor — §5.1 が今回の決着・後続課題と同期していない

**根拠:** [残作業表 #6](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:128) は依然「5回目の論点」を現在の残件として記載し、複数世代の「最終セット」が混在しています。F-cf10 は後続planに存在しますが、この表から辿れません。

**対案:** #5を今回の交差 restart の結果で閉じ、#6を第6回レビュー対応に更新してください。F-cf10へのポインタを追加し、過去世代の経緯は変更ログへ集約してください。

## 検証根拠

対象は `case/44.vitiated_air_wt/` 配下の [run索引・交差 restart 表](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:582) です。

| 比較 | 最終正規化最大差 | 確認結果 |
|---|---:|---|
| [`0177 ← 0172`](/home/sano/work/forge/case/44.vitiated_air_wt/run_0177_va3_M4.19_Lc8_noneq_inletTt_cell_lim1e_xr_cfl05_from0172/) | `6.2411e-4` | 全保存時点の最大も `9.6621e-4` |
| [`0178 ← 0175`](/home/sano/work/forge/case/44.vitiated_air_wt/run_0178_va3_M4.19_Lc8_noneq_inletTt_cell_lim1e_xr_cfl2_from0175/) | `6.0213e-4` | 6000 step以降、cell許容範囲内 |
| [`0179 ← 0170`](/home/sano/work/forge/case/44.vitiated_air_wt/run_0179_va3_M4.19_Lc8_noneq_inletTt_lim1e_xr_cfl05_from0170/) | `8.7069e-5` | node反復床内 |
| [`0180 ← 0171`](/home/sano/work/forge/case/44.vitiated_air_wt/run_0180_va3_M4.19_Lc8_noneq_inletTt_lim1e_xr_cfl2_from0171/) | `1.2692e-4` | 最終場はnode反復床内 |

cell許容上限は、独立再計算した最大反復床の2倍、**`1.1498e-3`**です。最終の `ro/P/T`・化学種・全モーメント保存量も床の2倍以内でした。

- restart入力の保存量11本は元の `res_24000.h5` と完全一致。メッシュ・BC・物性も一致しました。
- 全52スナップショットの `VALUE/*` に NaN/Infなし。最終凝縮域は `condLim=1`、補正G/Qとも0。
- 正式ツール再実行：4本とも **`NOT CONVERGED (stalled/plateau)`**、0.2%の系列判定は **`OVERALL: ALL STEADY`**。差が大きかった壁側領域の平均gも同じ閾値で `ALL STEADY`。
- 指定diffの瞬間速度形、floor前の同率制限、部分体積、平衡形除外、RK／dual-time降格を確認。新たな実装上の阻害事項は確認していません。
- [methods §4c](/home/sano/work/forge/methods/condensation.md:813) は復旧規則の条件・式・`r_min`依存・消滅への移行を記載しており、前回m1は解消しています。
- Wysłouzilの最終場差 `g_0=8.1963e-6`、`ro=1.2156e-6`を再現しました。GPUアクセスが遮断されており、CUDA単体試験・性能値の独立再実行はしていません。

## 推奨

**上記1→2→3の文書修正を行い、今回のレビューと採否を記録してから `accepted` へ移してください。** 受け入れる成果は、ソース残差のΔτ依存除去と、指定精度内での交差 restart 再現性です。低CFL側の停滞機構はF-cf10で追跡してください。

ファイルは変更していません。指摘・推奨は **plan未反映**です。

指摘数: Critical 0 / Major 0 / Minor 3
