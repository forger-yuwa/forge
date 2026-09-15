# codex レビュー: condensation-source-limiter-steady (result)

- **plan**: [`plans/active/condensation-source-limiter-steady.md`](../../plans/active/condensation-source-limiter-steady.md)
- **stage**: `result` (diff base `f1e40f41`)
- **date**: 2026-09-16
- **commit**: `fba3f529` (feature/sern-design)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M1/m1
- **focus**: 5 回目 (result-4 の M1–M3/m4/m5 反映後, 5 巡目バイナリ 355f4e40/cdb4fda5, 記録 fba3f529)。前回指摘の再指摘は不要。集中点: (1) §6 改定ゲート (a)(b)(c) に対する 5 巡目セット case/44 run_0170–0176 の実測 (README「5 巡目」表, plan §9 先頭) — 特に (b) cell の g_0 max 差 2.17e-3 が同一バイナリ反復ノイズ (4.4–5.8e-4) の 3.7–4.9 倍で不成立、壁側液相帯で cfl 0.5 が系統的に −0.41 %, 3 世代 warm 継続で頭打ち、node は同符号だが 4 世代で単調接近中 (1.55e-3)。これを『Δτ 依存の残り (accepted 不可)』とみなすか『§6 ノルム (g L1 1 % / 0.5 K / 5e-3) 内で積分量 0.2 % 一致・node 決定的 → accepted 可 (残件は followups)』とみなすか、根拠つきで判定してほしい。(2) Q0=0 減衰 (cond_evap_source_rate{,_f}) と単体 (k) の妥当性。(3) accepted に移してよいか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **NO-GO**

第5巡の実測値は再現しましたが、主目的に直結するゲート **§6(b) が未達**です。  
積分量の差が0.2%以内でも、液相場に残る系統差を反復ノイズとして扱えません。  
`Q0=0` の減衰処理と単体(k)の構成は妥当ですが、それだけでacceptedには進めません。

## 指摘一覧

### 1. Major — 許容ノルム内でも、Δτ非依存の完了条件を満たしていない

**根拠:** [plan §6(b)](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:137) は、CFL変更による場差を反復ノイズ床の2倍以内と定めています。[第5巡のrun索引・比較表](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:566) の最終HDF5から独立再計算しました。

| 比較 | `max|Δg_0|/max|g_ref|` |
|---|---:|
| cell反復 `0173−0172` | `4.3786e-4` |
| cell反復 `0174−0172` | `4.9928e-4` |
| cell反復 `0174−0173` | `5.7488e-4` |
| cell CFL変更 `0175−0172` | **`2.1741e-3`** |

最大の床を採用しても許容上限は `1.1498e-3`。実測差はその**1.89倍**です。差が `1e-3` を超える299セルは全て低CFL側の液相量が小さく、相対差中央値 **−0.4108%**、位置も記録の壁側液相帯を再現しました。

一方、以下も再現しています。

- node床の末尾／直前比：`0.972–1.008`、CFL間比：`0.865–1.053`。
- 第5巡7本：`check_quasisteady.py` の0.2%判定は **`OVERALL: ALL STEADY`**。
- 同7本：`check_convergence.py` は全て **`NOT CONVERGED (stalled/plateau)`**。
- 全49スナップショットの `VALUE/*` にNaN/Infなし。最終凝縮域の `condLim=1`、補正G/Q=0。

したがって、確認できたのは「報告量が指定窓で定常、クランプ非作動」です。**残差式からΔτを除去したことと、float32反復で同じ場に到達することは別**です。nodeの再現性も、異なるCFL間の固定点一致を証明しません。

ここから「ソース残差にΔτが再混入した」とは断定しません。しかし、限界サイクル平均のCFL依存であっても、現在の完了条件には不合格です。結果を見てノルム条件だけに切り替える根拠はありません。

**対案:** §5.1 #5を未完了のまま維持し、同一保存場からCFL 2／0.5の交差restartを行ってください。壁側液相帯の時系列、輸送・ソース・合計の全モーメント残差を比較し、初期場履歴と更新丸めの寄与を切り分け、§6(b)を満たしてから完了とするべきです。[推奨設定の「定常解が依存しなくなった」](/home/sano/work/forge/procedures/recommended-settings.md:120) も、現時点では「ソース残差の明示的Δτ依存を除去」に限定してください。

### 2. Minor — 新しい`Q0=0`復旧規則が現在仕様に記載されていない

**根拠:** [double実装](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSource_d.cuh:334) と [float実装](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:243) は、`q0<=1e-30` の場合、

`Sg = 3 ρg ṙ(r_min)/r_min`

で質量を減衰させます。しかし、[現在仕様](/home/sano/work/forge/methods/condensation.md:806) は通常の `Sg=4πρ_l q2 ṙ` と消滅条件のみを説明しています。`Q0=Q1=Q2=0, g>0` では、この二式の結果は異なります。

新分岐の符号・単位・Δτ非依存性は整合しています。ただし、これは不整合状態からの**復旧規則**であり、液滴分布から導いた蒸発速度ではありません。

**対案:** `methods/condensation.md` §4cへ、適用条件、減衰式、`r_min`への依存、`g<=g_rm`で消滅処理へ渡すことを明記してください。

## 実装・回帰の確認

- 指定diffを取得し、瞬間速度形、floor前の4増分同率制限、`volumePartial_d`、平衡形の除外、RK／dual-time降格を確認しました。既定変更はplanに記載されています。
- [単体(k)](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_limiter_steady.cu:323) は、2状態×double/floatについて実ソース→更新→消滅を通し、負のソースと全モーメント消滅を検査しています。対象の停止不具合を検出する構成として妥当です。GPUアクセスが遮断されており、**17／20ステップの独立再現はできていません**。
- `case/34.arthur_n2_nozzle/run_0105–0107` 配下の `lim1e*` は、各記録の参照・ノイズ条件で**場比較28/28 PASS**を再現しました。収束判定は別で、`NOT CONVERGED`です。
- `case/16.nozzle_wys/run_0470_limiter_regress_wys/lim1e,lim0e` は、`g`場差 `8.1963e-6`、series **ALL STEADY**を再現しました。全保存量の収束判定は両方 **NOT CONVERGED**です。
- `case/09.Taylor-Green/run_0064_periodic_seam_cond_source_final3/` の `Q0/dt` 比は `0.99999979–1.00000015`。
- §5.1にはcell未達が残り、F-cf8／F-cf9も登録済みです。`plans/README.md`のactive配置と`methods/index.md`の既存リンクは妥当です。

## 推奨

**`in_progress`を維持し、cellの系統差をこのplan内で切り分けてください。** 現在の主目的を未達のままfollowupsへ移し、acceptedにすることには反対します。

ファイルは変更していません。指摘・推奨は **plan未反映**です。

指摘数: Critical 0 / Major 1 / Minor 1
